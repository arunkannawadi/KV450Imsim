import numpy

def extract_header(pathname):
    header = ['#Garbage']
    with open(pathname,'r') as f:
        while header[-1][0]=='#':
            header.append(f.readline())

    ## Get rid of the first and last line
    fline = header.pop(0)
    lline = header.pop(-1)

    return header

def generate_header(current_psf_set, psfset_path, include_pointing_error=False, n_exposures=5, random_seed=None):
    header = [ ]
    # Random seed
    if random_seed==None:
        random_seed = numpy.random.randint(10000,1000000)*17-137
    else:
        if type(random_seed)==str:
            random_seed_from_file = open(random_seed,'r')
            random_seed = long(random_seed_from_file.readline().split()[-1])
            random_seed_from_file.close()
        else:
            pass
    header.append('Random seed: '+str(random_seed)+'\n')
    numpy.random.seed(random_seed)

    # Dithers and pointing errors
    dither_x=numpy.ceil(25/0.214)
    dither_y=numpy.ceil(85/0.214)

    # Obtain the pointing error regardless of whether we include them or not
    dither_pointing_error = numpy.round(numpy.random.normal(0,0.5,2*(n_exposures-1)),decimals=2)
    if not include_pointing_error:
        dither_pointing_error = numpy.zeros(8)

    header.append('Dither 0 pointing error:  0  0\n')
    for i in range(n_exposures-1):
            header.append('Dither %d + pointing error:  %f  %f \n' %(i+1,
                                            (i+1)*dither_x+dither_pointing_error[2*i],
                                            (i+1)*dither_y+dither_pointing_error[2*i+1]))

    # PSF parameters
    # PSF Moffat parameters measured by Reiko Nakajima using PSFEx on KiDS fields
    # Specifically chosen to be a nicely distributed collection of PSF FWHM
    psf_filename=psfset_path
    psf_field,psf_exp = numpy.loadtxt(psf_filename, usecols=(0,1),dtype=str,unpack=True)
    psf_chip = numpy.loadtxt(psf_filename, usecols=(2,),dtype=int,unpack=True)
    psf_fwhm,psf_ell1,psf_ell2,psf_beta = numpy.loadtxt(psf_filename, usecols=(3,4,5,6),dtype=float,unpack=True)

    # There are some catalog entries that have zeros everywhere (55 out of total).
    # These will crash the simmaker script, so need to be removed. If any of the 5 exposures is
    # going to have such a flawed line, try another random catalog line.
    for i in range(5):
            psf_ind = int(current_psf_set)*5+i
            header.append('PSF dither %d from field %s, exposure %s, chip %d, num %d\n'
                                            % (i,psf_field[psf_ind],psf_exp[psf_ind],psf_chip[psf_ind],0))
            header.append('PSF parameters, fwhm(arcsec) %f beta %f e1 %f e2 %f\n'
                                            % (psf_fwhm[psf_ind],psf_beta[psf_ind],psf_ell1[psf_ind],psf_ell2[psf_ind]))
    header.append('# (use g for PSF ellipticity in GalSim!)\n')

    # Parameter explanation
    header.append('# x  y   m   r   e1   e2   f    n')

    return header

def create_COSMOS_priors(current_psf_set, output_pathname, path_to_bsc_model, psfset_path, \
                         randomize_positions=True, randomize_orientations=False, scramble_mod_e=False, use_scrambled_e=False,
                         include_stars=True, include_pointing_error=False, random_seed=None, n_exposures=5, real_galaxy=False):

    header_list = generate_header(current_psf_set, psfset_path, include_pointing_error=include_pointing_error, random_seed=random_seed, n_exposures=n_exposures)
    header = ''
    for hd in header_list:
        header += hd

    ## ** Image specs ** $$
    image_size_x=16810
    image_size_y=16530
    dither_x=numpy.ceil(25/0.214)
    dither_y=numpy.ceil(85/0.214)
    buffer_pixels = 3
    all_image_size_x=(image_size_x+(n_exposures-1)*(dither_x+buffer_pixels))		# Rectangle containing all 5 exposures
    all_image_size_y=(image_size_y+(n_exposures-1)*(dither_y+buffer_pixels))		# with extra padding of 3 pixels
    
    ## *** Stars *** ##
    star_mag = numpy.loadtxt(path_to_bsc_model,unpack=True)
    n_stars = len(star_mag)
    star_x = numpy.random.uniform(-buffer_pixels,all_image_size_x, n_stars).T
    star_y = numpy.random.uniform(-buffer_pixels,all_image_size_y, n_stars).T
    star_e1 = numpy.zeros_like(star_mag)
    star_e2 = numpy.zeros_like(star_mag)
    star_r = numpy.zeros_like(star_mag)
    star_f = -numpy.ones_like(star_mag)
    star_n = numpy.zeros_like(star_mag)
    star_ZB4 = numpy.zeros_like(star_mag)
    star_ZB9 = numpy.zeros_like(star_mag)
    star_ID = -numpy.ones_like(star_mag)

    star_dat = numpy.array([star_x,star_y,star_mag,star_r,star_e1,star_e2,star_f,star_n,star_ZB4,star_ZB9,star_ID],)

    ## ** Galaxies ** ##
    if use_scrambled_e:
        input_truth_catalog = '/disks/shear15/KiDS/ImSim/pipeline/data/prior_reduced_all_scrambled_sersic_ZB9_OBJNO_q'
        if real_galaxy is True:
            print "'real_galaxy=True' and 'use_scrambled_e=True' are incompatible. Ignoring 'real_galaxy=True' condition."

    else:
        input_truth_catalog = '/disks/shear15/KiDS/ImSim/pipeline/data/prior_reduced_all_sersic_ZB9_OBJNO'
        #input_truth_catalog = '/disks/shear15/KiDS/ImSim/pipeline/data/reduced_KiDS_Griffith_iMS1_testing.cat'
        if real_galaxy is True:
            input_truth_catalog += '_real'

    gal_dat = numpy.loadtxt(input_truth_catalog, comments='#')

    ## Remove (mock) stars in the catalog
    cuts = numpy.abs(gal_dat[:,-4])>=0.3 ## TODO: Change it to a positive index!!!

    gal_dat = gal_dat[cuts]

    nobj = len(gal_dat)

    ## Randomize the positions
    if randomize_positions:
        new_x = numpy.random.uniform(-buffer_pixels,all_image_size_x,nobj)
        new_y = numpy.random.uniform(-buffer_pixels,all_image_size_y,nobj)

        gal_dat[:,0] = new_x
        gal_dat[:,1] = new_y

    if randomize_orientations:
        if randomize_positions:
            print "You have already opted to randomize positions. Randomizing orientation might not add much."
        else:
            print "You have chosen to randomize orientations. Any IA signal present will be lost."
        import cmath
        for gg in xrange(len(gal_dat)):
            cmplx_e = complex( gal_dat[gg,4], gal_dat[gg,5] )
            mod_e, true_phase = cmath.polar(cmplx_e)
            random_phase = numpy.random.uniform(0,numpy.pi) ## Phase must be in radians
            rndm_e = cmath.rect(mod_e, random_phase)
            gal_dat[gg,4], gal_dat[gg,5] = rndm_e.real, rndm_e.imag

    if scramble_mod_e:
        if use_scrambled_e:
            print "You are using scrambled ellipticity prior. Randomizing mod_e will only add noise."
        if real_galaxy:
            print "You are using 'real_galaxy=True'. I cannot scramble the ellipticities"
        else:
            new_indices = numpy.random.permutation(len(gal_dat))
            new_g1 = gal_dat[:,4][new_indices]
            new_g2 = gal_dat[:,5][new_indices]

            gal_dat[:,4] = new_g1
            gal_dat[:,5] = new_g2


    print star_dat.shape, gal_dat.shape
    if include_stars:
        dat = numpy.append(star_dat.T, gal_dat, axis=0)
    else:
        dat = gal_dat

    print star_dat.shape, gal_dat.shape, dat.shape
    numpy.savetxt(output_pathname, dat, header=header, fmt=['%.2f','%.2f','%.2f','%.3f','%.4f','%.4f','%.3f','%.3f','%.2f','%.2f','%d'])

def create_priorfile(current_psf_set, path_to_bsc_model, output_pathname, psfset_path, grid_positions=None, nr_of_stars=True, nr_of_bright_gals=True, nr_of_faint_gals=True,
        random_seed=None, randomize_positions=True, randomize_orientations=False, scramble_mod_e=False, use_scrambled_e=False, include_stars=True, include_pointing_error=False, n_exposures=5, realistic=True, real_galaxy=False):

        if realistic:
            print "Ignoring the keywords: grid_positions, nr_of_stars, nr_of_bright_gals, nr_of_faint_gals"
            create_COSMOS_priors(current_psf_set, output_pathname, path_to_bsc_model, psfset_path, randomize_positions, randomize_orientations,
                                    scramble_mod_e, use_scrambled_e, include_stars, include_pointing_error, random_seed, n_exposures, real_galaxy)

        else:
            print "Ignoring the keywords: randomize_positions, randomize_orientations, scramble_mod_e, use_scrambled_e, include_stars, include_pointing_error, n_exposures, real_galaxy"
            import priors_fixed_psf2_170324 as pfp

            pfp.create_priorfile(current_psf_set, path_to_bsc_model, output_pathname, psfset_path, grid_positions, nr_of_stars, nr_of_bright_gals, nr_of_faint_gals,
                                    random_seed=random_seed)

#if __name__=='__main__':
    ### For testing purpose only
    #create_COSMOS_priors(1, '/disks/shear15/KiDS/ImSim/pipeline/data/mock_galaxy_prior.dat', '/disks/shear15/KiDS/ImSim/pipeline/utils/bsc_model_175_0_m18_25',\
    #                    '/disks/shear15/KiDS/ImSim/pipeline/utils/65_psfs.txt')


