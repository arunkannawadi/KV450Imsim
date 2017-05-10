## Image simulations for KiDS
## Author: Arun Kannawadi

## To enter debugging mode, set the `debug' variable to be True
debug = False
diagnose = False
if debug:
    import pdb; pdb.set_trace()
import galsim
import math
import numpy as np
from astropy.io import fits
import pyfits
import time
import logging
import os, sys
import multiprocessing
from multiprocessing import Pool
from multiprocessing import Process, Queue

## Setup basic logging
#logging.basicConfig(format="%(message)s", level=logging.INFO, stream=sys.stdout) # if you want it on the screen
logging.basicConfig(filename="imsim.log",level=logging.INFO)
logger = logging.getLogger("imsim")

## Which bandpass are we simulating?
filter_name = 'r'

## Should we include shot noise?
include_shot_noise = False

## How many exposures do we simulate?
n_tot_exposures = 5

### How many rotations do we do to cancel shape noise?
#n_rotations = 4

#ditherArray = [galsim.PositionD(0.,0.), galsim.PositionD(5.3,-4.2), galsim.PositionD(-6.1,-4.9), galsim.PositionD(-3.6,5.2), galsim.PositionD(7.5,5.5)]
#ditherArray = [galsim.PositionD(0.,0.), galsim.PositionD(25/0.214,85/0.214), galsim.PositionD(2*25/0.214,2*85/0.214), galsim.PositionD(3*25/0.214,3*85/0.214), galsim.PositionD(4*25/0.214,4*85/0.214)]

## Load the combined input catalogue
input_cat_name = '/disks/shear14/arunkannawadi/catalogue_matching/KiDS_unmasked_Griffith_iMS1.cat'
input_catalogue = fits.open(input_cat_name)
input_data = input_catalogue[1].data

output_image_name = '/disks/shear14/arunkannawadi/imsim/canvas.fits'

## Should the invidual postage stamps be saved?
save_stamps = False
if save_stamps is True:
    print "Warning! save_stamps option is on. This would take a considerable amount of storage!"

## If saving the stamps, where should the postage stamps go?
postage_stamps_dir = '/disks/shear14/arunkannawadi/imsim/postage_stamps/{0}_band/trunc5/'.format(filter_name)

## What should the final image be called (exclude .fits extension)
#full_image_filename = 'short_expo{0}'.format(exp_id)

## What should the stamp records file be called as?
#stamp_records_filename = '/disks/shear14/arunkannawadi/imsim/short_expo{0}.fits'.format(exp_id)

## Setup the basic KiDS parameters
pix_scale = 0.214 ## arcsec per pixel
hst_cosmos_pixel_scale = 0.05
kids_noise_level = 0.0185095*0.9*np.sqrt(5) ## fudge factor 0.9 and account for co-adding
magAB_zeropoint = 24.66
#exp_time = 1.0/1800 #TNR
#gain = 3331.08077604 #TNR

## Random Number Generator
base_seed = 1234567

## Now some GalSim parameters
gsp = galsim.GSParams(maximum_fft_size=16384)

## Sensible values for the parameters
n_min, n_max = 0.3, 6.2 # limits set by GalSim
q_min, q_max = 0.05, 1.0
trunc_factor = 5 ## trunc_factor = 0 implies no truncation!

## Obtain the WCS
wcs_filename = "/disks/shear14/KiDS_simulations/Cosmos/Theli_image/KIDS_150p1_2p2_r_SDSS.V0.5.9A.swarp.cut.fits"
wcs = galsim.AstropyWCS(wcs_filename)
#wcs = galsim.AstropyWCS("/disks/shear14/KiDS_simulations/Cosmos/AW_images/KIDS450_150.1_2.2_r_sci.fits")
canvas_bounds = galsim.fits.read(wcs_filename).bounds


### Default cut
#default_cuts = (mask)&(rank>0)&(distance2d<0.5)&(weight>5)
#print "Default cuts have", default_cuts.sum(), " galaxies."

#if not 'KiDSCatalogue' in locals():
#    print "Loading KiDS catalogue ... "
#    #KiDSCatalogue = fits.open("/disks/shear10/KiDS/All_tiles/KIDS_150p1_2p2_r_SDSS.V0.5.9A_photoz_recalib_LF_mask_SG.cat.blind.cat")
#    #KiDSCatalogue = fits.open("/disks/shear14/KiDS_simulations/Cosmos/Theli_cat/KIDS_150p1_2p2_r_SDSS.V0.5.9A_photoz_recalib_LF_mask_SG.cat")
#    #KiDSCatalogue = fits.open('/disks/shear14/arunkannawadi/catalogue_matching/Griffith_COSMOS_noKiDS_iMS1.fits')
#    KiDSCatalogue = fits.open('/disks/shear14/KiDS_simulations/Cosmos/KIDS_HST_cat/KiDS_Griffith_iMS1.fits')
#    ## Apply mask before doing anything
#    try:
#        MASK_all = KiDSCatalogue[2].data.MASK
#        rank = KiDSCatalogue[2].data['rank']
#        mask = ~(np.array(MASK_all&0xfc3c,dtype=bool))
#    except:
#        mask = np.zeros(KiDSCatalogue[1].data.shape[0],dtype=bool)
#        print "Masked not obtained"
# 
#    KiDSData = KiDSCatalogue[2].data#[mask]
#
#    try:
#        RA2 = KiDSData.field("RA_THELI")
#    except KeyError:
#        RA2 = KiDSData.field('RA')
#
#    try:
#        DEC2 = KiDSData.field("DEC_THELI")
#    except KeyError:
#        DEC2 = KiDSData.field('DEC')
#
#    try:
#        MAG_AUTO2 = KiDSData.field("MAG_AUTO_THELI")
#    except KeyError:
#        MAG_AUTO2 = KiDSData.field('MAGR')
#    try:
#        LF_weights = KiDSData.field('weight')
#    except KeyError:
#        LF_weights = np.ones_like(MAG_AUTO2)
#
#    print "KiDS catalogue loaded"
#

## Routines
def getKiDSPSFs(psfset_id, exp_id, psf_params_pathname):
    """ Generate the Moffat PSF models

        @param psfset_id                    PSFSet ID, ranging from 0 to 12 (inclusive).
        @param exp_id                       Exposure ID, ranging from 0 to 5 (inclusive).
        @param psf_params_filename          The file where to read the Moffat parameters from. The
                                            format of the file should be such that 13 sets of 5
                                            exposures are group together, separated by an empty line.
                                            The lines must be in the following format:
                                            <blah blah blah> fwhm psf_e1 psf_e2 moffat_beta
                                            The 'e's are assumed to be 1 -q

        @returns                            A GalSim Moffat instance for the PSF
    """

    psf_params_file = open(psf_params_pathname,'r')
    psf_params_lines = psf_params_file.readlines()
    psf_params_file.close()

    psf_params = psf_params_lines[6*psfset_id+exp_id].split()
    ## 6 because the grouping is 5 PSF parameters + 1 empty line

    moffat_beta = float(psf_params[-1])
    psf_e2 = float(psf_params[-2])
    psf_e1 = float(psf_params[-3])
    fwhm = float(psf_params[-4]) #arcsec

    ## e = 1 - q
    psf_e = np.sqrt(psf_e1**2+psf_e2**2)
    ## g = (1-q)/(1+q)
    psf_g1,psf_g2 = psf_e1/(2-psf_e), psf_e2/(2-psf_e)

## HACK or REAL ?????????????????????????????????????????????????????????????????????????????
    psf_g1, psf_g2 = psf_e1, psf_e2

    psf = galsim.Moffat(beta=moffat_beta,fwhm=fwhm,trunc=4.5*fwhm).shear(g1=psf_g1,g2=psf_g2)

    return psf


## Get the PSF - very simplistic PSF for now
def getKiDSPSF(filter_name='r'):
    import warnings
    warnings.warn("This method is depreciated. Discontinue calling this routing ASAP.")

    ## Check if the bands are valid
    all_bands_set = set(['u','g','r','i'])
    if not set(filter_name).issubset(all_bands_set):
        raise ValueError("The 'bands' must be a subset of 'u','g','r','i' only. ")

    if filter_name=='r':
        moffat_beta = 2.41
        fwhm = 0.5 ## arcsec

        ## e = 1 - q
        psf_e1 = 0.08
        psf_e2 = -0.02
        psf_e = np.sqrt(psf_e1**2+psf_e2**2)
        psf_g1,psf_g2 = psf_e1/(2-psf_e), psf_e2/(2-psf_e)

        psf = galsim.Moffat(beta=moffat_beta,fwhm=fwhm).shear(g1=psf_g1,g2=psf_g2)
    
    return psf

def getStarImages(PSF,star_catalogue, ditherArray, exp_id):
    """ Obtain a canvas image containing stars.

        @param PSF          The PSF which appears as the star in the field
        @star_catalogue     The star catalogue containing X,Y positions and magnitudes.
                            This could either a FITS format and .npy file.
        @param ditherArray  Dither Array
        @param exp_id       Takes one of 0,1,2,3,4 [default: 0]

        @returns            A canvas image with the PSF images stitched.

    """

    ## Read in the star catalogue
    if isinstance(star_catalogue,str):
        path_star_catalogue = star_catalogue
        try:
            star_catalogue = fits.open(path_star_catalogue)
            star_dat = star_catalogue[1].data
        except:
            star_dat = np.load(path_star_catalogue)

    x_image = star_dat['X_IMAGE']
    y_image = star_dat['Y_IMAGE']
    star_mag = star_dat['MAG']
    star_flux = 10.**(-0.4*(star_mag-magAB_zeropoint))

    ## Set up an empty canvas
    star_field = galsim.Image(canvas_bounds,init_value=0.)

    for star_id in xrange(len(star_dat)):

	## Compute the integer and fractional offsets
	x_nominal = x_image[star_id] + 0.5 - ditherArray[exp_id,1]
	y_nominal = y_image[star_id] + 0.5 - ditherArray[exp_id,2]
	ix_nominal = int(math.floor(x_nominal+0.5))
	iy_nominal = int(math.floor(y_nominal+0.5))
	dx = x_nominal - ix_nominal
	dy = y_nominal - iy_nominal

	star_offset = galsim.PositionD(dx,dy) ## full offset

        star = PSF.copy()
        star = star.withFlux(star_flux[star_id])

        stamp = star.drawImage(scale=pix_scale,offset=star_offset)

        stamp.setCenter(ix_nominal, iy_nominal)

        stitch_flag = stitch_stamps(star_field, stamp, min_pixel_value=-1e-2, max_flux=np.inf)

    return star_field

def getPostageStamps(hst_indices, ditherArray, galaxy_dat, psfset_id, psf_params_pathname, exp_id, rot_id, g1, g2, n_rotations=4, sersic_only=True, fixed_stamp_size=None):
    """ Obtain the postage stamps for a given list of indices.

        @param hst_indices          A list of row indices that need to be simulated from the input catalogue.
        @param psfset_id            The ID of the CCD chip to determine the PSF.
        @param exp_id               Takes one of 0,1,2,3,4
        @param rot_id               Which of the n_rotations are we simulating?
        @param n_rotations          Number of rotations to use for shape noise cancellation.
        @param g1                   The true weak lensing shear g1 [default: 0.0]
        @param g2                   The true weak lensing shear g2 [default: 0.0]
        @param fixed_stamp_size     Whether or not to use a standard postage stamp size. Currently unused.

        @returns                    A canvas image with the postage stamps stitched and a NumPy record array
                                    with the logs.
    """

    print "Entering getPostageStamps ... "

    ## Get rid of the nan s in hst_indices
    hst_indices = hst_indices[~np.isnan(hst_indices)].astype(int)

    ## Get the PSF
    PSF = getKiDSPSFs(psfset_id=psfset_id, exp_id=exp_id,psf_params_pathname=psf_params_pathname) ##
    #PSF = getKiDSPSF(filter_name=filter_name) ## old

    ## Get the galaxy parameters
    ## Griffith parameters
    ## TODO: Move the reading operations in the outer level, preferably to the pipeline
    obj_no = galaxy_dat.field('OBJNO').astype(int) ## OBJNO

    n = galaxy_dat.field('N_GALFIT_HI').astype(float)  ## Sersic index n
    hlr = galaxy_dat.field('RE_GALFIT_HI').astype(float)  ## half-light radius / effective radius
    q = galaxy_dat.field('BA_GALFIT_HI')  ## axis ratio
    beta = galaxy_dat.field('PA_GALFIT_HI') ## position angle

    if sersic_only is False:
        ## Ref: http://great3.jb.man.ac.uk/leaderboard/data/public/COSMOS_25.2_training_sample_readme.txt

        ## Get the generic parameters
        use_bulgefit = galaxy_dat.field('use_bulgefit')
        btt = galaxy_dat.field('fit_dvc_btt')

        bulgefit = galaxy_dat.field('bulgefit')
        ## Get the disk parameters
        disk_hlr = bulgefit[:,1]
        disk_q = bulgefit[:,3]
        disk_beta = bulgefit[:,7] ## in radians

        ## Get the bulge parameters
        bulge_hlr = bulgefit[:,1+8]
        bulge_q = bulgefit[:,3+8]
        bulge_beta = bulgefit[:,7+8] ## in radians

    ## Matching parameters - used for cuts
    rank = galaxy_dat.field('rank') ## matching rank
    ## Other fields that might be used for the cuts
    weight = galaxy_dat.field('weight') ## lensfit weights

    try:
        MASK_all = galaxy_dat.MASK
        mask = ~(np.array(MASK_all&0xfc3c,dtype=bool))
    except:
        mask = np.zeros(galaxy_dat.shape[0],dtype=bool)
        logger.warning("Input catalogue does not have inbuilt masks")

    ## Positions (Theli astrometry when available, else use the Griffith ones)
    RA = galaxy_dat.field('RA_THELI').copy() ## RA
    DEC = galaxy_dat.field('DEC_THELI').copy() ## DEC
    RA[rank==0], DEC[rank==0] = galaxy_dat.field('RA')[rank==0], galaxy_dat.field('DEC')[rank==0]

    ## Use Theli magnitudes when available. Else, use MAGR
    mag = galaxy_dat.field('MAG_AUTO_THELI').copy()
    mag[rank==0] = galaxy_dat.field('MAGR')[rank==0]

    ## Avoid exponentiating magnitudes directly to minimize overflow/numerical errors
    #F0 =  gain*exp_time*(10**(0.4*magAB_zeropoint)) ## for FLUX_GAAP

    ## Flux is zero by default!, except for those with positive magnitudes (many)
    flux = np.zeros_like(mag,dtype=float)
    ##flux[mag>0] = 10**(-0.4*(mag[mag>0])-magAB_zeropoint)
    flux[mag>0] = 10**(-0.4*(mag[mag>0]-magAB_zeropoint))
    #assert flux[rank!=0]==galaxy_dat.field('FLUX_AUTO_THELI')[rank!=0]

#    for mag_id in xrange(len(mag)):
#      #try:
#      if mag[mag_id]>0:
#        flux[mag_id] = 10**(-0.4*(mag[mag_id]-magAB_zeropoint))
#      else:
#      #except:
#        flux[mag_id] = 0.
#
    ### WARNING: Do NOT name the unscaled_flux to flux and multiply itself with F0. I repeat, do NOT !!!
    #flux = (unscaled_flux*10**(0.4*magAB_zeropoint))*exp_time*gain
    #flux = flux.tolist()

    #postage_stamps = []
    stamp_record_list = []

    ## Initiate a canvas image
    canvas_images = galsim.Image(canvas_bounds)

    n_count, n_reset = 0, 500000 ## Report when n_count reaches n_reset
    print "Creating postage stamps now. I will report every {0} galaxies".format(n_reset)

    ## Beginning the loop over galaxies
    for hst_idx in hst_indices:
        t_p1 = time.time()
        if n_count%n_reset==0:
            t1 = time.time()

	failure = 0
        """ FAILURE FLAG BITS:
        Bit 0: Sersic n out of range
        Bit 1: Axis ratio out of range
        Bit 2: Error in generating the circular Sersic profile
        Bit 3: Error in incorporating intrinsic galaxy shape
        Bit 4: Error in applying the input lensing shear
        Bit 5: Error in drawing onto a postage stamp
        Bit 6: Don't care
        Bit 7 (Sign-bit): rank==0, ie, not found in KiDS catalogue
        """

        ## Get the galaxy ID
        gal_id = obj_no[hst_idx]

        ## Get the galaxy flux
        gal_flux = flux[hst_idx]

        ## Get the galaxy position
	gal_ra, gal_dec = RA[hst_idx]*galsim.degrees, DEC[hst_idx]*galsim.degrees
	gal_skypos = galsim.CelestialCoord(gal_ra,gal_dec)
	gal_imgpos = wcs.posToImage(gal_skypos)

	## Compute the integer and fractional offsets
	x_nominal = gal_imgpos.x + 0.5 - ditherArray[exp_id,1]
	y_nominal = gal_imgpos.y + 0.5 - ditherArray[exp_id,2]
	ix_nominal = int(math.floor(x_nominal+0.5))
	iy_nominal = int(math.floor(y_nominal+0.5))
	dx = x_nominal - ix_nominal
	dy = y_nominal - iy_nominal
	gal_offset = galsim.PositionD(dx,dy) ## fractional offset

        ## Get the structural parameters
        if sersic_only is True or (sersic_only is False and use_bulgefit[hst_idx]<=0):
            ## The order above matters, since use_bulgefit is undefined if sersic_only is True

            gal_n = n[hst_idx]
            gal_q = q[hst_idx]
            gal_hlr = np.sqrt(gal_q)*hlr[hst_idx]*hst_cosmos_pixel_scale ## in arcsec
            if n_rotations==1:
                gal_beta = (90+beta[hst_idx])
            else:
                ## Rotate in uniform steps from 0 to 90
                gal_beta = (90+beta[hst_idx]+rot_id*45.) ## in degrees
            ## Griffith PA are wrt y, so add a 90. It shouldn't matter for other catalogues if n_rotations>1

            ## See if the galaxy has sensible parameters and if not, dismiss them
            if (gal_n<n_min) or (gal_n>n_max):
              logger.error("Sersic n for {0} is out of range.".format(hst_idx),exc_info=0)
              failure += 2**0

            if (gal_q<q_min) or (gal_q>q_max):
              logger.error("Axis ratio for {0} is out of range.".format(hst_idx),exc_info=0)
              failure += 2**1

            ## Make the galaxy
            if not failure:
                try:
                    round_gal = galsim.Sersic(n=gal_n,half_light_radius=gal_hlr,flux=gal_flux,trunc=trunc_factor*gal_hlr,gsparams=gsp)
                except:
                    logger.exception("While creating a Sersic profile for {0}, an exception occurred".format(hst_idx),exc_info=0)
                    print "Exception occurred in creating the Sersic profile", hst_idx
                    print gal_n, gal_hlr, gal_flux
                    failure += 2**2

            if not failure:
                try:
                    gal = round_gal.shear(q=gal_q,beta=gal_beta*galsim.degrees)
                except:
                    logger.exception("While incorporating the intrinsic shape for {0}, an exception occurred".format(hst_idx),exc_info=0)
                    failure += 2**3

        else: ## if sersic_only is False and use_bulgefit[hst_idx] is True
            ## Assign so that stamp record can have these values
            gal_n, gal_beta, gal_q, gal_hlr = -9, -99, -9, -9

            ## Get the bulge component as DeVaucoulers
            gal_bulge_q = bulge_q[hst_idx]
            gal_bulge_hlr = np.sqrt(gal_bulge_q)*bulge_hlr[hst_idx]*hst_cosmos_pixel_scale #in arcsec
            if n_rotations==1:
                gal_bulge_beta = bulge_beta[hst_idx] ## in radians
            else:
                ## Rotate in uniform steps from 0 to pi/2
                gal_bulge_beta = bulge_beta[hst_idx]+rot_id*0.5*np.pi/(n_rotations-1) ## in radians
            try:
                round_bulge = galsim.DeVaucouleurs(half_light_radius=gal_bulge_hlr,trunc=trunc_factor*gal_bulge_hlr,gsparams=gsp)
            except:
                logger.exception("While creating a DeVaucouleurs profile for {0}, an exception occured.".format(hst_idx),exc_info=0)
                print gal_bulge_hlr
                failure = failure | 0x04

            if not failure:
                try:
                    bulge = round_bulge.shear(q=gal_bulge_q,beta=gal_bulge_beta*galsim.radians)
                except:
                    logger.exception("While incorporating the intrinsic shape to DeVaucouleurs for {0}, an exception occurred.".format(hst_idx),exc_info=0)
                    failure = failure | 0x08

            ## Get the disk component as Exponential
            gal_disk_q = disk_q[hst_idx]
            gal_disk_hlr = np.sqrt(gal_disk_q)*disk_hlr[hst_idx]*hst_cosmos_pixel_scale #in arcsec
            if n_rotations==1:
                gal_disk_beta = disk_beta[hst_idx] ## in radians
            else:
                ## Rotate in uniform steps from 0 to pi/2
                gal_disk_beta = disk_beta[hst_idx]+rot_id*0.5*np.pi/(n_rotations-1) ## in radians
            try:
                round_disk = galsim.Exponential(half_light_radius=gal_disk_hlr,gsparams=gsp)
            except:
                logger.exception("While creating an Exponential profile for {0}, an exception occurred.".format(hst_idx),exc_info=0)
                print gal_disk_hlr
                failure = failure | 0x04

            if not failure:
                try:
                    disk = round_disk.shear(q=gal_disk_q,beta=gal_disk_beta*galsim.radians)
                except:
                    logger.exception("While incorporating the instrinsic shape to Exponential for {0}, an excetion occurred.".format(hst_idx),exc_info=0)
                    failure = failure | 0x08

            if not failure:
                gal_btt = btt[hst_idx]
                # Add the bulge and the disk to get the galaxy profile
                gal = gal_flux*((gal_btt)*bulge+(1-gal_btt)*disk)


        if not failure:
            ## Apply the lensing shear to the galaxy
            try:
                gal = gal.shear(g1=g1,g2=g2)
            except:
                logger.exception("While applying the lensing shear for {0}, an exception occurred.".format(hst_idx),exc_info=0)
                failure += 2**4

        if not failure:
            ## Convolve with the PSF
            gal_conv = galsim.Convolve([gal,PSF])

        if not failure:
            try:
                t_g1 = time.time()
                draw_method = 'phot' if include_shot_noise else 'auto' ## potentially a flux-based choice. Keep it here.
                #stamp = gal_conv.drawImage(wcs=wcs.local(gal_imgpos),offset=gal_offset, method=draw_method)
                stamp = gal_conv.drawImage(scale=pix_scale, offset=gal_offset, method=draw_method)
                t_g2 = time.time()
                t_g = t_g2 - t_g1 ## time to draw the galaxy

                stamp_size = max(stamp.array.shape) ## in case, if stamp is rectangular
                #with open("drawImage_log.txt","a") as f:
                #    f.write("{0} {1}\n".format(hst_idx,t_g))
            ## let it choose a size automatically
            except (RuntimeError, MemoryError) as er:
                failure += 2**5
                logger.error("While drawing the image for {0}, an error occurred".format(hst_idx), exc_info=0)
                print type(er), " occurred while drawing ", hst_idx
                print "Continuing"
                sys.stdout.flush()
                #with open("drawImage_log.txt","a") as f:
                #    f.write("{0} -1\n".format(hst_idx))
                t_g = -1
                #continue

        if not failure:
            ## Calculate the flux captured
            flux_fraction = stamp.array.sum()/gal_flux
        else: 
            t_g = -1
            flux_fraction = -1
            stamp_size = -1

        ## Mark if the magnitudes correspond to the Theli or the HST catalogue
        if rank[hst_idx]==0:
            failure += 2**7

        if not failure:
            ## Recenter the stamp
            stamp.setCenter(ix_nominal,iy_nominal) ## integral offset

            ### Appending it to the 'postage_stamps' list
            #postage_stamps.append( stamp )

            ## Stitch the postage stamp to the canvas
            stitch_flag = stitch_stamps(canvas_images, stamp)
        else:
            stitch_flag = 0

        ## Append the count, if the galaxy was attempted to be drawn
        if failure&(0x10) | (not failure):
            n_count += 1 ## increment n_count to account for the time taken to attempt drawing the galaxy

        ## Save the stamps, if set
        if save_stamps and not failure:
            postage_stamp_filename = postage_stamps_dir+'{1}.fits'.format(filter_name,gal_id)
            ## It should be stored with gal_id/obj_no since the association might change with another catalogue matching
            stamp.write(postage_stamp_filename,overwrite=False)

        t_p2 = time.time()
        t_p = t_p2-t_p1 ## Time to process the galaxy

        ## Create a stamp record
        stamp_record = (gal_id, gal_n, gal_hlr, gal_q, gal_beta, gal_flux,
                        x_nominal, y_nominal, stamp_size, flux_fraction, t_g, t_p, failure, stitch_flag)

        ## Appending it to the list of stamp records
        stamp_record_list.append( stamp_record )

        if n_count%n_reset==0:
            t2 = time.time()
            logger.info("Progress report - Finished simulating {0} galaxies".format(hst_idx+1))
            print "Progress report - Finished simulating {0} galaxies".format(hst_idx+1)
            print "Total time for simulating {0} galaxies is {1}".format(n_reset,t2-t1)
            logger.info("Total time for simulating {0} galaxies is {1}".format(n_reset,t2-t1))
            sys.stdout.flush()

    print "Exiting from getPostage stamps."
    sys.stdout.flush()

    #canvas_image = canvas_from_postageStamps(wcs=wcs, canvas_filename='temp', postage_stamps=postage_stamps)

    return canvas_images, stamp_record_list

def workHorse(input_queue, output_queue):
    ## I have no idea how this parallel processing stuff works
    indices, kwargs = input_queue.get()
    results = getPostageStamps(indices, **kwargs)
    output_queue.put(results)

def imsim(exp_id, rot_id, psfset_id, g1, g2, n_rotations, ditherArray, dir_exp, psf_params_pathname, dir_psf=None, star_catalogue=None, galaxy_dat=None, sersic_only=True, cuts=None, parallelize=0):
    ## It is good to log the beginning time and some important input parameters
    logger.info("'imsim' called at: {0}".format(time.ctime()))

    kwds = {'psfset_id':psfset_id, 'exp_id':exp_id, 'rot_id':rot_id, 'n_rotations':n_rotations, 'sersic_only':sersic_only, \
            'ditherArray':ditherArray, 'g1':g1, 'g2':g2, 'fixed_stamp_size':None, 'galaxy_dat':galaxy_dat,'psf_params_pathname':psf_params_pathname}

    rng_seed = rng_generator(base_seed, **kwds)
    rng = galsim.BaseDeviate(rng_seed)

    ### Generate a star field
    if star_catalogue is not None:
        PSF = getKiDSPSFs(psfset_id=psfset_id, exp_id=exp_id,psf_params_pathname=psf_params_pathname)
        star_images = getStarImages(PSF=PSF,star_catalogue=star_catalogue,ditherArray=ditherArray,exp_id=exp_id)
    else:
        star_images = galsim.Image(canvas_bounds,init_value=0.)

    ## Get cuts on the catalogue. Only the galaxies that pass these cuts will be simulated
    ## Convert the cuts to row indices
    if cuts is None:
        indices = np.arange(0,len(galaxy_dat),1).astype(int)
    else:
        indices = np.where(cuts)[0]

    logger.info("'imsim' simulating {0} galaxies.".format(len(indices)))
    print "'imsim' simulating {0} galaxies.".format(len(indices))
    sys.stdout.flush()

    if parallelize:
        import warnings
        warnings.warn("Parallelisation within 'imsim' will mostly not work. Proceed with caution.")

        n_workers = 12 # because lensfit uses 12
        n_jobs = n_workers
        ## Pad hst_indices with nan, if needbe, and wrap it
        n_nanpad = (n_jobs - len(indices)%n_jobs) if len(indices)%n_jobs>0 else 0
        wrapped_indices = np.reshape(np.append(indices, np.nan*np.ones(n_nanpad)), (n_jobs, len(indices)/n_jobs+1), 'F')

        print "Starting with internal parallel processing (method: {0})...".format(parallelize); sys.stdout.flush()
        t1 = time.time()

        ## Queue method
        if parallelize==3:
            task_queue = Queue()
            done_queue = Queue()

            for proc_id in xrange(n_jobs):
                task_queue.put((wrapped_indices[proc_id], kwds))

            procs = [Process(target=workHorse,args=(task_queue, done_queue)) for proc_id in xrange(n_jobs)]
            for proc_id in xrange(n_jobs):
                procs[proc_id].start()
            sys.stdout.flush()

            results = [done_queue.get() for proc_id in xrange(n_jobs)]

        ## apply_async
        elif parallelize==2:
            p = Pool(processes=n_workers)
            results = [ ]
            def log_results(r):
                results.append(r)
            #workHorse = lambda x: getPostageStamps(x,psfset_id=0,exp_id=exp_id,n_rotations=n_rotations,g1=0.,g2=0.,fixed_stamp_size=None)
            def workHorse(indices):
                return getPostageStamps(indices,psfset_id=psfset_id, exp_id=exp_id, n_rotations=n_rotations, g1=g1, g2=g2, fixed_stamp_size=None)
            procs = [p.apply_async(func=getPostageStamps, args=(idx,),kwds=kwds,callback=log_results) for idx in wrapped_indices]
            p.close()
            p.join()

            # http://stackoverflow.com/questions/16224600/multiprocess-apply-async-how-do-i-wrap-args-and-kwargs
            # http://stackoverflow.com/questions/8533318/python-multiprocessing-pool-when-to-use-apply-apply-async-or-map
        
        ## map_async
        elif parallelize==1:
            p = Pool(processes=n_workers)
            def workHorse(idx):
                return_val = getPostageStamps(hst_indices=idx, **kwds)
                return return_val

            iterables = []
            for wi in wrapped_indices:
                item = kwds.copy
            map_result_get = p.map_async(func=workHorse, iterable=wrapped_indices)
            results = map_result.get(timeout=60*60*10) ## timeout = 10 hours

        print len(results)
        canvas_images, stamp_record_lists = zip(*results)

        t2 = time.time()
        print "Parallel processing finished in ", t2-t1, " seconds."
        sys.stdout.flush()

        ## Combine the canvas_images that each worker gives to obtain a full (noiseless) image
        gal_images = combine_canvas([canvas_image for canvas_image in canvas_images])

        ## Make the list of stamp records into a single list
        stamp_records = []
        for stamp_record_list in stamp_record_lists:
            stamp_records += stamp_record_list

    else: # if not parallelize
        gal_images, stamp_records = getPostageStamps(indices, **kwds)

    if save_stamps:
        logger.info("Individual postage stamps saved on the disk.")

    full_images = combine_canvas([gal_images, star_images])

    ## Store the PSF image
    if dir_psf is not None:
        psf_size = 32
        psf_image = galsim.Image(psf_size, psf_size)
        PSF = getKiDSPSFs(psfset_id=psfset_id, exp_id=exp_id,psf_params_pathname=psf_params_pathname)
        PSF_lf = PSF.shift(0.5*pix_scale, 0.5*pix_scale)
        psf_image = PSF_lf.drawImage(image=psf_image, scale=pix_scale)
        psf_filename = 'exp{0}.fits'.format(exp_id)
        psf_pathname = os.path.join(dir_psf,psf_filename)
        psf_image.write(psf_pathname)

    ## Convert the stamp records to a NumPy record array
    stamp_record_array = np.rec.array(stamp_records, \
                            formats='int64,float,float,float,float,float,float,float,int16,float32,float,float,int8,int8',\
                            names='OBJNO,SERSIC_n,HLR,q,PA,FLUX_INPUT,IMAGE_X,IMAGE_Y,STAMP_SIZE,FLUX_FRACTION,DRAWTIME,PROCTIME,FAILURE_FLAG,STITCH_FLAG')
        
    ## Save the stamp records in FITS format
    stamp_hdu = fits.BinTableHDU(data=stamp_record_array)
    ## Include the filter_name and exp_id in the header
    stamp_hdu.header['BANDPASS'] = filter_name
    stamp_hdu.header['EXPO_ID'] =  exp_id
    ## Include the PSF parameters in the header as well
    stamp_hdu.header['PSF_e1'] = -1
    stamp_hdu.header['PSF_e2'] = -1
    stamp_hdu.header['PSF_beta'] = -1 #PSF.getBeta()
    stamp_hdu.header['PSF_FWHM'] = -1 #PSF.getFWHM()
    stamp_hdu.header['rng_seed'] = rng_seed

    ## Write it on to the disk now
    stamp_records_filename = 'stamp_records{1}_rot0{0}_rng{2}.cat'.format(rot_id, exp_id,rng_seed-base_seed)
    stamp_records_pathname = os.path.join(dir_exp,stamp_records_filename)
    stamp_hdu.writeto(stamp_records_pathname, overwrite=True)

    ## Save the noiseless image
    if 0:
        full_images.write(full_image_filename+'.fits')
        logger.info("Noiseless Image written on to the file at {0}".format(time.ctime()))
        swap_pc_cd_header(full_image_filename+'.fits')

    ## Generate the noise field
    full_image_noise = galsim.GaussianNoise(rng, sigma=kids_noise_level)
    noise_field = galsim.Image(full_images.bounds)
    noise_field.addNoise(full_image_noise)
    #full_image.addNoise(full_image_noise)

    ## Add the noise field to the image
    full_images += noise_field

    ## Save the noisy image
    if 1:
        if rot_id == 0:
            full_image_noisy_filename = 'exp{0}.fits'.format(exp_id)
        else:
            full_image_noisy_filename = 'exp{0}_rot{1}.fits'.format(exp_id,str(rot_id).zfill(2))
        full_image_noisy_pathname = os.path.join(dir_exp,full_image_noisy_filename)
        #full_image_noisy_filename = full_image_filename+'_noisy.fits'
        ## Chop the image and save it
        chopped_image = chopImage(full_images, output_filename=full_image_noisy_pathname)
        #full_images[rot_id].write(full_image_noisy_pathname)
        logger.info("Noisy Image written on to the file at {0}".format(time.ctime()))
        swap_pc_cd_header(full_image_noisy_pathname)

    logger.info("'imsim' exited at: {0}".format(time.ctime()))

def chopImage(canvas_image, output_filename=None):
    import galsim
    if isinstance(canvas_image,str):
        canvas_image_filename = canvas_image
        canvas_image = galsim.fits.read(canvas_image)
        canvas_image_data = canvas_image.array
    elif isinstance(canvas_image,galsim.Image):
        canvas_image_data = canvas_image.array

    image_size_x = 16810
    image_size_y = 16530

    old_bounds = canvas_image.bounds
    center = old_bounds.center()
    xmin, xmax = center.x - image_size_x/2, center.x + image_size_x/2
    ymin, ymax = center.y - image_size_y/2, center.y + image_size_y/2

    new_bounds = galsim.BoundsI(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)

    new_image = canvas_image[new_bounds]

    if output_filename is not None:
        print "output_filename provided. Writing it on to the disk."
        new_image.write(output_filename, clobber=True)

    return new_image

def canvas_from_postageStamps(wcs=None, canvas_filename='canvas', postage_stamps='/disks/shear14/arunkannawadi/imsim/postage_stamps/r_band/', add_to_canvas=False, boost=False):
    logging.basicConfig(filename=canvas_filename+".log",level=logging.INFO)
    logger = logging.getLogger("canvas_from_postageStamps")
    logger.info("Starting 'canvas_from_postageStamps at {0} with".format(time.ctime()))
    for variables, values in locals().iteritems():
        logger.info("{0}: {1}".format(variables, values))
    cat_name = "/disks/shear14/KiDS_simulations/Cosmos/KIDS_HST_cat/KiDS_Griffith_iMS1.fits"
    #cat_name = "/disks/shear14/arunkannawadi/catalogue_matching/Griffith_COSMOS_noKiDS_iMS1.fits"
    catalogue = fits.open(cat_name)
    cat_data = catalogue[2].data
    obj_no = cat_data['OBJNO'].astype(int)
    if boost:
        CXX_IMAGE = cat_data['CXX_IMAGE_THELI']
	CYY_IMAGE = cat_data['CYY_IMAGE_THELI']
	CXY_IMAGE = cat_data['CXY_IMAGE_THELI']
	KRON_RADIUS = cat_data['KRON_RADIUS_THELI']
    #MAG_AUTO_THELI = cat_data['MAG_AUTO_THELI']
    #MAG_GAAP_r_CALIB = cat_data['MAG_GAAP_r_CALIB']
    idx3 = np.arange(0,cat_data.size,1).astype(int)

    #step_size = 20
    #selection_cut = (mask)&(LF_weights>=0) #&(np.arange(0,len(LF_weights),1).astype(int)%500==r) # % 500
    #additional_cut = (np.cumsum(selection_cut)>=step_size*0)&(np.cumsum(selection_cut)<step_size*(99+1))
    #selection_cut &= additional_cut

    ## Set up the canvas
    full_image = galsim.Image(canvas_bounds)
    #full_image.setOrigin(0,0) - this is a bug, since real images have 1,1 as the origin!
    if wcs is not None:
        full_image.wcs = wcs

    if isinstance(postage_stamps,str):
        postage_stamp_dir = postage_stamps
        postage_stamp_filenames = os.listdir(postage_stamp_dir)
        n_stamps = len(postage_stamp_filenames)
    else:
        n_stamps = len(postage_stamps)
    #postage_stamp_filenames = [str(l)+'.fits' for l in obj_no[idx3[selection_cut]]]

    t_d1 = time.time()
    n_nans, n_negs, n_bright, n_outofbounds = 0,0,0,0
    for postage_stamp_idx in xrange(n_stamps):
        ## Report progress
        if postage_stamp_idx%1000==0:
            t_d2 = time.time()
            logger.info("1000 galaxies read and incorporated in {0} seconds.".format(t_d2-t_d1))
            print "1000 galaxies read and incorporated in {0} seconds.".format(t_d2-t_d1)
            t_d1 = t_d2
        
        if isinstance(postage_stamps,str): 
            postage_stamp_filename = postage_stamp_filenames[postage_stamp_idx]
            try:    
                stamp = galsim.fits.read(os.path.join(postage_stamp_dir,postage_stamp_filename))
            except:
                logger.info(postage_stamp_filename+" not found.\n")
                print postage_stamp_filename, " not found." 
                continue
        else:
            postage_stamp_filename = postage_stamps[postage_stamp_idx][0]
            stamp = postage_stamps[postage_stamp_idx][1]

	if np.any(np.isnan(stamp.array)):
	    n_nans += 1
	    logger.info(postage_stamp_filename+" is full of nan s.\n")
	    print postage_stamp_filename, "is full of nan s"
	    continue

	if stamp.array.min()<-1e-3:    ## seems to be happening for stars
	    n_negs += 1
	    logger.info(postage_stamp_filename+" has negative values.\n")
	    print postage_stamp_filename, "has negative values"
	    continue

        if stamp.array.sum()>4000: ## way too bright galaxies
            n_bright += 1
            logger.info(postage_stamp_filename+" is too bright, with image_flux > 4000.\n")
            print postage_stamp_filename, " is too bright, with image_flux > 4000."
            continue

        gal_id = int(postage_stamp_filename.split('.')[-2])
        hst_idx = np.where(obj_no==gal_id)[0][0]
        if 0:
            tmp_gain = gain*exp_time*10**(0.4*magAB_zeropoint) ## temporarily only; if postage stamps dont have gain
        else:
            tmp_gain = 1.0 #10**(-0.4*(MAG_AUTO_THELI[hst_idx]-MAG_GAAP_r_CALIB[hst_idx]))
        if boost:
            ## Get the ellipse parameters
            cxx_image = CXX_IMAGE[hst_idx]
            cyy_image = CYY_IMAGE[hst_idx]
            cxy_image = CXY_IMAGE[hst_idx]
            kron_radius = KRON_RADIUS[hst_idx]

            center = stamp.center()
            x0, y0 = center.x, center.y
            x_lin = np.arange(stamp.xmin, stamp.xmax,1).astype(int)
            y_lin = np.arange(stamp.ymin, stamp.ymax,1).astype(int) ## which way is y oriented?
            xx, yy = np.meshgrid(x_lin,y_lin)
            
            stamp_mask = cxx_image*(xx-x0)**2 + cyy_image*(yy-y0)**2 + cxy_image*(xx-x0)*(yy-y0) < (2.5*kron_radius)**2
            flux_boost_factor = stamp.array.sum()/stamp.array[stamp_mask].sum()
	    mask_fraction = 1.*stamp_mask.sum()/stamp_mask.size
            print "flux boost factor = ", flux_boost_factor, ". Mask fraction = ", 1.*stamp_mask.sum()/stamp_mask.size
            logger.info("For gal_id = {0}, flux_boost_factor = {1} with mask fraction = {2}".format(gal_id, flux_boost_factor, mask_fraction))
	    if np.isinf(flux_boost_factor) or np.isnan(flux_boost_factor):
	       tmp_gain = 0.
	    else:
	       tmp_gain *= flux_boost_factor 

        stamp *= tmp_gain ## TEMPORARILY ONLY
        stitch_flag = stitch_stamps(full_image, stamp)
	if stitch_flag:
	    n_outofbounds += 1
	    logger.info("RunetimeError occured with gal_id = {0} because bounds = {1}".format(gal_id, bounds))
	    print "RuntimeError occured with gal_id = {0} because bounds = ".format(gal_id), bounds
    
    #full_image.write(canvas_filename+'.fits')
    #swap_pc_cd_header(canvas_filename+'.fits')

    #full_image_noise = galsim.GaussianNoise(rng,sigma=kids_noise_level)
    #full_image_noisy = full_image.copy()
    #full_image_noisy.addNoise(full_image_noise)
    #full_image_noisy.write(canvas_filename+'_noisy.fits')
    #swap_pc_cd_header(canvas_filename+'_noisy.fits')

    logger.info(" Total no. of postage stamps: {0}".format(len(postage_stamp_filenames)))
    logger.info(" No. of stamps with nan s: {0}".format(n_nans))
    logger.info(" No. of stamps with negative pixel values: {0}".format(n_negs))
    logger.info(" No. of stamps out of bounds: {0}".format(n_outofbounds))
    logger.info("Exiting 'canvas_from_postageStamps at {0}".format(time.ctime()))

    return full_image

def stitch_stamps(canvas_image, stamp, min_pixel_value=-1e-3, max_flux=4000):
    stitch_flag = 0
    if np.any(np.isnan(stamp.array)):
        stitch_flag += 2**0

    if stamp.array.min()<min_pixel_value:    ## seems to be happening for stars
	stitch_flag += 2**1

    if stamp.array.sum()>max_flux: ## way too bright galaxies
        stitch_flag += 2**2
 
    if not stitch_flag:
        bounds = stamp.bounds & canvas_image.bounds
        try:
            canvas_image[bounds] += stamp[bounds]
        except RuntimeError:
            stitch_flag = -1

    return stitch_flag

def generate_noisy_image(canvas_filename,noiseless_image,noise_field=None, f=None):
    if f is None or f==1.0:
        noisy_canvas_filename = canvas_filename+'_noisy.fits'
        f = 1.0
    else:
        noisy_canvas_filename = canvas_filename+'_noisy_f{0}.fits'.format(f)

    if noise_field is None:
        noise = galsim.GaussianNoise(rng,sigma=kids_noise_level)
        noise_field = galsim.Image(noiseless_image.xmax,noiseless_image.ymax)
        noise_field.addNoise(noise)

    noisy_image = noiseless_image + f*noise_field
    noisy_image.wcs = noiseless_image.wcs
    noisy_image.write(noisy_canvas_filename)
    swap_pc_cd_header(noisy_canvas_filename)

def swap_sersic_bulgedisc(canvas_filename, sersic_dir, bulgedisc_dir, output_filename, mode='sersic_to_bulgedisc'):
    if mode is 'sersic_to_bulgedisc': ## the standard mode of operation
        sgn = +1
    elif mode is 'bulgedisc_to_sersic':
        sgn = -1
    else:
    	print "mode = "+ str(mode) +" is unrecognized. Quitting ..."
	raise ValueError()

    ## Get the canvas
    full_image = galsim.fits.read(canvas_filename)

    ## Get the bulge+disc galaxies
    filenames = os.listdir(bulgedisc_dir)

    success_count = 0
    for filename in filenames:
      try:
        sersic_stamp = galsim.fits.read(os.path.join(sersic_dir,filename))
	bulgedisc_stamp = galsim.fits.read(os.path.join(bulgedisc_dir,filename))
      except IOError:
        print filename, "not found."
	continue 

	sersic_bounds = sersic_stamp.bounds & full_image.bounds
	bulgedisc_bounds = bulgedisc_stamp.bounds & full_image.bounds

	full_image[sersic_bounds] -= sgn*sersic_stamp[sersic_bounds]
	full_image[bulgedisc_bounds] += sgn*bulgedisc_stamp[bulgedisc_bounds]

        success_count += 1
	if success_count%100==0:
	    print success_count, "galaxies replaced so far ... "
    full_image.write(output_filename)
    swap_pc_cd_header(output_filename)
    
def swap_pc_cd_header(filename):
    hdulist = fits.open(filename)
    hdr = hdulist[0].header
    try:
        CD1_1, CD2_2 = hdr['CD1_1'], hdr['CD2_2']
    except KeyError:
        CD1_1, CD2_2 = None, None
    try:
        PC1_1, PC2_2 = hdr['PC1_1'], hdr['PC2_2']
    except KeyError:
        PC1_1, PC2_2 = None, None

    CD1_1, PC1_1 = PC1_1, CD1_1 ## swap
    CD2_2, PC2_2 = PC2_2, CD2_2 ## swap

    ## Assign the swapped values
    hdr['CD1_1'], hdr['CD2_2'] = CD1_1, CD2_2
    hdr['PC1_1'], hdr['PC2_2'] = PC1_1, PC2_2

    hdulist.writeto(filename,overwrite=True)

def combine_canvas(canvases, output_filename=None):
    """ Take a list of images or image filenames and return an image by adding all the images.

        @param canvases         list of images or image filenames to be added
        @param output_filename  The name with which the combined image must be saved.
                                The image will not be saved if `output_filename` is None.
                                [default: None]
    """
    ## TODO: Instead of simply adding the images, add by bounds, like done for postage stamps

    final_canvas = galsim.Image(canvas_bounds)
    for canvas in canvases:
        if isinstance(canvas,str):
            canvas_filename = canvas
            canvas = galsim.fits.read(canvas_filename)

        final_canvas += canvas

    wcs = canvas.wcs
    final_canvas.wcs = wcs
   
    if output_filename is not None:
        final_canvas.write(output_filename)
        swap_pc_cd_header(output_filename)

    return final_canvas
        
def subtract_images(image1,image2,output_filename,f=1.0):
    if isinstance(image1,str):
        image1 = galsim.fits.read(image1)
    if isinstance(image2,str):
        image2 = galsim.fits.read(image2)

    if image1.origin() == image2.origin():
        image3 = image1 - f*image2
    else:
        raise ValueError("Origins do not match")

    #image3.setOrigin(0,0)

    #if image1.wcs == image2.wcs:
    image3.wcs = image2.wcs

    image3.write(output_filename)
    swap_pc_cd_header(output_filename)
    

def combine_input_output_catalogue():
    ## Assume the catalogue has been created
    pass

def dumpDithers(path_dither=None, n_exposures=5):
    """ Generates the standard KiDS dither pattern and dumps it into a file for the pipeline to read. No pointing errors are included.

        @param path_dither      The path, if not None, where the 'ditherArray.txt' file needs to be stored. [default: None]
        @param n_exposures      No. of exposures or (1+No.of dither offsets). [default: 5]

        @returns                A n_exposures x 3 array, whose columns correspond to Exposure Id, X offset and Y offset.

    """
    ditherArray = np.dot(np.reshape(np.arange(0,n_exposures,1).astype(int),(n_exposures,1)), np.reshape(np.array([1, -25/0.214, -85/0.214]),(1,3)))
    if path_dither is not None:
        fname = 'ditherArray.txt'
        pathname = os.path.join(path_dither, fname)
        comments = 'EXP_ID \t +Dither x \t +Dither y'
        np.savetxt(fname=path_dither, X=ditherArray, header=comments, fmt=['%d', '%4.3f', '%4.3f'])

    return ditherArray

## Routine with Ricardo's signature
#def create_imsims(g1,g2,nproc,path_to_prior_file,\
#            path_to_psf_archive_folder, path_to_image_archive_folder,\
#            configFile, ditherArray, tmpDir, noise_sigma_in):
    
def create_imsims(psfSet,g1,g2,path_input_cat,path_star_cat,path_dither,dir_psf,dir_exp, rot_id, psf_params_pathname,
                    n_gal=True,stars=False,faint_gal=False, n_rotations=4, n_exposures=5,sersic_only=True,parallelize=True,internal_parallelize=2):

    ## Create a dither pattern over the multiple exposures and save it
    ditherArray =  dumpDithers(path_dither,n_exposures=n_exposures)

    ## Load the galaxy catalogue
    galaxy_cat = fits.open(path_input_cat)

    apply_cuts = False
    if apply_cuts is True:
        MASK_all = galaxy_cat[1].data.MASK
        mask = ~(np.array(MASK_all&0xfc3c,dtype=bool))
        rank = galaxy_cat[1].data['rank']
        distance2d = galaxy_cat[1].data['distance2d']
        weight = galaxy_cat[1].data['weight']
        mag1 = galaxy_cat[1].data['MAG_AUTO_THELI']
        size1 = galaxy_cat[1].data['FWHM_IMAGE_THELI'] # pix
        mag0 = galaxy_cat[1].data['MAGR']
        #size0 = galaxy_cat[1].data['FLUX_RADIUS_HI']
        easy_cuts = (mask)&(rank==1)&(distance2d<0.03)&(weight>12)
        default_cuts = (mask)&(rank>=0)&(distance2d<1)&(weight>=0)
        #star_cuts = (~(mag1<24)&(size1<2.8)&(rank>0) |~((rank==0)&(mag0<24)&(size0<2.8)&(mag0>0)&(size0>0)))
        handpicked_stars = galaxy_cat[1].data['handpicked_stars']
        star_cuts = (mag0>0)&(~handpicked_stars)
        #cuts = default_cuts&star_cuts
        #cuts = easy_cuts&star_cuts
        #cuts = cuts[np.cumsum(cuts)<2000]
        cuts = np.arange(0,len(mag0),1)<20

        galaxy_dat = galaxy_cat[1].data[cuts]
    else:
        galaxy_dat = galaxy_cat[1].data

    print "Input galaxy catalogue loaded ..."
    sys.stdout.flush()

    kwargs = {'g1':g1, 'g2':g2, 'n_rotations':n_rotations, 'dir_exp':dir_exp, 'dir_psf':dir_psf,
              'ditherArray':ditherArray,'parallelize':internal_parallelize, 'sersic_only':sersic_only,
              'star_catalogue':path_star_cat, 'galaxy_dat':galaxy_dat, 'psfset_id':psfSet,'psf_params_pathname':psf_params_pathname}

    if parallelize:
        n_workers = n_exposures
        p = Pool(processes=n_workers)
        procs = []
        for exp_id in xrange(n_exposures):
            proc = p.apply_async(func=imsim,args=(exp_id,rot_id),kwds=kwargs)
            procs.append( proc )
        p.close()
        p.join()

    else: ## if not parallelize
        for exp_id in xrange(n_exposures):
            imsim(exp_id=exp_id,rot_id=rot_id,**kwargs)


def create_imsim(psfSet, g1, g2,\
                path_input_cat, path_dither,\
                dir_psf, dir_exp, n_gal=True, \
                stars=False, faint_gal=False, rot=False):
    ## A Hack for now
    import subprocess
    import os

    for exp_id in range(5):
        subprocess.call(["cp", "/disks/shear14/arunkannawadi/imsim/exp{0}.fits".format(exp_id), dir_exp])
        subprocess.call(["cp", "/disks/shear14/arunkannawadi/imsim/psf_exp{0}.fits".format(exp_id),\
                        os.path.join(dir_psf,"exp{0}.fits".format(exp_id))])
        #subprocess.call(["cp","/disks/shear14/arunkannawadi/catalogue_matching/KiDS_unmasked_Griffith_iMS1_reducedforpipeline.cat", path_input_cat])
        subprocess.call(["cp", "/disks/shear14/arunkannawadi/imsim/ditherArray.txt", path_dither])

def rng_generator(base_seed, **kwargs):
    g1 = kwargs['g1']
    g2 = kwargs['g2']
    psfset_id = kwargs['psfset_id']
    exp_id = kwargs['exp_id']
    rot_id = kwargs['rot_id']

    ## All 5x416 square degrees get different rng_seed, so that the noise realisations don't repeat irrespective of parallelization
    rng_seed = base_seed + int(np.arctan(g2/g1)*4*2/np.pi) + 1000*psfset_id + 10*exp_id + 100*rot_id
    return rng_seed
    

#if __name__=='__main__':
#    ## Get cuts on the catalogue. Only the galaxies that pass these cuts will be simulated
#    ## Note: rank==-1 objects do not have the structural parameters and cannot be simulated ever.
#    cuts = (mask)&(rank>0)&(distance2d<1.2)&(weight>0)
#    print cuts.sum()
#
#    imsim(cuts=cuts, parallelize=True)
#    #imsim(0,500,start_afresh=True,parallelize=False,filter_name='r')
#    #imsim(1150,1500,start_afresh=False,parallelize=False,filter_name='r')
#    #imsim(1000,1500,start_afresh=False,parallelize=False,filter_name='r')
#    #imsim(1500,2000,start_afresh=False,parallelize=False,filter_name='r')
#    canvas_filename = 'canvas_all_deep_lt4000_noisy.fits'
#    output_filename = 'canvas_all_noisy_bd_tmp.fits'
#    bulgedisc_dir = 'postage_stamps/r_band/free_size/bd_images/'
#    sersic_dir = 'postage_stamps/r_band/free_size/'
#    #swap_sersic_bulgedisc(canvas_filename,sersic_dir,bulgedisc_dir,output_filename)
#    canvas_filename = output_filename
#    output_filename = 'canvas_all_noise_bd.fits'
#    sersic_dir = 'postage_stamps/r_band/free_size/noKiDS/'
#    #swap_sersic_bulgedisc(canvas_filename,sersic_dir,bulgedisc_dir,output_filename)
#
#    #canvas_from_postageStamps(wcs=wcs,canvas_filename='canvas_trial_boosted',postage_stamp_dir='/disks/shear14/arunkannawadi/imsim/postage_stamps/r_band//trunc5', boost=True)
#    #combine_canvas(canvas_filenames=['/disks/shear14/arunkannawadi/imsim/canvas_all.fits','/disks/shear14/arunkannawadi/imsim/canvas_noKiDS_lt4000.fits'], output_filename='/disks/shear14/arunkannawadi/imsim/canvas_all_deep_lt4000.fits')
#    f = 1.0
#    #subtract_images('/disks/shear14/KiDS_simulations/Cosmos/Theli_image/KIDS_150p1_2p2_r_SDSS.V0.5.9A.swarp.cut.fits','/disks/shear14/arunkannawadi/imsim/canvas_trial_boosted.fits',output_filename='/disks/shear14/arunkannawadi/imsim/canvas_trial_boosted_diff.fits'.format(f),f=f)
#    #subtract_images('/disks/shear14/KiDS_simulations/Cosmos/Theli_image/KIDS_150p1_2p2_r_SDSS.V0.5.9A.swarp.cut.fits','/disks/shear14/arunkannawadi/imsim/canvas_trial.fits',output_filename='/disks/shear14/arunkannawadi/imsim/canvas_trial_diff.fits'.format(f),f=f)
