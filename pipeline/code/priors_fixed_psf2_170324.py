"""
Create input catalog for KiDS simulations maker to read in as input for galaxy parameters
"""
import sys
#import os
import numpy
import numpy as np
#from matplotlib import pyplot as plt
#import pyfits as pf
#import time
#import random
#import scipy
#import scipy.optimize
from scipy import interpolate

#========================================== Functions ==============================================

# Bad old size prior
#def bad_old_size_prior(r,v606):
#	a = numpy.exp(-0.376174*(v606-23) -0.633215)
#	b = numpy.exp(-0.087752*(v606-23) +0.723772)
#	return numpy.exp(-(numpy.log10(r)- numpy.log10(a))**2/(2*numpy.log10(b)))

def size_prior(r,m):
	a = numpy.exp(-1.319698-0.277907*(m-23))/1.134
	return r*numpy.exp(-(r/a)**(4./3.))

def disk_ell_prior(e):
	return e*(1-numpy.exp((e-0.804)/0.2539))/((1+e)*numpy.sqrt(e**2+0.0256**2))

def bulge_ell_prior(e):
	return e*numpy.exp(-2.368*e-6.691*e**2)

def bulge_to_disk_ratio_prior(x):
	return numpy.exp(-x**2/(2*0.1**2))

def bulge_to_disk_ratio_prior2(nsamples):
	f=numpy.random.normal(scale=0.1,size=nsamples)
	for g in range(len(f)):
		while ((f[g]>=1) or (f[g]<0)):
			f[g]=numpy.random.normal(scale=0.1)
	a = numpy.random.uniform(size=nsamples)*10
	f[a>9]=1
	return f

# Bad old power law
#def bad_old_maghist_power_law(x,norm):
#	return norm*(1.30114292e-22)*x**(1.97479523e+01)

#def maghist_power_law(x,norm):
#	return norm* 1.05550893e-04 *10**(x*3.90781408e-01)

def mag_hist_power_law(x,norm):
	#[-12.46485528   0.98742738  -0.01346775]
	# run mag_power_law.py for these numbers in log10(N[per sq deg])
	# return norm* 10**(-12.46485528 + 0.98742738*x -0.01346775*x*x)
	# return norm * 10**(-8.04876893e+00 + x*6.64855006e-01 -7.51650932e-03*x*x)
	return norm* 10**(-8.85177133 +  0.716305260*x  -0.00832345561*x*x)
	

def draw_samples_specified_distribution(nsamples, fnc, z,*args):
	#zmin=1e-5
	#zmax=20
	#z = numpy.linspace(zmin,zmax) # Example indep variable
	fofz = fnc(z,*args)
	cdf_fofz = fofz.cumsum() # Next few steps to calculate CDF
	cdf_fofz -= cdf_fofz[0] # w/appropriate normalization
	cdf_fofz /= cdf_fofz[-1]
	max = numpy.where(cdf_fofz==1)[0].min()+1. #Nothing in CDF above 1
	cdf_fofz = cdf_fofz[:max] # cut off the trailing ones
	z = z[:max]
	model = interpolate.InterpolatedUnivariateSpline(cdf_fofz,z,k=4)
	samples = model(numpy.random.random(nsamples)) # Draw samples
	return samples


def get_random_gal_params(magn,nobj,r,ed,eb):
	scale_radius = np.array(
				[draw_samples_specified_distribution(1,size_prior,r,mm)[0] for mm in magn])
	# Something went wrong if the prior range is 0.01-20, so now manually truncate
	scale_radius[scale_radius<0.01]=0.01
	disk_ell = draw_samples_specified_distribution(nobj,disk_ell_prior,ed)
	bulge_ell = draw_samples_specified_distribution(nobj,bulge_ell_prior,eb)
	bd_ratio = bulge_to_disk_ratio_prior2(nobj)
	beta = numpy.random.uniform(0,180,size=nobj)

	# Change from    mag, scalelength, bulge ell, disk ell, pos angle, b/t ratio
	# 		   to    x, y, mag, scalelength, e1, e2, b/t ratio
	# Already account for when an object is bulge and when it's bulge+disk
	e1=numpy.cos(2*beta*(numpy.pi/180))
	e2=numpy.sin(2*beta*(numpy.pi/180))
	e1[bd_ratio==1]*=bulge_ell[bd_ratio==1]
	e2[bd_ratio==1]*=bulge_ell[bd_ratio==1]
	e1[bd_ratio<1]*=disk_ell[bd_ratio<1]
	e2[bd_ratio<1]*=disk_ell[bd_ratio<1]

	return scale_radius, e1, e2, bd_ratio

#================================================= Main ============================================


def create_priorfile(current_psf_set,path_to_bsc_model,
						prior_filename,psf_setpath,grid_positions=True,
						nr_of_stars=True,nr_of_bright_gals=True,nr_of_faint_gals=True,
						random_seed=None):
	"""
	Create a prior file to be read in by imsim.py
	Requires:
	A list of magnitudes measured in 109 KiDS DR2 fields from which magnitudes will be
	randomly sampled. A list with Moffat parameters for KiDS fields measured using PSFEx by Reiko
	Nakajima, to sample the PSF from. A list of star properties as given by the Besancon Galaxy
	model, used to sample the star magnitudes from. Name for the output prior file.
	Optional arguments:
	Determine the presence of stars, bright (m<25) galaxies and faint (25<m<28) galaxies
	All either a bool or an int or a float
	B Add objects to the images Y/N.
	I Add exactly this number of them to the image
	F Add any object below with mag<=this number

	WARNING: This script should not be run in parallel, as it uses numpy random number generator. It
	is not thread safe; i.e. different simultaneous runs will change the seed and the process will
	not be reproduceable.
	EDIT: Has been checked to function fine in Malta cluster
	"""

	#---------------------------------------- Set up -----------------------------------------------

	# Catalog format
	ofmt_s = '%.2f   %.2f   %.3f   %d   %d   %d   %d\n'			# Stars
	ofmt_g = '%.2f   %.2f   %.3f   %.3f   %.3f   %.3f   %.3f\n'	# Bright gals and 1000 faint gals
	ofmt_d = '%.2f   %.2f   %d   %d   %d   %d   %d\n'			# Faint duplicates


	# Image specifications
	# Single field size + chip gaps? or remove chip gaps again?
	# Plus the added size from four dithers
	# Plus 12 pixels for the varying error in the pointing which should not exceed 3 for each dither
	# In the simulation script there is already an if statement to differentiate between the
	# different dithers, so in the input file just add several redundant objects; I make a rectangle
	# which covers all the exposures, but also unobserved space.
	image_size_x=16810
	image_size_y=16530
	dither_x=numpy.ceil(25/0.214)
	dither_y=numpy.ceil(85/0.214)
	all_image_size_x=(image_size_x+4*dither_x+4*3)		# Rectangle containing all 5 exposures
	all_image_size_y=(image_size_y+4*dither_y+4*3)		# with extra padding of 3 pixels
	image_area_pix=all_image_size_x*all_image_size_y
	image_area_deg=image_area_pix*((0.214/3600)**2)

	# Range for galaxy parameters
	r=numpy.arange(0.,20,0.01)
	ed=numpy.arange(0,0.805,0.01)
	eb=numpy.arange(0,1.00,0.01)
	total_mag_range=numpy.arange(20,28,0.1)

	# Power law for n(m)
	mag_power_law = mag_hist_power_law(total_mag_range, image_area_deg)


	# Size of the sample of galaxies which will be duplicated from for the faint galaxies to
	# increase the speed of the simulation, but still have them in.
	nsample=1000


	#------------------------------------- Header comments -----------------------------------------

	# Input catalog name
	outfile = open(prior_filename,'w')

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
	outfile.write('# Random seed: '+str(random_seed)+'\n')
	numpy.random.seed(random_seed)

	# Dithers and pointing errors
	dither_pointing_error = numpy.round(numpy.random.normal(0,0.5,8),decimals=2)
	outfile.write('# Dither 0 pointing error:  0  0\n')
	for i in range(4):
		outfile.write('# Dither %d + pointing error:  %f  %f \n' %(i+1,
						(i+1)*dither_x+dither_pointing_error[2*i],
						(i+1)*dither_y+dither_pointing_error[2*i+1]))

	# PSF parameters
	# PSF Moffat parameters measured by Reiko Nakajima using PSFEx on KiDS fields
	# Specifically chosen to be a nicely distributed collection of PSF FWHM
	psf_filename=psf_setpath
	psf_field,psf_exp = numpy.loadtxt(psf_filename, usecols=(0,1),dtype=str,unpack=True)
	psf_chip = numpy.loadtxt(psf_filename, usecols=(2,),dtype=int,unpack=True)
	psf_fwhm,psf_ell1,psf_ell2,psf_beta = numpy.loadtxt(psf_filename, usecols=(3,4,5,6),dtype=float,unpack=True)

	# There are some catalog entries that have zeros everywhere (55 out of total).
	# These will crash the simmaker script, so need to be removed. If any of the 5 exposures is
	# going to have such a flawed line, try another random catalog line.
	for i in range(5):
		psf_ind = int(current_psf_set)*5+i
		outfile.write('# PSF dither %d from field %s, exposure %s, chip %d, num %d\n'
						% (i,psf_field[psf_ind],psf_exp[psf_ind],psf_chip[psf_ind],0))
		outfile.write('# PSF parameters, fwhm(arcsec) %f beta %f e1 %f e2 %f\n'
						% (psf_fwhm[psf_ind],psf_beta[psf_ind],psf_ell1[psf_ind],psf_ell2[psf_ind]))
	outfile.write('# (use g for PSF ellipticity in GalSim!)\n')

	# Parameter explanation
	outfile.write('# x  y   m   r   e1   e2   f\n')




	#-------------------------------- Stars --------------------------------------------------------


	if grid_positions:
		stamp_size=60
		pos = numpy.asarray([(xx*stamp_size+stamp_size/2,yy*stamp_size+stamp_size/2)
										for xx in range(int(all_image_size_x)/stamp_size)
										for yy in range(int(all_image_size_y)/stamp_size)],
										dtype=float)
		# Random position within the grid, to prevent correlation of position and magnitude
		numpy.random.shuffle(pos)
		x,y = zip(*pos)
		current_obj_nr=0



	# Take stars from the Besancon model in CFHT r filter
	# We don't need super bright stars as they would be masked in the real data, so no stars in the
	# catalog below mag=16
	star_mag = numpy.loadtxt(path_to_bsc_model,unpack=True,dtype=float)

	# Range for object positions
	if grid_positions:
		pass
	else:
		x = numpy.round( numpy.random.uniform(0,all_image_size_x,len(star_mag)),decimals=2)
		y = numpy.round( numpy.random.uniform(0,all_image_size_y,len(star_mag)),decimals=2)


	# All random processes are always done!	Only writing is omitted if the user wants.
	# This ensures that the seed progresses as though the user did want stars. This way the seed
	# will be kept the same for simulations with and without stars, if the user gives in the
	# same random seed.
	if nr_of_stars:

		if type(nr_of_stars)==int or type(nr_of_stars)==np.int64:
			for i in range(nr_of_stars):
				outfile.write(ofmt_s %(x[i],y[i],star_mag[i],0,0,0,-1))
			if grid_positions:
				current_obj_nr = i+1
		elif type(nr_of_stars)==float or type(nr_of_stars)==np.float64:
			for i in range(len(star_mag)):
				if star_mag[i]<=nr_of_stars:
					outfile.write(ofmt_s %(x[i],y[i],star_mag[i],0,0,0,-1))
					current_obj_nr+=1
		else:
			for i in range(len(star_mag)):
				outfile.write(ofmt_s %(x[i],y[i],star_mag[i],0,0,0,-1))
			if grid_positions:
				current_obj_nr = i+1



	#----------------------------------- 20<m<25 Galaxies ------------------------------------------


	#t0=time.time()

	# !!! All random processes are always done, this saves the change in RNG seed even if objects
	# !!! are not included in the catalog. So only writing is optional. This way the seed
	# !!! will be kept the same for simulations with and without stars, if the user gives in the
	# !!! same random seed.

	# All galaxies follow a power law behaviour and are randomly selected for each magnitude bin
	m=[]
	for	i in range(len(total_mag_range)):
		if total_mag_range[i]>24.91: break
		# Uniform randomly pick magnitudes from the bin
		m = numpy.append(m,
				numpy.random.uniform(total_mag_range[i],total_mag_range[i]+0.099,mag_power_law[i]))
	n=len(m)
	scale_radius, e1, e2, bd_ratio = get_random_gal_params(m,n,r,ed,eb)


	# Range for object positions
	if grid_positions:
		pass
	else:
		x = numpy.round( numpy.random.uniform(0,all_image_size_x,n),decimals=2)
		y = numpy.round( numpy.random.uniform(0,all_image_size_y,n),decimals=2)


	if nr_of_bright_gals:

		if type(nr_of_bright_gals)==int:
			if grid_positions:
				for i in range(nr_of_bright_gals):
					outfile.write(ofmt_g %(x[i+current_obj_nr],y[i+current_obj_nr],
									m[i],scale_radius[i],e1[i],e2[i],bd_ratio[i]))
				current_obj_nr+=i+1
			else:
				for i in range(nr_of_bright_gals):
					outfile.write(ofmt_g %(x[i],y[i],m[i],scale_radius[i],e1[i],e2[i],bd_ratio[i]))
		elif type(nr_of_bright_gals)==float:
			if grid_positions:
				current_obj_nr_counter=0
				for i in range(n):
					if m[i]<=nr_of_bright_gals:
						outfile.write(ofmt_g %(x[i+current_obj_nr],y[i+current_obj_nr],
									m[i],scale_radius[i],e1[i],e2[i],bd_ratio[i]))
						current_obj_nr_counter+=1
				current_obj_nr+=current_obj_nr_counter
			else:
				for i in range(n):
					if m[i]<=nr_of_bright_gals:
						outfile.write(ofmt_g %(x[i],y[i],m[i],scale_radius[i],e1[i],e2[i],
												bd_ratio[i]))
		else:
			if grid_positions:
				for i in range(n):
					outfile.write(ofmt_g %(x[i+current_obj_nr],y[i+current_obj_nr],
									m[i],scale_radius[i],e1[i],e2[i],bd_ratio[i]))
				current_obj_nr+=i+1
			else:
				for i in range(n):
					outfile.write(ofmt_g %(x[i],y[i],m[i],scale_radius[i],e1[i],e2[i],bd_ratio[i]))



	#------------------------------------ m>25 Galaxies --------------------------------------------

	# !!! Faint galaxies are created last. They are the most RNG intensive, so changing that would
	# !!! be a lot of work. Since after this it doesn't matter, I'll leave it as is...

	if nr_of_faint_gals:
		#t0=time.time()

		# Counts the number of faint galaxies, used to halt the for loop where user wants it
		counter=0

		# Galaxies with m>25 have random magnitudes sampled from the power law
		# For each mag bin of 0.1, pick 1000 random magnitudes and create all galaxy parameters for
		# those 1000. Then fill up to the required number by the power law by choosing only x and y
		# coordinates and a positive number below 1000. This number represents one of the first
		# thousand entries for the mag bin. This will duplicate objects in the simulation to
		# increase speed. Write out immediately for each mag bin. Should save them, though more I/O
		for	i in range(len(total_mag_range)):
			if total_mag_range[i]<24.99: continue

			if type(nr_of_faint_gals)==float:
				if total_mag_range[i]>nr_of_faint_gals: break

			npl=mag_power_law[i]

			# Uniform randomly pick magnitudes from the selection
			m = numpy.random.uniform(total_mag_range[i],total_mag_range[i]+0.099,nsample)
			scale_radius, e1, e2, bd_ratio = get_random_gal_params(m,nsample,r,ed,eb)


			# Range for object positions
			if grid_positions:
				pass
			else:
				x = numpy.round( numpy.random.uniform(0,all_image_size_x,nsample),decimals=2)
				y = numpy.round( numpy.random.uniform(0,all_image_size_y,nsample),decimals=2)



			if grid_positions:
				for j in range(nsample):
					if type(nr_of_faint_gals)==int:
						if counter==nr_of_faint_gals: break
					counter+=1
					outfile.write(ofmt_g %(x[j+current_obj_nr],y[j+current_obj_nr],
									m[j],scale_radius[j],e1[j],e2[j],bd_ratio[j]))
				current_obj_nr+=j+1
			else:
				for j in range(nsample):
					if type(nr_of_faint_gals)==int:
						if counter==nr_of_faint_gals: break
					counter+=1
					outfile.write(ofmt_g %(x[j],y[j],m[j],scale_radius[j],e1[j],e2[j],bd_ratio[j]))


			# Fill up the required number by the power law with coordinates and duplication number
			if grid_positions:
				pass
			else:
				x = numpy.round( numpy.random.uniform(0,all_image_size_x,npl-nsample),
								decimals=2)
				y = numpy.round( numpy.random.uniform(0,all_image_size_y,npl-nsample),
								decimals=2)
			duplication_number = numpy.random.randint(0,1000,npl-nsample)

			if grid_positions:
				for j in range(int(npl-nsample)):
					if type(nr_of_faint_gals)==int:
						if counter==nr_of_faint_gals: break
					counter+=1
					outfile.write(ofmt_d %(x[j+current_obj_nr],y[j+current_obj_nr],
											duplication_number[j],0,0,0,0))
				current_obj_nr+=j+1
			else:
				for j in range(int(npl-nsample)):
					if type(nr_of_faint_gals)==int:
						if counter==nr_of_faint_gals: break
					counter+=1
					outfile.write(ofmt_d %(x[j],y[j],duplication_number[j],0,0,0,0))

			# Extra break to break out of the for loop over magnitude bins
			if type(nr_of_faint_gals)==int:
				if counter==nr_of_faint_gals: break
			counter+=1

			#print time.time()-t0, total_mag_range[i], npl




	#------------------------------------- Done ----------------------------------------------------

	outfile.close()



if __name__ == "__main__":
	bla1=False
	bla2=True
	bla3=True
	create_priorfile(sys.argv[1],sys.argv[2],sys.argv[3],False,bla1,bla2,bla3,'priors/prior_m000m400_3')
