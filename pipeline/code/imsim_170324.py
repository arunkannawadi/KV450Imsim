"""
Image simulations for KiDS simplistic
"""
import sys
import os
import math
import numpy
#import logging
import time
#import shutil
import subprocess
from multiprocessing import Process, Queue, current_process, cpu_count
import pipeline_170324 as pipeline
import optparse
if os.path.isdir('/net/meije/data1/simulations/GalSim/'):
	sys.path.insert(0,'../../GalSim/')
import galsim
import pyfits
#---------------------------------- Functions ------------------------------------------------------



def create_gal_stamp(x,y,gal_mag,gal_size,gal_e1,gal_e2,gal_bt,
                        g1,g2,psf,pixel_scale,mag_zero):

    # Convert Lensfit scale radius to GalSim azimuthally averaged half light radius
    gal_e = numpy.sqrt(gal_e1**2+gal_e2**2)
    gal_size = gal_size*numpy.sqrt((1-gal_e)/(1+gal_e))         # arcsec

    # Create galaxy profile
    if gal_bt==1:
        #gal = galsim.DeVaucouleurs(flux=1., half_light_radius=gal_size)
        gal = galsim.Sersic(flux=1., n=4, half_light_radius=gal_size,
									trunc=4.5*gal_size,flux_untruncated=False)
    else:
        #bulge = galsim.DeVaucouleurs(flux=gal_bt, half_light_radius=gal_size)
        #disk = galsim.Exponential(flux=1.-gal_bt, scale_radius=gal_size)
        disk = galsim.Sersic(flux=1.-gal_bt, n=1, scale_radius=gal_size,
									trunc=4.5*gal_size,flux_untruncated=False)
        bulge = galsim.Sersic(flux=gal_bt, n=4, half_light_radius=gal_size,
									trunc=4.5*gal_size,flux_untruncated=False)
        gal = galsim.Add([bulge, disk])
    # Galaxy magnitude
    gal = gal.withFlux(10**(-0.4*(gal_mag-mag_zero)))
    # Galaxy ellipticity
    # No longer copy every object, but use new GS feature which works (see demo5.py)
    gal2 = gal.shear(g1=-gal_e2, g2=gal_e1)
    gal3 = gal.shear(g1=-gal_e1, g2=-gal_e2)
    gal4 = gal.shear(g1=gal_e2, g2=-gal_e1)
    gal = gal.shear(g1=gal_e1, g2=gal_e2)

    # Apply cosmological shear
    gal = gal.shear(g1=g1, g2=g2)
    gal2 = gal2.shear(g1=g1, g2=g2)
    gal3 = gal3.shear(g1=g1, g2=g2)
    gal4 = gal4.shear(g1=g1, g2=g2)

    # Convolve with the PSF and pixelize.
    final = galsim.Convolve(psf, gal)
    final2 = galsim.Convolve(psf, gal2)
    final3 = galsim.Convolve(psf, gal3)
    final4 = galsim.Convolve(psf, gal4)

    # Account for the fractional part of the position
    ix = int(math.floor(x+0.5))
    iy = int(math.floor(y+0.5))
    final = final.shift((x-ix)*pixel_scale,(y-iy)*pixel_scale)
    final2 = final2.shift((x-ix)*pixel_scale,(y-iy)*pixel_scale)
    final3 = final3.shift((x-ix)*pixel_scale,(y-iy)*pixel_scale)
    final4 = final4.shift((x-ix)*pixel_scale,(y-iy)*pixel_scale)

    # Draw it (with our desired stamp size)
    try:
        stamp = final.drawImage(scale=pixel_scale)
        stamp2 = final2.drawImage(scale=pixel_scale)
        stamp3 = final3.drawImage(scale=pixel_scale)
        stamp4 = final4.drawImage(scale=pixel_scale)
    except:
        stamp = galsim.ImageF(100,100,scale=pixel_scale)
        stamp2 = galsim.ImageF(100,100,scale=pixel_scale)
        stamp3 = galsim.ImageF(100,100,scale=pixel_scale)
        stamp4 = galsim.ImageF(100,100,scale=pixel_scale)
        final.drawImage(image=stamp)
        final2.drawImage(image=stamp2)
        final3.drawImage(image=stamp3)
        final4.drawImage(image=stamp4)
    return stamp, stamp2, stamp3, stamp4

#------------------------------------- Main --------------------------------------------------------



def create_imsims(g1,g2,nproc,
            path_to_prior_file,
            path_to_psf_archive_folder,
            path_to_image_archive_folder, configFile, ditherArray, tmpDir, noise_sigma_in):

    #---------------------------- Set up parameters for simulations --------------------------------

    # Image
    pixel_scale=0.214                   # arcsec/pixel
    image_size_x=16810                  # size (pixels) along y-axis of full image
    image_size_x_arcsec=image_size_x * pixel_scale
    image_size_x_arcmin=image_size_x_arcsec / 60.
    image_size_y=16530                  # size (pixels) along y-axis of full image
    image_size_y_arcsec=image_size_y * pixel_scale
    image_size_y_arcmin=image_size_y_arcsec / 60.
    #pix = galsim.Pixel(pixel_scale)


    # Noise
    mag_zero = 24.77+2.5*numpy.log10(360)   # From KiDS chip header
    # mag_zero = 24.79+2.5*numpy.log10(360)
    # mag_zero = 24.77+2.5*numpy.log10(360)
	# mag_zero = 24.67+2.5*numpy.log10(360)
    noise_sigma=noise_sigma_in          # From KiDS chip header
                                        # The background stdev in a single KiDS chip is around 20
                                        # But the magnitude zeropoint is given for a 1 second
                                        # exposure. So the choice is to change the noise level or
                                        # to change the magnitude zeropoint.

    print '      Magnitude Zeropoint = %f , Noise Sigma - %f' % (mag_zero, noise_sigma)

    #------------------------------ Image builder function -----------------------------------------

    def build_file(g1,g2, exposure_number, x_dither,y_dither, psf_params, prior_file,
                    gal_file_name, gal_file_name2,gal_file_name3, gal_file_name4, random_seed):

        """
        Does all the work building the two images
        """
        t0=time.time()

        # Duplication of faint galaxies
        duplicate_stamps=[]
        duplicate_stamps2=[]
        duplicate_stamps3=[]
        duplicate_stamps4=[]

        # PSF
        psf = galsim.Moffat(fwhm=psf_params[0], beta=psf_params[1], flux=1.,trunc=4.5*psf_params[0])
        psf = psf.shear(g1=psf_params[2],g2=psf_params[3])
        #final_psf = galsim.Convolve([psf, pix],real_space=True)
        # Write out the PSF for Lensfit, the centre has to be in pixel 16 (0 indexed)
        psf_lf = psf.shift(0.5*pixel_scale,0.5*pixel_scale)
        psf_image = galsim.ImageF(32,32,scale=pixel_scale)
        psf_lf.drawImage(image=psf_image)	# COMMENT IF YOU WANT TO REVERT BACK TO THE OLD WRONG WAY
        #psf.drawImage(image=psf_image)		# UNCOMMENT IF YOU WANT TO REVERT BACK TO THE OLD WRONG WAY
        psf_image.write(os.path.join(path_to_psf_archive_folder,
                                        'exp'+str(exposure_number)+'.fits'))
        # But use the unshifted PSF for the convolution. Because the galaxies are at their "unshifted"
        # positions, so convolving with the shifted PSFs would be wrong.


        # Setup the full image:
        full_image = galsim.ImageF(image_size_x, image_size_y,scale=pixel_scale)
        cenx = image_size_x/2 + 1
        ceny = image_size_y/2 + 1
        center = galsim.PositionD(cenx,ceny) * pixel_scale
        full_image2 = galsim.ImageF(image_size_x, image_size_y, scale=pixel_scale)
        full_image3 = galsim.ImageF(image_size_x, image_size_y, scale=pixel_scale)
        full_image4 = galsim.ImageF(image_size_x, image_size_y, scale=pixel_scale)

        # Loop over all objects that go into the image
        k=0; kk=0
        with open(path_to_prior_file) as fileobject:
            for line in fileobject:
                if line[0]=='#': continue
                k+=1
                (x,y,gal_mag,gal_size,
                gal_e1,gal_e2,gal_bt) = [float(lineobject) for lineobject in line.split()]

                #print k, x, y

                # Apply the dither + pointing error from prior file
                x-=x_dither
                y-=y_dither

                # Check whether this object is within the field of view after dithering
                if (x<0) or (x>image_size_x) or (y<0) or (y>image_size_y):
                    # Make sure the duplicate sample has exactly 1000 stamps
                    if gal_mag>=25 and gal_size>0:
                        stamp,stamp2,stamp3,stamp4 = create_gal_stamp(x,y,
                                                            gal_mag,gal_size,gal_e1,gal_e2,gal_bt,
                                                            g1,g2,psf,pixel_scale,mag_zero)
                        duplicate_stamps.append(stamp)
                        duplicate_stamps2.append(stamp2)
                        duplicate_stamps3.append(stamp3)
                        duplicate_stamps4.append(stamp4)
                    continue

                # Account for the fractional part of the position
                ix = int(math.floor(x+0.5))
                iy = int(math.floor(y+0.5))


                if gal_bt==-1:

                    # Create star
                    final_psf_model=psf.copy()
                    final_psf_model = final_psf_model.withFlux(10**(0.4*(mag_zero-gal_mag)))
                    final=final_psf_model
                    final2=final_psf_model.copy()
                    final3=final_psf_model.copy()
                    final4=final_psf_model.copy()
                    tijd1=time.time()

                    # Account for the fractional part of the position:
                    final = final.shift((x-ix)*pixel_scale,(y-iy)*pixel_scale)
                    final2 = final2.shift((x-ix)*pixel_scale,(y-iy)*pixel_scale)
                    final3 = final3.shift((x-ix)*pixel_scale,(y-iy)*pixel_scale)
                    final4 = final4.shift((x-ix)*pixel_scale,(y-iy)*pixel_scale)

                    # Draw objects (let GalSim try on it's own, otherwise enforce)
                    try:
                        stamp = final.drawImage(scale=pixel_scale)
                        stamp2 = final2.drawImage(scale=pixel_scale)
                        stamp3 = final3.drawImage(scale=pixel_scale)
                        stamp4 = final4.drawImage(scale=pixel_scale)
                    except:
                        stamp = galsim.ImageF(100,100, scale=pixel_scale)
                        stamp2 = galsim.ImageF(100,100, scale=pixel_scale)
                        stamp3 = galsim.ImageF(100,100, scale=pixel_scale)
                        stamp4 = galsim.ImageF(100,100, scale=pixel_scale)
                        final.drawImage(image=stamp)
                        final2.drawImage(image=stamp2)
                        final3.drawImage(image=stamp3)
                        final4.drawImage(image=stamp4)
                else:
                    # Create galaxy
                    kk+=1
                    if gal_size==0 and gal_bt==0:
                        #print len(duplicate_stamps),gal_mag
                        # Take stamp from duplicate sample
                        stamp = duplicate_stamps[int(gal_mag)]
                        stamp2 = duplicate_stamps2[int(gal_mag)]
                        stamp3 = duplicate_stamps3[int(gal_mag)]
                        stamp4 = duplicate_stamps4[int(gal_mag)]
                    else:
                        # Create an original stamp
                        stamp,stamp2,stamp3,stamp4 = create_gal_stamp(x,y,
                                                            gal_mag,gal_size,gal_e1,gal_e2,gal_bt,
                                                                g1,g2,psf,pixel_scale,mag_zero)
                        # Save every object from the duplication magnitude limit onwards
                        if gal_mag>=25:
                            duplicate_stamps.append(stamp)
                            duplicate_stamps2.append(stamp2)
                            duplicate_stamps3.append(stamp3)
                            duplicate_stamps4.append(stamp4)

                            # Due to rounding off of magnitudes to save storage space, some may be
                            # rounded up from 24.495 to 24.5 and so enter in the duplicate sample.
                            # To counteract this, use the knowledge that there should be only a
                            # thousand galaxies in the duplicate sample, so take the latest 1000
                            # and those should be the right ones
                            if len(duplicate_stamps)>1000:
                                duplicate_stamps=duplicate_stamps[-1000:]
                                duplicate_stamps2=duplicate_stamps2[-1000:]
                                duplicate_stamps3=duplicate_stamps3[-1000:]
                                duplicate_stamps4=duplicate_stamps4[-1000:]


                #---------------------------- Place object in full image ---------------------------

                # Recenter the stamp at the desired position:
                stamp.setCenter(ix,iy)
                stamp2.setCenter(ix,iy)
                stamp3.setCenter(ix,iy)
                stamp4.setCenter(ix,iy)


                # Find the overlapping bounds:
                bounds = stamp.bounds & full_image.bounds
                bounds2 = stamp2.bounds & full_image2.bounds
                bounds3 = stamp3.bounds & full_image3.bounds
                bounds4 = stamp4.bounds & full_image4.bounds


                # Add galaxy stamp to full image at specified position
                full_image[bounds] += stamp[bounds]
                full_image2[bounds2] += stamp2[bounds2]
                full_image3[bounds3] += stamp3[bounds3]
                full_image4[bounds4] += stamp4[bounds4]

                #if kk==1000:
                  #print k, time.time()-t0,time.time()-t1, t1-t0, '|', t4-t3,t3-t2,t2-t1
                #  break
        #---------------------------- Finishing touches --------------------------------------------

        # RNG seeds that vary for each step in for loop
        rng = galsim.BaseDeviate(long(random_seed))
        # rng = galsim.BaseDeviate(12345)
        #print rng, rng
        # Add Gaussian noise to the image with specified sigma
        full_image.addNoise(galsim.GaussianNoise(rng, sigma=noise_sigma))
        full_image2.addNoise(galsim.GaussianNoise(rng, sigma=noise_sigma))
        full_image3.addNoise(galsim.GaussianNoise(rng, sigma=noise_sigma))
        full_image4.addNoise(galsim.GaussianNoise(rng, sigma=noise_sigma))

        pipeline.chipImage(full_image.array,
                  int(configFile.get('chipping','chips_x')),
                  int(configFile.get('chipping','chips_y')),
                  int(configFile.get('chipping','chip_x_dim')),
                  int(configFile.get('chipping','chip_y_dim')),
                  float(configFile.get('chipping','chip_gap')),
                  ditherArray,
                  exposure_number,
                  configFile,
                  tmpDir,
                  isRotated=False,
                  rotationNumber=0)

        pipeline.chipImage(full_image2.array,
                  int(configFile.get('chipping','chips_x')),
                  int(configFile.get('chipping','chips_y')),
                  int(configFile.get('chipping','chip_x_dim')),
                  int(configFile.get('chipping','chip_y_dim')),
                  float(configFile.get('chipping','chip_gap')),
                  ditherArray,
                  exposure_number,
                  configFile,
                  tmpDir,
                  isRotated=True,
                  rotationNumber=1)

        pipeline.chipImage(full_image3.array,
                  int(configFile.get('chipping','chips_x')),
                  int(configFile.get('chipping','chips_y')),
                  int(configFile.get('chipping','chip_x_dim')),
                  int(configFile.get('chipping','chip_y_dim')),
                  float(configFile.get('chipping','chip_gap')),
                  ditherArray,
                  exposure_number,
                  configFile,
                  tmpDir,
                  isRotated=True,
                  rotationNumber=2)

        pipeline.chipImage(full_image4.array,
                  int(configFile.get('chipping','chips_x')),
                  int(configFile.get('chipping','chips_y')),
                  int(configFile.get('chipping','chip_x_dim')),
                  int(configFile.get('chipping','chip_y_dim')),
                  float(configFile.get('chipping','chip_gap')),
                  ditherArray,
                  exposure_number,
                  configFile,
                  tmpDir,
                  isRotated=True,
                  rotationNumber=3)

        t2 = time.time()
        return t2 - t0


    #---------------------------- Set up to run on different CPUs ----------------------------------


    def worker(input, output):
        """
        input is a queue with (args, info) tuples:
           args are the arguements to pass to build_file
           info is passed along to the output queue.
        output is a queue storing (result, info, proc) tuples:
           result is the return value of from build_file
           info is passed through from the input queue.
           proc is the process name.
        """
        for (args, info) in iter(input.get, 'STOP'):
            result = build_file(*args)
            output.put( (result, info, current_process().name) )


    #------------------------- Set up the task list ------------------------------------------------

    # Create task list
    task_queue = Queue()

    # Extract all information from the comment lines of the prior
    catheader=[]
    nobj=0
    with open(path_to_prior_file) as fileobject:
        for line in fileobject:
            if line[0]=='#': catheader.append(line.split())
            nobj+=1

    # Get RNG seed from prior file
    rng_seed=float(catheader[0][-1])

    # Loop over all dithers
    for exp_nr in range(5):

        # Dither values from input file
        x_dither = float(catheader[1+exp_nr][-2])
        y_dither = float(catheader[1+exp_nr][-1])
        # PSF parameters from input file
        psf_params = catheader[7+2*exp_nr]
        #print psf_params
        psf_params = [float(psf_params[i][:-1]) for i in [-7,-5,-3,-1]]
        #print x_dither, y_dither, psf_params

        # We put on the task queue the args to the build_file function and
        # some extra info to pass through to the output queue.
        # Our extra info is just the file name of the catalog that we use to write out which
        # file finished.

        task_queue.put( ((g1,g2, exp_nr, x_dither,y_dither, psf_params, path_to_prior_file,
                            os.path.join(path_to_image_archive_folder,
                                                'exp'+str(exp_nr)+'.fits'),
                            os.path.join(path_to_image_archive_folder,
                                                'exp'+str(exp_nr+5)+'.fits'),
                            os.path.join(path_to_image_archive_folder,
                                                'exp'+str(exp_nr+10)+'.fits'),
                            os.path.join(path_to_image_archive_folder,
                                                'exp'+str(exp_nr+15)+'.fits'),
                            rng_seed),
                            path_to_prior_file.split('/')[-1]) )

        # Dithers all use the same prior file, so only change the random seed for
        # the noise realization. This has to be the GalSim seed, since it will (hopefully) be
        # thread safe
        rng_seed+=nobj

    #print '\n tasks have been put on task queue \n'

    # Make a queue for processes that are done
    done_queue = Queue()

    # Run the tasks
    # Each Process command starts up a parallel process that will keep checking the queue
    # for a new task. If there is one there, it grabs it and does it. If not, it waits
    # until there is one to grab. When it finds a 'STOP', it shuts down.
    # The number of processors to be used is given by the master program
    for k in range(nproc):
        procThread=Process(target=worker, args=(task_queue, done_queue))
        procThread.daemon=True
        procThread.start()

    #print '\n tasks have been initiated \n'


    # In the meanwhile, the main process keeps going.
    # This loop is happening while the other processes are still working on their tasks.
    # You'll see that these logging statements get print out as the stamp images are still
    # being drawn.
    for i in range(5):
        result, info, proc = done_queue.get()
        file_name = info
        t = result
        #print proc,file_name, 'time:', t


    #print '\n tasks are drawn from done queue \n'

    # Stop the processes
    # The 'STOP's could have been put on the task list before starting the processes, or you
    # can wait.  In some cases it can be useful to clear out the done_queue (as we just did)
    # and then add on some more tasks.  We don't need that here, but it's perfectly fine to do.
    # Once you are done with the processes, putting nproc 'STOP's will stop them all.
    # This is important, because the program will keep running as long as there are running
    # processes, even if the main process gets to the end.  So you do want to make sure to
    # add those 'STOP's at some point!
    for k in range(nproc):
        task_queue.put('STOP')


    #print '\n STOPs have been added \n'

#------------------------------------ Actual script ------------------------------------------------

if __name__ == "__main__":
    create_imsims(-0.04, 0.00, 1,'../malta/prior_stars', 'psfs', 'ims','a','a','a',20.)
