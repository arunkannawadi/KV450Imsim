__author__ = 'Ian Fenech Conti'

"""
KiDS Image Simulations.
-----------------------

Version         0.2.7
Description     KiDS Image Simulations Pipeline.
Developed by    Ian Fenech Conti (University of Malta)
                Ricardo Herbonnet (Univeristy of Leiden)
Developed on    Monday 8th June 2015
Contact         ianfc89@gmail.com

"""

'''
Imports.
'''

import numpy as np
import math
import pyfits
import time
import warnings
import os
from subprocess import Popen, PIPE
from ConfigParser import SafeConfigParser
import os
import string
import random
import optparse
#import priors_fixed_psf2 as imsimpriors
import imsimpriors as imsimpriors
# import rewrite_prior
#import imsim as imsim
import imsim_RH_sersic as imsim
import time
import sys
import shutil
import tarfile
import ast
import datetime
# import pydrizzle.xytosky as xy_conv
from pyfits import Column
import socket
warnings.simplefilter("ignore")
from astropy.io import fits
from os.path import basename

'''
Global Params.
'''
CHIPNAMES = []
TMPDIR = []
ARCHDIR = []
RUNID = []

'''
Read in a large text-file and get a
line at a specific row index.

@path  : the path to which file we
need to open.
@index : the index of the string we
wish to extract from the file.
'''
def readLineAtIndex(path, index):
    objectsCount = 0
    with open(path) as fileobject:
        for count, line in enumerate(fileobject):
            if count == index:
                return line
                break

'''
Extracts the dither values from the cataloge
file. Some string manipulation is done to
be able to extract the two values.

@ditherArray : a 2D array storing x,y dither values
for each exposure.
@startIndex  : defines which line the dither values
are defined at in the catalog file.
'''
def extractDithers():
    print '      Extracting Dither Information';sys.stdout.flush()
    path = '%s/%s' % (ARCHDIR, parser.get('priors','prior_catalog'))
    exposures = int(parser.get('chipping','n_exposures'))
    ditherArray = np.zeros((exposures, 2))
    startIndex = 1
    for exposure in xrange(0, exposures):
        dither = readLineAtIndex(path, exposure + startIndex)
        ditherArray[exposure, 0] = float(dither.split(':')[1].split()[0])
        ditherArray[exposure, 1] = float(dither.split(':')[1].split()[1])
    return ditherArray

'''
Loads an image from file. (Mainly used for
testing)
'''
def loadData(path):
    hdulist = pyfits.open(path)
    data = hdulist[0].data
    hdulist.close()
    return data

'''
Saves an image to a file. We use this to dump the single
chips to file.
'''
def saveData(data, header, path):
    hdu = pyfits.PrimaryHDU(data)
    hdulist = pyfits.HDUList([hdu])
    hdulist[0].header = header
    hdulist.writeto(path, clobber=True)
    hdulist.close()

'''
Chip Image. Given an image and a set of dither values,
this function creates a set of M,N Chips defined on a size
with a specific chip gap.
'''
def chipImage(imageData, chipX, chipY, chipDimX, chipDimY, chipGap, ditherArray, exposure, parser, tmpDir, isRotated = False, rotationNumber = 0):
    # print '      Chipping Image ID = %d (isRotated : %s)' % (exposure, isRotated);sys.stdout.flush()

    # Set the initalise the chip id.
    chipID = 1

    # Set up the temporary directory.
    global TMPDIR
    TMPDIR= tmpDir

    # Begin the chipping process.
    for yIndex in xrange(0, chipY):
        for xIndex in xrange(0, chipX):
            # Extract the chip.
            chipData = extractChip(imageData, chipGap, chipDimX, chipDimY, xIndex, yIndex)

            # Set the chip header.
            chipHeader = setHeader(ditherArray[exposure,:], chipID, imageData.shape, xIndex, yIndex, chipDimX, chipDimY, chipGap, parser)

            # Define the chip name.
            chipName = 'exp%dchip_%dOFCS.sub' % (exposure, chipID)

            if isRotated == False:
                # Save chip data.
                saveData(chipData, chipHeader, '%s%s/%s.fits' % (TMPDIR, parser.get('directories','chip_directory') ,chipName))

                # Save chip weights.
                weights = chipData.copy()
                weights[:,:] = parser.get('chipping','default_weight')
                weightsName = '%s.weight' % (chipName)
                saveData(weights, chipHeader, '%s%s/%s.fits' % (TMPDIR, parser.get('directories','weight_directory'), weightsName))

                # Save chip .head file.
                headName = 'exp%dchip_%d.head' % (exposure, chipID)
                lensfitHeader('%s%s/%s' % (TMPDIR, parser.get('directories','head_directory'), headName), chipHeader)

                # Only append the chip names once. Avoid duplication for rotated images.
                CHIPNAMES.append(chipName)

            if isRotated == True:
                # Save chip data. Here we only need to save the actual chip data. The rest
                # is already done for the non rotated images.
                saveData(chipData, chipHeader, '%s%s%02d/%s.fits' % (TMPDIR, parser.get('directories','chiprot_directory'), rotationNumber ,chipName))

            # Increment the chipID.
            chipID += 1

'''
Prepare the list of chip names for lensift input.
'''
def prepareInputFile():

    # Begin the file naming process.
    for exposure in xrange(0, int(parser.get('chipping','n_exposures'))):
        # Set the initalise the chip id.
        chipID = 1
        for yIndex in xrange(0, int(parser.get('chipping','chips_y'))):
            for xIndex in xrange(0, int(parser.get('chipping','chips_x'))):
                # Define the chip name.
                chipName = 'exp%dchip_%dOFCS.sub' % (exposure, chipID)
                CHIPNAMES.append(chipName)
                # Increment the chipID.
                chipID += 1

'''
Extract Chip. Extract's the data of a given exposure into a
number of chips defined by layout, dimensions and chip gap.
'''
def extractChip(imageData, chipGap, chipDimX, chipDimY, xIndex, yIndex):
    # Determine the chip gap.
    spillX = chipGap * xIndex
    spillY = chipGap * yIndex

    # Get the start and ending indices to crop on x axis.
    xStart = (chipDimX * xIndex) + spillX
    xEnd = (xStart + chipDimX)

    # Get the start and ending indices to crop on y axis.
    yStart = (chipDimY * yIndex) + spillY
    yEnd = (yStart + chipDimY)

    # Crop the data.
    chipData = imageData[yStart:yEnd, xStart:xEnd]

    return chipData

'''
Set Chip Header. Set the chip header to incorporate
the dither, gap and ImageID in prep for lensfit.
'''
def setHeader(ditherValue, chipID, imageDimensions, xIndex, yIndex, chipDimX, chipDimY, chipGap, parser):

    fitsHDU = pyfits.open(parser.get('chipping_io','sample_header'))
    fitsHeader = fitsHDU[0].header

    fitsHeader.set('CTYPE1', 'RA---TAN')
    fitsHeader.set('CTYPE2', 'DEC--TAN')

    try: ## THIS WAS THE ORIGINAL CODE THAT STOPPED WORKING ON JULY 22, 2017
        fitsHeader.update('EXPTIME', 360)
        fitsHeader.update('PHOTPLAM', 0.000000)
        fitsHeader.update('PHOTZPT', 0.000000)
        fitsHeader.update('PHOTFLAM', 0.000000)
    except ValueError: ## set the same instead of update, by AKJ
        fitsHeader.set('EXPTIME', 360)
        fitsHeader.set('PHOTPLAM', 0.000000)
        fitsHeader.set('PHOTZPT', 0.000000)
        fitsHeader.set('PHOTFLAM', 0.000000)

    fitsHeader.set('RADECSYS', 'FK5     ')

    CRVAL1 = 3.450000000E+01
    CRVAL2 = -7.000000000E+00

    fitsHeader.set('CRVAL1', CRVAL1)
    fitsHeader.set('CRVAL2', CRVAL2)

    fitsHeader.set('GAIN', 1)

    naxis1 = imageDimensions[1]
    naxis2 = imageDimensions[0]

    ditherX = -ditherValue[0]
    ditherY = -ditherValue[1]

    spillX = chipGap * xIndex
    spillY = chipGap * yIndex

    xStart = (chipDimX * xIndex) + spillX
    yStart = (chipDimY * yIndex) + spillY

    fitsHeader.set('IMAGEID', chipID)
    fitsHeader.set('CRPIX1', ((naxis1 / 2.0) + ditherX) - xStart)
    fitsHeader.set('CRPIX2', ((naxis2 / 2.0) + ditherY) - yStart)

    # GET THE GS_SCALE AND CONVERT TO RA/DEC BASED ON THE PIXEL SCALE
    fitsHeader.set('GS_SCALE', 0.214)
    gs_scale = fitsHeader['GS_SCALE']

    cd_value = gs_scale / 3600
    fitsHeader.set('CD1_1', -cd_value)
    fitsHeader.set('CD1_2', 0.0000000000000)
    fitsHeader.set('CD2_1', 0.0000000000000)
    fitsHeader.set('CD2_2', +cd_value)

    fitsHeader.set('PV1_0', 0.0000000000000)
    fitsHeader.set('PV1_1', 1.0000000000000)
    fitsHeader.set('PV1_2', 0.0000000000000)
    fitsHeader.set('PV1_4', 0.0000000000000)
    fitsHeader.set('PV1_5', 0.0000000000000)
    fitsHeader.set('PV1_6', 0.0000000000000)
    fitsHeader.set('PV1_7', 0.0000000000000)
    fitsHeader.set('PV1_8', 0.0000000000000)
    fitsHeader.set('PV1_9', 0.0000000000000)
    fitsHeader.set('PV1_10', 0.0000000000000)

    fitsHeader.set('PV2_0', 0.0000000000000)
    fitsHeader.set('PV2_1', 1.0000000000000)
    fitsHeader.set('PV2_2', 0.0000000000000)
    fitsHeader.set('PV2_4', 0.0000000000000)
    fitsHeader.set('PV2_5', 0.0000000000000)
    fitsHeader.set('PV2_6', 0.0000000000000)
    fitsHeader.set('PV2_7', 0.0000000000000)
    fitsHeader.set('PV2_8', 0.0000000000000)
    fitsHeader.set('PV2_9', 0.0000000000000)
    fitsHeader.set('PV2_10', 0.0000000000000)

    try:
        fitsHeader.update('CUNIT1', 'deg     ')
        fitsHeader.update('CUNIT2', 'deg     ')
    except ValueError:
        fitsHeader.set('CUNIT1', 'deg     ')
        fitsHeader.set('CUNIT2', 'deg     ')

    fitsHDU.close()

    return fitsHeader

'''
Create Headfile. Prepare a LensFIT .head file
from the FITS output header which is used to read
in the Astrometry.

@chipHeader : a PYFITS object header with all the chip
'''
def lensfitHeader(path, chipHeader):
    f = open(path, 'w')
    for line in chipHeader.cards:
        ## line has cardimage attribute until JULY 22, 2017
        #f.write(line.cardimage + "\n")
        f.write(line.image + "\n") ## change by AKJ
    f.close()

'''
Create Input File. Prepare a LensFIT input asci file
containing a list of all the chip names.

@chipNames : a list with each of the elements corresponding
to the LensfitFormat name specification.
'''
def createInputFile():
    print '      Creating LensFIT Input File';sys.stdout.flush()
    f = open('%s%s' % (TMPDIR, parser.get('lensfit','input_file')), 'w')
    for chipName in CHIPNAMES:
        f.write(chipName + "\n")
    f.close()

'''
Prepares the LensFIT type asci catalog. With rows
corresponding to.
'''
def createLensFITCatalog(rot_number):
    print '      Creating LensFIT Catalog for rotation: %s' % rot_number
    if rot_number==0:
        sextractorCatalog = '%s%s' % (TMPDIR, parser.get('sextractor','cataloge_path'))
        outputPath = '%s%s' % (TMPDIR, parser.get('lensfit','input_catalog'))
    else:
        sextractorCatalog = '%s%s%02d%s' % (TMPDIR, 'sexrot', rot_number, '.cat')
        outputPath = TMPDIR+'catalog_rot0{0}.asc'.format(rot_number) # A HACK FOR NOW
    data = np.loadtxt(sextractorCatalog,comments='#')
    f = open(outputPath, 'w')
    for row in data:
        f.write('%f %f %f %f \n' % (row[2], row[3], row[6], row[13]))
    f.close()

'''
Prepares the LensFIT type asci catalog. Generated
from the priors file. RH.
'''
def createFitsTable(pathPriorFile, pathSaveFile):
    c1 = Column(name='x', format='D')
    c2 = Column(name='y', format='D')
    c3 = Column(name='magnitude', format='D')
    col_definitions = pyfits.ColDefs([c1, c2, c3])

    x,y,m,r,f = np.loadtxt(pathPriorFile,usecols=(0,1,2,3,6),unpack=True)

    sel = np.invert([((f==0)&(r==0))|(m>24.5)][0])

    xsel=x[sel]
    ysel=y[sel]
    msel=m[sel]

    image_hdu = pyfits.new_table(col_definitions, nrows=len(xsel))

    # for c, tiled_position in enumerate(zip(xsel, ysel, msel)):
    for index, (value1, value2, value3) in enumerate(zip(xsel, ysel, msel)):
        image_hdu.data[index] = [((value1)),
                             ((value2)),
                             (value3)]

    sys.stdout = open(os.devnull, "w")
    image_hdu.writeto(pathSaveFile, clobber=True)
    sys.stdout = sys.__stdout__

'''
Prepares the LensFIT type asci catalog. Generated
from the priors file. RH.
'''
def createLensFITCatalogXYSKY(pathPriorFitsFile, pathCatalog, pathMatchFile):
    sys.stdout = open(os.devnull, "w")
    output = xy_conv.XYtoSky_pars('%s[0]' % parser.get('priors','reference_header'),
                                  None, None,
                                  pathPriorFitsFile,
                                  'x,y', xy_conv.yes,
                                  'IDCTAB', xy_conv.no, None,
                                  xy_conv.no)
    sys.stdout = sys.__stdout__


    fx = pyfits.open(pathPriorFitsFile, memmap=False)
    d = fx[1].data

    f = open(pathCatalog, 'w')
    f2 = open(pathMatchFile, 'w')

    for index, (RA, DEC) in enumerate(zip(output[0], output[1])):
        f.write("%s %s %s %s\n" % (RA, DEC, d[index][2], index+1))
        f2.write("%s %s %s %s %s %s\n" % (d[index][0], d[index][1],RA, DEC, d[index][2], index+1))

    f.close()
    f2.close()

'''
Call LensFIT using subprocess command, set all the
relevant environement variables and pass the relevant args.
'''
def flensfit(isRotated = False, rotationNumber = 0):
    print '      Running LensFIT (isRotated : %s, rotationNumber : %02d)' % (isRotated, rotationNumber), ;sys.stdout.flush()
    lensfitPath = parser.get('lensfit','lensfit_path')
    startTime = time.time()
    envVars = os.environ.copy()
    envVars['SWARP_CONFIG'] = parser.get('lensfit','swarp_path')
    if sexonrot==0 or (sexonrot==1)&(rotationNumber==0):
        ## Take the same SExtractor -> input catalogue irrespective of the rotationNumber
        CATALOGUE_GALAXIES = '%s%s' % (TMPDIR, parser.get('lensfit','input_catalog'))
    else:
        CATALOGUE_GALAXIES = TMPDIR+'catalog_rot0{0}.asc'.format(rotationNumber) # A HACK FOR NOW
    envVars['CATALOGUE_GALAXIES'] = CATALOGUE_GALAXIES

    envVars['PSF_DIR'] = '%s%s/' % (TMPDIR, parser.get('directories','psf_directory'))
    envVars['HEAD_DIR'] = '%s%s/' % (TMPDIR, parser.get('directories','head_directory'))
    envVars['PSF_OVERSAMPLING'] = parser.get('lensfit','psf_oversampling')
    envVars['PRIOR_PARAMETERS'] = parser.get('lensfit','prior_path')
    envVars['PECUT'] = '0.02'
    envVars['PRCUT'] = '0.02'
    envVars['LCUT'] = '0.05'
    envVars['WAVEBAND'] = 'R'
    envVars['CAMERA'] = 'KIDS'
    # envVars['SKIP'] = '1'

    if isRotated == False:
        outputPath = '%s/%s' % (ARCHDIR, parser.get('lensfit','output_file'))
        stdoutPath = '%s/%s.log' % (ARCHDIR, parser.get('lensfit','output_file'))
        stderrPath = '%s/%s.err.log' % (ARCHDIR, parser.get('lensfit','output_file'))
        envVars['DATA_DIR'] = '%s%s/' % (TMPDIR, parser.get('directories','chip_directory'))
    else:
        outputPath = '%s/%02d.%s' % (ARCHDIR, rotationNumber, parser.get('lensfit','outputrot_file'))
        stdoutPath = '%s/%02d.%s.log' % (ARCHDIR, rotationNumber, parser.get('lensfit','outputrot_file'))
        stderrPath = '%s/%02d.%s.err.log' % (ARCHDIR, rotationNumber, parser.get('lensfit','outputrot_file'))
        envVars['DATA_DIR'] = '%s%s%02d/' % (TMPDIR, parser.get('directories','chiprot_directory'), rotationNumber)

    outLog = open(stdoutPath, "w")
    errLog = open(stderrPath, "w")

    process = Popen([lensfitPath,
                     '%s%s' % (TMPDIR, parser.get('lensfit','input_file')),
                     outputPath,
                     parser.get('lensfit','postage_size'),
                     parser.get('lensfit','start_exposure'),
                     parser.get('lensfit','end_exposure'),
                     parser.get('lensfit','start_mag'),
                     parser.get('lensfit','end_mag')],
                     stdout=outLog, stderr=errLog, env=envVars, cwd=TMPDIR)

    stdout, stderr = process.communicate()

    outLog.flush()
    outLog.close()

    errLog.flush()
    errLog.close()
    print ' [%s]' % (str(datetime.timedelta(seconds=(time.time()-startTime))))

'''
Convert GalSIM PSF to a quadrant swapped centered
PSF in the format.
'''
def convertPSF(exposure):
    converterPath = parser.get('psf_convert','convert_path')
    galsimPSFPath = '%s%s/%s%d.fits' % (TMPDIR,
                                     parser.get('directories','galsim_psf_directory'),
                                     parser.get('psf_convert','galsim_name'),
                                     exposure)
    outputPSFPath = '%s%s/exp%dchip.psfcoeffs.fits' % (TMPDIR,
                                                       parser.get('directories','psf_directory'),
                                                       exposure)
    process = Popen([converterPath, galsimPSFPath, outputPSFPath], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()

'''
Call Swarp using subprocess command, configure the correct
params file depending on a weight file coadd or a chip file coadd.
'''
def swarpChips(isWeights=False, chip_rot_name='chip', rot_number=0):
    print '      Running SWarp (isWeights : %s)' % isWeights
    print '      Running SWarp on rotation: %s' % rot_number
    sys.stdout.flush()
    startTime = time.time()
    processPath = parser.get('swarp','swarp_path')
    if isWeights:
      if rot_number==0: ## Do this only for one rotation, since the weights are all the same anyway
        workingDirectory = '%s%s/' % (TMPDIR, parser.get('directories','weight_directory'))
        paramPath = parser.get('swarp','weights_param')
        imagePath = '%s%s/%s' % (TMPDIR, parser.get('directories','weight_directory'), parser.get('swarp','weights_wildcard'))
        savePath = '%s%s/%s' % (TMPDIR, parser.get('directories','weight_directory'), parser.get('swarp','weights_save'))
        ''' Saves the useless output weights file '''
        swarpWeightPath = '%s%s/%s.weight' % (TMPDIR, parser.get('directories','weight_directory'), parser.get('swarp','weights_save'))
      else:
          ## We called this by mistake. Just exit
          return None
    else:
        if rot_number>0:
            chip_rot_name += 'rot0' + str(rot_number)
        workingDirectory = '%s%s/' % (TMPDIR, chip_rot_name) 
        paramPath = parser.get('swarp','chips_param')
        imagePath = '%s%s/%s' % (TMPDIR, chip_rot_name, parser.get('swarp','chips_wildcard'))
        savePath = '%s%s/%s' % (TMPDIR, chip_rot_name, parser.get('swarp','chips_save'))
        ''' Saves the useless output weights file '''
        swarpWeightPath = '%s%s/%s.weight' % (TMPDIR, chip_rot_name, parser.get('swarp','chips_save'))

    swarpRun = '%s %s -c %s -IMAGEOUT_NAME %s -WEIGHTOUT_NAME %s' % (processPath, imagePath, paramPath, savePath, swarpWeightPath)
    print 'The swarp command is '
    print swarpRun

    process = Popen(swarpRun, stdout=PIPE, stderr=PIPE, shell=True, cwd=workingDirectory)
    stdout, stderr = process.communicate()
    print ' [%s]' % (str(datetime.timedelta(seconds=(time.time()-startTime))))

'''
Call SExtractor using subprocess command, configure
the right paramaters file and input args.
'''
def sextractorPositions(chip_rot_name='chip', rot_number=0):
    print '      Running SEXtractor on rotation:    %s' % (rot_number)
    sys.stdout.flush()
    startTime = time.time()
    if rot_number > 0:
        chip_rot_name = chip_rot_name + 'rot0' + str(rot_number)
    processPath = parser.get('sextractor','sextractor_path')
    imagePath = '%s%s/%s' % (TMPDIR, chip_rot_name, parser.get('swarp','chips_save'))
    ## The weights are same regardless of the rotation
    weightingPath ='%s%s/%s' % (TMPDIR, parser.get('directories','weight_directory'), parser.get('swarp','weights_save'))
    paramsPath = parser.get('sextractor','sex_param')
    configPath = parser.get('sextractor','sex_config')
    #catalogPath = '%s%s' % (TMPDIR, parser.get('sextractor','cataloge_path'))
    if rot_number==0:
        catalogPath = '%s%s' % (TMPDIR, 'sex.cat')
        catalogPathArchive = ARCHDIR+'/sex.cat'
    else:
        catalogPath = TMPDIR+'sexrot0'+str(rot_number)+'.cat'
        catalogPathArchive = ARCHDIR+'/sexrot0'+str(rot_number)+'.cat'
    catalogPathArchive = '%s/%s' % (ARCHDIR, parser.get('sextractor','cataloge_path'))
    workingDirectory = parser.get('directories','config_directory')

    sexRun = '%s %s -WEIGHT_IMAGE %s -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s ' % \
             (processPath,
              imagePath,
              weightingPath,
              configPath,
              paramsPath,
              catalogPath)
    print ''
    print sexRun ## Uncomment this if you want to know what command it ran
    process = Popen(sexRun, stdout=PIPE, stderr=PIPE, shell=True, cwd=workingDirectory)
    stdout, stderr = process.communicate()
    
    # Copy the sextractor catalog to the archive.
    shutil.copy(catalogPath, catalogPathArchive)
    print ' [%s]' % (str(datetime.timedelta(seconds=(time.time()-startTime))))

def createRandomKeyFolder(randomKey):
    # For the archive folder
    randomKeyDirectory = '%s%s' % (parser.get('directories','archive_directory'), randomKey)
    configDirectory = '%s/%s' % (randomKeyDirectory, '.misc')
    try:
        os.mkdir(randomKeyDirectory)
        os.mkdir(configDirectory)
    except OSError:
        print "This random key already exists."

    # Copy the config files
    configPath = parser.get('sextractor','sex_config')
    shutil.copy(configPath, configDirectory)

    # Copy the source code for reference
    pipeline_abspath = os.path.abspath(__file__) ## get self absolute path
    shutil.copy(pipeline_abspath, configDirectory) ## save itself first
    imsim_abspath = pipeline_abspath.replace(__file__, imsim.__file__)
    shutil.copy(imsim_abspath, configDirectory) ## save the image simulation code as well

    # For the temp folder
    randomKeyDir = '%s%s' % (parser.get('directories','root_directory'), randomKey)
    try:
        os.mkdir(randomKeyDir)
    except OSError:
        print "This random key already exists."

    return randomKeyDirectory

'''
Temp Directory Structure. The IMSIM Pipeline will create
a set of directories and subdirectories in order to handle
with the data creation. These will be removed once done.
'''
def createWorkFolder(g1, g2, psfSet):
    print '      Creating Temp. Directories';sys.stdout.flush()

    # Load the Root Directory from the parser.
    rootDirectory = parser.get('directories','root_directory')

    # Add shear information.
    g1Part = ''
    if g1 > 0:
        g1Part += 'p'
    else:
        g1Part += 'm'
    g1Part += '%03d' % abs((g1*10000.1))
    g2Part = ''
    if g2 > 0:
        g2Part += 'p'
    else:
        g2Part += 'm'
    g2Part += '%03d' % abs((g2*10000.1))

    # Set up the RUNID.
    runID = '%s%s_%s_%s' % (g1Part, g2Part, psfSet, randomKey)

    # Create the sub-directories.
    TMPDIR = '%s%s/%s/' % (rootDirectory, randomKey, runID)
    chipDir = '%s%s' % (TMPDIR, parser.get('directories','chip_directory'))
    headDir = '%s%s' % (TMPDIR, parser.get('directories','head_directory'))
    psfDir = '%s%s' % (TMPDIR, parser.get('directories','psf_directory'))
    weightsDir = '%s%s' % (TMPDIR, parser.get('directories','weight_directory'))
    galsimPSFDir = '%s%s' % (TMPDIR, parser.get('directories','galsim_psf_directory'))

    try:
        os.mkdir(TMPDIR)
        os.mkdir(chipDir)
        os.mkdir(headDir)
        os.mkdir(psfDir)
        os.mkdir(weightsDir)
        os.mkdir(galsimPSFDir)
    except:
        print "The working directories already exist."

    # Make the rotated folders.
    nRot = int(parser.get('imsim','n_rot'))
    for ii in range(0, nRot):
        chipDirRot = '%s%s%02d' % (TMPDIR, parser.get('directories','chiprot_directory'), ii+1)
        try:
            os.mkdir(chipDirRot)
        except:
            print "%s already exists." % (chipDirRot)

    return TMPDIR, runID

'''
Create a folder in the data archive. This will only store
the priors file and the LensFIT output files and these will
also be compressed to save space.
'''
def createArchiveFolder(runID):
    print '      Creating Archive Directory';sys.stdout.flush()
    archiveDirectory = '%s%s/%s' % (parser.get('directories','archive_directory'), randomKey, runID)
    try:
        os.mkdir(archiveDirectory)
    except:
        print "The archive directory already exists."
    return archiveDirectory

'''
Compresses an archive folder into tar.gz. No compression levels
are passed as it was found default compression level is the optimal.
'''
def compressArchive():
    print '      Compressing Archive';sys.stdout.flush()
    with tarfile.open('%s%s.tar.gz' % (parser.get('directories','archive_directory'), RUNID), "w:gz") as tar:
        tar.add(ARCHDIR, arcname=os.path.basename(ARCHDIR))

'''
Creates a ciruclar range of g1, g2 values with a given |g|
and a number of data points, n
'''
def createG1G2Range(modG,n=32):
    return np.round(np.array([(math.cos(2*np.pi/n*x)*modG,math.sin(2*np.pi/n*x)*modG) for x in xrange(0,n+1)]), 4)

'''
Prints out a g1, g2 range to be used for multiple runs
args.
'''
def printG1G2Range():
    g1g2 = createG1G2Range(0.04)
    print g1g2[:, 0].tolist()
    print g1g2[:, 1].tolist()
    exit()

def cut_priorfile_on_area(path_to_priorfile,path_to_save_location,xlimit=4000,ylimit=4000):
	header_inpriorfile = open(path_to_priorfile,'r')
	header_info=''; faint_mag_ind=-1
	for line in header_inpriorfile:
		if line[0]=='#':
			# Remove the hashtag otherwise double
			header_info +=line[2:]
		else:
			faint_mag_ind+=1
			if line.split()[3]=='0' and line.split()[6]=='0':
				break
	the_priorfile = np.loadtxt(path_to_priorfile,unpack=True)
	reduced_priorfile = the_priorfile[:,:(faint_mag_ind-1000)]
	area_selection = (reduced_priorfile[0]<xlimit) & (reduced_priorfile[1]<ylimit)
	reduced_priorfile = reduced_priorfile[:,area_selection]
	new_priorfile = np.concatenate((reduced_priorfile,the_priorfile[:,(faint_mag_ind-1000):]),axis=1)
	# Remove last end-of-line in header, otherwise empty hashtag line
	np.savetxt(path_to_save_location,new_priorfile.T,fmt='%.3f',header=header_info[:-1])

'''
Convert ASCII Output to FITS Table
'''
def AsciiToFits(catalogue_name):
	
	data = np.loadtxt('%s/%s' % (ARCHDIR, catalogue_name))
	c1 = fits.Column(name='wcsx', format='D', array=data[:,0])
	c2 = fits.Column(name='wcsy', format='D', array=data[:,1])
	c3 = fits.Column(name='bias-corrected <e1>', format='E', array=data[:,2])
	c4 = fits.Column(name='bias-corrected <e2>', format='E', array=data[:,3])
	c5 = fits.Column(name='weight', format='E', array=data[:,4])
	c6 = fits.Column(name='fitclass', format='J', array=data[:,5])
	c7 = fits.Column(name='bias-corrected scalelength /pixels', format='E', array=data[:,6])
	c8 = fits.Column(name='bulge-fraction', format='E', array=data[:,7])
	c9 = fits.Column(name='model-flux', format='E', array=data[:,8])
	c10 = fits.Column(name='pixel SNratio', format='E', array=data[:,9])
	c11 = fits.Column(name='model SNratio', format='D', array=data[:,10])
	c12 = fits.Column(name='contamination radius', format='E', array=data[:,11])
	c13 = fits.Column(name='PSF-e1', format='E', array=data[:,12])
	c14 = fits.Column(name='PSF-e2', format='E', array=data[:,13])
	c15 = fits.Column(name='PSF-Strehl-ratio', format='E', array=data[:,14])
	c16 = fits.Column(name='PSF-Q11', format='E', array=data[:,15])
	c17 = fits.Column(name='PSF-Q22', format='E', array=data[:,16])
	c18 = fits.Column(name='PSF-Q12', format='E', array=data[:,17])
	c19 = fits.Column(name='star-galaxy f-probability', format='E', array=data[:,18])
	c20 = fits.Column(name='r correction', format='E', array=data[:,19])
	c21 = fits.Column(name='2D measurement variance', format='E', array=data[:,20])
	c22 = fits.Column(name='mean-likelihood |e|', format='E', array=data[:,21])
	c23 = fits.Column(name='e1 correction', format='E', array=data[:,22])
	c24 = fits.Column(name='e2 correction', format='E', array=data[:,23])
	c25 = fits.Column(name='neighbour mag', format='E', array=data[:,24])
	c26 = fits.Column(name='neighbour distance', format='E', array=data[:,25])
	c27 = fits.Column(name='catmag', format='E', array=data[:,26])
	c28 = fits.Column(name='n-exposures-used', format='J', array=data[:,27])
	c29 = fits.Column(name='cat-ID', format='J', array=data[:,28])
	c30 = fits.Column(name='old weight', format='E', array=data[:,29])

	coldefs = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10,
				c11, c12, c13, c14, c15, c16, c17, c18,
				c19, c20, c21, c22, c23, c24, c25, c26, 
				c27, c28, c29, c30])

	tbhdu = fits.BinTableHDU.from_columns(coldefs)

	tbhdu.writeto('%s/%s.fits' % (ARCHDIR, catalogue_name), overwrite=True)

def append_n_column(prior):
    """ Takes in a prior file that imsimpriors generated and add two extra columns - for n and z.
    """
    newlines = []
    with open(prior,'r') as f:
        lines = f.readlines()

    for line in lines:
        if not line[0]=='#':
            newline = line[:-1]+"   -1   -1\n" ## strip the newline char and add it later
        else: ## it is the header
            if line[1]=='x':
                line = line[:-1]+'    n    z\n'
            newline = line
        newlines.append(newline)

    with open(prior,'w') as f:
        f.writelines(newlines)

def make_it_sersic(prior):
    """ Takes in a prior file that imsimpriors generated and see B/T==9 while filling in n
    """
    newlines = []
    with open(prior,'r') as f:
        lines = f.readlines()

    for line in lines:
        if not line[0]=='#':
            words = line.split()
            words[6] = '9' ## B/T==9
            ## TO DO: Replace by random values
            words[7] = '5.0' ## Sersic n
            newline = ' '.join(words)+'\n'
        else:
            newline = line

        newlines.append(newline)

    with open(prior,'w') as f:
        f.writelines(newlines)


def fetch_random_seed(RUNID):
        prefix = RUNID.rpartition('_')[0]
        parent_directory = '/disks/shear12/herbonnet/kids/data/vault_03052016_newmagdistr/'
        subdirectory = [subdir for subdir in os.listdir(parent_directory) if prefix in subdir and not '.tar.gz' in subdir][0]

        ## The following functionality is available within imsimpriors as well, so might as well return the full path of the subdirectory
        with open(os.path.join(parent_directory,subdirectory,'prior'),'r') as f:
            line = f.readline()
        random_seed = line.split()[-1]
        print "Fetch %s as the random seed from " % random_seed, os.path.join(parent_directory, subdirectory, 'prior')
        return long(random_seed)

def rng_generator(base_seed, g1, g2, psfset_id):
    import cmath

    g = complex(g1,g2)
    ## Phase of g takes [-3,-2,-1,0,1,2,3,4]*np.pi/4

    ## All 5x416 square degrees get different rng_seed, so that the noise realisations don't repeat irrespective of parallelization
    rng_seed = base_seed + int(3+cmath.phase(g)*4/np.pi) + 10*psfset_id
    return rng_seed
    
"""
Main Method.
--------------------------------------
"""
if __name__ == '__main__':

    print ''
    print '  ------------------------------------------------------------------------- '
    print '  IMSIM Pipeline Started '
    print '  Version 0.3.1'
    print '  Current Version (Full Run, Sextractor + Stars + Faint Gals + 5 Fixed PSFs)'
    print '  Started at : %s ' % time.strftime("%c")
    print '  Running on : %s' % socket.gethostname()
    print '  ------------------------------------------------------------------------- '

    # Set up run-time args.
    argopts = optparse.OptionParser("usage: %prog [options] arg1")
    argopts.add_option("-a", "--g1", dest="g1range",
                       default="0.037",
                       type="string",
                       help="Specify [Start, Stop, Step] for the g1 values")
    argopts.add_option("-b", "--g2", dest="g2range",
                       default="0.0153",
                       type="string",
                       help="Specify [Start, Stop, Step] for the g2 values")
    argopts.add_option("-n", "--noise", dest="noiseSigma",
                       default="17.03",
                       type="string",
                       help="Specify Noise Sigma Level")
    argopts.add_option("-c", "--config", dest="config",
                       default="/home/ian/Documents/KIDS/config/config.ini",
                       type="string",
                       help="Specify Config File for IMSIM Pipeline.")
    argopts.add_option("-r", "--rmtemp", dest="rmtemp",
                       default="True",
                       help="Specify if Temp Folder is to be deleted.")
    argopts.add_option("-d", "--ramdisk", dest="ramdisk",
                       default="False",
                       help="Specify if Use of RAMDISK")
    argopts.add_option("-p", "--psfrange", dest="psfrange",
                       default="0,1,2,3,4",
                       help="Specify the range of PSF runs")
    argopts.add_option("-k", "--randomKey", dest="randomKey",
                        default="testing",
                        help="Specify the Job ID through a randomKey")
    argopts.add_option("-s", "--skipblock", dest="skipblock",
                        default="-1",
                        help="Specify the blocks you want to skip.")

    (options, args) = argopts.parse_args()
    configPath = options.config

    # Load Configuration Settings.
    print '  Loading Config from %s' % configPath
    parser = SafeConfigParser()
    parser.read(configPath)

    # Define a unique and global random key for a run.
    global randomKey
    randomKey = options.randomKey

    ## RandomKey must not contain an underscore, otherwise, only a part will be extracted in later operations.
    if '_' in randomKey:
        print "The value for random key contained an underscore. Removing it."
        randomKey = randomKey.replace('_','')
        print "Updated random key is", randomKey

    # Create the randomKey directories
    randomKeyDirectory = createRandomKeyFolder(randomKey)

    # Check if user specifices to remove the temporary directory.
    removeTemp = ast.literal_eval(options.rmtemp)
    print '  Removing TMP Files = %s' % removeTemp

    # Check if user specifices to use ram disk.
    ramDisk = ast.literal_eval(options.ramdisk)
    print '  Using RAM Disk = %s' % ramDisk

    # Start the timer
    t0 = time.time()

    # Loop through g1, g2 range and start sims.
    g1Range = options.g1range.split(',')
    g1Range = [float(i) for i in g1Range]
    g2Range = options.g2range.split(',')
    g2Range = [float(i) for i in g2Range]

    # Get noise sigma level.
    noise_sigma = float(options.noiseSigma)

    # Get the PSF Range.
    psfRange = options.psfrange.split(',')
    psfRange = [int(i) for i in psfRange]

    # Get the skipblocks
    skipBlock = options.skipblock.split(',')
    skipBlock = [int(i) for i in skipBlock]

    print '  Random key          = %s' %randomKey
    print '  Running Range of g1 = %s' % g1Range
    print '  Running Range of g2 = %s' % g2Range
    print '  Running PSF Sets    = %s' % psfRange
    print '  Noise Sigma         = %f' % noise_sigma

    # Get some prior parameters
    realistic = bool(int(parser.get('priors', 'realistic')))
    if realistic:
        randomize_positions = bool(int(parser.get('priors', 'randomize_positions')))
        randomize_orientations = bool(int(parser.get('priors', 'randomize_orientations')))
    else:
        randomize_positions, randomize_orientations = True, True ## does not matter

    print '  Realistic priors? :     ', realistic
    print '  Randomizing positions? :   ', randomize_positions
    print '  Randomizing_orientations? :    ', randomize_orientations

    # Check if user specifices to run SExtractor on all rotations
    sexonrot = int(parser.get('sextractor', 'sexonrot'))
    if sexonrot==1:
        print '   Running SExtractor on all the rotations.'
    elif sexonrot==0:
        print '   Running SExtractor only on rotation 0.'

    # Check if user specifices to use only Sersic or mixed models
    sersic_only = bool(int(parser.get('imsim', 'sersic_only')))
    if sersic_only:
        print '   Generating only Sersic models.'
    else:
        print '   Generating Sersic and B+D models.'

    for psfSet in psfRange:

        for g1g2 in zip(g1Range, g2Range):

            print '    Running for g1, g2 = %.4f \t %.4f [PSF Set : %d]' % (g1g2[0], g1g2[1], psfSet)

            # Create unique temp directory
            TMPDIR, RUNID = createWorkFolder(g1g2[0], g1g2[1], psfSet)

            # Create data archive entry.
            ARCHDIR = createArchiveFolder(RUNID)
            print '      Run ID = %s' % RUNID

            g1g2ShearRange = RUNID.split('_')[0]
            g1g2ShearPSF = RUNID.split('_')[1]

            # priorFilePrevious = '/users/ianfc/kids/priors4/prior_%s_%s' % (g1g2ShearRange, g1g2ShearPSF)
            # priorFilePrevious = '/users/ianfc/priors/bright_small_gals_prior_psf2'
            # Call priors code to generate a set of priors.
            print '      Generating Priors'
            startTime = time.time()

            ## Fetch the random seed only if you want to reproduce the previous runs. Else, see it to None
            if psfSet<5 and (g1g2[0]**2+g1g2[1]**2 > 0.03**2): ## either the PSFs in the vault or if we are running 0 shear tests
                random_seed = fetch_random_seed(RUNID)
            else:
                base_seed = 1234567
                random_seed = rng_generator(base_seed, g1g2[0], g1g2[1], psfSet)
            #random_seed = None

##AKJ:  Copy the prior file from REP030516 to force a comparison
#            input_prior = '/disks/shear15/KiDS/ImSim/pipeline/archive/REP030516/'+RUNID.rpartition('_')[0]+'_REP030516/prior'
#            prior_prefix = '/disks/shear12/herbonnet/kids/data/vault_03052016_newmagdistr/'
#            input_prior_list = [prior_prefix+ip+'/prior' for ip in os.listdir(prior_prefix) if RUNID.rpartition('_')[0] in ip and not '.tar.gz' in ip#]
#            assert len(input_prior_list)==1
#            input_prior = input_prior_list[0]
## AKJ: Turn off creating priorfile to force a comparison

            if not 1 in skipBlock: ## Functional block 1 -  Generate the prior file
                imsimpriors.create_priorfile(psfSet,
                                             parser.get('priors', 'besancon_path'),
                                             '%s/%s' % (ARCHDIR, parser.get('priors', 'prior_catalog')),
                                             parser.get('priors', 'psfset_path'),
                                             grid_positions=False,
                                             nr_of_stars=True,
                                             nr_of_bright_gals=True,
                                             nr_of_faint_gals=True,
                                             random_seed=random_seed,
                                             randomize_positions=randomize_positions,
                                             randomize_orientations=randomize_orientations,
                                             realistic=realistic
                                             )
                if not realistic:
                    ## Append a dummy column, since the Sersic n column will be missing, which the imsim routine expects
                    append_n_column('%s/%s' %(ARCHDIR, parser.get('priors', 'prior_catalog')))

                if realistic & sersic_only:
                    ## Randomly plug in Sersic N values and set B/T == 9
                    make_it_sersic('%s/%s' %(ARCHDIR, parser.get('priors', 'prior_catalog')))

    #            # rewrite_prior.rewrite_old_prior(priorFilePrevious,'%s/%s' % (ARCHDIR, parser.get('priors', 'prior_catalog')))
                # shutil.copyfile(priorFilePrevious, '%s/%s' % (ARCHDIR, parser.get('priors', 'prior_catalog')))
                # cut_priorfile_on_area(priorFilePrevious, '%s/%s' % (ARCHDIR, parser.get('priors', 'prior_catalog')))
    #            print ' [%s]' % (str(datetime.timedelta(seconds=(time.time() - startTime))))

                # Generate the input catalogue for Lensfit from the priors file.
                '''
                print '      Generating Lensfit Input File'
                createFitsTable('%s/%s' % (ARCHDIR, parser.get('priors', 'prior_catalog')),
                                '%s/%s.fits' % (ARCHDIR, parser.get('priors', 'prior_catalog')))

                createLensFITCatalogXYSKY('%s/%s.fits' % (ARCHDIR, parser.get('priors', 'prior_catalog')),
                                          '%s%s' % (TMPDIR, parser.get('lensfit','input_catalog')),
                                          '%s/%s.match.asci' % (ARCHDIR, parser.get('priors', 'prior_catalog')))
                
                '''

                # Flush the stdoutput
                sys.stdout.flush()

            # Extract the dither information
            ditherArray = extractDithers()

            if not 2 in skipBlock: ## Functional block 2 - Generate the images
                # Call the image generator code. Generates 10 images,
                # 5 at nominal position and 5 at 90deg rotatation.
                print '      Rendering Images ** USING NO IMSIM REMORPH **',
                startTime = time.time()
                imsim.create_imsims(g1g2[0], g1g2[1],
                                    int(parser.get('imsim', 'n_threads')),
                                    '%s/%s' % (ARCHDIR, parser.get('priors', 'prior_catalog')),
                                    '%s%s/' % (TMPDIR, parser.get('directories', 'galsim_psf_directory')),
                                    'REMOVE THIS PARAMATER',
                                    parser, ditherArray, TMPDIR, noise_sigma)
                print ' [%s]' % (str(datetime.timedelta(seconds=(time.time() - startTime))))

                # Flush the stdoutput
                sys.stdout.flush()

                # Convert the PSFs
                for exposure in xrange(0, int(parser.get('chipping', 'n_exposures'))):
                    convertPSF(exposure)

                # Create a LensFIT Input file, clear the list of chip names.
                prepareInputFile()
                createInputFile()
                CHIPNAMES = []

                '''
                Removed the Swarping and Sextractor for now. We now create the input
                catalogue directly from the truth priors catalogue.
                '''

            if not 3 in skipBlock: ## Functional block 3 - Run SWarp & SExtractor
                # Swarp Chips/#Weights
                swarpChips(isWeights=True)
                if sexonrot==0:
                    swarpChips(isWeights=False, chip_rot_name='chip')
                else:
                    for nRot in range(0, 1+int(parser.get('imsim', 'n_rot'))):
                        swarpChips(isWeights=False, chip_rot_name='chip', rot_number=nRot)

                # Sextractor Positions
                if sexonrot==1:
                    for nRot in range(0, 1+int(parser.get('imsim', 'n_rot'))):
                        sextractorPositions(chip_rot_name='chip', rot_number=nRot)
                else:
                    sextractorPositions(chip_rot_name='chip', rot_number=0)

                # Create LensFIT Catalogue
                if sexonrot==1:
                    for nRot in range(0, 1+int(parser.get('imsim', 'n_rot'))):
                        createLensFITCatalog(nRot)
                else:
                    createLensFITCatalog(0)

            ## If it passed all the 3 functional blocks, run lensfit
            if not 4 in skipBlock: ## Functional block 4 - Run lensfit

                # Run LensFIT for Nominal and Rotated Images
                flensfit(isRotated=False, rotationNumber=0)
                for nRot in range(0, int(parser.get('imsim', 'n_rot'))):
                    flensfit(isRotated=True, rotationNumber=nRot + 1)

                # Apply the weight recalibration.
                weight_recal_script = parser.get('lensfit', 'weights_recal_path')

                command = 'python %s --input=%s/output.fits.asc' % (weight_recal_script, ARCHDIR)
                process = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
                stdout, stderr = process.communicate()

                for nRot in range(0, int(parser.get('imsim', 'n_rot'))):
                    command = 'python %s --input=%s/%02d.output.rot.fits.asc' % (weight_recal_script, ARCHDIR, nRot + 1)
                    process = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
                    stdout, stderr = process.communicate()

                    # Convert ASCII to FITS Table
    #            AsciiToFits('output.fits.asc.scheme2b_corr')

    #            for nRot in range(0, int(parser.get('imsim', 'n_rot'))):
    #                AsciiToFits('%02d.output.rot.fits.asc.scheme2b_corr' % (nRot+1))

                # Compress the data into an archive
                compressArchive()
        
                # Remove the temp folders if specificed to do so - PERMANENTLY DISABLED FOR NOW
    #            if removeTemp:
    #                shutil.rmtree(TMPDIR)
        
                # Remove uncompressed archive folder - PERMANENTLY DISABLED FOR NOW
    #            shutil.rmtree(ARCHDIR)
    
            # Copy over archive to vault
            if ramDisk:
                archivePathRamDisk = '%s%s.tar.gz' % (parser.get('directories', 'archive_directory'), RUNID)
                archivePathVault = '%s%s.tar.gz' % (parser.get('directories', 'vault_directory'), RUNID)
                shutil.move(archivePathRamDisk, archivePathVault)
                os.remove(archivePathRamDisk)

    # Stop timer
    t1 = time.time()
    print ''
    print '  ------------------------------------------------------------------------- '
    print '  IMSIM Pipeline Completed in (%.4fs) ' % (t1 - t0)
    print '  Finished at : %s ' % time.strftime("%c")
    print '  ------------------------------------------------------------------------- '
