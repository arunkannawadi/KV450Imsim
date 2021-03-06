__author__ = 'Ian Fenech Conti, Arun Kannawadi'

"""
KiDS Image Simulations.
-----------------------

Version         0.2.8
Description     KiDS Image Simulations Pipeline.
Developed by    Ian Fenech Conti (University of Malta)
                Ricardo Herbonnet (Univeristy of Leiden)
                Arun Kannawadi (University of Leiden)
Developed on    Monday 8th June 2015
Modified on     Wednesday 26th April 2017
Contact         ianfc89@gmail.com, arunkannawdi@strw.leidenuniv.nl

"""

'''
Imports.
'''
#import pdb; pdb.set_trace() # Do not have parallelization on. Do not run from a script
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
import priors_fixed_psf2 as imsimpriors
# import rewrite_prior
import imsim as imsim
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
from multiprocessing import Pool, Process, Queue

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
    path = '%s/%s' % (ARCHDIR, 'ditherArray.txt')
    ditherArray = np.loadtxt(path, usecols=(1,2))
    #exposures = int(parser.get('chipping','n_exposures'))
    #ditherArray = np.zeros((exposures, 2))
    #startIndex = 1
    #for exposure in xrange(0, exposures):
        #dither = readLineAtIndex(path, exposure + startIndex)
        #ditherArray[exposure, 0] = float(dither.split(':')[1].split()[0])
        #ditherArray[exposure, 1] = float(dither.split(':')[1].split()[1])
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
Chip Exposures
'''
def ChipExposures(TMPDIR, ditherArray, rot_id):
    if rot_id==0:
        # Chip nominal the nominal image.
        for nExp in range(0, int(parser.get('chipping', 'n_exposures'))):
            image = loadData('%s%s/exp%d.fits' % (TMPDIR, parser.get('directories', 'exp_directory'), nExp))

            chipImage(image,
                      int(parser.get('chipping', 'chips_x')),
                      int(parser.get('chipping', 'chips_y')),
                      int(parser.get('chipping', 'chip_x_dim')),
                      int(parser.get('chipping', 'chip_y_dim')),
                      float(parser.get('chipping', 'chip_gap')),
                      ditherArray,
                      nExp,
                      isRotated=False,
                      rotationNumber=0,
                      TMPDIR=TMPDIR)

    else:
        # Chip rotated images.         
        for nExp in range(0, int(parser.get('chipping', 'n_exposures'))):
            # Load the nominal image.
            image = loadData('%s%s/exp%d_rot%02d.fits' % (TMPDIR, parser.get('directories', 'exp_directory'), nExp, rot_id))
            chipImage(image,
                      int(parser.get('chipping', 'chips_x')),
                      int(parser.get('chipping', 'chips_y')),
                      int(parser.get('chipping', 'chip_x_dim')),
                      int(parser.get('chipping', 'chip_y_dim')),
                      float(parser.get('chipping', 'chip_gap')),
                      ditherArray,
                      nExp,
                      isRotated=True,
                      rotationNumber=rot_id,
                      TMPDIR=TMPDIR)

'''
Chip Image. Given an image and a set of dither values,
this function creates a set of M,N Chips defined on a size
with a specific chip gap.
'''
def chipImage(imageData, chipX, chipY, chipDimX, chipDimY, chipGap, ditherArray, exposure, isRotated = False, rotationNumber = 0, TMPDIR=None):
    # print '      Chipping Image ID = %d (isRotated : %s)' % (exposure, isRotated);sys.stdout.flush()

    # Set the initalise the chip id.
    chipID = 1

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
def prepareInputFile(CHIPNAMES):

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

    return CHIPNAMES

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

    fitsHeader.update('EXPTIME', 360)
    fitsHeader.update('PHOTPLAM', 0.000000)
    fitsHeader.update('PHOTZPT', 0.000000)
    fitsHeader.update('PHOTFLAM', 0.000000)

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

    fitsHeader.update('CUNIT1', 'deg     ')
    fitsHeader.update('CUNIT2', 'deg     ')

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
        f.write(line.cardimage + "\n")
    f.close()

'''
Create Input File. Prepare a LensFIT input asci file
containing a list of all the chip names.

@chipNames : a list with each of the elements corresponding
to the LensfitFormat name specification.
'''
def createInputFile(TMPDIR, CHIPNAMES):
    print '      Creating LensFIT Input File';sys.stdout.flush()
    f = open('%s%s' % (TMPDIR, parser.get('lensfit','input_file')), 'w')
    for chipName in CHIPNAMES:
        f.write(chipName + "\n")
    f.close()

'''
Prepares the LensFIT type asci catalog. With rows
corresponding to.
'''
def createLensFITCatalog(TMPDIR, rot_id=0):
    print '      Creating LensFIT Catalog'

    if rot_id==0:
        sextractorCatalog = '%s%s' % (TMPDIR, parser.get('sextractor','cataloge_path'))
        outputPath = '%s%s' % (TMPDIR, parser.get('lensfit','input_catalog'))
    else:
        sextractorCatalog = '%s%s%02d%s' % (TMPDIR, 'sexrot', rot_id, '.cat')
        outputPath = TMPDIR+'catalog_rot0{0}.asc'.format(rot_id) # A HACK FOR NOW
    data = pyfits.getdata(sextractorCatalog)
    reduced_data = np.array([data['X_WORLD'], data['Y_WORLD'], data['MAG_AUTO'], data['NUMBER']]).T
    np.savetxt(outputPath, reduced_data)

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
def flensfit(TMPDIR, ARCHDIR, isRotated = False, rotationNumber = 0):
    print '      Running LensFIT (isRotated : %s, rotationNumber : %02d)' % (isRotated, rotationNumber), ;sys.stdout.flush()
    lensfitPath = parser.get('lensfit','lensfit_path')
    startTime = time.time()
    envVars = os.environ.copy()
    envVars['SWARP_CONFIG'] = parser.get('lensfit','swarp_path')
    if rotationNumber==0:
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
def swarpChips(TMPDIR, isWeights = False, chip_rot_name='chip', rot_number=0):
    print '      Running SWarp (isWeights : %s)' % isWeights, ;sys.stdout.flush()
    print ' TMPDIR = ', TMPDIR; sys.stdout.flush()
    startTime = time.time()
    processPath = parser.get('swarp','swarp_path')
    if isWeights:
      if rot_number==0:
        workingDirectory = '%s%s/' % (TMPDIR, parser.get('directories','weight_directory'))
        paramPath = parser.get('swarp','weights_param')
        imagePath = '%s%s/%s' % (TMPDIR, parser.get('directories','weight_directory'), parser.get('swarp','weights_wildcard'))
        savePath = '%s%s/%s' % (TMPDIR, parser.get('directories','weight_directory'), parser.get('swarp','weights_save'))
        ''' Saves the useless output weights file '''
        swarpWeightPath = '%s%s/%s.weight' % (TMPDIR, parser.get('directories','weight_directory'), parser.get('swarp','weights_save'))
    else:
        if rot_number > 0:
            chip_rot_name = chip_rot_name + 'rot0' + str(rot_number)
        workingDirectory = '%s%s/' % (TMPDIR, chip_rot_name)
        paramPath = parser.get('swarp','chips_param')
        imagePath = '%s%s/%s' % (TMPDIR, chip_rot_name, parser.get('swarp','chips_wildcard'))
        savePath = '%s%s/%s' % (TMPDIR, chip_rot_name, parser.get('swarp','chips_save'))
        ''' Saves the useless output weights file '''
        swarpWeightPath = '%s%s/%s.weight' % (TMPDIR, chip_rot_name, parser.get('swarp','chips_save'))

    swarpRun = '%s %s -c %s -IMAGEOUT_NAME %s -WEIGHTOUT_NAME %s' % (processPath, imagePath, paramPath, savePath, swarpWeightPath)
    #print ''
    #print swarpRun
    process = Popen(swarpRun, stdout=PIPE, stderr=PIPE, shell=True, cwd=workingDirectory)
    stdout, stderr = process.communicate()
    print ' [%s]' % (str(datetime.timedelta(seconds=(time.time()-startTime))))

'''
Call SExtractor using subprocess command, configure
the right paramaters file and input args.
'''
def sextractorPositions(TMPDIR, ARCHDIR, chip_rot_name='chip', rot_number=0):
    if rot_number > 0:
        chip_rot_name = chip_rot_name + 'rot0' + str(rot_number)
    print '      Running SEXtractor',;sys.stdout.flush()
    startTime = time.time()
    processPath = parser.get('sextractor','sextractor_path')
    imagePath = '%s%s/%s' % (TMPDIR, chip_rot_name, parser.get('swarp','chips_save'))
    weightingPath ='%s%s/%s' % (TMPDIR, parser.get('directories','weight_directory'), parser.get('swarp','weights_save'))
    paramsPath = parser.get('sextractor','sex_param')
    configPath = parser.get('sextractor','sex_config')
    if rot_number==0:
        catalogPath = '%s%s' % (TMPDIR, parser.get('sextractor','cataloge_path'))
    else:
        catalogPath = '%s%s%02d%s' % (TMPDIR, 'sexrot',rot_number,'.cat') # A HACK FOR NOW
    catalogPathArchive = '%s/%s' % (ARCHDIR, parser.get('sextractor','cataloge_path'))
    workingDirectory = parser.get('directories','config_directory')

    sexRun = '%s %s -WEIGHT_IMAGE %s -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s ' % \
             (processPath,
              imagePath,
              weightingPath,
              configPath,
              paramsPath,
              catalogPath)
    #print ''
    #print sexRun
    process = Popen(sexRun, stdout=PIPE, stderr=PIPE, shell=True, cwd=workingDirectory)
    stdout, stderr = process.communicate()
    
    # Copy the sextractor catalog to the archive.
    shutil.copy(catalogPath, catalogPathArchive)
    print ' [%s]' % (str(datetime.timedelta(seconds=(time.time()-startTime))))

'''
Temp Directory Structure. The IMSIM Pipeline will create
a set of directories and subdirectories in order to handle
with the data creation. These will be removed once done.
'''
def createWorkFolder(g1, g2, psfSet, randomKey=None):
    print '      Creating Temp. Directories';sys.stdout.flush()

    # Load the Root Directory from the parser.
    rootDirectory = parser.get('directories','root_directory')

    if randomKey is None:
        # Define a unique random key for a run.
        randomKey = ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(10))

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
    expDir = '%s%s' % (TMPDIR, parser.get('directories','exp_directory'))

    try:
        os.mkdir(TMPDIR)
        os.mkdir(chipDir)
        os.mkdir(headDir)
        os.mkdir(psfDir)
        os.mkdir(weightsDir)
        os.mkdir(galsimPSFDir)
        os.mkdir(expDir)

        # Make the rotated folders.
        nRot = int(parser.get('imsim','n_rot'))
        for ii in range(0, nRot):
            chipDirRot = '%s%s%02d' % (TMPDIR, parser.get('directories','chiprot_directory'), ii+1)
            os.mkdir(chipDirRot)

    except OSError:
        pass

    return TMPDIR, runID
'''
Create a folder to hold the data archive. This will also store the config files used for the run.
'''
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
Create a folder in the data archive. This will only store
the priors file and the LensFIT output files and these will
also be compressed to save space.
'''
def createArchiveFolder(runID):
    print '      Creating Archive Directory';sys.stdout.flush()
    randomKey = runID.split('_')[-1]
    archiveDirectory = '%s%s/%s' % (parser.get('directories','archive_directory'), randomKey, runID)
    try:
        os.mkdir(archiveDirectory)
    except OSError:
        pass
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

def combine_galaxy_star_catalogues(path_input_cat, path_star_cat,path_output_cat):
    galaxy_cat = fits.open(path_input_cat)
    star_cat = fits.open(path_star_cat)

    galaxy_data = galaxy_cat[1].data
    star_data = star_cat[1].data

    rank = galaxy_data['rank']

    x_offset, y_offset = 1627.71, 646.21

    X = np.append(star_data['X_IMAGE']-x_offset, galaxy_data['Xpos_THELI']-x_offset)
    Y = np.append(star_data['Y_IMAGE']-y_offset, galaxy_data['Ypos_THELI']-y_offset)
    MAG = np.append(star_data['MAG'], galaxy_data['MAG_AUTO_THELI']*(rank>0) + galaxy_data['MAGR']*(rank==0))
    SERSICN = np.append(-np.ones(len(star_data)), galaxy_data['N_GALFIT_HI'])
    RE = np.append(-np.ones(len(star_data)), galaxy_data['RE_GALFIT_HI'])
    q = np.append(np.ones(len(star_data)), galaxy_data['BA_GALFIT_HI'])
    phi = np.append(np.zeros(len(star_data)), (galaxy_data['PA_GALFIT_HI']+90.)*np.pi/180)
    photoz = np.append(np.zeros(len(star_data)), galaxy_data['Z_B'])
    SeqNr = np.append(-np.ones(len(star_data)), galaxy_data['SeqNr'])

    flag = np.array( [-1]*len(star_data)+[2]*(rank>0).sum()+[0]*(rank==0).sum() )

    c1 = fits.Column(name='X_IMAGE',format='D',array=X,unit='pix')
    c2 = fits.Column(name='Y_IMAGE',format='D',array=Y,unit='pix')
    c3 = fits.Column(name='MAG_AUTO',format='D',array=MAG)
    c4 = fits.Column(name='N_GALFIT_HI',format='D',array=SERSICN)
    c5 = fits.Column(name='RE_GALFIT_HI',format='D',array=RE,unit='pix')
    c6 = fits.Column(name='BA_GALFIT_HI',format='D',array=q)
    c7 = fits.Column(name='PA_GALFIT_HI',format='D',array=phi,unit='radians')
    c8 = fits.Column(name='Z_B',format='D',array=photoz)
    c9 = fits.Column(name='FLAG',format='J',array=flag)
    c10 = fits.Column(name='SeqNr',format='J',array=SeqNr)

    cols = fits.ColDefs([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10])
    tbhdu = fits.BinTableHDU.from_columns(cols)

    tbhdu.writeto(path_output_cat)


def pipeline1(g1g2, psfSet, rot_id, randomKey, q=None):

    # Create unique temp directory
    TMPDIR, RUNID = createWorkFolder(g1g2[0], g1g2[1], psfSet, randomKey)

    # Create data archive entry.
    ARCHDIR = createArchiveFolder(RUNID)

    #g1g2ShearRange = RUNID.split('_')[0] # UNKNOWN FUNCTION
    #g1g2ShearPSF = RUNID.split('_')[1] # UNKNOWN FUNCTION

    print '      Rendering Images for rotation = ', rot_id

    # Flush the stdoutput
    sys.stdout.flush()
    startTime = time.time()

    # Call the image generator code.
    imsim.create_imsims(psfSet, g1g2[0], g1g2[1],
                        galCatPath, starCatPath,
                        '%s/%s' % (ARCHDIR, 'ditherArray.txt'),
                        '%s%s/' % (TMPDIR, parser.get('directories', 'galsim_psf_directory')),
                        '%s%s/' % (TMPDIR, parser.get('directories', 'exp_directory')),
                        n_gal=True, stars=stars, faint_gal=False,
                        rot_id=rot_id, n_rotations=n_rot, sersic_only=sersic_only)

    print ' [%s]' % (str(datetime.timedelta(seconds=(time.time() - startTime))))

    # Flush the stdoutput
    sys.stdout.flush()

    # Extract the dither information
    ditherArray = extractDithers()

    # Flush the stdoutput
    sys.stdout.flush()

    # Chip the exposures
    print '      Chipping Exposures',
    startTime = time.time()
    ChipExposures(ditherArray, rot_id)
    print ' [%s]' % (str(datetime.timedelta(seconds=(time.time() - startTime))))

    # Flush the stdoutput
    sys.stdout.flush()

    # Convert the PSFs
    print '      Converting PSFs'
    for exposure in xrange(0, int(parser.get('chipping', 'n_exposures'))):
        convertPSF(exposure)

    # Create a LensFIT Input file, clear the list of chip names.
    CHIPNAMES = []
    prepareInputFile()
    createInputFile()
    CHIPNAMES = []

    if rot_id==0:
        ## To do this for all rot_id are not is something to be discussed.

        # Swarp Chips/Weights
        swarpChips(isWeights=False, chip_rot_name='chip')
        swarpChips(isWeights=True)

        # Sextractor Positions
        sextractorPositions(chip_rot_name='chip', rot_number='00')

        # Create LensFIT Catalogue
        createLensFITCatalog()

    if q is not None:
        q.put(1)

def pipeline2(g1g2, psfSet, rot_id, randomKey):

    # Create unique temp directory
    TMPDIR, RUNID = createWorkFolder(g1g2[0], g1g2[1], psfSet, randomKey)

    # Create data archive entry.
    ARCHDIR = createArchiveFolder(RUNID)

    # Run LensFIT
    flensfit(isRotated=bool(rot_id), rotationNumber=rot_id)

    # Apply the weight recalibration.
    weight_recal_script = parser.get('lensfit', 'weights_recal_path')

    if rot_id == 0:
        command = 'python %s --input=%s/output.fits.asc' % (weight_recal_script, ARCHDIR)
    else:
        command = 'python %s --input=%s/%02d.output.rot.fits.asc' % (weight_recal_script, ARCHDIR, rot_id)
    process = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = process.communicate()

    # Convert ASCII to FITS Table
    if rot_id == 0:
        AsciiToFits('output.fits.asc.scheme2b_corr')
    else:
        AsciiToFits('%02d.output.rot.fits.asc.scheme2b_corr' % (rot_id))

def pipeline(g1g2, psfSet, rot_id, randomKey, q=None):

    # Create unique temp directory
    TMPDIR, RUNID = createWorkFolder(g1g2[0], g1g2[1], psfSet, randomKey)

    # Create data archive entry.
    ARCHDIR = createArchiveFolder(RUNID)

    #g1g2ShearRange = RUNID.split('_')[0] # UNKNOWN FUNCTION
    #g1g2ShearPSF = RUNID.split('_')[1] # UNKNOWN FUNCTION

    print '      Rendering Images for rotation = ', rot_id

    # Flush the stdoutput
    sys.stdout.flush()
    startTime = time.time()

    # Call the image generator code.
    imsim.create_imsims(psfSet, g1g2[0], g1g2[1],
                        galCatPath, starCatPath,
                        '%s/%s' % (ARCHDIR, 'ditherArray.txt'),
                        '%s%s/' % (TMPDIR, parser.get('directories', 'galsim_psf_directory')),
                        '%s%s/' % (TMPDIR, parser.get('directories', 'exp_directory')),
                        n_gal=True, stars=stars, faint_gal=False,
                        psf_params_pathname=psf_params_pathname,
                        rot_id=rot_id, n_rotations=n_rot, sersic_only=sersic_only,
                        parallelize=False, internal_parallelize=0)

    print ' [%s]' % (str(datetime.timedelta(seconds=(time.time() - startTime))))

    # Flush the stdoutput
    sys.stdout.flush()

    # Extract the dither information
    ditherArray = extractDithers()

    # Flush the stdoutput
    sys.stdout.flush()

    # Chip the exposures
    print '      Chipping Exposures',
    startTime = time.time()
    ChipExposures(TMPDIR, ditherArray, rot_id)
    print ' [%s]' % (str(datetime.timedelta(seconds=(time.time() - startTime))))

    # Flush the stdoutput
    sys.stdout.flush()

    # Convert the PSFs
    print '      Converting PSFs'
    for exposure in xrange(0, int(parser.get('chipping', 'n_exposures'))):
        convertPSF(exposure)

    # Create a LensFIT Input file, clear the list of chip names.
    CHIPNAMES = []
    CHIPNAMES = prepareInputFile(CHIPNAMES)
    createInputFile(TMPDIR, CHIPNAMES)
    CHIPNAMES = []

    # Swarp Chips/Weights (Weights only for rot_id==0)
    swarpChips(TMPDIR, isWeights=False, chip_rot_name='chip',rot_number=rot_id)
    if rot_id==0:
        swarpChips(TMPDIR, isWeights=True)

    # Sextractor Positions
    sextractorPositions(TMPDIR, ARCHDIR, chip_rot_name='chip', rot_number=rot_id)

    # Create LensFIT Catalogue
    createLensFITCatalog(TMPDIR, rot_id)

    # Run LensFIT
    flensfit(TMPDIR, ARCHDIR, isRotated=bool(rot_id), rotationNumber=rot_id)

    # Apply the weight recalibration.
    weight_recal_script = parser.get('lensfit', 'weights_recal_path')

    if rot_id == 0:
        command = 'python %s --input=%s/output.fits.asc' % (weight_recal_script, ARCHDIR)
    else:
        command = 'python %s --input=%s/%02d.output.rot.fits.asc' % (weight_recal_script, ARCHDIR, rot_id)
    process = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = process.communicate()

    # Convert ASCII to FITS Table
    if rot_id == 0:
        AsciiToFits('output.fits.asc.scheme2b_corr')
    else:
        AsciiToFits('%02d.output.rot.fits.asc.scheme2b_corr' % (rot_id))

    if q is not None:
        q.put(1)


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

    (options, args) = argopts.parse_args()
    configPath = options.config

    # Load Configuration Settings.
    print '  Loading Config from %s' % configPath
    parser = SafeConfigParser()
    parser.read(configPath)

    # Define a unique random key for a run.
    #randomKey = ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(10))
    #randomKey = 'LC0U9WKPPF' # for testing
    #randomKey = 'JVHX2VWORO' # for Massimo
    randomKey = options.randomKey

    ## RandomKey must not contain an underscore, otherwise, only a part will be extracted in later operations.
    if '_' in randomKey:
        print "The value for random key contained an underscore. Removing it."
        randomKey = randomKey.replace('_','')
        print "Updated random key is", randomKey

    # Create the randomKey directories
    randomKeyDirectory = createRandomKeyFolder(randomKey)

    # Copy the config.ini to the randomKeyDirectory
    shutil.copy(configPath, '%s/%s' % (randomKeyDirectory,'.misc'))

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

    # Get noise sigma leve.
    noise_sigma = float(options.noiseSigma)

    # Get the PSF Range.
    psfRange = options.psfrange.split(',')
    psfRange = [int(i) for i in psfRange]

    # Check if rotations are on.
    rot_on = int(parser.get('imsim', 'rot'))
    n_rot = int(parser.get('imsim', 'n_rot'))

    # Check if we are simulating B+D models along with Sersic
    sersic_only = bool(int(parser.get('imsim', 'sersic_only')))
    print "  Are we using only Sersic profiles?: ", sersic_only
    if sersic_only is True:
        galCatPath = parser.get('imsim', 'sersic_catalogue')
    else:
        galCatPath = parser.get('imsim', 'BD_catalogue')

    print "  Using {0} as the input galaxy catalogue ".format(galCatPath)

    stars = bool(int(parser.get('imsim', 'stars')))
    print "  Are we including stars? :", stars

    if stars:
        starCatPath = parser.get('imsim', 'star_catalogue')
        print "  Using {0} as the input star catalogue ".format(starCatPath)

    else:
        starCatPath = None

    ## Set whether to parallelize the pipeline (recommended) or not
    parallelize_pipeline = False

    print '  Random key          = %s' %randomKey
    print '  Running Range of g1 = %s' % g1Range
    print '  Running Range of g2 = %s' % g2Range
    print '  Running PSF Sets    = %s' % psfRange
    print '  Noise Sigma         = %f' % noise_sigma
    print '  Rotations Flag      = %d' % rot_on

    for g1g2 in zip(g1Range, g2Range):

        for psfSet in psfRange:

            print '    Running for g1, g2 = %.4f \t %.4f [PSF Set : %d]' % (g1g2[0], g1g2[1], psfSet)

            # Create unique temp directory
            TMPDIR, RUNID = createWorkFolder(g1g2[0], g1g2[1], psfSet, randomKey)

            # Create data archive entry.
            ARCHDIR = createArchiveFolder(RUNID)
            print '      Run ID = %s' % RUNID

            g1g2ShearRange = RUNID.split('_')[0] #UNKNOWN FUNCTION
            g1g2ShearPSF = RUNID.split('_')[1] #UNKNOWN FUNCTION

            # Flush the stdoutput
            sys.stdout.flush()

            if parallelize_pipeline is False:
#                for rot_id in xrange(0):
                for rot_id in xrange(n_rot+1):
                    # Call the image generator code.
                    print '      Rendering Images'; sys.stdout.flush()

                    startTime = time.time()
                    imsim.create_imsims(psfSet, g1g2[0], g1g2[1],
                                        galCatPath, starCatPath,
                                        '%s/%s' % (ARCHDIR, 'ditherArray.txt'),
                                        '%s%s/' % (TMPDIR, parser.get('directories', 'galsim_psf_directory')),
                                        '%s%s/' % (TMPDIR, parser.get('directories', 'exp_directory')),
                                        psf_params_pathname = parser.get('priors', 'psfset_path'),
                                        n_gal=True, stars=stars, faint_gal=False,
                                        rot_id=rot_id, sersic_only=sersic_only, parallelize=False)

                    print ' [%s]' % (str(datetime.timedelta(seconds=(time.time() - startTime))))

                    # Flush the stdoutput
                    sys.stdout.flush()

                    # Extract the dither information
                    ditherArray = extractDithers()

                    # Flush the stdoutput
                    sys.stdout.flush()

                    # Chip the exposures
                    print '      Chipping Exposures',
                    startTime = time.time()
                    ChipExposures(TMPDIR, ditherArray, rot_id)
                    print ' [%s]' % (str(datetime.timedelta(seconds=(time.time() - startTime))))

                    # Flush the stdoutput
                    sys.stdout.flush()

                    # Convert the PSFs
                    print '      Converting PSFs'
                    for exposure in xrange(0, int(parser.get('chipping', 'n_exposures'))):
                        convertPSF(exposure)

                    # Create a LensFIT Input file, clear the list of chip names.
                    CHIPNAMES = []
                    prepareInputFile(CHIPNAMES)
                    createInputFile(TMPDIR, CHIPNAMES)
                    CHIPNAMES = []

                    '''
                    Removed the Swarping and Sextractor for now. We now create the input
                    catalogue directly from the truth priors catalogue.
                    '''

                    # Swarp Chips/Weights
                    swarpChips(TMPDIR, isWeights=False, chip_rot_name='chip', rot_number=rot_id)
                    if rot_id==0:
                        swarpChips(TMPDIR, isWeights=True)

                    # Sextractor Positions
                    sextractorPositions(TMPDIR, ARCHDIR, chip_rot_name='chip', rot_number=rot_id)

                    # Create LensFIT Catalogue
                    createLensFITCatalog(TMPDIR, rot_id)

                # Run LensFIT for Nominal and Rotated Images
                flensfit(TMPDIR, ARCHDIR, isRotated=False, rotationNumber=0)
                if rot_on == 1:
                    for nRot in range(0, int(parser.get('imsim', 'n_rot'))):
                        flensfit(TMPDIR, ARCHDIR, isRotated=True, rotationNumber=nRot + 1)

                # Apply the weight recalibration.
                weight_recal_script = parser.get('lensfit', 'weights_recal_path')

                command = 'python %s --input=%s/output.fits.asc' % (weight_recal_script, ARCHDIR)
                process = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
                stdout, stderr = process.communicate()

                if rot_on == 1:
                    for nRot in range(0, int(parser.get('imsim', 'n_rot'))):
                        command = 'python %s --input=%s/%02d.output.rot.fits.asc' % (weight_recal_script, ARCHDIR, nRot + 1)
                        process = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
                        stdout, stderr = process.communicate()

                # Convert ASCII to FITS Table
                AsciiToFits('output.fits.asc.scheme2b_corr')

                if rot_on == 1:
                    for nRot in range(0, int(parser.get('imsim', 'n_rot'))):
                        AsciiToFits('%02d.output.rot.fits.asc.scheme2b_corr' % (nRot+1))

                # Compress the data into an archive
                compressArchive()

                # Remove the temp folders if specificed to do so
                if removeTemp:
                    shutil.rmtree(TMPDIR)

                # Remove uncompressed archive folder
#                shutil.rmtree(ARCHDIR)

                # Copy over archive to vault
                if ramDisk:
                    archivePathRamDisk = '%s%s.tar.gz' % (parser.get('directories', 'archive_directory'), RUNID)
                    archivePathVault = '%s%s.tar.gz' % (parser.get('directories', 'vault_directory'), RUNID)
                    shutil.move(archivePathRamDisk, archivePathVault)
                    os.remove(archivePathRamDisk)

    if parallelize_pipeline is True:
        ## WARNING: UNCHECKED
        procs_all = []
        n_workers = len(g1Range)*len(psfRange)*(n_rot+1)
        p = Pool(processes=n_workers)
        for psfSet in psfRange:
            for g1g2 in zip(g1Range, g2Range):
                for rot_id in xrange(n_rot+1):
                    proc = p.apply_async(func=pipeline, args=(g1g2, psfSet, rot_id, randomKey))
                    procs_all.append(proc)

        p.close()
        p.join()

        #procs_all = []
        ##p = Pool(processes=64)
        #for psfSet in psfRange:
        #    for g1g2 in zip(g1Range, g2Range):
        #        for rot_id in xrange(n_rot+1):
        #            proc = Process(target=pipeline1, args=(g1g2, psfSet, rot_id, randomKey))
        #            #proc = p.apply_async(func=pipeline1, args=(g1g2, psfSet, rot_id, randomKey))
        #            procs_all.append(proc)

        ### Wait for 'all' image simulations to finish before running lensfit on 'all' fields

        ##useless_results = [q.get() for proc in procs_all]

        #procs2_all = []
        ##p2 = Pool(processes=64)
        #for psfSet in psfRange:
        #    for g1g2 in zip(g1Range, g2Range):
        #        for rot_id in xrange(n_rot+1):
        #            proc2 = Process(target=pipeline2, args=(g1g2, psfSet, rot_id, randomKey))
        #            #proc2 = p2.apply_async(func=pipeline2, args=(g1g2, psfSet, rot_id, randomKey))
        #            procs2_all.append(proc2)

    # Stop timer
    t1 = time.time()
    print ''
    print '  ------------------------------------------------------------------------- '
    print '  IMSIM Pipeline Completed in (%.4fs) ' % (t1 - t0)
    print '  Finished at : %s ' % time.strftime("%c")
    print '  ------------------------------------------------------------------------- '
