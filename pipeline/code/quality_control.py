import math
import numpy as np
import pyfits as pf
import galsim
import os, sys
import time
import scipy
from scipy.spatial import cKDTree
## Astropy modules
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
## Matplotlib modules
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
#debug = bool(int(sys.argv[2]))
#if debug:
#import pdb; pdb.set_trace()


def check_consistency(randomKey, psfIDs=[0,1,2,3,4]):
    gRange = ['p400m000','m400m000','m000p400', 'm000m400',
              'm283m283','p283m283','m283p283', 'p283p283' ]

    ## Load the input/true WCS
    input_imgname = '/disks/shear14/KiDS_simulations/Cosmos/Theli_image/KIDS_150p1_2p2_r_SDSS.V0.5.9A.swarp.cut.fits'
    wcs_input = galsim.AstropyWCS(input_imgname)

    ## Load the input/truth catalogue
    input_catname = '/disks/shear14/KiDS_simulations/Cosmos/KIDS_HST_cat/KiDS_Griffith_iMS1_handpicked_stars.cat'
    input_catalogue = fits.open(input_catname)
    input_data = input_catalogue[1].data
    print "Loaded the input data"

    ## Obtain the cuts on the input/truth catalogue
    MASK_all = input_data.MASK
    mask = ~np.array(MASK_all&0xfc3c,dtype=bool)
    handpicked_stars = input_data['handpicked_stars']
    rank = input_data['rank']
    distance2d = input_data['distance2d']
    assert handpicked_stars.dtype==bool

    cuts = mask&(rank>=0)&(distance2d<1)&(~handpicked_stars)

    OBJNO = input_data['OBJNO'][cuts]
    RA = input_data['RA'][cuts]
    DEC = input_data['DEC'][cuts]

    X,Y = [],[]
    x_offset, y_offset = 2500, 2500
    ## Converting the sky position to image positions (takes about a minute)...
    ## This is needed because Xpos_THELI and Ypos_THELI aren't filled for the faint galaxies
    for gg in xrange(cuts.sum()):
        pos = wcs_input.posToImage(galsim.CelestialCoord(RA[gg]*galsim.degrees, DEC[gg]*galsim.degrees))
        x, y = pos.x - x_offset, pos.y - y_offset

        X.append(x)
        Y.append(y)

    X = np.array(X)
    Y = np.array(Y)
    ## "Building a kd-tree with {0} galaxies using their input positions...".format(cuts.sum())
    tree = cKDTree(np.vstack([X,Y]).T)

    for psfID in psfIDs:
      for g_id in xrange(len(gRange)):
        runID = gRange[g_id]+'_'+str(psfID)+'_'+randomKey
        print "Comparing ", runID
        #shearID, psfID, randomKey = runID.split('_')
        ARCHDIR = os.path.join('/disks/shear15/KiDS/ImSim/pipeline/archive/',randomKey,runID)
        TMPDIR = os.path.join('/disks/shear15/KiDS/ImSim/temp',randomKey,runID)

        prior_catname = 'prior'
        prior_pathname = os.path.join(ARCHDIR,prior_catname)
        prior_dat = np.loadtxt(prior_pathname)

        sex_arrs, lf_arrs = [ ], [ ]
        indices = [ ]
        for rot_id in xrange(4):
            sex_catname = 'sexrot0{0}.cat'.format(rot_id)
            lf_catname = '0{0}.output.rot.fits.asc.scheme2b_corr'.format(rot_id)
            if rot_id==0:
                    sex_catname = 'sex.cat'
            #        lf_catname = 'output.fits.asc.scheme2b_corr'
            sex_pathname = os.path.join(TMPDIR,sex_catname)
            lf_pathname = os.path.join(ARCHDIR,lf_catname)

            sex_params_filename = 'kidssims.param'
            sex_params_pathname = os.path.join('/disks/shear15/KiDS/ImSim/pipeline/backup/pipeline/config/',sex_params_filename)
            with open(sex_params_pathname,'r') as f:
                sex_fieldnames = f.readlines()
            ## Remove the empty lines
            n_emptylines = sex_fieldnames.count('\n')
            for ii in xrange(n_emptylines):
                sex_fieldnames.remove('\n')
            ## Strip of the newline character from the rest
            for ii in xrange(len(sex_fieldnames)):
                sex_fieldnames[ii] = sex_fieldnames[ii][:-1]

            lf_fieldnames = []
            with open(lf_pathname,'r') as f:
                for lineno in xrange(31):
                    line = f.readline()
                    words = line.split()
                    ## Omit the first line. It is not a column name
                    if lineno>0:
                        ## Because somebody thought giving a space in between is legible
                        lf_fieldname = ' '.join(words[2:])
                        lf_fieldnames.append(lf_fieldname)

            sex_arr = np.loadtxt(sex_pathname)
            lf_arr = np.loadtxt(lf_pathname)

            assert len(sex_arr)==len(lf_arr)
            assert sex_arr.shape[1]==len(sex_fieldnames)
            assert lf_arr.shape[1]==len(lf_fieldnames)

            d2d, idx = tree.query(np.array([ sex_arr[:,sex_fieldnames.index('X_IMAGE')], sex_arr[:,sex_fieldnames.index('Y_IMAGE')] ]).T)

            sex_arrs.append( sex_arr )
            lf_arrs.append( lf_arr )
            indices.append( idx )

        ## Assuming that the columns mean the same for all rotations, append them
        sex_dat = np.vstack(tuple(sex_arrs))
        lf_dat = np.vstack(tuple(lf_arrs))

        ## Make the QC directory, if it doesn't exist already
        QC_dirname = os.path.join(ARCHDIR,'QC')
        if not 'QC' in os.listdir(ARCHDIR):
            os.mkdir(QC_dirname)

        ## Make the overall distributions

        ## Magnitude plots
        fig, ax = plt.subplots()
        bins = np.arange(16,27,0.05)
        prior_mag_col_id = 2
        _n, _bins, _patches = ax.hist(prior_dat[:,prior_mag_col_id], bins=bins, histtype='step', color='k', label='Input magnitude')
        _n, _bins, _patches = ax.hist(sex_dat[:,sex_fieldnames.index('MAG_AUTO')], bins=bins, histtype='step', weights=0.25*np.ones(len(sex_dat)), color='r', label='Output magnitude')
        ax.set_yscale('log')
        _lgnd = ax.legend(loc='best')
        fig.suptitle('Magnitude distributions')
        fig_filename = 'magnitudes.png'
        fig_pathname = os.path.join(QC_dirname,fig_filename)
        fig.savefig(fig_pathname)

    #    ## SExtractor SNR plots
    #    fig, ax = plt.subplots()
    #    bins = np.logspace(-2,4,60)
    #    _n, _bins, _patches = ax.hist(input_data[cuts]['FLUX_AUTO_THELI']/input_data[cuts]['FLUXERR_AUTO_THELI'], bins=bins, histtype='step', color='k', label='SNR in data')
    #    _n, _bins, _patches = ax.hist(sex_dat[:,sex_fieldnames.index('FLUX_AUTO')]/sex_dat[:,sex_fieldnames.index('FLUXERR_AUTO')], bins=bins, histtype='step', weights=0.25*np.ones(len(sex_dat)), color='r', label='SNR in sims')
    #    ax.set_xscale('log')
    #    _lgnd = ax.legend(loc='best')
    #    fig.suptitle('SNR from SExtractor')
    #    fig_filename = 'snr.png'
    #    fig_pathname = os.path.join(QC_dirname,fig_filename)
    #    fig.savefig(fig_pathname)

    #    ## FWHM plots
    #    fig, ax = plt.subplots()
    #    bins = np.logspace(-2,2,20)
    #    _n, _bins, _patches = ax.hist(input_data[cuts]['FWHM_IMAGE_THELI'], bins=bins, histtype='step', color='k', label='FWHM_IMAGE in data')
    #    _n, _bins, _patches = ax.hist(sex_dat[:,sex_fieldnames.index('FWHM_IMAGE')], bins=bins, histtype='step', weights=0.25*np.ones(len(sex_dat)), color='r', label='FWHM_IMAGE in sims')
    #    ax.set_xscale('log')
    #    _lgnd = ax.legend(loc='best')
    #    fig.suptitle('FWHM_IMAGE from SExtractor')
    #    fig_filename = 'fwhm.png'
    #    fig_pathname = os.path.join(QC_dirname,fig_filename)
    #    fig.savefig(fig_pathname)

        ## Scalelength
        fig, ax = plt.subplots()
        bins = np.logspace(-2,2,20)
        _n, _bins, _patches = ax.hist(input_data[cuts]['bias_corrected_scalelength_pixels'], bins=bins, histtype='step', color='k', label='LF scalength in data')
        _n, _bins, _patches = ax.hist(lf_dat[:,lf_fieldnames.index('bias-corrected scalelength /pixels')], bins=bins, histtype='step', weights=0.25*np.ones(len(lf_dat)), color='r', label='LF scalelength in sims')
        ax.set_xscale('log')
        _lgnd = ax.legend(loc='best')
        fig.suptitle('Scalelength from LF')
        fig_filename = 'scalelength.png'
        fig_pathname = os.path.join(QC_dirname,fig_filename)
        fig.savefig(fig_pathname)

        ## model SNR
        fig, ax = plt.subplots()
        bins = np.logspace(-2,4,60)
        _n, _bins, _patches = ax.hist(input_data[cuts]['model_SNratio'], bins=bins, histtype='step', color='k', label='Model SNR in data')
        _n, _bins, _patches = ax.hist(lf_dat[:,lf_fieldnames.index('model SNratio')], bins=bins, histtype='step', weights=[0.25]*len(lf_dat), color='r', label='Model SNR in sims')
        ax.set_xscale('log')
        _lgnd = ax.legend(loc='best')
        fig.suptitle('Model SNR from LF')
        fig_filename = 'snr_model.png'
        fig_pathname = os.path.join(QC_dirname,fig_filename)
        fig.savefig(fig_pathname)

        ## pixel SNR
        fig, ax = plt.subplots()
        bins = np.logspace(-2,4,60)
        _n, _bins, _patches = ax.hist(input_data[cuts]['pixel_SNratio'], bins=bins, histtype='step', color='k', label='Pixel SNR in data')
        _n, _bins, _patches = ax.hist(lf_dat[:,lf_fieldnames.index('pixel SNratio')], bins=bins, histtype='step', weights=[0.25]*len(lf_dat), color='r', label='Pixel SNR in sims')
        ax.set_xscale('log')
        _lgnd = ax.legend(loc='best')
        fig.suptitle('Pixel SNR from LF')
        fig_filename = 'snr_pixel.png'
        fig_pathname = os.path.join(QC_dirname,fig_filename)
        fig.savefig(fig_pathname)

        plt.close('all')

def check_psf(runID, refID='FC17SB'):
    TMPDIR = '/disks/shear15/KiDS/ImSim/temp/'
    PSFDIRS = ['psf','galpsf']

    gRange = ['p400m000','m400m000','m000p400', 'm000m400',
              'm283m283','p283m283','m283p283', 'p283p283' ]

    for g_id in gRange:
        for psfDir in PSFDIRS:
            for psfID in xrange(5):
                refPsfFiles = os.listdir(os.path.join(TMPDIR, refID, g_id+'_'+str(psfID)+'_'+refID, psfDir))
                psfFiles = os.listdir(os.path.join(TMPDIR, runID, g_id+'_'+str(psfID)+'_'+runID, psfDir))

                for psfFile in psfFiles:
                    refpsf = pf.getdata(os.path.join(TMPDIR, refID, g_id+'_'+str(psfID)+'_'+refID, psfDir, psfFile))
                    psf = pf.getdata(os.path.join(TMPDIR, refID, g_id+'_'+str(psfID)+'_'+refID, psfDir, psfFile))

                    np.testing.assert_array_equal(refpsf, psf)

    print "psf and galpsf for {0} are consistent with those of {1}".format(runID, refID)

def check_prior_headers(runID, refID,psfIDs=[0]):
    ARCHDIR = '/disks/shear15/KiDS/ImSim/pipeline/archive/'
    gRange = ['p400m000','m400m000','m000p400', 'm000m400',
              'm283m283','p283m283','m283p283', 'p283p283' ]
    
    n_mismatches = 0
    for psfID in psfIDs:
        for g_id in gRange:
            refPriorPath = os.path.join(ARCHDIR, refID, g_id+'_'+str(psfID)+'_'+refID, 'prior')
            priorPath = os.path.join(ARCHDIR, runID, g_id+'_'+str(psfID)+'_'+runID, 'prior')

            if not os.path.exists(priorPath):
                print priorPath, "does not exist. Moving on ... "
                continue

            refPrior, prior = [ ], [ ]
            with open(refPriorPath,'r') as f:
                while True:
                    line = f.readline()
                    if line[0]=='#':
                        refPrior.append(line)
                    else:
                        break

            with open(priorPath,'r') as f:
                while True:
                    line = f.readline()
                    if line[0]=='#':
                        prior.append(line)
                    else:
                        break

            try:
                assert len(prior)>0
            except AssertionError:
                print priorPath, " has no header"

            try:
                assert refPrior[:-2] == prior[:-2]
            except AssertionError:
                n_mismatches += 1
                print priorPath, " does not have a matching header with ", refPriorPath

    if n_mismatches==0:
        print "'check_prior_headers' passed successfully"
    else:
        print "'check_prior_headers' found ", n_mismatches, " mismatch(es)."

def check_prior_data(runID, refID, psfIDs=[0]):
    ARCHDIR = '/disks/shear15/KiDS/ImSim/pipeline/archive/'
    gRange = ['p400m000','m400m000','m000p400', 'm000m400',
              'm283m283','p283m283','m283p283', 'p283p283' ]
    
    n_mismatches = 0
    for psfID in psfIDs:
        for g_id in gRange:
            refPriorPath = os.path.join(ARCHDIR, refID, g_id+'_'+str(psfID)+'_'+refID, 'prior')
            priorPath = os.path.join(ARCHDIR, runID, g_id+'_'+str(psfID)+'_'+runID, 'prior')

            if not os.path.exists(priorPath):
                print priorPath, "does not exist. Moving on ... "
                continue

            refPrior = np.loadtxt(refPriorPath, comments='#')
            prior = np.loadtxt(priorPath, comments='#')

            try:
                np.testing.assert_array_equal(refPrior, prior)
            except AssertionError:
                n_mismatches += 1
                print "Data of ", priorPath, " does not match with that of ", refPriorPath

    if n_mismatches==0:
        print "'check_prior_data' passed successfully"
    else:
        print "'check_prior_data' found ", n_mismatches, " mismatch(es)."
    
def check_exposure(runID, refID, expoIDs=[0], psfIDs=[0,1,2]):
    TMPDIR = '/disks/shear15/KiDS/ImSim/temp'
    gRange = ['p400m000','m400m000','m000p400', 'm000m400',
              'm283m283','p283m283','m283p283', 'p283p283' ]
    
    for psfID in psfIDs:
      for expoID in expoIDs:
        for g_id in gRange:
          for rot_id in xrange(4):

            chipdir = 'chip/' if rot_id==0 else 'chiprot0{0}/'.format(rot_id)

            refImagePath = os.path.join(TMPDIR, refID, g_id+'_'+str(psfID)+'_'+refID, chipdir)
            imagePath = os.path.join(TMPDIR, runID, g_id+'_'+str(psfID)+'_'+runID, chipdir)

            if not os.path.exists(imagePath):
                print imagePath, " does not exist. Moving on ..."
                break

            refChipnames = os.listdir(refImagePath)
            chipnames = os.listdir(imagePath)
            chipname = 'exp{0}chip_10OFCS.sub.fits'.format(expoID)

            t1 = time.time()
            refImage = pf.getdata(refImagePath+chipname)
            image = pf.getdata(imagePath+chipname)

            try:
                np.testing.assert_array_equal(refImage, image)
            except AssertionError:
                print imagePath+chipname, " does not match with ", refImagePath+chipname

            t2 = time.time()
            print imagePath, " verified in only ", t2-t1, " seconds."

def check_prior_priornew(runID, psfIDs=[0,1,2,3,4]):
    ARCHDIR = '/disks/shear15/KiDS/ImSim/pipeline/archive/'
    gRange = ['p400m000','m400m000','m000p400', 'm000m400',
              'm283m283','p283m283','m283p283', 'p283p283' ]

    for psfID in psfIDs:
        for g_id in gRange:
            prior_filename = os.path.join(ARCHDIR,runID,g_id+'_'+str(psfID)+'_'+runID,'prior_test')
            priornew_filename = prior_filename.replace('prior_test','prior_new')

            prior = np.loadtxt(prior_filename)
            priornew = np.loadtxt(priornew_filename)

            try:
                np.testing.assert_array_equal(prior[:,:7],priornew[:,:7])
                np.testing.assert_array_equal(prior[:,8:],priornew[:,8:])
                np.testing.assert_array_equal(np.abs(prior[:,7]),np.abs(priornew[:,7]))
            except AssertionError:
                print "One of the non-N columns in ", prior_filename, " don't match."

            try:
                np.testing.assert_array_equal(prior[:,7],priornew[:,7])
                print "Sign doesn't appear to have changed."
            except AssertionError:
                print "Sign appears to have changed."
            
if __name__=='__main__':
    runID = sys.argv[1]
    if len(sys.argv)>2:
        refID = sys.argv[2]
    #check_consistency(runID)
    #check_psf(runID)
    #check_prior_headers(runID, refID='FC17SB',psfIDs=[0,1,2,3,4])
    #check_prior_data(runID, refID='TSTGr0615truepos',psfIDs=[0,1,2,3,4])
    check_exposure(runID, refID=refID,expoIDs=[4],psfIDs=[0,1,2,3,4])

    #check_exposure('TST77LF0727','TSTGr0727',expoIDs=[0,1,2,3,4],psfIDs=[0,1,2,3,4])
    #check_exposure('TST77LF0729','TSTGr0729truepos',expoIDs=[0,1,2,3,4],psfIDs=[0,1,2,3,4])
    #check_prior_priornew(runID,range(1,5))
