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
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
#debug = bool(int(sys.argv[2]))
#if debug:
#    import pdb; pdb.set_trace()

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
weight = input_data['weight'][cuts]

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

runID = 'm400m000_0_TSTGr0729truepos'
shearID, psfID, randomKey = runID.split('_')
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
#        if rot_id==0:
#                sex_catname = 'sex.cat'
#                lf_catname = 'output.fits.asc.scheme2b_corr'
    sex_pathname = os.path.join(ARCHDIR,sex_catname)
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
if not 'QC2' in os.listdir(ARCHDIR):
    os.mkdir(os.path.join(ARCHDIR,'QC2'))

lf_weight = lf_dat[:,lf_fieldnames.index('corrected weight')]
print lf_weight.shape, prior_dat.shape, weight.shape

len(weight[idx]), len(lf_dat)

## Make the overall weight distribution

weight1 = weight[idx]
weight2 = lf_dat[:,lf_fieldnames.index('corrected weight')]

fig, ax = plt.subplots()
bins = np.arange(0,16,0.1)
_n, _bins, _patches = ax.hist(weight1, weights=4*np.ones_like(weight1), bins=bins, histtype='step', color='k', label='True weights')
_n, _bins, _patches = ax.hist(weight2, bins=bins, histtype='step', color='r', label='Sim weights')
ax.set_yscale('log')

## Make the overall distributions

## Magnitude plots
fig, ax = plt.subplots()
bins = np.arange(16,27,0.05)
prior_mag_col_id = 2
_n, _bins, _patches = ax.hist(prior_dat[:,prior_mag_col_id], bins=bins, weights=4*np.ones(len(prior_dat)), histtype='step', color='k', label='Input magnitude')
_n, _bins, _patches = ax.hist(sex_dat[:,sex_fieldnames.index('MAG_AUTO')], weights=np.ones_like(weight2), bins=bins, histtype='step', color='r', label='Output magnitude')
ax.set_yscale('log')
_lgnd = ax.legend(loc=3)
fig.suptitle('Magnitude distributions')
fig_filename = 'input_magnitudes.png'
fig_pathname = os.path.join(ARCHDIR,'QC2',fig_filename)
fig.savefig(fig_pathname)

weighted = True
if weighted:
    weights1 = 4.*weight1
    weights2 = weight2
else:
    weights1 = 4.*np.ones_like(weight1)
    weights2 = np.ones_like(weight2)
    
fig, ax = plt.subplots()
bins = np.arange(19,25,0.05)
_n, _bins, _patches = ax.hist(input_data[cuts]['MAG_AUTO_THELI'][idx], weights=weights1, bins=bins, histtype='step', color='k', label='Input magnitude')
_n, _bins, _patches = ax.hist(sex_dat[:,sex_fieldnames.index('MAG_AUTO')], weights=weights2, bins=bins, histtype='step', color='r', label='Output magnitude')
ax.set_yscale('log')
ax.set_ylim(1e3,1e5)
_lgnd = ax.legend(loc=8)
fig.suptitle('Magnitude distributions')
fig_filename = 'magnitudes.png'
fig_pathname = os.path.join(ARCHDIR,'QC2',fig_filename)
fig.savefig(fig_pathname)


## SExtractor SNR plots
weighted=False
if weighted:
    weights1 = 4*weight1
    weights2 = weight2
else:
    weights1 = 4.*np.ones_like(weight1)
    weights2 = np.ones_like(weight2)
    
fig, ax = plt.subplots()
bins = np.logspace(0,4,60)
sex_snr  = input_data[cuts]['FLUX_AUTO_THELI'][idx]/input_data[cuts]['FLUXERR_AUTO_THELI'][idx]
plot_cuts = ~(np.isnan(sex_snr)|np.isinf(sex_snr))
                  
_n, _bins, _patches = ax.hist(sex_snr[plot_cuts], weights=weights1[plot_cuts], bins=bins, histtype='step', color='k', label='SNR in data')
_n, _bins, _patches = ax.hist(sex_dat[:,sex_fieldnames.index('FLUX_AUTO')]/sex_dat[:,sex_fieldnames.index('FLUXERR_AUTO')], weights=weights2, bins=bins, histtype='step', color='r', label='SNR in sims')
ax.set_xscale('log')
_lgnd = ax.legend(loc='best')
fig.suptitle('SNR from SExtractor')
fig_filename = 'snr.png'
fig_pathname = os.path.join(ARCHDIR,'QC2',fig_filename)
fig.savefig(fig_pathname)

## Scalelength
weighted = True
if weighted:
    weights1 = 4*weight1
    weights2 = weight2
else:
    weights1 = 4*np.ones_like(weight1)
    weights2 = np.ones_like(weight2)
fig, ax = plt.subplots()
bins = np.logspace(-1,2,60)
_n, _bins, _patches = ax.hist(input_data[cuts]['bias_corrected_scalelength_pixels'][idx], weights=weights1, bins=bins, histtype='step', color='k', label='LF scalength in data')
_n, _bins, _patches = ax.hist(lf_dat[:,lf_fieldnames.index('bias-corrected scalelength /pixels')], weights=weights2, bins=bins, histtype='step', color='r', label='LF scalelength in sims')
ax.set_xscale('log')
_lgnd = ax.legend(loc='best')
fig.suptitle('Scalelength from LF')
fig_filename = 'scalelength.png'
fig_pathname = os.path.join(ARCHDIR,'QC2',fig_filename)
fig.savefig(fig_pathname)

## Ellipticity
weighted = True
if weighted:
    weights1 = 4*weight1
    weights2 = weight2
else:
    weights1 = 4*np.ones_like(weight1)
    weights2 = np.ones_like(weight2)
fig, ax = plt.subplots()
bins = np.logspace(-1,2,60)
_n, _bins, _patches = ax.hist(np.sqrt(input_data[cuts]['e1']**2+input_data[cuts]['e2']**2)[idx], weights=weights1, bins=bins, histtype='step', color='k', label='LF ellip in data')
_n, _bins, _patches = ax.hist(np.sqrt(lf_dat[:,lf_fieldnames.index('e1')]**2+lf_dat[:,lf_fieldnames.index('e2')]**2), weights=weights2, bins=bins, histtype='step', color='r', label='LF ellip in sims')
_lgnd = ax.legend(loc='best')
fig.suptitle('Ellipticity from LF')
fig_filename = 'ellipticity.png'
fig_pathname = os.path.join(ARCHDIR,'QC2',fig_filename)
fig.savefig(fig_pathname)

## model SNR
fig, ax = plt.subplots()
bins = np.logspace(0,4,60)
_n, _bins, _patches = ax.hist(input_data[cuts]['model_SNratio'][idx], weights=4*weight[idx], bins=bins, histtype='step', color='k', label='Model SNR in data')
_n, _bins, _patches = ax.hist(lf_dat[:,lf_fieldnames.index('model SNratio')], weights=weight2, bins=bins, histtype='step', color='r', label='Model SNR in sims')
ax.set_xscale('log')
_lgnd = ax.legend(loc=1)
fig.suptitle('Model SNR from LF')
fig_filename = 'snr_model.png'
fig_pathname = os.path.join(ARCHDIR,'QC2',fig_filename)
fig.savefig(fig_pathname)

## pixel SNR
weighted = True
if weighted:
    weights1 = 4*weight1
    weights2 = weight2
else:
    weights1 = 4*np.ones_like(weight1)
    weights2 = np.ones_like(weight2)
    
fig, ax = plt.subplots()
bins = np.logspace(0,3,60)
_n, _bins, _patches = ax.hist(input_data[cuts]['pixel_SNratio'][idx], weights=weights1, bins=bins, histtype='step', color='k', label='Pixel SNR in data')
_n, _bins, _patches = ax.hist(lf_dat[:,lf_fieldnames.index('pixel SNratio')], weights=weights2, bins=bins, histtype='step', color='r', label='Pixel SNR in sims')
ax.set_xscale('log')
_lgnd = ax.legend(loc='best')
fig.suptitle('Pixel SNR from LF')
fig_filename = 'snr_pixel.png'
fig_pathname = os.path.join(ARCHDIR,'QC2',fig_filename)
fig.savefig(fig_pathname)

## FWHM plots
weighted = True
if weighted:
    weights1 = 4*weight1
    weights2 = weight2
else:
    weights1 = 4*np.ones_like(weight1)
    weights2 = np.ones_like(weight2)
    
fig, ax = plt.subplots()
bins = np.logspace(0,2,60)
_n, _bins, _patches = ax.hist(input_data[cuts]['FWHM_IMAGE_THELI'][idx], weights=weights1, bins=bins, histtype='step', color='k', label='Data')
_n, _bins, _patches = ax.hist(sex_dat[:,sex_fieldnames.index('FWHM_IMAGE')], weights=weights2, bins=bins, histtype='step', color='r', label='Sims')
ax.set_xscale('log')
_lgnd = ax.legend(loc=2)
fig.suptitle('FWHM_IMAGE from SExtractor')
fig_filename = 'fwhm.png'
fig_pathname = os.path.join(ARCHDIR,'QC2',fig_filename)
fig.savefig(fig_pathname)

