import numpy as np
from astropy.io import fits
import os, sys
sys.path.append('/disks/shear14/arunkannawadi/catalogue_matching/')
from match_catalogues import match_catalogues, make_joint_catalogue 
import pdb; pdb.set_trace()

def append_Griffith_catalogue():
    kids_pathname = '/disks/shear14/KiDS_simulations/Cosmos/COSMOS_photoz.cat'
    griffith_pathname = '/disks/shear14/arunkannawadi/catalogue_matching/Griffith_COSMOS_catalogue.fits'

    kids_cat = fits.open(kids_pathname)
    griffith_cat = fits.open(griffith_pathname)

    kids_dat = kids_cat[1].data
    griffith_dat = griffith_cat[1].data

    mask_all = kids_dat.MASK
    kids_cuts = ~np.array(mask_all&0xfc3c,dtype=bool)

    ## Query points
    RA1 = kids_dat['ALPHA_J2000']
    DEC1 = kids_dat['DELTA_J2000']
    MAG1 = kids_dat['MAG_AUTO']

    ## kd-tree nodes
    RA2 = griffith_dat['RA']
    DEC2 = griffith_dat['DEC']
    MAG2 = griffith_dat['MAGR']

    matching_indices = match_catalogues(use_kdtree=True, inverse_MAG_AUTO_SCALE=1, greedy_level=4, RA1=RA1, RA2=RA2, DEC1=DEC1, DEC2=DEC2,
                                        MAG_AUTO1=MAG1, MAG_AUTO2=MAG2, mask1=kids_cuts)

    output_filename = '/disks/shear15/KiDS/ImSim/pipeline/data/KiDS_Griffith_iMS1_testing.cat'

    make_joint_catalogue(cat1=kids_cat,cat2=griffith_cat,matching_indices=matching_indices,suffix1='_KIDS',suffix2='_HST',mask1=kids_cuts,output_filename=output_filename)

def append_GalSim_catalogue():
    galsim_cat1_pathname = '/disks/shear15/KiDS/ImSim/pipeline/data/COSMOS_25.2_training_sample/real_galaxy_catalog_25.2.fits'
    galsim_cat2_pathname = '/disks/shear15/KiDS/ImSim/pipeline/data/COSMOS_25.2_training_sample/real_galaxy_catalog_25.2_fits.fits'

    galsim_cat1 = fits.open(galsim_cat1_pathname)
    galsim_cat2 = fits.open(galsim_cat2_pathname)

    galsim_dat1 = galsim_cat1[1].data
    galsim_dat2 = galsim_cat2[1].data

    np.testing.assert_array_equal(galsim_dat1['IDENT'], galsim_dat2['IDENT'])

    ## Skip joining the two catalogues for now

    RA_galsim = galsim_dat1['RA']
    DEC_galsim = galsim_dat1['DEC']
    MAG_galsim = galsim_dat1['MAG']

    ## Load the KiDS+Griffith catalogue
    #main_pathname = '/disks/shear14/KiDS_simulations/Cosmos/KIDS_HST_cat/KiDS_Griffith_iMS1_handpicked_stars.cat'
    main_pathname = '/disks/shear15/KiDS/ImSim/pipeline/data/KiDS_Griffith_iMS1_testing.cat'
    main_cat = fits.open(main_pathname)
    main_dat = main_cat[1].data

    main_cuts = main_dat['rank']>=0
    mask1 = main_cuts

    RA_main = main_dat['RA']
    DEC_main = main_dat['DEC']
    MAG_main = main_dat['MAGR']

    matching_indices = match_catalogues(use_kdtree=True, inverse_MAG_AUTO_SCALE=0, greedy_level=4, RA1=RA_galsim, RA2=RA_main, DEC1=DEC_galsim, DEC2=DEC_main,
                                        MAG_AUTO1=MAG_galsim, MAG_AUTO2=MAG_main, mask2=main_cuts)

    output_filename = '/disks/shear15/KiDS/ImSim/pipeline/data/KiDS_Griffith_iMS1_handpicked_stars_GalSim_iMS0.cat'
    make_joint_catalogue(cat1=galsim_cat1,cat2=main_cat,matching_indices=matching_indices,suffix1='_galsim',suffix2='_Griffith',mask2=main_cuts,output_filename=output_filename)

if __name__=='__main__':
    append_Griffith_catalogue()
    #append_GalSim_catalogue()
