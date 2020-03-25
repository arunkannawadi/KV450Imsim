import numpy as np
from astropy.io import fits
import pyfits as pf
import os,sys
from match import match
from ConfigParser import SafeConfigParser
import optparse

def assignment(objno, catID, ZB9, sort=True):
    superset = set(objno)
    subset = set(catID)
    try:
        assert subset.issubset(superset)
    except AssertionError:
        print "The prior catalogue has OBJNO that does not belong in the HST catalogues. You should look into it"
    if set(catID).issubset(superset):
        subset_cuts = True
    else:
        difference_set = set(catID).difference(superset)
        print "Size of difference set = ", len(difference_set)
        subset_cuts = np.ones_like(catID,dtype=bool)
        for ii in difference_set:
            subset_cuts[catID==ii] = False

    objno_argsort = np.argsort(objno)
    matching_indices = np.searchsorted(objno, catID[subset_cuts][0], sorter=objno_argsort)
    return ZB9[subset_cuts][objno_argsort[matching_indices]], subset_cuts
    
runID = sys.argv[1]
psfID = sys.argv[2]
configPath = sys.argv[3]

# Load Configuration Settings.
print '  Loading Config from %s' % configPath
parser = SafeConfigParser()
parser.read(configPath)

ARCHDIR = parser.get('assign_ZB', 'ARCHDIR')
prior_filename = 'prior'
#shearIDs = [sID for sID in os.listdir(os.path.join(ARCHDIR,runID.split('/')[0].replace('globalRecal',''))) if sID.count('_')==2 and '.fits' not in sID]
#print "No. of sub directories = ", len(shearIDs)

all_galaxies_catname = '/disks/shear14/KiDS_simulations/Cosmos/KIDS_HST_cat/KiDS_Griffith_iMS1_handpicked_stars_GOLDonly.cat'
all_galaxies_cat = fits.open(all_galaxies_catname)
all_galaxies_dat = all_galaxies_cat[1].data
all_cuts =  (all_galaxies_dat['rank']==1)&(all_galaxies_dat['distance2d']<1.)&(all_galaxies_dat['OBJNO']>0) ## Z_B can be assigned only if its there in KiDS and in Griffith
all_OBJNO = all_galaxies_dat['OBJNO'][all_cuts]
all_SeqNr = all_galaxies_dat['SeqNr'][all_cuts]

cosmos_photoz_catname = parser.get('assign_ZB','COSMOSCatalogue')
cosmos_photoz_cat = fits.open(cosmos_photoz_catname)
cosmos_photoz_dat = cosmos_photoz_cat[1].data

cosmos_idx, all_idx = match(R=cosmos_photoz_dat['SeqNr'], Q=all_SeqNr)

## Test the consistency between 4-band photo-z without having to provide all_idx as a sanity check
if len(cosmos_photoz_dat)>=143297:
    np.testing.assert_array_equal(all_SeqNr[all_idx], cosmos_photoz_dat['SeqNr'][cosmos_idx])
    np.testing.assert_array_equal(all_galaxies_dat[all_cuts]['Z_B'], cosmos_photoz_dat['Z_B_ugri'][cosmos_idx])

ZB4 = cosmos_photoz_dat['Z_B_ugri'][cosmos_idx]
ZB9 = cosmos_photoz_dat['Z_B'][cosmos_idx]
TB4 = cosmos_photoz_dat['T_B_ugri'][cosmos_idx]
TB9 = cosmos_photoz_dat['T_B'][cosmos_idx]
SEQNR = cosmos_photoz_dat['SeqNr'][cosmos_idx]

SOM_flags = [ 'Flag_SOM_Fid_NONE',
 'Flag_SOM_multispec3_NONE',
 'Flag_SOM_noDEEP2_NONE',
 'Flag_SOM_noVVDS_NONE',
 'Flag_SOM_nozCOSMOS_NONE',
 'Flag_SOM_speczquality4_NONE',
 'GroupFactor_SOM_Fid_NONE',
 'GroupFactor_SOM_multispec3_NONE',
 'GroupFactor_SOM_noDEEP2_NONE',
 'GroupFactor_SOM_noVVDS_NONE',
 'GroupFactor_SOM_nozCOSMOS_NONE',
 'GroupFactor_SOM_speczquality4_NONE',
 'RedshiftGoldClass_NONE',
 'SurveyGoldClass_NONE',
 'SurveyGoldFlag_NONE']


## Set the photo-z to -2 if the GAaP flag is set
#ZB4[cosmos_photoz_dat['GAAP_Flag_ugriZYJHKs'][cosmos_idx]!=0] = -9
ZB9[cosmos_photoz_dat['GAAP_Flag_ugriZYJHKs'][cosmos_idx]!=0] = -9
TB4[cosmos_photoz_dat['GAAP_Flag_ugriZYJHKs'][cosmos_idx]!=0] = -9
TB9[cosmos_photoz_dat['GAAP_Flag_ugriZYJHKs'][cosmos_idx]!=0] = -9

## We have to end with [0] for dimension matching since we won't be providing all_idx anywhere

#for sID in #shearIDs:
#    print "On ", sID
#    prior_pathname = os.path.join(ARCHDIR,runID,sID,prior_filename)
#    prior = np.loadtxt(prior_pathname)
#
#    assert prior.shape[1]==11
#    print "Using last column and last but one for OBJNO and ZB9"
#    objno = prior[:,-1]
#    prior[:,-2] = assignment(objno=all_OBJNO, catID=objno, ZB9=all_ZB9)
#    np.savetxt(os.path.join(ARCHDIR,runID,sID,prior_pathname+'_new'),prior)

if int(psfID)<13:
    masterCat_filename = os.path.join(ARCHDIR,runID, 'MasterCat_'+runID+'_set_'+psfID+'.fits')
else:
    masterCat_filename = os.path.join(ARCHDIR,runID, 'MasterCat_'+runID+'_all_'+psfID+'_PSF.fits')

masterCat = fits.open(masterCat_filename)
masterCat_data = masterCat[1].data

masterCat_cuts = (masterCat_data['Cat_ID']>-1)
objno = masterCat_data['Cat_ID'][masterCat_cuts]
all_OBJNO = all_OBJNO.astype(int)
objno = objno.astype(int)
#subset_cuts = (objno>=all_OBJNO.min())&(objno<=all_OBJNO.max())
matching_indices, subset_cuts = match(R=all_OBJNO, Q=objno)
try:
   np.testing.assert_array_equal(all_OBJNO[matching_indices],objno[subset_cuts])
except AssertionError:
    print "The function 'assignment' doesn't seem to be doing the right thing. Refraining from assigning ZB9"

## Somehow, assignment after double 'view' doesn't work. So doing it one by one
sub_array = -5.0*np.ones(masterCat_cuts.sum())
sub_array[subset_cuts] = ZB9[matching_indices]
new_array = -5.0*np.ones(len(masterCat_data))
new_array[masterCat_cuts] = sub_array

sub_array_Z4 = -5.0*np.ones(masterCat_cuts.sum())
sub_array_Z4[subset_cuts] = ZB4[matching_indices]
new_array_Z4 = -5.0*np.ones(len(masterCat_data))
new_array_Z4[masterCat_cuts] = sub_array_Z4

sub_array_TB = -5.0*np.ones(masterCat_cuts.sum())
sub_array_TB[subset_cuts] = TB9[matching_indices]
new_array_TB = -5.0*np.ones(len(masterCat_data))
new_array_TB[masterCat_cuts] = sub_array_TB

sub_array_SEQNR = -5*np.ones(masterCat_cuts.sum(),dtype=int)
sub_array_SEQNR[subset_cuts] = SEQNR[matching_indices]
new_array_SEQNR = -5*np.ones(len(masterCat_data),dtype=int)
new_array_SEQNR[masterCat_cuts] = sub_array_SEQNR

new_array_dict = dict.fromkeys(SOM_flags)
for SOM_colname in SOM_flags:
    sub_array = -1*np.ones(masterCat_cuts.sum(),dtype=int)
    SOM_col = cosmos_photoz_dat[SOM_colname][cosmos_idx]
    sub_array[subset_cuts] = SOM_col[matching_indices]
    new_array_dict[SOM_colname] = -1*np.ones(len(masterCat_data),dtype=int)
    new_array_dict[SOM_colname][masterCat_cuts] = sub_array

try:
    np.testing.assert_array_equal(masterCat_data['ZB4_in'][masterCat_cuts][subset_cuts], ZB4[matching_indices])
except AssertionError as ae:
    print "The 4-band ZBs do not seem to be in agreement."

masterCat_data['ZB4_in'] = new_array_Z4
new_hdu = fits.BinTableHDU(data=masterCat_data)

if not 'ZB9_in' in masterCat_data.names:
    new_col = fits.Column(name='ZB9_in',format='D',array=new_array)
    new_coldefs = masterCat_data.columns + new_col
    new_hdu = fits.BinTableHDU.from_columns(new_coldefs)
else:
    #masterCat_data['ZB9_in'][masterCat_cuts][subset_cuts] = ZB9[matching_indices] ## equiv. to all_galaxies_dat['ZB9'][all_cuts][matching_indices]
    masterCat_data['ZB9_in'] = new_array
    new_hdu = fits.BinTableHDU(data=masterCat_data)

if not 'TB9_in' in masterCat_data.names:
    new_col = fits.Column(name='TB9_in',format='D',array=new_array_TB)
    new_coldefs = new_hdu.data.columns + new_col
    new_hdu = fits.BinTableHDU.from_columns(new_coldefs)
else:
    #masterCat_data['ZB9_in'][masterCat_cuts][subset_cuts] = ZB9[matching_indices] ## equiv. to all_galaxies_dat['ZB9'][all_cuts][matching_indices]
    new_hdu.data['TB9_in'] = new_array_TB
    new_hdu = fits.BinTableHDU(data=new_hdu.data)

if not 'SeqNr' in masterCat_data.names:
    new_col = fits.Column(name='SeqNr',format='J',array=new_array_SEQNR)
    new_coldefs = new_hdu.data.columns + new_col
    new_hdu = fits.BinTableHDU.from_columns(new_coldefs)

for SOM_colname in SOM_flags:
    if not SOM_colname in masterCat_data.names:
        new_col = fits.Column(name=SOM_colname, format='J',array=new_array_dict[SOM_colname])
        new_coldefs = new_hdu.data.columns + new_col
        new_hdu = fits.BinTableHDU.from_columns(new_coldefs)

## HACK ALERT: Insert high-quality 30-band photo-z / spec-z for 4-band Z_B
if False:
    sub_array_Z = -5.0*np.ones(masterCat_cuts.sum())
    sub_array_Z[subset_cuts] = all_galaxies_dat['Z'][all_cuts][matching_indices]
    new_array_Z = -5.0*np.ones(len(masterCat_data))
    new_array_Z[masterCat_cuts] = sub_array_Zi

    new_hdu.data['ZB4_in'] = new_array_Z
    new_hdu = fits.BinTableHDU(data=new_hdu.data)
if False:
    griffith_laigle_catname = parser.get('assign_ZB','Griffith_Laigle')
    griffith_laigle_cat = fits.open(griffith_laigle_catname)
    griffith_laigle_dat = griffith_laigle_cat[1].data
    idx, cuts = match(R=griffith_laigle_dat['OBJNO'], Q=objno)
    sub_array_Z = -5.0*np.ones(masterCat_cuts.sum())
    sub_array_Z[cuts] = griffith_laigle_dat['PHOTOZ_Laigle'][idx]
    new_array_Z = -5.0*np.ones(len(masterCat_data))
    new_array_Z[masterCat_cuts] = sub_array_Z
    new_hdu.data['ZB4_in'] = new_array_Z
    new_hdu = fits.BinTableHDU(data=new_hdu.data)

## Uncomment this to test the assignment here, in an ineffecient but safe way
if False:
    print "Testing the assignment now. This can take a while ... "
    for row_id in xrange((all_cuts.sum())):
        ref_seqnr, ref_objno = all_galaxies_dat['SeqNr'][all_cuts][row_id], all_galaxies_dat['OBJNO'][all_cuts][row_id]
        ref_zb9 = cosmos_photoz_dat['Z_B'][np.where(cosmos_photoz_dat['SeqNr']==ref_seqnr)[0][0]]
        try:
            test_zb9 = masterCat_data['ZB9_in'][masterCat_data['Cat_ID']==ref_objno]
            assert (test_zb9==ref_zb9).all()
        except:
            if (test_zb9>=0).any():
                print "Assertion failure in ", row_id, ref_seqnr, ref_objno, ref_zb9, all_galaxies_dat['rank'][all_cuts][row_id], test_zb9
            pass
    print "Assignment test complete!"
output_pathname = masterCat_filename
new_hdu.writeto(output_pathname, overwrite=True )
