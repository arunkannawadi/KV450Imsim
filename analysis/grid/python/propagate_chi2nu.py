from astropy.io import fits
from astropy.table import Table
import numpy as np
import os, sys
import time
#import pdb; pdb.set_trace()

def propagate_field(path_to_mastercat, fieldnames=['CHI2NU_HI']):

    t_start = time.time()
    ## Read the master catalogue
    masterDat = fits.open(path_to_mastercat)[1].data

    existing_fieldnames = [ ]
    for fieldname in fieldnames:
        if fieldname in masterDat.names:
            print fieldname, " already exists in the master catalogue. Skipping ..."
            existing_fieldnames.append(fieldname)
    for fn in existing_fieldnames:
        fieldnames.remove(fn)

    if len(fieldnames)==0:
        print "All the fieldnames already exist."
        return 0

    ## Read the input file
    path_to_inputCat = '/disks/shear14/KiDS_simulations/Cosmos/KIDS_HST_cat/KiDS_Griffith_iMS1_handpicked_stars.cat'
    inputCat = fits.open(path_to_inputCat)
    inputDat = inputCat[1].data
    objno = inputDat['OBJNO']
    objno_argsort = np.argsort(objno)

    print "Input catalogue loaded"

    t1 = time.time()
    print "Time to read the two catalogues: ", np.round(t1-t_start, 2), " seconds."
    sys.stdout.flush()

    ## Need to remove the objects which have no match in the prior
    if 'prior_matched' in masterDat.names:
        prior_matched = masterDat['prior_matched']
    else:
        prior_matched = np.ones(len(masterDat),dtype=int)

    catID = masterDat['Cat_ID']

    new_array = np.empty((len(masterDat),len(fieldnames)))
    if False:
        ## Brute force: Bad
        for lineno in xrange(len(masterDat)):
            t1 = time.time()
            if (prior_matched[lineno]==1)& (~(catID[lineno]==-1)):
                try:
                    row_id = np.where(objno==catID[lineno])[0][0]
                except:
                    pass
                for field_id, fieldname in enumerate(fieldnames):
                    new_array[lineno,field_id] = inputDat[fieldname][row_id]
            else:
                for field_id, fieldname in enumerate(fieldnames):
                    new_array[lineno,field_id] = -99
            t2 = time.time()
            print "Time for lineno: ", np.round(t2-t1,2), " seconds."
            sys.stdout.flush()
    ## More optimal

    if True:
        cuts = (prior_matched==1)&(catID>-1)
        matching_indices = np.searchsorted(objno[objno_argsort],catID[cuts])
        for field_id, fieldname in enumerate(fieldnames):
            new_array[cuts,field_id] = inputDat[fieldname][objno_argsort[matching_indices]]
            new_array[~cuts,field_id] = -99*np.ones( (~cuts).sum() )

    formats = 'float,'*len(fieldnames)
    new_recarray = np.rec.array(new_array, formats='float',names=fieldnames)

    new_tbhdu = fits.BinTableHDU(data=new_recarray)

    combined_columns = masterDat.columns + new_tbhdu.columns

    combined_hdu = fits.BinTableHDU.from_columns(combined_columns)

#    output_filename = 'testing.fits'
    output_filename = path_to_mastercat
    combined_hdu.writeto(output_filename, overwrite=True)

if __name__=='__main__':
    propagate_field(sys.argv[1])
            

    





