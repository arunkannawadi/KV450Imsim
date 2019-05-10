#import pdb; pdb.set_trace()
import os, sys
import numpy as np
import scipy, scipy.optimize
from astropy.table import Table
from astropy.io import fits
from collections import Counter
import time

def bias_function(x,m,c):
    return (1.+m)*x+c

def rearrange(master_cat_pathname):
    rot_ids = np.linspace(0,3,1).astype(int)
    g1Range = np.array([-0.04,0.00,0.04,0.00,-0.0283,+0.0283,+0.0283,-0.0283])
    g2Range = np.array([0.00,0.04,0.00,-0.04,+0.0283,+0.0283,-0.0283,-0.0283])
    gRange = zip(g1Range,g2Range)

    t_pre1 = time.time()
    ### Table.read takes time to load and stuff, but allows fast access of fields
    t_load1 = time.time()
    masterDat = Table.read(master_cat_pathname) ## slow loading, fast access
#    print "Choosing only HI snr galaxies"
#    masterDat = masterDat[(masterDat['snr_model']>100)&(masterDat['snr_model']<150)]

#    masterCat = fits.open(master_cat_pathname)
#    masterDat = masterCat[1].data ## fast loading, slow access
    t_load2 = time.time()
    print "Time to load the catalogue data: ", np.round(t_load2-t_load1,2), " seconds."

    ## Varying output quantities
    LFweight = masterDat['LFweight']
    oldLFweight = masterDat['oldLFweight']
    e1 = masterDat['e1']
    e2 = masterDat['e2']
    mag_out = masterDat['mag_out']
    snr_model = masterDat['snr_model']
    size_out = masterDat['size_out']
    fitclass= masterDat['fitclass']
    nd = masterDat['nd']
    prior_matched = masterDat['prior_matched']
    try:
        chi2nu = masterDat['CHI2NU_HI']
    except:
        pass

    ## Constant input quantities
    N_in = masterDat['N_in']
    mag_in = masterDat['mag_in']
    size_in = masterDat['size_in']
    ZB5_in = masterDat['ZB4_in']
    rotation = masterDat['rotation']
    e1_in = masterDat['e1_in']
    e2_in = masterDat['e2_in']
    g1 = masterDat['g1']
    g2 = masterDat['g2']

    ## e1_in, e2_in are unsheared, but rotated intrinsic ellipticities
    ## Obtain the sheared ellipticities
    e_complex_in = e1_in + 1j*e2_in
    g_complex = g1 + 1j*g2
    sheared_e_in = (e_complex_in+g_complex)/(1.+g_complex.conj()*e_complex_in)
    sheared_e1_in = sheared_e_in.real
    sheared_e2_in = sheared_e_in.imag

    catID = masterDat['Cat_ID']
    if 'FC17' in master_cat_pathname:
        ## For FC simulations, the prior catalogue changes with shear. So assign them different IDs
        g_index = np.zeros_like(catID)
        for g_idx in xrange(8):
            g_index[ (g1==g1Range[g_idx])&(g2==g2Range[g_idx]) ] = g_idx
        catID += g_index * 1e6

    counter = Counter(catID)
    row_indices = np.array(counter.keys())
    counts = np.array(counter.values())

    ## Kick out -1 and 20170193 (dummy)
    row_indices = row_indices[counts<100]
    counts = counts[counts<100]
    ## Note that the order of the two lines is important

    nRows = len(row_indices)
    nCols = 32
#    nRows = 10
    print "Number of galaxies: ", nRows

    ## Initialize the empty arrays
    e1_in_array = np.empty(nRows)
    e2_in_array = np.empty(nRows)
    mag_in_array = np.empty(nRows)
    size_in_array = np.empty(nRows)
    N_in_array = np.empty(nRows)
    ZB5_in_array = np.empty(nRows)

    mag_out_avg_array = np.empty(nRows)
    size_out_avg_array = np.empty(nRows)
    e1_avg_array = np.empty(nRows)
    e2_avg_array = np.empty(nRows)
    sheared_e1_in_avg_array = np.empty(nRows)
    sheared_e2_in_avg_array = np.empty(nRows)
    snr_avg_array = np.empty(nRows)
    rotation_avg_array = np.empty(nRows)
    LFweight_avg_array = np.empty(nRows)
    LFweight_tot_array = np.empty(nRows)
    oldLFweight_avg_array = np.empty(nRows)
    oldLFweight_tot_array = np.empty(nRows)

#    mag_out_array = -99*np.ones((nRows,nCols))
#    size_out_array = -99*np.ones((nRows,nCols))
#    snr_array = -99*np.ones((nRows,nCols))
#    e1_array = -99*np.ones((nRows,nCols))
#    e2_array = -99*np.ones((nRows,nCols))
#    LFweight_array = -99*np.ones((nRows,nCols))
#    SEweight_array = np.zeros((nRows,nCols))
#    fitclass_array = -99*np.ones((nRows,nCols))


    ## Shear estimates
    g_cuts = [ (g1==g1Range[g_index])&(g2==g2Range[g_index]) for g_index in xrange(8) ]
    g1_meas_array = -99*np.ones((nRows,nCols))
    g2_meas_array = -99*np.ones((nRows,nCols))
    m1_array = np.empty(nRows)
    m2_array = np.empty(nRows)
    c1_array = np.empty(nRows)
    c2_array = np.empty(nRows)
    m1err_array = np.empty(nRows)
    m2err_array = np.empty(nRows)
    c1err_array = np.empty(nRows)
    c2err_array = np.empty(nRows)

    mask_array = np.empty(nRows)
    n_occurences = np.empty(nRows)
    n_good = np.empty(nRows)

    ## Specify the cuts you want to impose here. False will get 0 weight, True will get LF weights
    elimination_cuts = ( ~((fitclass==1)|(fitclass==-3)|(fitclass==-7)|(fitclass==-1)|(fitclass==-4)|(fitclass==-10)) )&(nd>4.25)&(size_out>0.5)
    #elimination_cuts = (fitclass==0)&(nd>4.25)&(size_out>0.5)

    t_pre2 = time.time()
    print "Preprocessing time: ", np.round(t_pre2-t_pre1,2), " seconds."
    t_start = time.time()
    for row_id in xrange(nRows):
        if row_id%500==0: t1 = time.time()
        cID = row_indices[row_id]
        id_cuts = (catID==cID)

        xdata1, ydata1 = [], []
        xdata2, ydata2 = [], []
        for g_index in xrange(8):
            row_cuts =  g_cuts[g_index] & id_cuts

            if row_cuts.sum()>0:
                weights = LFweight[row_cuts]*elimination_cuts[row_cuts]

                if weights.sum()>0:
                    g1_meas_array[row_id,g_index] = np.average(e1[row_cuts],weights=weights)
                    g2_meas_array[row_id,g_index] = np.average(e2[row_cuts],weights=weights)

                    assert np.abs(g1_meas_array[row_id,g_index])<=1.
                    assert np.abs(g2_meas_array[row_id,g_index])<=1.
                else:
                    #print "For cID: {0} and g_id:{1}, total LF weight is zero.".format(cID,g_index)
                    g1_meas_array[row_id,g_index] = 0.0
                    g2_meas_array[row_id,g_index] = 0.0
                mask_array[row_id] += 2**g_index

                ## Collect data for regression
                xdata1.append( g1Range[g_index])
                xdata2.append( g2Range[g_index])
                ydata1.append( g1_meas_array[row_id, g_index] )
                ydata2.append( g2_meas_array[row_id, g_index] )

            else:
                g1_meas_array[row_id,g_index] = -99
                g2_meas_array[row_id,g_index] = -99

        ## More than one data point needed to perform regression
        if len(np.unique(xdata1))>1:
            m1c1, err1 = scipy.optimize.curve_fit(bias_function, xdata=xdata1, ydata=ydata1)
            if (np.abs(m1c1[0])>171) :
                pass
        else:
            #print "For cID:", cID, "only one data point found. Not calculating m1,c1."
            m1c1 = [-999.,-999.]
            err1 = np.array([[9801.,9801.],[9801.,9801.]])

        if len(np.unique(xdata2))>1:
            m2c2, err2 = scipy.optimize.curve_fit(bias_function, xdata=xdata2, ydata=ydata2)
            #print m1c1, m2c2
            #print xdata1
            #print ydata1
            #print xdata2
            #print ydata2
            if (np.abs(m2c2[0])>171) :
                pass
        else:
            #print "For cID:", cID, "only one data point found. Not calculating m1,c1."
            m2c2 = [-999.,-999.]
            err2 = np.array([[9801.,9801.],[9801.,9801.]])

        m1_array[row_id] = m1c1[0]
        c1_array[row_id] = m1c1[1]

        m2_array[row_id] = m2c2[0]
        c2_array[row_id] = m2c2[1]

        m1err_array[row_id] = np.sqrt(err1[0,0])
        c1err_array[row_id] = np.sqrt(err1[1,1])

        m2err_array[row_id] = np.sqrt(err2[0,0])
        c2err_array[row_id] = np.sqrt(err2[1,1])

        ## Fill in the input array & average array
        ## Calculate the weights
        weights = LFweight[id_cuts]*elimination_cuts[id_cuts]

        ## How many times does the galaxy contribute to the bias?
        n_good[row_id] = np.sum(weights>0)
        ## How many times has the galaxy been detected?
        n_occurences[row_id] = counter[cID]

        ## Keep the input quantities as it is
        e1_in_array[row_id] = e1_in[id_cuts][0]
        e2_in_array[row_id] = e2_in[id_cuts][0]
        mag_in_array[row_id] = mag_in[id_cuts][0]
        size_in_array[row_id] = size_in[id_cuts][0]
        N_in_array[row_id] = N_in[id_cuts][0]
        ZB5_in_array[row_id] = ZB5_in[id_cuts][0]

        ## Take the weighted average of output quantities
        mag_out_avg_array[row_id] = -99 if not weights.sum() else np.average(mag_out[id_cuts], weights=weights)
        size_out_avg_array[row_id] = -99 if not weights.sum() else np.average(size_out[id_cuts], weights=weights)
        e1_avg_array[row_id] = -9 if not weights.sum() else np.average(e1[id_cuts], weights=weights)
        e2_avg_array[row_id] = -9 if not weights.sum() else np.average(e2[id_cuts], weights=weights)
        snr_avg_array[row_id] = -9 if not weights.sum() else np.average(snr_model[id_cuts], weights=weights)

        ## Take the weighted average of input sheared ellipticities - LF selection bias
        sheared_e1_in_avg_array[row_id] = -9 if not weights.sum() else np.average(sheared_e1_in[id_cuts], weights=weights)
        sheared_e2_in_avg_array[row_id] = -9 if not weights.sum() else np.average(sheared_e2_in[id_cuts], weights=weights)

        ## Unweighted rotations, average and total lensfit weights
        rotation_avg_array[row_id] = rotation[id_cuts].mean()
        LFweight_avg_array[row_id] = 0 if weights.sum()==0 else np.average(LFweight[id_cuts], weights=elimination_cuts[id_cuts])
        LFweight_tot_array[row_id] = np.sum( LFweight[id_cuts]*(elimination_cuts[id_cuts]) )

        if (row_id+1)%500==0: 
            t2 = time.time()
            print "Time for ", row_id, np.round(t2-t1,2), "seconds."
            sys.stdout.flush()

    ## Fill in the output array
#    for lineno in xrange(len(catID)):
#        cID = catID[lineno]
#        row_id = row_indices.index(cID)
#
#        g_index = gRange.index( (g1[lineno], g2[lineno]) )
#        r_index = rotation[lineno]
#        ggr_index = 4*g_index + r_index
#
#        snr_array[row_id, ggr_index] = snr_model[lineno]
#        e1_array[row_id, ggr_index] = e1[lineno]
#        e2_array[row_id, ggr_index] = e2[lineno]
#        mag_out_array[row_id, ggr_index] = mag_out[lineno]
#        size_out_array[row_id, ggr_index] = size_out[lineno]
#        fitclass_array[row_id, ggr_index] = fitclass[lineno]
#        LFweight_array[row_id, ggr_index] = LFweight[lineno]
#        SEweight_array[row_id, ggr_index] = 1

    t_finish = time.time()

    print "Total time for ", row_id, " galaxies: " , np.round(t_finish-t_start,2), "seconds"
    print "Average time per galaxy: ", np.round((t_finish-t_start)/row_id,2), " seconds."
    sys.stdout.flush()
        
    ## Create collapsed arrays
    arr1 = np.rec.array([row_indices, n_occurences, n_good, mask_array, mag_in_array, size_in_array, e1_in_array, e2_in_array, N_in_array, ZB5_in_array,\
                         rotation_avg_array, LFweight_avg_array, LFweight_tot_array, oldLFweight_avg_array, oldLFweight_tot_array, snr_avg_array,\
                         mag_out_avg_array, size_out_avg_array, e1_avg_array, e2_avg_array, sheared_e1_in_avg_array, sheared_e2_in_avg_array,\
                         m1_array, m1err_array, m2_array, m2err_array, c1_array, c1err_array, c2_array, c2err_array],\
                         formats='int,int,int,int,float,float,float,float,float,float, float,float,float,float,float,float,float,float,float,float,float,float,\
                         float,float,float,float,float,float,float,float',
                         names='Cat_ID,n_occur,n_good,MASK,mag_in,size_in,e1_in,e2_in,N_in,ZB5_in,rotations_avg,weight_avg,weight_tot,oldWeight_avg,oldWeight_tot,snr_avg,\
                                mag_out_avg,size_out_avg,e1_avg,e2_avg,sheared_e1_in_avg,sheared_e2_in_avg,m1,m1err,m2,m2err,c1,c1err,c2,c2err')

#    arr2 = np.rec.array([mag_out_array,size_out_array,snr_array,weight_array,fitclass_array,e1_array,e2_array],
#                            formats='float,float,float,float,int,float,float', names='mag_out,size_out,snr_model,LFweight,fitclass,e1,e2')

#    arr3 = np.rec.array([g1_meas_array,g2_meas_array],
#                            formats='float,float', names='g1_meas,g2_meas')

    tbhdu1 = fits.BinTableHDU(data=arr1)
    output_filename = master_cat_pathname.replace('MasterCat','rearranged_MasterCat_fitcuts')
    tbhdu1.writeto(output_filename, overwrite=True)

#    tbhdu2= fits.BinTableHDU(data=arr2)
#    output_filename = master_cat_pathname.replace('MasterCat','rearrange_MasterCat_fitcuts_full')
#    tbhdu2.writeto(output_filename, overwrite=True)

if __name__=='__main__':
    rearrange(sys.argv[1])
    
