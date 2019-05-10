import numpy
from numpy import *
import scipy
from scipy import *
import matplotlib
if __name__=='__main__':
    matplotlib.use('Agg')
    matplotlib.rcParams['ps.useafm']=True
    matplotlib.rcParams['pdf.use14corefonts']=True
    matplotlib.rcParams['text.usetex']=True
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pyfits
import sys
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import os.path
#IFC17 calibration

def make_money_plot(infile, TomoFlag, nband):
    if False: ## data lost with veersemeer
        DirOld_JM='/disks/shear10/viola/KIDS/CosmicShear/ImageSim/src/kids_sim_analysis/python/ReSampling/Correction_cat/Final'
        field=['G9', 'G12', 'G15', 'G23', 'GS']
        m_1_new_JM=[]
        m_2_new_JM=[]
        m_3_new_JM=[]
        m_4_new_JM=[]
        m_5_new_JM=[]
        for xx in range(0,5):
            dataNew_JM=genfromtxt(DirOld_JM+'/Calib_30520_SignalToNoise_'+field[xx]+'.contRadius.dat')
            m_1_new_JM.append((dataNew_JM[0,2]+dataNew_JM[0,3])/2.0)
            m_2_new_JM.append((dataNew_JM[1,2]+dataNew_JM[1,3])/2.0)
            m_3_new_JM.append((dataNew_JM[2,2]+dataNew_JM[2,3])/2.0)
            m_4_new_JM.append((dataNew_JM[3,2]+dataNew_JM[3,3])/2.0)
            m_5_new_JM.append((dataNew_JM[4,2]+dataNew_JM[4,3])/2.0)


        m_all_new_JM=[average(m_1_new_JM),average(m_2_new_JM),average(m_3_new_JM),average(m_4_new_JM),average(m_5_new_JM) ]
        err=array([0.002,0.002,0.002,0.002,0.002])
    #Import table

    my_file = os.path.exists(infile+'/'+'Summary_multiplicative.dat')
    if(TomoFlag=='Yes'):
        my_file_tomo = os.path.exists(infile+'/'+'Summary_multiplicative_tomo'+nband+'.dat')
    elif(TomoFlag=='No'):
        my_file_tomo = os.path.exists(infile+'/'+'Summary_multiplicative.dat')
    else:
        raise ValueError("TomoFlag has to be Yes or No")

    if((my_file==True)&(my_file_tomo==True)): # This is buggy if TomoFlag is 'No'
        data=genfromtxt(infile+'/'+'Summary_multiplicative.dat')
        bin_mean=data[:,0]
        m_bin=data[:,1]
        m_err=data[:,2]
        data_tomo=genfromtxt(infile+'/'+'Summary_multiplicative_tomo'+nband+'.dat')
        bin_mean_tomo=data_tomo[:,0]
        m_bin_tomo=data_tomo[:,1]
        m_err_tomo=data_tomo[:,2]
    else:
        #Name of the KiDS patches
        fields=['G9', 'G12', 'G15', 'G23', 'GS']

        weight=[]
        m=[]
        ZB=[]
        
        for xx in range(0,len(fields)):
            fits_data = pyfits.open(infile+'Multiplicative_'+fields[xx]+'_tomo_'+nband+'band.cat')
            #fits_data = pyfits.open(infile+'Multiplicative_'+fields[xx]+'_'+nband+'band.cat')
            weight.extend(fits_data[1].data['w'])
            m.extend(fits_data[1].data['m'])
            ZB.extend(fits_data[1].data['ZB'+nband])

        weight=array(weight)
        m=array(m)
        ZB=array(ZB)

        weight_tomo=[]
        m_tomo=[]
        ZB_tomo=[]
        
        for xx in range(0,len(fields)):
            if(TomoFlag=='Yes'):
                fits_data_tomo = pyfits.open(infile+'Multiplicative_'+fields[xx]+'_tomo_'+nband+'band.cat')
            elif(TomoFlag=='No'):
                fits_data_tomo = pyfits.open(infile+'Multiplicative_'+fields[xx]+'_'+nband+'band.cat')
            else:
                raise ValueError("TomoFlag must be Yes or No")

            weight_tomo.extend(fits_data_tomo[1].data['w'])
            m_tomo.extend(fits_data_tomo[1].data['m'])
            ZB_tomo.extend(fits_data_tomo[1].data['ZB'+nband])

        weight_tomo=array(weight_tomo)
        m_tomo=array(m_tomo)
        ZB_tomo=array(ZB_tomo)

        #Import surface of doom (for errors)

#        SDfile=sys.argv[3]
#        SD=genfromtxt(SDfile, comments='#')

#        err_SD=array(SD[:,8])
#        weight_SD=array(SD[:,23])
        stat_err = -1
#        stat_err=sqrt(sum(err_SD**2))/len(err_SD)
#    #    stat_err_tomo=stat_err*sqrt(5.)
        
        ZB_bins=[0.1, 0.301,0.501,0.701,0.901,1.201] ## changed 10 to 1.201

        m_bin=[]
        m_err=[]
        bin_mean=[]
        m_bin_tomo=[]
        m_err_tomo=[]
        bin_mean_tomo=[]
     
        for kk in range(0, len(ZB_bins)-1):
            mask_tomo=(ZB_tomo>=ZB_bins[kk])&(ZB_tomo<ZB_bins[kk+1])&(m_tomo>-9999)
            SDfile_tomo = infile+'/2bin/MV_Tomo{0}{1}_100_SignalToNoise_ResolutionAltMV_binning_global.txt'.format(nband,kk+1)
            SD_tomo = genfromtxt(SDfile_tomo, comments='#')
            err_SD_tomo = array(SD_tomo[:,8])
            weight_SD_tomo = array(SD_tomo[:,23])
            stat_err_tomo = sqrt(sum(err_SD_tomo**2))/len(err_SD_tomo)
            m_err_tomo.append(stat_err_tomo)

        ## m>-9999 condition added by AKJ
        for kk in range(0, len(ZB_bins)-1):
            mask=(ZB>=ZB_bins[kk])&(ZB<ZB_bins[kk+1])&(m>-9999)
            m_bin.append((average(m[mask], weights=weight[mask])))
            m_err.append(stat_err)
            bin_mean.append(0.5*(ZB_bins[kk]+ZB_bins[kk+1]))

        for kk in range(0, len(ZB_bins)-1):
            mask_tomo=(ZB_tomo>=ZB_bins[kk])&(ZB_tomo<ZB_bins[kk+1])&(m_tomo>-9999)
            m_bin_tomo.append((average(m_tomo[mask_tomo], weights=weight_tomo[mask_tomo])))
    #        m_err_tomo.append(stat_err_tomo)
            bin_mean_tomo.append(0.5*(ZB_bins[kk]+ZB_bins[kk+1]))

#        fileOut=open(infile+'/'+'Summary_multiplicative.dat', 'wb')
#        for kk in range(0, 5):
#            print >> fileOut,  bin_mean[kk], m_bin[kk], m_err[kk]
#        fileOut.close()
        fileOut_tomo=open(infile+'/'+'Summary_multiplicative_tomo'+nband+'.dat', 'wb')
        for kk in range(0, 5):
            print >> fileOut_tomo,  bin_mean_tomo[kk], m_bin_tomo[kk], m_err_tomo[kk]
        fileOut_tomo.close()

    ## Added by AKJ begins
    #with open(infile+'/'+'FC17.dat','w') as fileOutFC17:
    #    for kk in xrange(5):
    #        print >> fileOutFC17, bin_mean_tomo[kk], m_all_new_JM[kk], err[kk]
    FC17values = numpy.loadtxt(infile+'/FC17.dat')
    m_all_new_JM = FC17values[:,1]
    err = FC17values[:,2]
        
    ## Added by AKJ ends

    xpos=array([1.25, 2.25, 3.25, 4.25, 5.25])
    TomoLabel=['0.1<ZB<=0.3','0.3<ZB<=0.5',  '0.5<ZB<=0.7', '0.7<ZB<=0.9','0.9<ZB<=1.2']
    TomoLabel=[r'$0.1\leq z_B<0.3$', r'$0.3\leq z_B < 0.5$', r'$0.5\leq z_B < 0.7$', r'$0.7\leq z_B <0.9$', r'$0.9 \leq z_B < 1.2$']

    fig=plt.figure(figsize = (15,8))
    gs1 = gridspec.GridSpec(12, 8)
    gs1.update(wspace=0.025, hspace=0.025)
    ax1=fig.add_subplot(1,1,1)
    #ax1.add_patch(patches.Rectangle((0,m_all_new_JM[0]-0.01), 1.8, 0.02,hatch='\\',fill=False, label=r'$\mathrm{KiDS-450 \, requirement}$'))
    #ax1.add_patch(patches.Rectangle((1.8,m_all_new_JM[1]-0.01), 1.0, 0.02,hatch='\\',fill=False))
    #ax1.add_patch(patches.Rectangle((2.8,m_all_new_JM[2]-0.01), 1.0, 0.02,hatch='\\',fill=False))
    #ax1.add_patch(patches.Rectangle((3.8,m_all_new_JM[3]-0.01), 1.0, 0.02,hatch='\\',fill=False))
    #ax1.add_patch(patches.Rectangle((4.8,m_all_new_JM[4]-0.01), 1.0, 0.02,hatch='\\',fill=False))
    ax1.add_patch(patches.Rectangle((0,m_bin_tomo[0]-0.02), 1.8, 0.04, hatch='\\', fill=False, label=r'$\mathrm{KV450 \, requirement \, (\pm 2\%) }$'))
    ax1.add_patch(patches.Rectangle((1.8,m_bin_tomo[1]-0.02), 1.0, 0.04, hatch='\\', fill=False))
    ax1.add_patch(patches.Rectangle((2.8,m_bin_tomo[2]-0.02), 1.0, 0.04, hatch='\\', fill=False))
    ax1.add_patch(patches.Rectangle((3.8,m_bin_tomo[3]-0.02), 1.0, 0.04, hatch='\\', fill=False))
    ax1.add_patch(patches.Rectangle((4.8,m_bin_tomo[4]-0.02), 1.0, 0.04, hatch='\\', fill=False))
    eb1=plt.errorbar(xpos+0.2,m_all_new_JM, yerr=err,fmt='8', label=r'$\mathrm{FC17}$', linewidth=2, color='red', markersize=8, fillstyle='none')

    #ax1.add_patch(patches.Rectangle((0.0, -0.01), 6, 0.02,hatch='\\',fill=False))
    plt.errorbar(xpos, m_bin, yerr=m_err, fmt='s',label='$\mathrm{New \, calibration \, (2D)}$', color='blue')
    plt.errorbar(xpos+0.1, m_bin_tomo, yerr=m_err_tomo, fmt='o',label='$\mathrm{Final \, calibration \, (3D)}$', color='black',linewidth=3)

    plt.axhline(0, linestyle='-.', linewidth=3, color='black')
    plt.axvline(0.8, color='black')
    plt.axvline(1.8, color='black')
    plt.axvline(2.8, color='black')
    plt.axvline(3.8, color='black')
    plt.axvline(4.8, color='black')
    plt.axvline(5.8, color='black')
    #plt.ylim(-0.08, 0.005)
    plt.ylim(-0.08,0.05)
    plt.xlim(0.8,5.8)
    plt.xticks(xpos,TomoLabel, rotation='horizontal', fontsize=22)
    plt.yticks(fontsize=22)
    plt.ylabel(r'$\mathrm{m}$', fontsize=28)
    plt.legend(loc=0, prop={'size':18}, shadow=True)
    plt.tight_layout()
    #plt.savefig(infile+'Money_plot.png')
    plt.savefig(infile+'Money_plot_'+nband+'band.pdf')
    print "Saved the plot in ", infile
    plt.show()

if __name__=='__main__':
    infile=sys.argv[1]
    TomoFlag=sys.argv[2]
    nband=sys.argv[4] ## must be one of 4 or 9
    if not(nband is '4' or nband is '9'):
        print "'nband' has to be 4 or 9. Exiting"
        exit()

    make_money_plot(infile, TomoFlag, nband)
