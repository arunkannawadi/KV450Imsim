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
from astropy.io import fits
import sys
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import os
#IFC17 calibration

def make_money_plot(infile, gold_flag, blind):
    infileGold = infile+'_'+gold_flag+'/'
    outputName = 'Summary_multiplicative_tomo9'+blind+'.dat'
    my_file = os.path.exists(os.path.join(infileGold,outputName))

    if(my_file==True):
        data=genfromtxt(os.path.join(infileGold, outputName))
        bin_mean=data[:,0]
        m_bin=data[:,1]
        m_err=data[:,2]
    else:
        #Name of the KiDS patches
        fields=['N', 'S']

        weight=[]
        m=[]
        ZB=[]

        for xx in range(0,len(fields)):
            fits_data = fits.open(infileGold+'Multiplicative_'+fields[xx]+blind+'_tomo_9band.cat')
            weight.extend(fits_data[1].data['w'])
            m.extend(fits_data[1].data['m'])
            ZB.extend(fits_data[1].data['ZB9'])

        weight=array(weight)
        m=array(m)
        ZB=array(ZB)

        m_bin=[]
        m_err=[]
        bin_mean=[]

        for kk in range(0, len(ZB_bins)-1):
            mask=(ZB>=ZB_bins[kk])&(ZB<ZB_bins[kk+1])&(m>-9999)
            SDfile = infileGold+'/2bin/MV_Tomo9{0}_100_SignalToNoise_ResolutionAltMV_binning_global.txt'.format(kk+1)
            if not os.path.exists(SDfile):
                bin_mean.append(0.5*(ZB_bins[kk]+ZB_bins[kk+1]))
                m_err.append(-99)
                m_bin.append(-1)
                continue
            SD = genfromtxt(SDfile, comments='#')
            err_SD = array(SD[:,8])
            weight_SD = array(SD[:,23])
            stat_err = sqrt(sum(err_SD**2))/len(err_SD)
            m_err.append(stat_err)
            m_bin.append((average(m[mask], weights=weight[mask])))
            bin_mean.append(0.5*(ZB_bins[kk]+ZB_bins[kk+1]))

    fileOut=open(os.path.join(infileGold,'Summary_multiplicative_tomo9'+blind+'.dat'), 'wb')
    for kk in range(len(m_bin)):
        print >> fileOut,  bin_mean[kk], m_bin[kk], m_err[kk]
    fileOut.close() 

    xpos=array([1.25, 2.25, 3.25, 4.25, 5.25, 6.25])
    TomoLabel=['0.1<ZB<=0.3','0.3<ZB<=0.5',  '0.5<ZB<=0.7', '0.7<ZB<=0.9','0.9<ZB<=1.2', '1.2<ZB<=2.0']
    TomoLabel=[r'$0.1\leq z_B<0.3$', r'$0.3\leq z_B < 0.5$', r'$0.5\leq z_B < 0.7$', r'$0.7\leq z_B <0.9$', r'$0.9 \leq z_B < 1.2$', r'1.2 \leq z_B \leq 2.0']

    fig=plt.figure(figsize = (15,8))
    gs1 = gridspec.GridSpec(12, 8)
    gs1.update(wspace=0.025, hspace=0.025)
    ax1=fig.add_subplot(1,1,1)
    #ax1.add_patch(patches.Rectangle((0,m_all_new_JM[0]-0.01), 1.8, 0.02,hatch='\\',fill=False, label=r'$\mathrm{KiDS-450 \, requirement}$'))
    #ax1.add_patch(patches.Rectangle((1.8,m_all_new_JM[1]-0.01), 1.0, 0.02,hatch='\\',fill=False))
    #ax1.add_patch(patches.Rectangle((2.8,m_all_new_JM[2]-0.01), 1.0, 0.02,hatch='\\',fill=False))
    #ax1.add_patch(patches.Rectangle((3.8,m_all_new_JM[3]-0.01), 1.0, 0.02,hatch='\\',fill=False))
    #ax1.add_patch(patches.Rectangle((4.8,m_all_new_JM[4]-0.01), 1.0, 0.02,hatch='\\',fill=False))
    ax1.add_patch(patches.Rectangle((0,m_bin[0]-0.02), 1.8, 0.04, hatch='\\', fill=False, label=r'$\mathrm{KV450 \, requirement \, (\pm 2\%) }$'))
    ax1.add_patch(patches.Rectangle((1.8,m_bin[1]-0.02), 1.0, 0.04, hatch='\\', fill=False))
    ax1.add_patch(patches.Rectangle((2.8,m_bin[2]-0.02), 1.0, 0.04, hatch='\\', fill=False))
    ax1.add_patch(patches.Rectangle((3.8,m_bin[3]-0.02), 1.0, 0.04, hatch='\\', fill=False))
    ax1.add_patch(patches.Rectangle((4.8,m_bin[4]-0.02), 1.0, 0.04, hatch='\\', fill=False))

    #ax1.add_patch(patches.Rectangle((0.0, -0.01), 6, 0.02,hatch='\\',fill=False))
    plt.errorbar(xpos, m_bin, yerr=m_err, fmt='o',label='$\mathrm{Final \, calibration \, (3D)}$', color='black',linewidth=3)

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
    plt.savefig(infileGold+'Money_plot'+blind+'.pdf')
    print "Saved the plot in ", infileGold
    plt.show()

if __name__=='__main__':
    infile=sys.argv[1]

    ZB_bins=[0.1, 0.301,0.501,0.701,0.901,1.201, 2.001]
    gold_flags = ['Fid'] #['nogold','Fid','noDEEP2','noVVDS','nozCOSMOS','multispec3','speczquality4']
    blinds = ['_A','_B','_C'][2:]
    for gold_flag in gold_flags:
        for blind in blinds:
            make_money_plot(infile, gold_flag, blind)
