import numpy
from numpy import *
import pyfits
import sys
import scipy 
from scipy.interpolate import griddata
#import pdb; pdb.set_trace()

#Code to assign m-correction to each galaxy
#The correction is calculated on a 20x20 grid Signal-to-noise & Resolution
#where the Resolution is defined as:
#
#R=size_psf/((size_ab)**2+size_psf)
#
#where
#
#size_ab=LF_scalength*((1.-|e|)/(1.+|e|))**0.5
#
#and
#
#size_psf=(psf_Q11*psf_Q22-psf_Q12**2)**0.5
#
#Massimo Viola & Julian Merten & Arun Kannawadi(August 2017)

#################################################################################

def SizeFromMom(psf_Q11_sim,psf_Q22_sim,psf_Q12_sim):
    return (psf_Q11_sim*psf_Q22_sim-psf_Q12_sim**2)**0.5

def bias_2D_bin_query(obs1, obs2,data_table):
    query_index = -1
    for index in range(0,len(data_table)):
        if((obs1>=data_table[index][1]) & (obs1 < data_table[index][2]) & (obs2 >= data_table[index][4]) & (obs2 < data_table[index][5])):
            query_index = index
            break
      
  
    if(query_index != -1):
        output = (data_table[query_index][7],data_table[query_index][11],data_table[query_index][9],data_table[query_index][13])
    else:
        output = (0.0,0.0,0.0,0.0)
        #output = (-99,-99,0.0,0.0)

    output=array(output)

    return output

#####################################################################################

## Compare the KiDS-450 and KV450 catalogues
def compare_catalogues():
    #Name of the KiDS patches
    patches=['G9', 'G12', 'G15', 'G23', 'GS']
    #Name of the fields to check
    fieldnames = ['e1_A','e2_A','weight_A','PSF_Q11','PSF_Q12','PSF_Q22','model_SNratio','bias_corrected_scalelength_pixels']

    
    #Location of the KiDS-450 catalogues
    Dir_KiDS450='/disks/shear10/KiDS/KiDS-450/SHEAR-PZ-CATALOGUES/'
    #Name of the LF reweighting scheme
    data_reweight='reweight_5x5x5'

    #Location of the KV450 catalogues
    Dir_KV450='/disks/shear10/hendrik/KV450_CATALOGUES_PATCH_V0.5.9/'
    Dir_KV450='/disks/shear14/KIDS/KV450_CATALOGUES_PATCH_V0.5.9/'

    for xx in range(0,len(patches)):
        print 'Comparing ', patches[xx]

        old_cat = pyfits.open(Dir_KiDS450+'KiDS_'+patches[xx]+'_'+data_reweight+'_BLIND_PF.cat')
        new_cat = pyfits.open(Dir_KV450+'KV450_'+patches[xx]+'.cat')

        ## The data is stored at different indices
        old_dat = old_cat[2].data
        new_dat = new_cat[1].data

        for fieldname in fieldnames:
            try:
                numpy.testing.assert_array_equal(old_dat[fieldname],new_dat[fieldname])
            except AssertionError:
                print "Mismatch found in ", patches[xx], " field for ", fieldname

        ## Compare the 4-band photo-z
        try:
            numpy.testing.assert_array_equal(old_dat['Z_B'],new_dat['Z_B_ugri'])
        except AssertionError:
            print "Mismatch found in ", patches[xx], " field for 4-band photo-z"

#####################################################################################

#Read in the 20x20 m-bias surface
SD=sys.argv[1]
DirOut=sys.argv[2]
#Tomo flag
TomoFlag=int(sys.argv[3])
nband=sys.argv[4]

## The following catalogue location is for KiDS-450:
#####Location of the KiDS  catalogues
####Dir_KiDS='/disks/shear10/KiDS/KiDS-450/SHEAR-PZ-CATALOGUES/'
#####Name of the LF reweighting scheme
####data_reweight='reweight_5x5x5'

use_old_data = False ## set True to use the catalogues with wrong LF weights and wrong 9-bband ZB
## The following catalogue location is for KV450:
if use_old_data:
    Dir_KV450='/disks/shear10/hendrik/KV450_CATALOGUES_PATCH_V0.5.9/'  ## these have wrong LF weights and wrong 9-band ZB
else:
    Dir_KV450='/disks/shear14/KIDS/KV450_CATALOGUES_PATCH_V0.5.9_GOLD/' ## @Angus: Enter the path to the KV450 GOLD catalogues here
#Name of the KiDS patches
patches=['G9', 'G12', 'G15', 'G23', 'GS']

## Toggle comment to compare the catalogues
#compare_catalogues()

#Looping through the patches
for xx in range(0,len(patches)):
    print 'Patch: ', patches[xx]
    binary_file_name='Multiplicative_'+patches[xx]+'_'+nband+'band.cat'
    binary_file_name_tomo='Multiplicative_'+patches[xx]+'_tomo_'+nband+'band.cat'
        
    #Open the IMSIM catalogue
    #hdulist2 = pyfits.open(Dir_KiDS+'/KiDS_'+patches[xx]+'_'+data_reweight+'_BLIND_PF.cat')
    if use_old_data:
        hdulist2 = pyfits.open(Dir_KV450+'KV450_'+patches[xx]+'.cat') ## this is to open the catalogues with wrong LF weights and wrong 9-band ZB
    else:
        hdulist2 = pyfits.open(Dir_KV450+'KV450_'+patches[xx]+'_reweight_3x4x4_v2_good_goldclasses.cat')

    #tbdata = hdulist2[2].data ## data is in location 2 for KiDS-450
    tbdata = hdulist2[1].data ## data is in location 1 for KV450
    ## Old data had all the cuts already imposed. For new data, we have to do it.
    if not use_old_data:
        ## Really shouldn't have to use the binary star cuts again
        data_cuts = ((tbdata['bias_corrected_scalelength_pixels']<tbdata['binary_cut'])|(numpy.hypot(tbdata['bias_corrected_e1'],tbdata['bias_corrected_e2'])<0.8))&(tbdata['GAAP_Flag_ugriZYJHKs']==0)
        #data_cuts &= (tbdata['T_B']>1.9)
        tbdata = tbdata[data_cuts]
    else:
        data_cuts = tbdata['Flag_SOM_speczquality4_NONE']==1 ## @Angus: Enter the cuts you want to impose on the KV450 GOLD catalogues here
        tbdata = tbdata[data_cuts]

    ID=tbdata.SeqNr
    tile=tbdata.THELI_INT
    e1_data_all_A = tbdata.bias_corrected_e1 #.e1_A
    e2_data_all_A = tbdata.bias_corrected_e2 #.e2_A
    ZB4_data_all_A = tbdata.Z_B_ugri
    ZB9_data_all_A = tbdata.Z_B
    if nband=='4':
        ZB_data_all_A = ZB4_data_all_A
    elif nband=='9':
        ZB_data_all_A = ZB9_data_all_A
    else:
        print "'nband' has to be 4 or 9. Exiting..."
        exit()

    size_data_all = tbdata.bias_corrected_scalelength_pixels
    SNR_data_all = tbdata.model_SNratio
    Q11_psf_data_all = tbdata.PSF_Q11
    Q22_psf_data_all = tbdata.PSF_Q22
    Q12_psf_data_all = tbdata.PSF_Q12
    #weight_data_all_A = tbdata.weight_A
    weight_data_all_A = tbdata.recal_weight

    #Define PSF size
    size_psf_data_all= SizeFromMom(Q11_psf_data_all,Q22_psf_data_all,Q12_psf_data_all)

    #Define |e| for the 3 blindings
    mode_data_all_A=(e1_data_all_A**2+e2_data_all_A**2)**0.5

    #Define circularised galaxy size
    size_ab_data_all_A=size_data_all*((1.-mode_data_all_A)/(1.+mode_data_all_A))**0.5

    #Define galaxy 'resolution'
    res_ab_data_all_A=size_psf_data_all/((size_ab_data_all_A)**2+size_psf_data_all)
    
    if(TomoFlag==0):
        print 'standard analysis'
        #Read in surface of doom 
        data=genfromtxt(SD, comments='#')

        #m1_A=zeros(len(SNR_data_all))
        #m2_A=zeros(len(SNR_data_all))

        m1_A=empty(len(SNR_data_all))
        m2_A=empty(len(SNR_data_all))
        m1_A.fill(-9999)
        m2_A.fill(-9999)
#        m1_A.fill(0.)
#        m2_A.fill(0.)
####### Changed by AKJ and reverted back

        #Apply calibration

        SNR_min=data[:,1]
        SNR_max=data[:,2]
        R_min=data[:,4]
        R_max=data[:,5]
        m1_grid=data[:,7]
        m2_grid=data[:,11]

        dim=len(SNR_min)
        
        for kk in range(0, dim):
            mask=(SNR_data_all>=SNR_min[kk]) & (SNR_data_all < SNR_max[kk]) & (res_ab_data_all_A >= R_min[kk]) & (res_ab_data_all_A < R_max[kk])
            m1_A[mask]=m1_grid[kk]
            m2_A[mask]=m2_grid[kk]

        #Compute average m bias
    
        m1_A=array(m1_A)
        m2_A=array(m2_A)
        m_A=(m1_A+m2_A)/2.0
   
        #Write output on a fits file
        c1=pyfits.Column(name='m',format='D',array=m_A)
        c2=pyfits.Column(name='SeqNr',format='D',array=ID)
        c3=pyfits.Column(name='THELI_INT',format='D',array=tile)
        c4=pyfits.Column(name='w',format='D',array=weight_data_all_A)
        c5=pyfits.Column(name='SNR',format='D',array=SNR_data_all)
        c6=pyfits.Column(name='R',format='D',array=res_ab_data_all_A)
        c7=pyfits.Column(name='ZB4',format='D',array=ZB4_data_all_A)
        c8=pyfits.Column(name='ZB9',format='D',array=ZB9_data_all_A)
        
        tbhdu = pyfits.new_table([c1, c2, c3, c4, c5, c6,c7,c8])
        tbhdu.writeto(DirOut+'/'+binary_file_name, clobber=True)
    else:
        #Tomo=[0.101,0.301,0.501,0.701,0.901,10.]
        Tomo=[0.101,0.301,0.501,0.701,0.901,1.201]
        m_A=[]
        SNR_full=[]
        R_full=[]
        w_full=[]
        SeqNr_full=[]
        THELI_INT_full=[]
        ZB_full=[]
        
        for tt in range(0,5):
            maskBin=(ZB_data_all_A>=Tomo[tt])&(ZB_data_all_A<Tomo[tt+1])
            SNR_bin=SNR_data_all[maskBin]
            R_bin=res_ab_data_all_A[maskBin]
            w_bin=weight_data_all_A[maskBin]
            SeqNr_bin=ID[maskBin]
            THELI_INT_bin=tile[maskBin]
            ZB_bin=ZB_data_all_A[maskBin]
            
            #m1_A=zeros(len(SNR_bin))
            #m2_A=zeros(len(R_bin))
            m1_A=empty(len(SNR_bin))
            m2_A=empty(len(SNR_bin))
            m1_A.fill(-9999)
            m2_A.fill(-9999)
            #m1_A.fill(0.)
            #m2_A.fill(0.)
            ### changed by AKJ and reverted back
            #Read in surface of doom 
           
            try:
                data=genfromtxt(DirOut+'/2bin/MV_Tomo'+nband+str(tt+1)+'_100_SignalToNoise_ResolutionAltMV_binning_global.txt', comments='#')
            except:
                print "Reading through genfromtxt failed for tt={0}. Using loadtxt".format(tt)
                def line_to_array(line):
                    return array( [float(num) for num in line[:-1].split(' ')] )
                with open(DirOut+'/2bin/MV_Tomo'+nband+str(tt+1)+'_100_SignalToNoise_ResolutionAltMV_binning_global.txt','r') as fp:
                    datalines = fp.readlines()
                data = stack([line_to_array(dline) for dline in datalines[1:]])

            SNR_min=data[:,1]
            SNR_max=data[:,2]
            R_min=data[:,4]
            R_max=data[:,5]
            m1_grid=data[:,7]
            m2_grid=data[:,11]

            dim=len(SNR_min)
            for kk in range(0, dim):
                mask=(SNR_bin>=SNR_min[kk]) & (SNR_bin < SNR_max[kk]) & (R_bin >= R_min[kk]) & (R_bin < R_max[kk])
                m1_A[mask]=m1_grid[kk]
                m2_A[mask]=m2_grid[kk]
                #print mask
            #Compute average m bias

            print 'tomo bin..', tt, len(m1_A), len(ZB_bin)
            m1_A=array(m1_A)
            m2_A=array(m2_A)
            m_bin_A=(m1_A+m2_A)/2.0   
            
            m_A.extend(m_bin_A)
            SNR_full.extend(SNR_bin)
            R_full.extend(R_bin)
            w_full.extend(w_bin)
            SeqNr_full.extend(SeqNr_bin)
            THELI_INT_full.extend(THELI_INT_bin)
            ZB_full.extend(ZB_bin)
        m_A=array(m_A)
        SNR_full=array(SNR_full)
        R_full=array(R_full)
        w_full=array(w_full)
        SeqNr_full=array(SeqNr_full)
        THELI_INT_full=array(THELI_INT_full)
        ZB_full=array(ZB_full)
        #Write output on a fits file
        c1=pyfits.Column(name='m',format='D',array=m_A)
        c2=pyfits.Column(name='SeqNr',format='D',array=SeqNr_full)
        c3=pyfits.Column(name='THELI_INT',format='D',array=THELI_INT_full)
        c4=pyfits.Column(name='w',format='D',array=w_full)
        c5=pyfits.Column(name='SNR',format='D',array=SNR_full)
        c6=pyfits.Column(name='R',format='D',array=R_full)
        c7name = 'ZB4' if nband=='4' else 'ZB9'
        c7=pyfits.Column(name=c7name,format='D',array=ZB_full)
        tbhdu = pyfits.new_table([c1, c2, c3, c4, c5, c6,c7])
        tbhdu.writeto(DirOut+'/'+binary_file_name_tomo, clobber=True)
        #tbhdu.writeto(binary_file_name_tomo, clobber=True)

