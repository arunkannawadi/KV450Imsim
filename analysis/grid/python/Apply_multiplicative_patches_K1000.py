import numpy
from numpy import *
import pyfits
from astropy.io import fits
import os, sys
import optparse
import scipy
from scipy.interpolate import griddata

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
def line_to_array(line):
    return array( [float(num) for num in line[:-1].split(' ')] )

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

# Set up run-time args.
argopts = optparse.OptionParser("usage: %prog [options] arg1")
#argopts.add_option("-g", "--gold", dest="gold_flag",
#                   default="Fid",
#                   type="string",
#                   help="Specify the gold flag (Fid, noDEEP2 etc.)")
argopts.add_option("-p", "--path", dest="Dir_KV450",
                    default="/disks/shear14/KiDS/K1000/",
                    type="string",
                    help="Specify the path to the K-1000 catalogues")

defaultDirOut = os.path.join(*(SD.split('/')[:-2]))
argopts.add_option("-o", "--result", dest="DirOut",
                    default=defaultDirOut,
                    type="string",
                    help="Specify the path to the results directory")
argopts.add_option("-r", "--reuse", dest="reuse",
                    default="False",
                    type="string",
                    help="Set to True to skip generating existing files again")
(options,args) = argopts.parse_args()
#gold_flag = options.gold_flag
reuse = options.reuse
reuse = True if reuse=='True' else False
DirOut = options.DirOut
Dir_KV450 = options.Dir_KV450
#resultsDir = options.resultsDir

print("DirOut = ", DirOut)
DirOut = '/'+DirOut.rstrip('/') ## strip out any trailing / and add a leading /

## The following catalogue location is for KiDS-450:
#####Location of the KiDS  catalogues
####Dir_KiDS='/disks/shear10/KiDS/KiDS-450/SHEAR-PZ-CATALOGUES/'
#####Name of the LF reweighting scheme
####data_reweight='reweight_5x5x5'

data_type = 'K1000' # KV450/ K1000
Tomo=[0.101,0.301,0.501,0.701,0.901,1.201, 2.001]
use_old_data = False ## set True to use the catalogues with wrong LF weights and wrong 9-bband ZB
## The following catalogue location is for KV450:
if data_type=='KV450':
    if use_old_data:
        Dir_KV450='/disks/shear10/hendrik/KV450_CATALOGUES_PATCH_V0.5.9/'  ## these have wrong LF weights and wrong 9-band ZB
    else:
        Dir_KV450='/disks/shear14/KIDS/KV450_CATALOGUES_PATCH_V0.5.9_GOLD/' ## @Angus: Enter the path to the KV450 GOLD catalogues here
    #Name of the KiDS patches
    patches=['G9', 'G12', 'G15', 'G23', 'GS']
    blinds = ['']

    ## Toggle comment to compare the catalogues
    #compare_catalogues()
else: ## K1000
    patches = ['N','S']
    blinds = ['_A','_B','_C'][2:]
    gold_flags = ['Fid'] #['nogold', 'Fid','noDEEP2', 'noVVDS', 'nozCOSMOS', 'multispec3', 'speczquality4']

#Looping through the patches
for pid, patch in enumerate(patches):
    print 'Patch: ', patch
    #binary_file_name='Multiplicative_'+patch+'_'+nband+'band.cat' #DELETEIT
        
    #Open the K-1000 catalogue
    hdulist2 = fits.open(os.path.join(Dir_KV450, 'K1000_{}_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_THELI_INT.cat'.format(patch)))
    tbdata = hdulist2[1].data ## data is in location 1 for KV450 and K-1000

    #tbdata = hdulist2[2].data ## data is in location 2 for KiDS-450
    
    ## Really shouldn't have to use the binary star cuts again
    if numpy.sum(tbdata['FLAG_GAAP_ugriZYJHKs']):
        data_cuts = tbdata['FLAG_GAAP_ugriZYJHKs']==0 ## This is likely to hold for all records
        data = tbdata[data_cuts]

    for blind in blinds:
        for gold_flag in gold_flags:
            DirOutGold = DirOut+'_'+gold_flag
            binary_file_name_tomo='Multiplicative_'+patch+blind+'_tomo_9band.cat'
            if os.path.exists(os.path.join(DirOutGold,binary_file_name_tomo)) and reuse: continue
            
            if gold_flag=='nogold':
                data = tbdata
            else:
                som_cuts = tbdata['Flag_SOM_{0}{1}'.format(gold_flag, blind)]==1 ## @Angus: Enter the cuts you want to impose on the KV450 GOLD catalogues here
                data = tbdata[som_cuts]
            
            ID=data.SeqNr
            tile=data.THELI_INT
            weight_data_all_A = data['recal_weight'+blind]
            e1_data_all_A = data['autocal_e1'+blind] #.e1_A
            e2_data_all_A = data['autocal_e2'+blind] #.e2_A
            ZB9_data_all_A = data.Z_B
            ZB_data_all_A = ZB9_data_all_A

            size_data_all = data.bias_corrected_scalelength_pixels
            SNR_data_all = data.model_SNratio
            Q11_psf_data_all = data.PSF_Q11
            Q22_psf_data_all = data.PSF_Q22
            Q12_psf_data_all = data.PSF_Q12
            #weight_data_all_A = data.weight_A

            #Define PSF size
            size_psf_data_all= SizeFromMom(Q11_psf_data_all,Q22_psf_data_all,Q12_psf_data_all)
            try:
                numpy.testing.assert_array_equal(size_psf_data_all, data.PSFsize)
            except AssertionError:
                print("The PSF size definitions do not agree")

            #Define |e| for the 3 blindings
            mode_data_all_A= numpy.hypot(e1_data_all_A, e2_data_all_A)

            #Define circularised galaxy size
            size_ab_data_all_A=size_data_all*((1.-mode_data_all_A)/(1.+mode_data_all_A))**0.5

            #Define galaxy 'resolution'
            res_ab_data_all_A=size_psf_data_all/((size_ab_data_all_A)**2+size_psf_data_all)
    
            m_blind=[]
            SNR_full=[]
            R_full=[]
            w_blind=[]
            SeqNr_full=[]
            THELI_INT_full=[]
            ZB_full=[]
            
            for tt in range(0,len(Tomo)-1):
                maskBin=(ZB_data_all_A>=Tomo[tt])&(ZB_data_all_A<Tomo[tt+1])
                SNR_bin=SNR_data_all[maskBin]
                R_bin=res_ab_data_all_A[maskBin]
                w_bin=weight_data_all_A[maskBin]
                SeqNr_bin=ID[maskBin]
                THELI_INT_bin=tile[maskBin]
                ZB_bin=ZB_data_all_A[maskBin]
                
                m1_blind=empty(len(SNR_bin))
                m2_blind=empty(len(SNR_bin))
                m1_blind.fill(-9999)
                m2_blind.fill(-9999)
                
                #Read in surface of doom
                sod_pathname = DirOutGold+'/2bin/MV_Tomo9'+str(tt+1)+'_100_SignalToNoise_ResolutionAltMV_binning_global.txt'
                if not os.path.exists(sod_pathname): break ## happens if we are not processing the current gold_flag
                try:
                    sod=genfromtxt(sod_pathname, comments='#')
                except:
                    ## This block shouldn't get executed normally, unless there are some horrible formatting issues with the surface-of-doom file
                    print "WARNING: Reading through genfromtxt failed for tt={0}. Using loadtxt".format(tt)
                    with open(sod_pathname,'r') as fp:
                        datalines = fp.readlines()
                    sod = stack([line_to_array(dline) for dline in datalines[1:]])

                SNR_min=sod[:,1]
                SNR_max=sod[:,2]
                R_min=sod[:,4]
                R_max=sod[:,5]
                m1_grid=sod[:,7]
                m2_grid=sod[:,11]

                dim=len(SNR_min)
                for kk in range(0, dim):
                    mask=(SNR_bin>=SNR_min[kk]) & (SNR_bin < SNR_max[kk]) & (R_bin >= R_min[kk]) & (R_bin < R_max[kk])
                    m1_blind[mask]=m1_grid[kk]
                    m2_blind[mask]=m2_grid[kk]

                #Compute average m bias
                print 'tomo bin..', tt, len(m1_blind), len(ZB_bin)
                m1_blind=array(m1_blind)
                m2_blind=array(m2_blind)
                m_bin_blind=0.5*(m1_blind+m2_blind)
                
                m_blind.extend(m_bin_blind)
                SNR_full.extend(SNR_bin)
                R_full.extend(R_bin)
                w_blind.extend(w_bin)
                SeqNr_full.extend(SeqNr_bin)
                THELI_INT_full.extend(THELI_INT_bin)
                ZB_full.extend(ZB_bin)

            if len(m_blind)==0: continue ## happens if we are not processing the current gold_flag

            m_blind=array(m_blind)
            SNR_full=array(SNR_full)
            R_full=array(R_full)
            w_blind=array(w_blind)
            SeqNr_full=array(SeqNr_full)
            THELI_INT_full=array(THELI_INT_full)
            ZB_full=array(ZB_full)
            
            #Write output on a fits file
            c1=pyfits.Column(name='m',format='D',array=m_blind)
            c2=pyfits.Column(name='SeqNr',format='D',array=SeqNr_full)
            c3=pyfits.Column(name='THELI_INT',format='D',array=THELI_INT_full)
            c4=pyfits.Column(name='w',format='D',array=w_blind)
            c5=pyfits.Column(name='SNR',format='D',array=SNR_full)
            c6=pyfits.Column(name='R',format='D',array=R_full)
            c7=pyfits.Column(name='ZB9',format='D',array=ZB_full)
            tbhdu = pyfits.new_table([c1, c2, c3, c4, c5, c6,c7])
            print("Writing out the file ", os.path.join(DirOutGold, binary_file_name_tomo))
            tbhdu.writeto(os.path.join(DirOutGold, binary_file_name_tomo), clobber=True)
