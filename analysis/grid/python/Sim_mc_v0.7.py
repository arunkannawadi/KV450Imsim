import sys
from sys import * 
from scipy import integrate
import scipy
import re 
import io
import os
import numpy
import weighted
from numpy import *
from numpy import average
#The line below are for runs on systems without an output display system
import matplotlib
matplotlib.use('Agg')
from subprocess import Popen, PIPE
from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from fitting import * 
import collections
from scipy import spatial
from numpy import arccos
from numpy import arcsin
from numpy import arctan2
from numpy import cos
from numpy import sin
from os.path import isfile 
from MyFunction import *
from fitting_functions import *
import pyfits
#import pdb; pdb.set_trace()
###########################
#Command line inputs

#Directory where the result folder will be.
#Master catalogue should be in there.
Dir=sys.argv[1]
#Master file
file=str(sys.argv[2])
#Yes/No: if yes analysis in bins of strehl
sep_strehl_flag=str(sys.argv[3])
#Yes/No: do want to do the analysis in the PSF frame ? 
rotation_flag=str(sys.argv[4])
#Number of bins for the main parameter used to compute the bias. 
Nbin=int(sys.argv[5])
#Number of bins for the second parameter used to compute the bias. 
Nbin2=int(sys.argv[6])

#None: no calibration, just row analysy. This needs to be done first!
#After you have run the pipleine the first time yopu can use the output to
#compute residuals. You can then use one of the options below.
#
#Spline/Spline2/Multi_spline: different spline fits (not used in any publication)
#Bin: calibration type in IFC17. You need the surface of doom 
#File: reads the (m,c) calibration from the master file.
#Table: fine sampled calibration (100x100 or so). Calibration is based on
#a nearest neighbour search in the table.
cal_type=str(sys.argv[7])

#Which quantity to use to define the calibration
#Fiducial: snr_res_alt_mv (SNR+R calibration)
#Other options below
cal_select=str(sys.argv[8])

#Auxilliary quantity.
#If cal_type=Bin/Table --> m_func is the path to the file containing the
#surface of doom
#If cal_type=file --> column base name of the calibration (e.g. mycalibration) 
m_func=sys.argv[9]

#Set it to an empty string (you do not need it)
c_func=sys.argv[10]

#Integer ID you can give to the run (it just write the identifier to the output
#file)
Run_id= sys.argv[11]

##HISTORIC REASONS:
Psf_set = Run_id

##HISTORIC REASONS:
Flag="MV"


##HISTORIC REASONS:

file_ending=""

if(not(os.path.isdir(Dir+"./Results"))):
    os.mkdir(Dir+"./Results")

if(not(os.path.isdir(Dir+"./Results/1bin"))):
    os.mkdir(Dir+"./Results/1bin")

if(not(os.path.isdir(Dir+"./Results/2bin"))):
    os.mkdir(Dir+"./Results/2bin")


##ALL FUNCTIONALITY REGARDING CALIBRATION WITH FUNCTIONS IS CURRENTLY COMMENTED OUT, I LEFT IT IN FOR LATER REVIVAL

#m1_params=(float(sys.argv[11]),float(sys.argv[12]),float(sys.argv[13]),float(sys.argv[14]))
#m2_params=(float(sys.argv[15]),float(sys.argv[16]),float(sys.argv[17]),float(sys.argv[18]))
#c1_params=(float(sys.argv[19]),float(sys.argv[20]),float(sys.argv[21]),float(sys.argv[22]))
#c2_params=(float(sys.argv[22]),float(sys.argv[23]),float(sys.argv[24]),float(sys.argv[25]))


####INITIAL PARAMETER CHECKS, DEFINES FUNCTIONALITY LATER ON##########


#function=False
spline=False
spline2=False
spline_multi=False
table_cal = False
bin_cal = False
file_cal = False

if(cal_type=='spline'):
    spline=True
elif(cal_type=='spline2'):
    spline=True
    spline2=True
elif(cal_type=='multi_spline'):
    spline=True
    spline_multi=True
#elif(cal_type=='func'):
#    function=True
elif(cal_type=='table'):
    table_cal=True
elif(cal_type=='bin'):
    bin_cal=True
elif(cal_type=='file'):
    file_cal=True

#if(function==True):
#    func_m=function_set[str(m_func)]
#    func_c=function_set[str(c_func)]


alt_file_ending=""


if(rotation_flag=='True' or rotation_flag=='true' or rotation_flag=='Yes' or rotation_flag=='yes' or rotation_flag=='Y' or rotation_flag=='y' or  rotation_flag=='1'):
    psf_rotation=True
    alt_file_ending = alt_file_ending+"_rot"
else:
    psf_rotation=False
    alt_file_ending = alt_file_ending

##Trying to read FITS master

if(os.path.exists(file)):
    print "Master FITS binary exists, reading..."
    fits_data = pyfits.open(file)
    #TB_sim_mask=fits_data[1].data['TB9_in']
    catID_sim_mask=fits_data[1].data['Cat_ID']
    w_sim_mask=fits_data[1].data['LFweight']
    e1_sim_mask=fits_data[1].data['e1']
    e2_sim_mask=fits_data[1].data['e2']
    mag_sim_weight=fits_data[1].data['mag_out']
    g1_in_mask=fits_data[1].data['g1']
    g2_in_mask=fits_data[1].data['g2']
    strehl_sim_weight=fits_data[1].data['strehl']
    snr_sim_weight=fits_data[1].data['snr_model']
    size_sim_weight=fits_data[1].data['size_out']
    mag_in_sim_weight=fits_data[1].data['mag_in']
    size_in_sim_weight=fits_data[1].data['size_in']
    e1_in_sim_weight=fits_data[1].data['e1_in']
    e2_in_sim_weight=fits_data[1].data['e2_in']
    psf_size_sim_weight=fits_data[1].data['psf_size_in']
    psf_e1_sim_weight=fits_data[1].data['psf_e1_in']
    psf_e2_sim_weight=fits_data[1].data['psf_e2_in']
    r_corr_sim_weight=fits_data[1].data['size_corr']
    e1_corr_sim_weight=fits_data[1].data['e1_corr']
    e2_corr_sim_weight=fits_data[1].data['e2_corr']
    nd_sim_weight=fits_data[1].data['nd']
    ls_var_sim_weight=fits_data[1].data['LS variance']
    fitClass_sim_weight=fits_data[1].data['fitclass']
    contamination_radius_sim=fits_data[1].data['contamination_radius']
    if(file_cal):
        m1_sim_weight=fits_data[1].data['m1_'+m_func]
        m2_sim_weight=fits_data[1].data['m2_'+m_func]
        c1_sim_weight=fits_data[1].data['c1_'+m_func]
        c2_sim_weight=fits_data[1].data['c2_'+m_func]

    del fits_data
    print "FITS data read...proceeding with analysis."
    
else:
    print 'Given FITS file does not exist. Exiting'
    sys.exit()

size_raw_sim = size_sim_weight-r_corr_sim_weight


###############################################################################################################
#Masking objects with w==0 and fitclass.neq.-6, -1 or -10 AKJ: or 1 or 2

#ask_weight=(w_sim_mask>0)&(fitClass_sim_weight!=-6)&(fitClass_sim_weight!=-1)&(fitClass_sim_weight!=-10)&(size_raw_sim > 0.5)&(fitClass_sim_weight!=1)&(fitClass_sim_weight!=2)
mask_weight = (np.hypot(e1_sim_mask,e2_sim_mask)<=0.8)|(size_sim_weight>=(0.5 * exp(0.65788*(24.2 - mag_sim_weight))))
mask_weight &=(w_sim_mask>0)&(fitClass_sim_weight!=-1)&(fitClass_sim_weight!=-10)&(size_raw_sim > 0.5)&(fitClass_sim_weight!=1)&(fitClass_sim_weight!=2)&(contamination_radius_sim>4.25)
#mask_weight &= (TB_sim_mask>1.9) ## Benjamin Joachimi's blue/red split

## Impose a stricter star cut
#strict_star_cuts = True
#star_objnos = np.loadtxt('/disks/shear15/KiDS/ImSim/pipeline/data/smallgals.asc')
#for star_objno in star_objnos:
#    strict_star_cuts &= (catID_sim_mask==star_objno)
strict_star_cuts = (size_in_sim_weight>0.1-0.01*(mag_in_sim_weight-16.)) | (size_in_sim_weight==0) ## Filter out stars in HST masquerading as galaxies
#mask_weight &= (strict_star_cuts)

w_sim_mask=w_sim_mask[mask_weight]
e1_sim_mask=e1_sim_mask[mask_weight]
e2_sim_mask=e2_sim_mask[mask_weight]
mag_sim_weight=mag_sim_weight[mask_weight]
g1_in_mask=g1_in_mask[mask_weight]
g2_in_mask=g2_in_mask[mask_weight]
strehl_sim_weight=strehl_sim_weight[mask_weight]
snr_sim_weight=snr_sim_weight[mask_weight]
size_sim_weight=size_sim_weight[mask_weight]
mag_in_sim_weight=mag_in_sim_weight[mask_weight]
size_in_sim_weight=size_in_sim_weight[mask_weight]
e1_in_sim_weight=e1_in_sim_weight[mask_weight]
e2_in_sim_weight=e2_in_sim_weight[mask_weight]
psf_size_sim_weight=psf_size_sim_weight[mask_weight]
psf_e1_sim_weight=psf_e1_sim_weight[mask_weight]
psf_e2_sim_weight=psf_e2_sim_weight[mask_weight]
r_corr_sim_weight=r_corr_sim_weight[mask_weight]
e1_corr_sim_weight=e1_corr_sim_weight[mask_weight]
e2_corr_sim_weight=e2_corr_sim_weight[mask_weight]
nd_sim_weight=nd_sim_weight[mask_weight]
ls_var_sim_weight=ls_var_sim_weight[mask_weight]



if(psf_rotation==True):
    g1_in_mask, g2_in_mask=rotation(g1_in_mask,g2_in_mask,psf_e1_sim_weight,psf_e2_sim_weight)
    e1_sim_mask, e2_sim_mask=rotation(e1_sim_mask,e2_sim_mask,psf_e1_sim_weight,psf_e2_sim_weight)
    e1_in_sim_weight, e2_in_sim_weight=rotation(e1_in_sim_weight,e2_in_sim_weight,psf_e1_sim_weight,psf_e2_sim_weight)
    e1_corr_sim_weight, e2_corr_sim_weight=rotation(e1_corr_sim_weight,e2_corr_sim_weight,psf_e1_sim_weight,psf_e2_sim_weight)



#Calculating the true ellipticity of objects

##THIS IS A RELIC FROM A FORMER ANALYSIS, I LEAVE IT IN, BUT COMMENTED IF NEEDED LATER

#print "Calculating true ellipticities."

#e1_true=[]
#e2_true=[]
#for i in range(0,len(e1_sim_mask)):
#    g = complex(g1_in_mask[i],g2_in_mask[i])
#    e_int = complex(e1_in_sim_weight[i],e2_in_sim_weight[i])
#    e_true = (e_int+g)/(1.+g.conjugate()*e_int)
#    e1_true.append(e_true.real)
#    e2_true.append(e_true.imag)
#    #w_sim_mask[i] = 1.
#e1_true=array(e1_true)
#e2_true=array(e2_true)
#print "Done."


####FROM HERE ON IS ACTUAL ANALYSIS, NOT FILE_READING AND SELECTION######

#Computing the total weight

TotWeight=sum(w_sim_mask)

#Minimum number of shear values to perform the linear regression (this is relevant for psf properties)

NumShearMin=5.0

###############################################################################################################

#ObsType=["mag","mag_out", "size_in", "size_in_alt", "snr", "size_out","size_out_raw","size_out_alt","size_out_alt_raw","resolution", "resolution_alt", "nn", "lsv","resolution_alt_raw"]
#ObsType=["mag_out"]

#ObsType=["resolution_alt", "resolution_alt_raw"]

#ObsType=["size_out","size_out_raw","size_out_alt","size_out_alt_raw","resolution_alt_mv", "nn"]

#ObsType=["mag","size_in_alt","snr","resolution_alt_mv"]

#ObsType=["psf_e1","psf_e2","psf_size","psf_strehl"] 

#ObsType=["size_in","size_in_alt","size_in_alt_sheared","size_out_raw","size_out_alt_raw"]

ObsType=["snr","size_out","resolution_alt_mv","size_out_alt","size_out_alt_raw"]

#ObsType=["snr"]


#Calculating some helper quantities, needed later on
e_in=sqrt(e1_in_sim_weight*e1_in_sim_weight + e2_in_sim_weight*e2_in_sim_weight)
#e_in_sheared=sqrt(e1_true*e1_true + e2_true*e2_true)
e_out=sqrt(e1_sim_mask*e1_sim_mask + e2_sim_mask*e2_sim_mask)
e_out_raw=sqrt((e1_sim_mask-e1_corr_sim_weight)**2 + (e2_sim_mask-e2_corr_sim_weight)**2)
size_sim_alt_weight=size_sim_weight*sqrt((1-e_out)/(1+e_out))
resolution=psf_size_sim_weight/(size_sim_weight**2)
resolution_alt=psf_size_sim_weight/(size_sim_alt_weight**2)
size_in_sim_alt_weight=size_in_sim_weight*sqrt((1-e_in)/(1+e_in))
#size_in_sim_alt_sheared_weight=size_in_sim_weight*sqrt((1-e_in_sheared)/(1+e_in_sheared))
size_sim_alt_weight_raw=(size_sim_weight-r_corr_sim_weight)*sqrt((1-e_out_raw)/(1+e_out_raw))
resolution_alt_raw=psf_size_sim_weight/(size_sim_alt_weight_raw**2)
resolution_alt_mv=psf_size_sim_weight/((size_sim_alt_weight**2)+psf_size_sim_weight)



if(cal_select=='mag'):
    cal1=mag_in_sim_weight
    cal2=mag_in_sim_weight
    spline2=False
    table_cal=False
    bin_cal=False
elif(cal_select=='mag_out'):
    cal1=mag_sim_weight
    cal2=mag_sim_weight
    spline2=False
    table_cal=False
    bin_cal = False
elif(cal_select=='size_alt'):
    cal1=size_in_sim_weight*sqrt((1-e_in)/(1+e_in))
    cal2=cal1
    spline2=False
    table_cal=False
    bin_cal = False
elif(cal_select=='snr'):
    cal1=snr_sim_weight
    cal2=snr_sim_weight
    spline2=False
    table_cal=False
    bin_cal= False
elif(cal_select=='res'):
    cal1=resolution
    cal2=resolution
    spline2=False
    table_cal=False
    bin_cal = False
elif(cal_select=='res_alt'):
    cal1=resolution_alt
    cal2=resolution_alt
    spline2=False
    table_cal=False
    bin_cal = False
elif(cal_select=='snr_res'):
    cal1=snr_sim_weight
    cal2=resolution
elif(cal_select=='snr_res_alt'):
    cal1=snr_sim_weight
    cal2=resolution_alt
elif(cal_select=='snr_res_alt_mv'):
    cal1=snr_sim_weight
    cal2=resolution_alt_mv
elif(cal_select=='snr_res_alt_raw'):
    cal1=snr_sim_weight
    cal2=resolution_alt_raw
elif(cal_select=='snr_size_alt'):
    cal1=snr_sim_weight
    cal2=size_sim_alt_weight

elif(cal_select=='mag_size_alt'):
    cal1=mag_in_sim_weight
    cal2=size_in_sim_weight*sqrt((1-e_in)/(1+e_in))
elif(cal_select=='mag_res_alt'):
    cal1=mag_sim_weight
    cal2=resolution_alt

else:
#    function=False
    spline=False
    spline2=False
    table_cal=False
    bin_cal = False

#if(function==True):
#    alt_file_ending=alt_file_ending+"func_cal"


if(spline==True):
    m1_splines=[]
    m2_splines=[]
    c1_splines=[]
    c2_splines=[]
    if(spline2==True):
        if(m_func =='none'):
            string_start=Dir+'/Results/1bin/'+Flag+"_"+str(Psf_set)+'_100_'
            if(cal_select=='snr_res'):
                string_mid="SignalToNoise"
                string_mid2="Resolution"
            else:
                string_mid="SignalToNoise"
                string_mid2="ResolutionAlt"
                
            string_end='_binning_global'
            spline_data=genfromtxt(string_start+string_mid+string_end+alt_file_ending+".txt",comments='#')
            spline_data2=genfromtxt(string_start+string_mid2+string_end+alt_file_ending+".txt",comments='#')
        else:
            spline_data=genfromtxt(Dir+'/'+str(m_func),comments='#')
            spline_data2=genfromtxt(Dir+'/'+str(c_func),comments='#')
        
        spline_x=spline_data[:,0]
        spline_m1=spline_data[:,1]
        spline_c1=spline_data[:,3]
        spline_m2=spline_data[:,5]
        spline_c2=spline_data[:,7]
        spline2_x=spline_data2[:,0]
        spline2_m1=spline_data2[:,1]
        spline2_c1=spline_data2[:,3]
        spline2_m2=spline_data2[:,5]
        spline2_c2=spline_data2[:,7]
        alt_file_ending=alt_file_ending+"_2spline_cal"
        cal_spline_m1=interpolate.interp1d(spline_x,spline_m1,kind='slinear',bounds_error=False)
        cal_spline_m2=interpolate.interp1d(spline_x,spline_m2,kind='slinear',bounds_error=False)
        cal_spline_c1=interpolate.interp1d(spline_x,spline_c1,kind='slinear',bounds_error=False)
        cal_spline_c2=interpolate.interp1d(spline_x,spline_c2,kind='slinear',bounds_error=False)
        cal_spline2_m1=interpolate.interp1d(spline2_x,spline2_m1,kind='slinear',bounds_error=False)
        cal_spline2_m2=interpolate.interp1d(spline2_x,spline2_m2,kind='slinear',bounds_error=False)
        cal_spline2_c1=interpolate.interp1d(spline2_x,spline2_c1,kind='slinear',bounds_error=False)
        cal_spline2_c2=interpolate.interp1d(spline2_x,spline2_c2,kind='slinear',bounds_error=False)
        m1_splines.append(cal_spline_m1)
        m1_splines.append(cal_spline2_m1)
        m2_splines.append(cal_spline_m2)
        m2_splines.append(cal_spline2_m2)
        c1_splines.append(cal_spline_c1)
        c1_splines.append(cal_spline2_c1)
        c2_splines.append(cal_spline_c2)
        c2_splines.append(cal_spline2_c2)
        m1_splines=array(m1_splines)
        m2_splines=array(m2_splines)
        c1_splines=array(c1_splines)
        c2_splines=array(c2_splines)

    elif(spline_multi==True):
        m1_splines, m2_splines, c1_splines, c2_splines, spline_thresholds = create_binned_splines(Dir+'/'+str(m_func),int(c_func))
        alt_file_ending=alt_file_ending+"_multi_spline_cal"
        #mag_data=genfromtxt(Dir+"/MV_23.0_100_MagnitudeOut_binning_globalscheme2b_corr_multi_spline_cal.txt",comments='#')
        #mag_x=mag_data[:,0]
        #mag_m1=mag_data[:,1]
        #mag_m2=mag_data[:,5]
        #mag_spline1=interpolate.interp1d(mag_x,mag_m1,kind='slinear',bounds_error=False)
        #mag_spline2=interpolate.interp1d(mag_x,mag_m2,kind='slinear',bounds_error=False)

    else:
        if(m_func =='none'):
            string_start=Dir+'/Results/1bin/'+Flag+"_"+str(Psf_set)+'_100_'
            if(cal_select=='mag'):
                string_mid="Magnitude"
            elif(cal_select=='snr'):
                string_mid="SignalToNoise"
            else:
                string_mid="inputSize"
                
                string_end='_binning_global'
                spline_data=genfromtxt(string_start+string_mid+string_end+alt_file_ending+".txt",comments='#')
        else:
            spline_data=genfromtxt(Dir+'/'+str(m_func),comments='#')
        
        spline_x=spline_data[:,0]
        spline_m1=spline_data[:,1]
        spline_c1=spline_data[:,3]
        spline_m2=spline_data[:,5]
        spline_c2=spline_data[:,7]
        alt_file_ending=alt_file_ending+"_spline_cal"
        cal_spline_m1=interpolate.interp1d(spline_x,spline_m1,kind='slinear',bounds_error=False)
        cal_spline_m2=interpolate.interp1d(spline_x,spline_m2,kind='slinear',bounds_error=False)
        cal_spline_c1=interpolate.interp1d(spline_x,spline_m1,kind='slinear',bounds_error=False)
        cal_spline_c2=interpolate.interp1d(spline_x,spline_m2,kind='slinear',bounds_error=False)
        m1_splines.append(cal_spline_m1)
        m2_splines.append(cal_spline_m2)
        c1_splines.append(cal_spline_c1)
        c2_splines.append(cal_spline_c2)
        m1_splines=array(m1_splines)
        m2_splines=array(m2_splines)
        c1_splines=array(c1_splines)
        c2_splines=array(c2_splines)

elif(bin_cal==True):
    alt_file_ending=alt_file_ending+"_bin_cal"
    table_data=genfromtxt(Dir+'/'+str(m_func),comments='#')

elif(table_cal==True):
    alt_file_ending=alt_file_ending+"_table_cal"
    table_data=genfromtxt(Dir+'/'+str(m_func),comments='#')
    x_table=table_data[:,0]
    y_table=table_data[:,3]
    m1_table=table_data[:,7]
    c1_table=table_data[:,9]
    m2_table=table_data[:,11]
    c2_table=table_data[:,13]

    shift_x = min(x_table)
    shift_y = min(y_table)
    x_table = x_table - shift_x
    y_table = y_table - shift_y
    scale_x = max(x_table)
    scale_y = max(y_table)

    x_table=x_table/scale_x
    y_table=y_table/scale_y

    tree_pos = numpy.dstack([x_table,y_table])[0]
    Tree=scipy.spatial.cKDTree(tree_pos,leafsize=128)

elif(file_cal==True):
    alt_file_ending=alt_file_ending+"_file_cal"

#Creating the unique shears in the total data set

g1_theo=[]
g2_theo=[]
shear_pairs=column_stack((g1_in_mask,g2_in_mask)) 
unique_indices=unique_rows(shear_pairs)
for index in range(0,len(unique_indices)):
    g1_theo.append(g1_in_mask[unique_indices[index]])
    g2_theo.append(g2_in_mask[unique_indices[index]])

g1_theo=array(g1_theo)
g2_theo=array(g2_theo)
    


for kk in range(len(ObsType)):
    print 'I am doing...', ObsType[kk]
    if(ObsType[kk]=="snr"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight=snr_sim_weight
        fileOut="SignalToNoise_binning"+psf_ending+alt_file_ending+".txt"
        fileOutAll="SignalToNoise_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotName="Bias_signalToNoise."+psf_ending+alt_file_ending+".png"
        PlotNameAll="Bias_signalToNoise_global"+psf_ending+alt_file_ending+".png"
        ObsLabel='snr'
    if(ObsType[kk]=="mag"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight=mag_in_sim_weight
        fileOut="Magnitude_binning"+psf_ending+alt_file_ending+".txt"
        fileOutAll="Magnitude_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotName="Bias_magnitude."+psf_ending+alt_file_ending+".png"
        PlotNameAll="Bias_magnitude_global"+psf_ending+alt_file_ending+".png"
        ObsLabel="mag"
    if(ObsType[kk]=="mag_out"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight=mag_sim_weight
        fileOut="MagnitudeOut_binning"+psf_ending+alt_file_ending+".txt"
        fileOutAll="MagnitudeOut_binning_global"+psf_ending+alt_file_ending+".txt\
"
        PlotName="Bias_magnitude_out."+psf_ending+alt_file_ending+".png"
        PlotNameAll="Bias_magnitude_out_global"+psf_ending+alt_file_ending+".png"
        ObsLabel="mag_out"
    if(ObsType[kk]=="size_in"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight=size_in_sim_weight
        fileOut="InputSize_binning"+psf_ending+alt_file_ending+".txt"
        fileOutAll="InputSize_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotName="Bias_inputSize."+psf_ending+alt_file_ending+".png"
        PlotNameAll="Bias_inputSize_global"+psf_ending+alt_file_ending+".png"
        ObsLabel='size (in)'
    if(ObsType[kk]=="size_in_alt"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight=size_in_sim_alt_weight
        fileOut="InputSizeAlt_binning"+psf_ending+alt_file_ending+".txt"
        fileOutAll="InputSizeAlt_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotName="Bias_inputSizeAlt."+psf_ending+alt_file_ending+".png"
        PlotNameAll="Bias_inputSizeAlt_global"+psf_ending+alt_file_ending+".png"
        ObsLabel='alt size (in)'
#    if(ObsType[kk]=="size_in_alt_sheared"):
#        if(psf_rotation==True):
#            ObsType[kk]=ObsType[kk]+"PsfRot"
#            psf_ending="PsfRot."
#        else:
#            psf_ending=""
#        obs_sim_weight=size_in_sim_alt_sheared_weight
#        fileOut="InputSizeAltSheared_binning"+psf_ending+alt_file_ending+".txt"
#        fileOutAll="InputSizeAltSheared_binning_global"+psf_ending+alt_file_ending+".txt"
#        PlotName="Bias_inputSizeAltSheared."+psf_ending+alt_file_ending+".png"
#        PlotNameAll="Bias_inputSizeAltSheared_global"+psf_ending+alt_file_ending+".png"
#        ObsLabel='alt size sheared (in)'

    if(ObsType[kk]=="size_out"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight=size_sim_weight
        fileOut="OutputSize_binning"+psf_ending+alt_file_ending+".txt"
        fileOutAll="OutputSize_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotName="Bias_outputSize."+psf_ending+alt_file_ending+".png"
        PlotNameAll="Bias_outputSize_global"+psf_ending+alt_file_ending+".png"
        ObsLabel='size (out)'
    if(ObsType[kk]=="size_out_raw"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight=size_sim_weight-r_corr_sim_weight
        fileOut="OutputSizeRaw_binning"+psf_ending+alt_file_ending+".txt"
        fileOutAll="OutputSizeRaw_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotName="Bias_outputSizeRaw."+psf_ending+alt_file_ending+".png"
        PlotNameAll="Bias_outputSizeRaw_global"+psf_ending+alt_file_ending+".png"
        ObsLabel='size (out no selfcal)'
    if(ObsType[kk]=="size_out_alt"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight=size_sim_alt_weight
        fileOut="OutputSizeAlt_binning"+psf_ending+alt_file_ending+".txt"
        fileOutAll="OutputSizeAlt_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotName="Bias_outputSizeAlt."+psf_ending+alt_file_ending+".png"
        PlotNameAll="Bias_outputSizeAlt_global"+psf_ending+alt_file_ending+".png"
        ObsLabel='size (sqrt(ab) selfcal)'
    if(ObsType[kk]=="size_out_alt_raw"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight=(size_sim_weight-r_corr_sim_weight)*sqrt((1-e_out_raw)/(1+e_out_raw))
        fileOut="OutputSizeAltRaw_binning"+psf_ending+alt_file_ending+".txt"
        fileOutAll="OutputSizeAltRaw_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotName="Bias_outputSizeAltRaw."+psf_ending+alt_file_ending+".png"
        PlotNameAll="Bias_outputSizeAltRaw_global"+psf_ending+alt_file_ending+".png"
        ObsLabel='size (sqrt(ab) no selfcal)'
    if(ObsType[kk]=="resolution"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight=resolution
        fileOut="Resolution_binning"+psf_ending+alt_file_ending+".txt"
        fileOutAll="Resolution_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotName="Bias_Resolution."+psf_ending+alt_file_ending+".png"
        PlotNameAll="Bias_Resolution_global"+psf_ending+alt_file_ending+".png"
        ObsLabel='R'
    if(ObsType[kk]=="resolution_alt"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight=resolution_alt
        fileOut="ResolutionAlt_binning"+psf_ending+alt_file_ending+".txt"
        fileOutAll="ResolutionAlt_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotName="Bias_ResolutionAlt."+psf_ending+alt_file_ending+".png"
        PlotNameAll="Bias_ResolutionAlt_global"+psf_ending+alt_file_ending+".png"
        ObsLabel='R (alt)'
    if(ObsType[kk]=="resolution_alt_mv"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight=resolution_alt_mv
        fileOut="ResolutionAltMV_binning"+psf_ending+alt_file_ending+".txt"
        fileOutAll="ResolutionAltMV_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotName="Bias_ResolutionAltMV."+psf_ending+alt_file_ending+".png"
        PlotNameAll="Bias_ResolutionAltMV_global"+psf_ending+alt_file_ending+".png"
        ObsLabel='R (alt MV)'
    if(ObsType[kk]=="resolution_alt_raw"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight=resolution_alt_raw
        fileOut="ResolutionAltRaw_binning"+psf_ending+alt_file_ending+".txt"
        fileOutAll="ResolutionAltRaw_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotName="Bias_ResolutionAltRaw."+psf_ending+alt_file_ending+".png"
        PlotNameAll="Bias_ResolutionAltRaw_global"+psf_ending+alt_file_ending+".png"
        ObsLabel='R (alt raw)'
    if(ObsType[kk]=="nn"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight=nd_sim_weight
        fileOut="NearestNeighbour_binning"+psf_ending+alt_file_ending+".txt"
        fileOutAll="NearestNeighbour_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotName="Bias_NearestNeighbour."+psf_ending+alt_file_ending+".png"
        PlotNameAll="Bias_NearestNeighbour_global"+psf_ending+alt_file_ending+".png"
        ObsLabel="nn"
    if(ObsType[kk]=="psf_e"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight=sqrt(psf_e1_sim_weight*psf_e1_sim_weight + psf_e2_sim_weight*psf_e2_sim_weight)
        fileOut="PsfE_binning"+psf_ending+alt_file_ending+".txt"
        fileOutAll="PsfE_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotName="Bias_PsfE."+psf_ending+alt_file_ending+".png"
        PlotNameAll="Bias_PsfE_global"+psf_ending+alt_file_ending+".png"
        ObsLabel="psf_e"
    if(ObsType[kk]=="psf_e1"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight=psf_e1_sim_weight
        fileOut="PsfE1_binning"+psf_ending+alt_file_ending+".txt"
        fileOutAll="PsfE1_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotName="Bias_PsfE1."+psf_ending+alt_file_ending+".png"
        PlotNameAll="Bias_PsfE1_global"+psf_ending+alt_file_ending+".png"
        ObsLabel="psf_e1"
    if(ObsType[kk]=="psf_e2"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight=psf_e2_sim_weight
        fileOut="PsfE2_binning"+psf_ending+alt_file_ending+".txt"
        fileOutAll="PsfE2_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotName="Bias_PsfE2."+psf_ending+alt_file_ending+".png"
        PlotNameAll="Bias_PsfE2_global"+psf_ending+alt_file_ending+".png"
        ObsLabel="psf_e2"
    if(ObsType[kk]=="psf_size"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight=psf_size_sim_weight
        fileOut="PsfSize_binning"+psf_ending+alt_file_ending+".txt"
        fileOutAll="PsfSize_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotName="Bias_PsfSize."+psf_ending+alt_file_ending+".png"
        PlotNameAll="Bias_PsfSize_global"+psf_ending+alt_file_ending+".png"
        ObsLabel="psf_size"
    if(ObsType[kk]=="psf_strehl"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight=strehl_sim_weight
        fileOut="PsfStrehl_binning"+psf_ending+alt_file_ending+".txt"
        fileOutAll="PsfStrehl_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotName="Bias_PsfStrehl."+psf_ending+alt_file_ending+".png"
        PlotNameAll="Bias_PsfStrehl_global"+psf_ending+alt_file_ending+".png"
        ObsLabel="psf_strehl"

    if(ObsType[kk]=="lsv"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight=ls_var_sim_weight
        fileOut="LsVar_binning"+psf_ending+alt_file_ending+".txt"
        fileOutAll="LsVar_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotName="Bias_Ls_Var."+psf_ending+alt_file_ending+".png"
        PlotNameAll="Bias_LsVar_global"+psf_ending+alt_file_ending+".png"
        ObsLabel="ls_var"

    # This switch is only relevant if you want to analyse true shears instead of measured ones
#    e1_sim_mask = e1_true
#    e2_sim_mask = e2_true

        
    ObsBin=ComputeBins_1D_weighted(obs_sim_weight, Nbin, w_sim_mask)
    ObsBin=unique(ObsBin)
    print ObsType[kk],ObsBin
    if(len(ObsBin)>1):
    #if(1==2):
        #Run with Nbins bins. Each psf-set is analysed globally. 
        if(spline==True):
            if(spline2==True):
                DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100),-99,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines)
            elif(spline_multi==True):
                #DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100),-99,file_ending,cal1,cal2,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines, spline_thresholds,mag_sim_weight,mag_spline1,mag_spline2)
                DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100),-99,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines, spline_thresholds)
            else:
                DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100),-99,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines)
#        elif(function==True):
#            DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100),-99,file_ending,cal1,cal2,None,None,None,None,func_m,m1_params,m2_params,func_c,c1_params,c2_params)
        elif(table_cal==True):
            print "Doing this."
            DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100),-99,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,m1_table,m2_table,c1_table,c2_table,shift_x,shift_y,scale_x, scale_y,Tree)
        elif(bin_cal==True):
            print "Doing this."
            DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100),-99,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,table_data)
        elif(file_cal==True):
            DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100),-99,file_ending,cal1,cal2,m1_sim_weight,m2_sim_weight,c1_sim_weight,c2_sim_weight,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None)
        else:
            DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100),-99,file_ending)
            

            
    else:
        continue

    if(sep_strehl_flag=='Yes' or sep_strehl_flag=='Y' or sep_strehl_flag=='True' or  sep_strehl_flag=='yes' or sep_strehl_flag=='y' or sep_strehl_flag=='true'):
        #Run with Nbins bins. Each psf-set is further splits in 5 as the average PSF properties are different for galaxies having different number of exposures. The binning is done as a function of strehl ratio. Note that defining the binning using the quantiles doesn't work optimally here as the distribution is not smooth. After computing the quantiles I am filtering out the non-unique values and use the resulting array to define the bins. This can be improved.

        Psf_bins=ComputeBins_1D_weighted(strehl_sim_weight,10, w_sim_mask)    
        Psf_bins=unique(Psf_bins)
        print Psf_bins
        Nunique=len(Psf_bins)
        print 'Unique psf..', Nunique

        for gg in range(0, Nunique-1):
            name_add = str(gg) + "_strehl"
            mask_gg=(strehl_sim_weight>=Psf_bins[gg])&(strehl_sim_weight<Psf_bins[gg+1])
            strehl_set=numpy.median(strehl_sim_weight[mask_gg])
            ObsBin_set=ComputeBins_1D_weighted(obs_sim_weight[mask_gg], Nbin, w_sim_mask[mask_gg])
            ObsBin_set=unique(ObsBin_set)
            print Psf_bins[gg], Psf_bins[gg+1], ObsBin_set
            if(len(ObsBin_set)>1):
                #DO ALL MUST BE ADJUESTED TO TAKE CALIBRATION QUANTITIES, CALIBRATION FUNCTION AND ITS PARAMETERS.  
                if(spline==True):
                    if(spline2==True):
                        DoAll(ObsType[kk],ObsBin_set, obs_sim_weight[mask_gg], w_sim_mask[mask_gg], e1_sim_mask[maskgg], e2_sim_mask[mask_gg], e1_in_sim_weight[mask_gg], e2_in_sim_weight[mask_gg], e1_corr_sim_weight[mask_gg],e2_corr_sim_weight[mask_gg], g1_in_mask[mask_gg],g2_in_mask[mask_gg], g1_theo, g2_theo, fileOut, Flag, PlotName, ObsLabel,psf_e1_sim_weight[mask_gg],psf_e2_sim_weight[mask_gg],psf_size_sim_weight[mask_gg],Dir,TotWeight,Psf_set, name_add,strehl_set,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines)
                    elif(spline_multi==True):
                        #DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100),-99,file_ending,cal1,cal2,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines, spline_thresholds,mag_sim_weight,mag_spline1,mag_spline2)
                        DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines, spline_thresholds)
                    else:
                        DoAll(ObsType[kk],ObsBin_set, obs_sim_weight[mask_gg], w_sim_mask[mask_gg], e1_sim_mask[mask_gg], e2_sim_mask[mask_gg], e1_in_sim_weight[mask_gg], e2_in_sim_weight[mask_gg], e1_corr_sim_weight[mask_gg],e2_corr_sim_weight[mask_gg], g1_in_mask[mask_gg],g2_in_mask[mask_gg], g1_theo, g2_theo, fileOut, Flag, PlotName, ObsLabel,psf_e1_sim_weight[mask_gg],psf_e2_sim_weight[mask_gg],psf_size_sim_weight[mask_gg],Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines)
#                elif(function==True):
#                    DoAll(ObsType[kk],ObsBin_set, obs_sim_weight[mask_gg], w_sim_mask[mask_gg], e1_sim_mask[mask_gg], e2_sim_mask[mask_gg], e1_in_sim_weight[mask_gg], e2_in_sim_weight[mask_gg], e1_corr_sim_weight[mask_gg],e2_corr_sim_weight[mask_gg], g1_in_mask[mask_gg],g2_in_mask[mask_gg], g1_theo, g2_theo, fileOut, Flag, PlotName, ObsLabel,psf_e1_sim_weight[mask_gg],psf_e2_sim_weight[mask_gg],psf_size_sim_weight[mask_gg],Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending,cal1,cal2,None,None,None,None,func_m,m1_params,m2_params,func_c,c1_params,c2_params)
                elif(table_cal==True):
                    DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,m1_table,m2_table,c1_table,c2_table,shift_x,shift_y,scale_x, scale_y,Tree)
                elif(bin_cal==True):
                    DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,table_data)
                elif(file_cal==True):
                    DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending,cal1,cal2,m1_sim_weight,m2_sim_weight,c1_sim_weight,c2_sim_weight,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None)
                else:
                    DoAll(ObsType[kk],ObsBin_set, obs_sim_weight[mask_gg], w_sim_mask[mask_gg], e1_sim_mask[mask_gg], e2_sim_mask[mask_gg], e1_in_sim_weight[mask_gg], e2_in_sim_weight[mask_gg], e1_corr_sim_weight[mask_gg],e2_corr_sim_weight[mask_gg], g1_in_mask[mask_gg],g2_in_mask[mask_gg], g1_theo, g2_theo, fileOut, Flag, PlotName, ObsLabel,psf_e1_sim_weight[mask_gg],psf_e2_sim_weight[mask_gg],psf_size_sim_weight[mask_gg],Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending)

            else:
                continue


        Psf_bins=ComputeBins_1D_weighted(psf_e1_sim_weight,10, w_sim_mask)    
        Psf_bins=unique(Psf_bins)
        print Psf_bins
        Nunique=len(Psf_bins)
        print 'Unique psf..', Nunique

        for gg in range(0, Nunique-1):
            name_add = str(gg) + "_e1"
            mask_gg=(psf_e1_sim_weight>=Psf_bins[gg])&(psf_e1_sim_weight<Psf_bins[gg+1])
            strehl_set=numpy.median(psf_e1_sim_weight[mask_gg])
            ObsBin_set=ComputeBins_1D_weighted(obs_sim_weight[mask_gg], Nbin, w_sim_mask[mask_gg])
            ObsBin_set=unique(ObsBin_set)
            print Psf_bins[gg], Psf_bins[gg+1], ObsBin_set
            if(len(ObsBin_set)>1):
                #DO ALL MUST BE ADJUESTED TO TAKE CALIBRATION QUANTITIES, CALIBRATION FUNCTION AND ITS PARAMETERS.  
                if(spline==True):
                    if(spline2==True):
                        DoAll(ObsType[kk],ObsBin_set, obs_sim_weight[mask_gg], w_sim_mask[mask_gg], e1_sim_mask[maskgg], e2_sim_mask[mask_gg], e1_in_sim_weight[mask_gg], e2_in_sim_weight[mask_gg], e1_corr_sim_weight[mask_gg],e2_corr_sim_weight[mask_gg], g1_in_mask[mask_gg],g2_in_mask[mask_gg], g1_theo, g2_theo, fileOut, Flag, PlotName, ObsLabel,psf_e1_sim_weight[mask_gg],psf_e2_sim_weight[mask_gg],psf_size_sim_weight[mask_gg],Dir,TotWeight,Psf_set, name_add,strehl_set,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines)
                    elif(spline_multi==True):
                        #DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100),-99,file_ending,cal1,cal2,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines, spline_thresholds,mag_sim_weight,mag_spline1,mag_spline2)
                        DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines, spline_thresholds)
                    else:
                        DoAll(ObsType[kk],ObsBin_set, obs_sim_weight[mask_gg], w_sim_mask[mask_gg], e1_sim_mask[mask_gg], e2_sim_mask[mask_gg], e1_in_sim_weight[mask_gg], e2_in_sim_weight[mask_gg], e1_corr_sim_weight[mask_gg],e2_corr_sim_weight[mask_gg], g1_in_mask[mask_gg],g2_in_mask[mask_gg], g1_theo, g2_theo, fileOut, Flag, PlotName, ObsLabel,psf_e1_sim_weight[mask_gg],psf_e2_sim_weight[mask_gg],psf_size_sim_weight[mask_gg],Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines)
#                elif(function==True):
#                    DoAll(ObsType[kk],ObsBin_set, obs_sim_weight[mask_gg], w_sim_mask[mask_gg], e1_sim_mask[mask_gg], e2_sim_mask[mask_gg], e1_in_sim_weight[mask_gg], e2_in_sim_weight[mask_gg], e1_corr_sim_weight[mask_gg],e2_corr_sim_weight[mask_gg], g1_in_mask[mask_gg],g2_in_mask[mask_gg], g1_theo, g2_theo, fileOut, Flag, PlotName, ObsLabel,psf_e1_sim_weight[mask_gg],psf_e2_sim_weight[mask_gg],psf_size_sim_weight[mask_gg],Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending,cal1,cal2,None,None,None,None,func_m,m1_params,m2_params,func_c,c1_params,c2_params)
                elif(table_cal==True):
                    DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,m1_table,m2_table,c1_table,c2_table,shift_x,shift_y,scale_x, scale_y,Tree)
                elif(bin_cal==True):
                    DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,table_data)
                elif(file_cal==True):
                    DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending,cal1,cal2,m1_sim_weight,m2_sim_weight,c1_sim_weight,c2_sim_weight,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None)
                else:
                    DoAll(ObsType[kk],ObsBin_set, obs_sim_weight[mask_gg], w_sim_mask[mask_gg], e1_sim_mask[mask_gg], e2_sim_mask[mask_gg], e1_in_sim_weight[mask_gg], e2_in_sim_weight[mask_gg], e1_corr_sim_weight[mask_gg],e2_corr_sim_weight[mask_gg], g1_in_mask[mask_gg],g2_in_mask[mask_gg], g1_theo, g2_theo, fileOut, Flag, PlotName, ObsLabel,psf_e1_sim_weight[mask_gg],psf_e2_sim_weight[mask_gg],psf_size_sim_weight[mask_gg],Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending)

        Psf_bins=ComputeBins_1D_weighted(psf_e2_sim_weight,10, w_sim_mask)    
        Psf_bins=unique(Psf_bins)
        print Psf_bins
        Nunique=len(Psf_bins)
        print 'Unique psf..', Nunique

        for gg in range(0, Nunique-1):
            name_add = str(gg) + "_e2"
            mask_gg=(psf_e2_sim_weight>=Psf_bins[gg])&(psf_e2_sim_weight<Psf_bins[gg+1])
            strehl_set=numpy.median(psf_e2_sim_weight[mask_gg])
            ObsBin_set=ComputeBins_1D_weighted(obs_sim_weight[mask_gg], Nbin, w_sim_mask[mask_gg])
            ObsBin_set=unique(ObsBin_set)
            print Psf_bins[gg], Psf_bins[gg+1], ObsBin_set
            if(len(ObsBin_set)>1):
                #DO ALL MUST BE ADJUESTED TO TAKE CALIBRATION QUANTITIES, CALIBRATION FUNCTION AND ITS PARAMETERS.  
                if(spline==True):
                    if(spline2==True):
                        DoAll(ObsType[kk],ObsBin_set, obs_sim_weight[mask_gg], w_sim_mask[mask_gg], e1_sim_mask[maskgg], e2_sim_mask[mask_gg], e1_in_sim_weight[mask_gg], e2_in_sim_weight[mask_gg], e1_corr_sim_weight[mask_gg],e2_corr_sim_weight[mask_gg], g1_in_mask[mask_gg],g2_in_mask[mask_gg], g1_theo, g2_theo, fileOut, Flag, PlotName, ObsLabel,psf_e1_sim_weight[mask_gg],psf_e2_sim_weight[mask_gg],psf_size_sim_weight[mask_gg],Dir,TotWeight,Psf_set, name_add,strehl_set,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines)
                    elif(spline_multi==True):
                        #DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100),-99,file_ending,cal1,cal2,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines, spline_thresholds,mag_sim_weight,mag_spline1,mag_spline2)
                        DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines, spline_thresholds)
                    else:
                        DoAll(ObsType[kk],ObsBin_set, obs_sim_weight[mask_gg], w_sim_mask[mask_gg], e1_sim_mask[mask_gg], e2_sim_mask[mask_gg], e1_in_sim_weight[mask_gg], e2_in_sim_weight[mask_gg], e1_corr_sim_weight[mask_gg],e2_corr_sim_weight[mask_gg], g1_in_mask[mask_gg],g2_in_mask[mask_gg], g1_theo, g2_theo, fileOut, Flag, PlotName, ObsLabel,psf_e1_sim_weight[mask_gg],psf_e2_sim_weight[mask_gg],psf_size_sim_weight[mask_gg],Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines)
#                elif(function==True):
#                    DoAll(ObsType[kk],ObsBin_set, obs_sim_weight[mask_gg], w_sim_mask[mask_gg], e1_sim_mask[mask_gg], e2_sim_mask[mask_gg], e1_in_sim_weight[mask_gg], e2_in_sim_weight[mask_gg], e1_corr_sim_weight[mask_gg],e2_corr_sim_weight[mask_gg], g1_in_mask[mask_gg],g2_in_mask[mask_gg], g1_theo, g2_theo, fileOut, Flag, PlotName, ObsLabel,psf_e1_sim_weight[mask_gg],psf_e2_sim_weight[mask_gg],psf_size_sim_weight[mask_gg],Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending,cal1,cal2,None,None,None,None,func_m,m1_params,m2_params,func_c,c1_params,c2_params)
                elif(table_cal==True):
                    DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,m1_table,m2_table,c1_table,c2_table,shift_x,shift_y,scale_x, scale_y,Tree)
                elif(bin_cal==True):
                    DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,table_data)
                elif(file_cal==True):
                    DoAll(ObsType[kk],ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending,cal1,cal2,m1_sim_weight,m2_sim_weight,c1_sim_weight,c2_sim_weight,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None)
                else:
                    DoAll(ObsType[kk],ObsBin_set, obs_sim_weight[mask_gg], w_sim_mask[mask_gg], e1_sim_mask[mask_gg], e2_sim_mask[mask_gg], e1_in_sim_weight[mask_gg], e2_in_sim_weight[mask_gg], e1_corr_sim_weight[mask_gg],e2_corr_sim_weight[mask_gg], g1_in_mask[mask_gg],g2_in_mask[mask_gg], g1_theo, g2_theo, fileOut, Flag, PlotName, ObsLabel,psf_e1_sim_weight[mask_gg],psf_e2_sim_weight[mask_gg],psf_size_sim_weight[mask_gg],Dir,TotWeight,Psf_set,name_add,strehl_set,file_ending)

            else:
                continue


        
###############################################################################################################
#Here I am doing the 2D binning


#ObsType=["observed","observed_alt", "true", "mag_nn", "size_nn"]
ObsType=["observed_alt_mv"]
for kk in range(len(ObsType)):
    print 'I am doing...', ObsType[kk]
    if(ObsType[kk]=="observed"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight_1=snr_sim_weight
        obs_sim_weight_2=resolution
        fileOutAll="SignalToNoise_Resolution_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotNameAll="Bias_signaToNoise_Resolution_global"+psf_ending+alt_file_ending+".png"
        ObsLabel_1='r$\mathrm{S/N}$'
        ObsLabel_2='r$\mathrm{resolution (a)}$'
        namePlotBinning="Binning_resolution_sn."+psf_ending+alt_file_ending+".png"
    if(ObsType[kk]=="observed_alt"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight_1=snr_sim_weight
        obs_sim_weight_2=resolution_alt
        fileOutAll="SignalToNoise_ResolutionAlt_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotNameAll="Bias_signaToNoise_ResolutionAlt_global"+psf_ending+alt_file_ending+".png"
        ObsLabel_1='r$\mathrm{S/N}$'
        ObsLabel_2='r$\mathrm{resolution (sqrt(ab))}$'
        namePlotBinning="Binning_resolution_alt__sn."+psf_ending+alt_file_ending+".png"
    if(ObsType[kk]=="observed_alt_mv"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight_1=snr_sim_weight
        obs_sim_weight_2=resolution_alt_mv
        fileOutAll="SignalToNoise_ResolutionAltMV_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotNameAll="Bias_signaToNoise_ResolutionAltMV_global"+psf_ending+alt_file_ending+".png"
        ObsLabel_1='r$\mathrm{S/N}$'
        ObsLabel_2='r$\mathrm{resolution (sqrt(ab) MV)}$'
        namePlotBinning="Binning_resolution_alt_sn."+psf_ending+alt_file_ending+".png"
    if(ObsType[kk]=="observed_alt_raw"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight_1=snr_sim_weight
        obs_sim_weight_2=resolution_alt_raw
        fileOutAll="SignalToNoise_ResolutionAltRaw_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotNameAll="Bias_signaToNoise_ResolutionAltRaw_global"+psf_ending+alt_file_ending+".png"
        ObsLabel_1='r$\mathrm{S/N}$'
        ObsLabel_2='r$\mathrm{resolution (sqrt(ab))}$ (uncorr)'
        namePlotBinning="Binning_resolution_alt__sn."+psf_ending+alt_file_ending+".png"
    if(ObsType[kk]=="observed_size_alt_raw"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight_1=snr_sim_weight
        obs_sim_weight_2=size_sim_alt_weight
        fileOutAll="SignalToNoise_SizeAlt_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotNameAll="Bias_SignaToNoise_SizeAlt_global"+psf_ending+alt_file_ending+".png"
        ObsLabel_1='r$\mathrm{S/N}$'
        ObsLabel_2='r$\mathrm{size (sqrt(ab))}$'
        namePlotBinning="Binning_resolution_alt__sn."+psf_ending+alt_file_ending+".png"
    if(ObsType[kk]=="true"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight_1=mag_in_sim_weight
        obs_sim_weight_2=size_in_sim_weight
        fileOutAll="Magnitude_InputSize_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotNameAll="Bias_magnitude_inputSize_global"+psf_ending+alt_file_ending+".png"
        ObsLabel_1='r$\mathrm{mag}$'
        ObsLabel_2='r$\mathrm{size_{in}}$'
        namePlotBinning="Binning_size_mag."+psf_ending+alt_file_ending+".png"
    if(ObsType[kk]=="mag_nn"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight_1=mag_sim_weight
        obs_sim_weight_2=nd_sim_weight
        fileOutAll="Magnitude_nearestNeighbour_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotNameAll="Bias_magnitude_nearestNeighbour_global"+psf_ending+alt_file_ending+".png"
        ObsLabel_1='r$\mathrm{mag}$'
        ObsLabel_2='r$\mathrm{nn}$'
        namePlotBinning="Binning_Magnitude_nearestNeighbour."+psf_ending+alt_file_ending+".png"
    if(ObsType[kk]=="size_nn"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight_1=size_in_sim_weight
        obs_sim_weight_2=nd_sim_weight
        fileOutAll="InputSize_nearestNeighbour_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotNameAll="Bias_InputSize_nearestNeighbour_global"+psf_ending+alt_file_ending+".png"
        ObsLabel_1='r$\mathrm{size_{in}}$'
        ObsLabel_2='r$\mathrm{nn}$'
        namePlotBinning="Binning_InputSize_nearestNeighbour."+psf_ending+alt_file_ending+".png"

    if(ObsType[kk]=="observed_mag_alt"):
        if(psf_rotation==True):
            ObsType[kk]=ObsType[kk]+"PsfRot"
            psf_ending="PsfRot."
        else:
            psf_ending=""
        obs_sim_weight_1=mag_sim_weight
        obs_sim_weight_2=resolution_alt
        fileOutAll="MagnitudeOut_ResolutionAlt_binning_global"+psf_ending+alt_file_ending+".txt"
        PlotNameAll="Bias_MagnitudeOut_ResolutionAlt_global"+psf_ending+alt_file_ending+".png"
        ObsLabel_1='r$\mathrm{mag}$'
        ObsLabel_2='r$\mathrm{resolution (sqrt(ab))}$'
        namePlotBinning="Binning_resolution_alt_mag."+psf_ending+alt_file_ending+".png"



    #Calculate the bins such that each bin contains the same number of points.
    #bins2D=ComputeBins_2D_emp_weighted_plus(obs_sim_weight_1,obs_sim_weight_2,w_sim_mask, Nbin,Nbin2,200,1000,1,50,1,1,Dir, ObsType, namePlotBinning)
    bins2D=ComputeBins_2D_emp_weighted(obs_sim_weight_1,obs_sim_weight_2,w_sim_mask, Nbin,Nbin2,Dir, ObsType, namePlotBinning)
    #bins2D=ComputeBins_2D_emp_linear(obs_sim_weight_1,obs_sim_weight_2, Nbin,Nbin2,2.,180.,0.0,25.,Dir, ObsType, namePlotBinning)
    
    #Run     
    #DO ALL MUST BE ADJUESTED TO TAKE CALIBRATION QUANTITIES, CALIBRATION FUNCTION AND ITS PARAMETERS. 
    if(spline==True):
        if(spline2==True):
            DoAll2Bins(ObsType[kk],bins2D,obs_sim_weight_1,obs_sim_weight_2,mag_in_sim_weight, size_in_sim_alt_weight, mag_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel_1, ObsLabel_2,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100), Nbin,-99,file_ending,cal1,cal2,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines)
        elif(spline_multi==True):
            #DoAll2Bins(ObsType[kk],bins2D,obs_sim_weight_1,obs_sim_weight_2, mag_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel_1, ObsLabel_2,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100), Nbin,-99,file_ending,cal1,cal2,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines,spline_thresholds,mag_spline1,mag_spline2)
            DoAll2Bins(ObsType[kk],bins2D,obs_sim_weight_1,obs_sim_weight_2, mag_in_sim_weight, size_in_sim_alt_weight, mag_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel_1, ObsLabel_2,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100), Nbin,-99,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines,spline_thresholds)
        else:
            DoAll2Bins(ObsType[kk],bins2D,obs_sim_weight_1,obs_sim_weight_2,mag_in_sim_weight, size_in_sim_alt_weight, mag_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel_1, ObsLabel_2,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100), Nbin,-99,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines)
#    elif(function==True):
#        DoAll2Bins(ObsType[kk],bins2D,obs_sim_weight_1,obs_sim_weight_2,mag_in_sim_weight, mag_sim_weight, size_in_sim_alt_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel_1, ObsLabel_2,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100), Nbin,-99,file_ending,cal1,cal2,None,None,None,None,func_m,m1_params,m2_params,func_c,c1_params,c2_params)
    elif(table_cal):
        print "Doing This"
        DoAll2Bins(ObsType[kk],bins2D,obs_sim_weight_1,obs_sim_weight_2,mag_in_sim_weight, size_in_sim_alt_weight, mag_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel_1, ObsLabel_2,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100), Nbin,-99,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,m1_table,m2_table,c1_table,c2_table,shift_x,shift_y,scale_x, scale_y,Tree)
    elif(bin_cal):
        print "Doing This"
        DoAll2Bins(ObsType[kk],bins2D,obs_sim_weight_1,obs_sim_weight_2,mag_in_sim_weight, size_in_sim_alt_weight, mag_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel_1, ObsLabel_2,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100), Nbin,-99,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,table_data)
    elif(file_cal):
        print "Doing This"
        DoAll2Bins(ObsType[kk],bins2D,obs_sim_weight_1,obs_sim_weight_2,mag_in_sim_weight, size_in_sim_alt_weight, mag_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel_1, ObsLabel_2,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100), Nbin,-99,file_ending,cal1,cal2,m1_sim_weight,m2_sim_weight,c1_sim_weight,c2_sim_weight,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None)
    else:
        DoAll2Bins(ObsType[kk],bins2D,obs_sim_weight_1,obs_sim_weight_2,mag_in_sim_weight, size_in_sim_alt_weight,mag_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel_1, ObsLabel_2,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100), Nbin,-99,file_ending)
            

    if(sep_strehl_flag=='Yes' or sep_strehl_flag=='Y' or sep_strehl_flag=='True' or  sep_strehl_flag=='yes' or sep_strehl_flag=='y' or sep_strehl_flag=='true'):

        #Run with Nbins bins. Each psf-set is further splits in 5 as the average PSF properties are different for galaxies having different number of exposures. The binning is done as a function of strehl ratio. Note that defining the binning using the quantiles doesn't work optimally here as the distribution is not smooth. After computing the quantiles I am filtering out the non-unique values and use the resulting array to define the bins. This can be improved.

        Psf_bins=ComputeBins_1D(strehl_sim_weight,5, w_sim_mask)    
        Psf_bins=unique(Psf_bins)
        print Psf_bins
        Nunique=len(Psf_bins)
        print 'Unique psf..', Nunique

        for gg in range(0, Nunique-1):
            mask_gg=(strehl_sim_weight>=Psf_bins[gg])&(strehl_sim_weight<Psf_bins[gg+1])
            strehl_set=numpy.median(strehl_sim_weight[mask_gg])
            bins2D=ComputeBins_2D(obs_sim_weight_1[mask_gg],obs_sim_weight_2[mask_gg], Nbin,Dir, ObsType, namePlotBinning)
            #DO ALL MUST BE ADJUESTED TO TAKE CALIBRATION QUANTITIES, CALIBRATION FUNCTION AND ITS PARAMETERS.

            if(spline==True):
                if(spline2==True):
                    DoAll2Bins(ObsType[kk],bins2D,obs_sim_weight_1,obs_sim_weight_2,mag_in_sim_weight, size_in_sim_alt_weight, mag_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel_1, ObsLabel_2,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100), Nbin,-99,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines)
                elif(spline2==True):
                    DoAll2Bins(ObsType[kk],bins2D,obs_sim_weight_1,obs_sim_weight_2,mag_in_sim_weight, size_in_sim_alt_weight, mag_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel_1, ObsLabel_2,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100), Nbin,-99,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines,spline_thresholds)
                else:
                    DoAll2Bins(ObsType[kk],bins2D,obs_sim_weight_1,obs_sim_weight_2,mag_in_sim_weight, size_in_sim_alt_weight,mag_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel_1, ObsLabel_2,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100), Nbin,-99,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,m1_splines, m2_splines,c1_splines, c2_splines)
#            elif(function==True):
#                DoAll2Bins(ObsType[kk],bins2D,obs_sim_weight_1,obs_sim_weight_2,mag_in_sim_weight, size_in_sim_alt_weight, mag_sim-weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel_1, ObsLabel_2,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100), Nbin,-99,file_ending,cal1,cal2,None,None,None,None,func_m,m1_params,m2_params,func_c,c1_params,c2_params)
            elif(table_cal):
                DoAll2Bins(ObsType[kk],bins2D,obs_sim_weight_1,obs_sim_weight_2,mag_in_sim_weight, size_in_sim_alt_weight, mag_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel_1, ObsLabel_2,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100), Nbin,-99,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,m1_table,m2_table,c1_table,c2_table,shift_x,shift_y,scale_x, scale_y,Tree)
            elif(bin_cal):
                print "Doing This"
                DoAll2Bins(ObsType[kk],bins2D,obs_sim_weight_1,obs_sim_weight_2,mag_in_sim_weight, size_in_sim_alt_weight, mag_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel_1, ObsLabel_2,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100), Nbin,-99,file_ending,cal1,cal2,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,table_data)
            elif(file_cal):
                print "Doing This"
                DoAll2Bins(ObsType[kk],bins2D,obs_sim_weight_1,obs_sim_weight_2,mag_in_sim_weight, size_in_sim_alt_weight, mag_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel_1, ObsLabel_2,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100), Nbin,-99,file_ending,cal1,cal2,m1_sim_weight,m2_sim_weight,c1_sim_weight,c2_sim_weight,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None)
            else:
                DoAll2Bins(ObsType[kk],bins2D,obs_sim_weight_1,obs_sim_weight_2,mag_in_sim_weight, size_in_sim_alt_weight, mag_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOutAll, Flag, PlotNameAll, ObsLabel_1, ObsLabel_2,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight,Dir,TotWeight,Psf_set, str(100), Nbin,-99,file_ending)


        
