from numpy import cos
from numpy import sin
from numpy import arccos
from numpy import arcsin
from numpy import arctan2
import sys
from sys import * 
from scipy import integrate
import scipy
import re 
import io
import os
import numpy
from numpy import *
from numpy import average
from subprocess import Popen, PIPE
from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from fitting import *
from weighted import *
import collections
from scipy import spatial
from scipy.stats.mstats import mquantiles
import matplotlib.patches as patches
from scipy.interpolate import *

def create_splines(file):
    data=genfromtxt(file,comments='#')
    x_1d=data[:,0]
    m1_1d=data[:,1]
    m2_1d=data[:,5]
    spline1_1d=interpolate.interp1d(x_1d,m1_1d,kind='slinear',bounds_error=False)
    spline2_1d=interpolate.interp1d(x_1d,m2_1d,kind='slinear',bounds_error=False)

    return spline1_1d, spline2_1d

def query_multi_bin_splines(snr,R,spline_array,thresholds,mag=None,mag_spline=None):
    bins=len(spline_array)
    spline_index=-1
    for index in range(0,bins):
        if(R>thresholds[index] and R <= thresholds[index+1]):
            spline_index=index
            break
    if(spline_index != -1):
        output=query_spline((snr,R),spline_array[spline_index])
    else:
        output=array(0.0).ravel()

    if(mag!=None):
        output2=query_spline((mag,mag),mag_spline)
        output=output+output2

    return output

def create_binned_splines(file_name, num_splines):
    m1_splines=[]
    m2_splines=[]
    c1_splines=[]
    c2_splines=[]
    thresholds=[]

    data=genfromtxt(file_name,comments='#')

    thresholds.append(data[0][4])
    thresholds.append(data[0][5])
    for index in range(1,num_splines):
        thresholds.append(data[index][5])
    thresholds=array(thresholds)
    
    for index in range(0,num_splines):
        data_m1=data[index:len(data[7]):num_splines,7]
        data_m2=data[index:len(data[11]):num_splines,11]
        data_c1=data[index:len(data[9]):num_splines,9]
        data_c2=data[index:len(data[13]):num_splines,13]
        data_snr=data[index:len(data[0]):num_splines,0]

        m1_spline=interpolate.interp1d(data_snr,data_m1,kind='slinear',bounds_error=False)      
        m1_splines.append(m1_spline)
        m2_spline=interpolate.interp1d(data_snr,data_m2,kind='slinear',bounds_error=False)      
        m2_splines.append(m2_spline)
        c1_spline=interpolate.interp1d(data_snr,data_c1,kind='slinear',bounds_error=False)      
        c1_splines.append(c1_spline)
        c2_spline=interpolate.interp1d(data_snr,data_c2,kind='slinear',bounds_error=False)      
        c2_splines.append(c2_spline)

    m1_splines=array(m1_splines)
    m2_splines=array(m2_splines)
    c1_splines=array(c1_splines)
    c2_splines=array(c2_splines)

    return m1_splines, m2_splines, c1_splines, c2_splines, thresholds

        


def query_spline(xy_tuple, spline):
    (value, dummy) = xy_tuple
    query = spline(value)

    if(math.isnan(query)):
        spline_max=max(spline.x)
        spline_min=min(spline.x)
        if(value >= spline_max):
            output=spline.y[argmax(spline.x)]
        else:
            output=spline.y[argmin(spline.x)]
    else:
        output = query

    return output.ravel()

def query_spline_emp(xy_tuple, spline):
    (value1, value2) = xy_tuple
    if(value1 > 80. and value2 < 4):
        output=array(0.0)
    else:
        query = spline(value1)
        
        if(math.isnan(query)):
            spline_max=max(spline.x)
            spline_min=min(spline.x)
            if(value1 >= spline_max):
                output=spline.y[argmax(spline.x)]
            else:
                output=spline.y[argmin(spline.x)]
        else:
            output = query
                
    return output.ravel()

def query_spline2(xy_tuple, spline1, spline2):
    (value1, value2) = xy_tuple
    query1 = spline1(value1)
    query2 = spline2(value2)

    if(math.isnan(query1)):
        spline_max=max(spline1.x)
        spline_min=min(spline1.x)
        if(value1 >= spline_max):
            output1=spline1.y[argmax(spline1.x)]
        else:
            output1=spline1.y[argmin(spline1.x)]
    else:
        output1 = query1

    if(math.isnan(query2)):
        spline_max=max(spline2.x)
        spline_min=min(spline2.x)
        if(value2 >= spline_max):
            output2=spline2.y[argmax(spline2.x)]
        else:
            output2=spline2.y[argmin(spline2.x)]
    else:
        output2 = query2

    return (output1+output2).ravel()

def query_spline2_emp(xy_tuple, spline1, spline2):
    (value1, value2) = xy_tuple
    if(value1 > 80. and value2 < 4.):
        output1=array(0.0)
        output2=array(0.0)
    else:
        query1 = spline1(value1)
        query2 = spline2(value2)

        if(math.isnan(query1)):
            spline_max=max(spline1.x)
            spline_min=min(spline1.x)
            if(value1 >= spline_max):
                output1=spline1.y[argmax(spline1.x)]
            else:
                output1=spline1.y[argmin(spline1.x)]
        else:
            output1 = query1

        if(math.isnan(query2)):
            spline_max=max(spline2.x)
            spline_min=min(spline2.x)
            if(value2 >= spline_max):
                output2=spline2.y[argmax(spline2.x)]
            else:
                output2=spline2.y[argmin(spline2.x)]
        else:
            output2 = query2

    return (output1+output2).ravel()



def unique_rows(data):
    """ convenience function to find unique elements in an array; there seem to be issues with float precision though"""
    b = ascontiguousarray(data).view(dtype((void, data.dtype.itemsize*data.shape[1])))
    _, idx = unique(b, return_index=True)

    return idx

def SizeFromMom(psf_Q11_sim,psf_Q22_sim,psf_Q12_sim):
    return (psf_Q11_sim*psf_Q22_sim-psf_Q12_sim**2)**0.5

def ComputeBins_1D(sn, Nbin, weight):
    #Remove objects with weight==0 
    mask=(weight>0)
    snMask=sn[mask]
    weightMask=weight[mask]

    #Define the quantiles limits based on the number of bins

    Nbin_quanta=[]
    for x in range(0,Nbin+1):
        Nbin_quanta.append(x/float(Nbin))

    #Calculate quantiles for s/n
    sn_quanta=mquantiles(snMask, prob=Nbin_quanta)

    return (sn_quanta)

def ComputeBins_1D_weighted(sn, Nbin, weight):
    #Remove objects with weight==0 
    mask=(weight>0)
    snMask=sn[mask]
    weightMask=weight[mask]

    #Define the quantiles limits based on the number of bins

    Nbin_quanta=[]
    for x in range(0,Nbin+1):
        Nbin_quanta.append(x/float(Nbin))

    #Calculate quantiles for s/n
    sn_quanta=[]
    for x in range(0,len(Nbin_quanta)):
        sn_quanta.append(quantile_1D(sn, weight,Nbin_quanta[x]))
    sn_quanta=array(sn_quanta)
    sn_quanta=unique(sn_quanta)



    return (sn_quanta)

def ComputeBins_2D(sn, size, Nbin,Dir,ObsType,namePlot):
    
    #Define the quantiles limits based on the number of bins

    Nbin_quanta=[]
    for x in range(0,Nbin+1):
        Nbin_quanta.append(x/float(Nbin))

    #Calculate quantiles for s/n
    sn_quanta=mquantiles(sn, prob=Nbin_quanta)
    sn_quanta=unique(sn_quanta)
    #Now computing quantiles for size in each sn bin 

    y_bin=[]
    x_bin=[]
    for xx in range(0,len(sn_quanta)-1):
        maskTmp=((sn>=sn_quanta[xx])&(sn<sn_quanta[xx+1]))
        size_quanta=mquantiles(size[maskTmp], prob=Nbin_quanta)
        size_quanta=unique(size_quanta)
        y_bin.append(size_quanta)
        x_bin.append(sn_quanta)

    x_bin=array(x_bin)
    y_bin=array(y_bin)
    sn=array(sn)
    size=array(size)

    #Make a nice plot
    plt.figure()
    if(ObsType=='observed'):
        plt.semilogx()
        plt.semilogy()
        plt.xlim(1,300)
        plt.ylim(0.2,8)
        plt.xlabel(r'$\mathrm{S/N}$', fontsize=20)
        plt.ylabel(r'$\mathrm{size}$', fontsize=20)
    if(ObsType=='true'):
        plt.semilogy()
        plt.xlim(19,25.5)
        plt.ylim(0.2,8)
        plt.xlabel(r'$\mathrm{mag}$', fontsize=20)
        plt.ylabel(r'$\mathrm{size}$', fontsize=20)  
    plt.scatter(sn,size, s=1)
    currentAxis = plt.gca()
    for x in range(0,len(x_bin[0])-1):
        for y in range(0,len(y_bin[x])-1):
            len_x=x_bin[0][x+1]-x_bin[0][x]
            len_y=y_bin[x][y+1]-y_bin[x][y]
            maskBin=((sn>=x_bin[0][x])&(sn<x_bin[0][x+1])&(size>=y_bin[x][y])&(size<y_bin[x][y+1]))
            currentAxis.add_patch(patches.Rectangle((x_bin[0][x], y_bin[x][y]),len_x,len_y,fill=None, alpha=1, color='red', linewidth=3))

    plt.savefig(Dir+'/Results/2bin/'+namePlot)
    return (x_bin, y_bin)

def ComputeBins_2D_emp(sn, size, Nbin1,Nbin2,Dir,ObsType,namePlot):
    
    #Define the quantiles limits based on the number of bins

    Nbin1_quanta=[]
    Nbin2_quanta=[]
    for x in range(0,Nbin1+1):
        Nbin1_quanta.append(x/float(Nbin1))
    for x in range(0,Nbin2+1):
        Nbin2_quanta.append(x/float(Nbin2))

    #Calculate quantiles for s/n
    sn_quanta=mquantiles(sn, prob=Nbin1_quanta)
    sn_quanta=unique(sn_quanta)
    #Now computing quantiles for size in each sn bin 

    y_bin=[]
    x_bin=[]
    #dummy=(0.0,2.5,7.0,1.e6)
    for xx in range(0,len(sn_quanta)-1):
        maskTmp=((sn>=sn_quanta[xx])&(sn<sn_quanta[xx+1]))
        size_quanta=mquantiles(size[maskTmp], prob=Nbin2_quanta)
        size_quanta=unique(size_quanta)
        y_bin.append(size_quanta)
        x_bin.append(sn_quanta)

    x_bin=array(x_bin)
    y_bin=array(y_bin)
    sn=array(sn)
    size=array(size)

    return (x_bin, y_bin)

def ComputeBins_2D_emp_weighted(sn, size,weight, Nbin1,Nbin2,Dir,ObsType,namePlot):
    
    #Define the quantiles limits based on the number of bins

    Nbin1_quanta=[]
    Nbin2_quanta=[]
    for x in range(0,Nbin1+1):
        Nbin1_quanta.append(x/float(Nbin1))
    for x in range(0,Nbin2+1):
        Nbin2_quanta.append(x/float(Nbin2))

    #Calculate quantiles for s/n
    sn_quanta=[]
    for x in range(0,len(Nbin1_quanta)):
        sn_quanta.append(quantile_1D(sn, weight,Nbin1_quanta[x]))
    sn_quanta=unique(sn_quanta)
    #Now computing quantiles for size in each sn bin 
    
    y_bin=[]
    x_bin=[]
    #dummy=(0.0,2.5,7.0,1.e6)
    for xx in range(0,len(sn_quanta)-1):
        maskTmp=((sn>=sn_quanta[xx])&(sn<sn_quanta[xx+1]))
        size_quanta=[]
        for kk in range(0,len(Nbin2_quanta)):
            size_quanta.append(quantile_1D(size[maskTmp], weight[maskTmp],Nbin2_quanta[kk]))
        size_quanta=unique(size_quanta)
        y_bin.append(size_quanta)
        x_bin.append(sn_quanta)

    x_bin=array(x_bin)
    y_bin=array(y_bin)
    sn=array(sn)
    size=array(size)

    return (x_bin, y_bin)

def ComputeBins_2D_emp_linear(sn, size, Nbin1,Nbin2,lower1,upper1,lower2,upper2,Dir,ObsType,namePlot):
    
    #Define the quantiles limits based on the number of bins

    interval1 = upper1-lower1
    step1 = interval1/Nbin1
    interval2 = upper2-lower2
    step2 = interval2/Nbin2

    y_bin=[]
    x_bin=[]
    current_x=[]
    current_y=[]
    for x in range(0,Nbin1+1):
        if(x==Nbin1):
            current_x.append(1000.)
        else:
            current_x.append(lower1+(x+1)*step1)
        
    for x in range(0,Nbin2+1):
        if(x==Nbin2):
            current_y.append(1000.)
        else:
            current_y.append(lower2+(x+1)*step2)       
    current_x=array(current_x)
    current_y=array(current_y) 

    for x in range(0,Nbin1+1):
        x_bin.append(current_x)
    for x in range(0,Nbin2+1):
        y_bin.append(current_y)
        
    x_bin=array(x_bin)
    y_bin=array(y_bin)


    return (x_bin, y_bin)


def ComputeBins_2D_emp_weighted_plus(sn, size,weight, Nbin1,Nbin2,thresh1, thresh1_2, thresh2, thresh2_2, Nbin1_2,Nbin2_2, Dir,ObsType,namePlot):
    
    sn=array(sn)
    size=array(size)


    #Creating additional bins

    snr_field = thresh1_2 - thresh1
    size_field = thresh2_2 - thresh2
    snr_step = snr_field / Nbin1_2
    size_step = size_field / Nbin2_2

    add_x=[]
    for i in range(1,Nbin1_2+2):
        if(i==(Nbin1_2+1)):
           add_x.append(600.)
        else:
           add_x.append(thresh1+i*snr_step)


    add_y=[]
    for i in range(1,Nbin2_2+2):
        if(i==(Nbin2_2+1)):
           add_y.append(1000.)
        else:
           add_y.append(thresh2+i*size_step)
        
    add_x=array(add_x)
    add_y=array(add_y)

    #First limit the data within the thresholds
    snr_mask=(sn<thresh1)
    sn_one=sn[snr_mask]
    size_one=size[snr_mask]
    weight_one=weight[snr_mask]
    size_mask=size_one<thresh2
    sn_limited=sn_one[size_mask]
    size_limited=size_one[size_mask]
    weight_limited=weight_one[size_mask]



    #Define the quantiles limits based on the number of bins
    Nbin1_quanta=[]
    Nbin2_quanta=[]
    for x in range(0,Nbin1+1):
        Nbin1_quanta.append(x/float(Nbin1))
    for x in range(0,Nbin2+1):
        Nbin2_quanta.append(x/float(Nbin2))

    #Calculate quantiles for s/n
    sn_quanta=[]
    for x in range(0,len(Nbin1_quanta)):
        sn_quanta.append(quantile_1D(sn_limited, weight_limited,Nbin1_quanta[x]))
    sn_quanta=unique(sn_quanta)
    #Now computing quantiles for size in each sn bin 
    
    y_bin=[]
    x_bin=[]
    #dummy=(0.0,2.5,7.0,1.e6)
    for xx in range(0,len(sn_quanta)-1):
        #x_bin_aux=[]
        #y_bin_aux=[]
        maskTmp=((sn_limited>=sn_quanta[xx])&(sn_limited<sn_quanta[xx+1]))
        size_quanta=[]
        for kk in range(0,len(Nbin2_quanta)):
            size_quanta.append(quantile_1D(size_limited[maskTmp], weight_limited[maskTmp],Nbin2_quanta[kk]))
        size_quanta=unique(size_quanta)
        #for i in range(0,len(size_quanta)):
        #    y_bin_aux.append(size_quanta[i])
        #for i in range(0,len(add_y)):
        #    y_bin_aux.append(add_y[i])
        #for i in range(0,len(sn_quanta)):
        #    x_bin_aux.append(sn_quanta[i])
        #for i in range(0,len(add_x)):
        #    x_bin_aux.append(add_x[i])
        #x_bin_aux=array(x_bin_aux)
        #y_bin_aux=array(y_bin_aux)
        x_bin.append(sn_quanta)
        y_bin.append(size_quanta)

    x_bin=array(x_bin)
    y_bin=array(y_bin)

    x_bin_out=[]
    for i in range(0,len(x_bin)+len(add_x)+1):
        current_x = append(x_bin[0],add_x)
        x_bin_out.append(current_x)

    y_bin_out=[]
    for i in range(0,len(y_bin)+len(add_y)+1):
        current_y=append(y_bin[0],add_y)
        y_bin_out.append(current_y)

    x_bin_out = array(x_bin_out)
    y_bin_out = array(y_bin_out)

    print shape(x_bin_out), shape(y_bin_out)

    return (x_bin_out, y_bin_out)



def func(x,a,b):
    return (1.0+a)*x+b

def rotation(e1,e2,e1psf,e2psf):
    "This function rotates galaxies with respect to a reference frame in which the psf has e2_psf=0."
    angle_gal=0.5*arctan2(e2,e1)
    angle_psf=0.5*arctan2(e2psf,e1psf)

    #Rotation angle

    angle_rot=-angle_psf#-angle_gal

    e1_rot=e1*cos(2.*angle_rot)-e2*sin(2.*angle_rot)
    e2_rot=e1*sin(2.*angle_rot)+e2*cos(2.*angle_rot)

    return [e1_rot,e2_rot]

def DoAll(Observable,ObsBin, obs_sim_weight, w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOut, Flag, PlotName, ObsLabel,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight, Dir,TotWeight,Psf_set,subSet,strehl_set,ending_flag,snr=None,res=None,m1=None,m2=None,c1=None,c2=None,m_func=None,m1_params=None,m2_params=None,c_func=None,c1_params=None,c2_params=None,m1_splines=None,m2_splines=None,c1_splines=None,c2_splines=None,spline_thresholds=None, m1_table=None, m2_table=None, c1_table=None, c2_table=None, x_shift=None, y_shift=None, x_scale=None, y_scale=None, tree=None,table_data=None):

    cal=False
    spline_cal=False
    spline_cal2=False
    spline_cal_multi=False
    table_cal=False
    bin_cal=False
    file_cal=False
    if(snr!=None):
        cal=True
        if(m1_splines!=None):
            spline_cal=True
            if(spline_thresholds!=None):
                spline_cal_multi=True
            elif(len(m1_splines)>1):
                spline_cal2=True
        elif(tree!=None):
            table_cal=True
        elif(table_data!=None):
            bin_cal=True
        elif(m1!=None):
            file_cal=True
        
    
    if(ending_flag=='nocorr'):
        if(cal==False):
            fLog=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'.log','wb')
            fOutShear=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear.dat','wb')
        elif(spline_cal==False):
            if(table_cal==True):
                fLog=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_table_cal.log','wb')
                fOutShear=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'shear_table_cal.dat','wb')
            elif(bin_cal==True):
                fLog=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_bin_cal.log','wb')
                fOutShear=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'shear_bin_cal.dat','wb')
            elif(file_cal==True):
                fLog=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_file_cal.log','wb')
                fOutShear=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'shear_file_cal.dat','wb')
            else:
                fLog=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_cal.log','wb')
                fOutShear=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear_cal.dat','wb')
        elif(spline_cal==True):
            if(spline_cal2==True):
                fLog=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_2spline_cal.log','wb')
                fOutShear=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear_2spline_cal.dat','wb')
            elif(spline_cal_multi==True):
                fLog=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_multi_spline_cal.log','wb')
                fOutShear=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear_multi_spline_cal.dat','wb')
            else:
                fLog=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_spline_cal.log','wb')
                fOutShear=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear_spline_cal.dat','wb')

    else:
        if(cal==False):
            fLog=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'.'+ending_flag+'.log','wb')
            fOutShear=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear.'+ending_flag+'.dat','wb')
        elif(spline_cal==False):
            if(table_cal==True):
                fLog=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'.'+ending_flag+'_table_cal.log','wb')
                fOutShear=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear.'+ending_flag+'_table_cal.dat','wb')
            elif(bin_cal==True):
                fLog=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'.'+ending_flag+'_bin_cal.log','wb')
                fOutShear=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear.'+ending_flag+'_bin_cal.dat','wb')
            elif(file_cal==True):
                fLog=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'.'+ending_flag+'_file_cal.log','wb')
                fOutShear=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear.'+ending_flag+'_file_cal.dat','wb')
            else:
                fLog=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'.'+ending_flag+'_cal.log','wb')
                fOutShear=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear.'+ending_flag+'_cal.dat','wb')
        elif(spline_cal==True):
            if(spline_cal2==True):
                fLog=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'.'+ending_flag+'_2spline_cal.log','wb')
                fOutShear=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear.'+ending_flag+'_2spline_cal.dat','wb')
            elif(spline_cal2==True):
                fLog=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'.'+ending_flag+'_multi_spline_cal.log','wb')
                fOutShear=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear.'+ending_flag+'_multi_spline_cal.dat','wb')
            else:
                fLog=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'.'+ending_flag+'_spline_cal.log','wb')
                fOutShear=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear.'+ending_flag+'_spline_cal.dat','wb')

    print >> fOutShear, "#Median bin", "Median strehl", "g1_in","g2_in","g1_out","g2_out","g1_out_nc","g2_out_nc","g1_out_err","g2_out_err","g1_out_err_nc","g2_out_err_nc", "m1", "m2","c1","c2"
    ObsBin_mean=[]
    BinWeight=[]
    BinCount=[]
    p0=[0,0]

    m1_m=[]
    m1_m_err=[]
    m1_m_err_BS=[]
    c1_m=[]
    c1_m_err=[]
    c1_m_err_BS=[]
    m2_m=[]
    m2_m_err=[]
    m2_m_err_BS=[]
    c2_m=[]
    c2_m_err=[]
    c2_m_err_BS=[]

    m1_m_nc=[]
    m1_m_err_nc=[]
    m1_m_err_BS_nc=[]
    c1_m_nc=[]
    c1_m_err_nc=[]
    c1_m_err_BS_nc=[]
    m2_m_nc=[]
    m2_m_err_nc=[]
    m2_m_err_BS_nc=[]
    c2_m_nc=[]
    c2_m_err_nc=[]
    c2_m_err_BS_nc=[]

    average_m1_bin=[]
    average_m2_bin=[]
    average_c1_bin=[]
    average_c2_bin=[]
    
    print >> fLog, '############################################'
    print >> fLog, 'total number of used objects : ', len(w_sim_mask)
    print >> fLog, '############################################'

    for x in range(0,len(ObsBin)-1):
        if(Observable=="nn"):
            maskBinObs=(obs_sim_weight>=ObsBin[x])&(obs_sim_weight<ObsBin[x+1])
            #ObsBin_mean.append(ObsBin[x])
            w_sim_mask_obs=w_sim_mask[maskBinObs]
            ObsBin_mean.append(average(obs_sim_weight[maskBinObs],weights=w_sim_mask_obs))

        else:
            maskBinObs=(obs_sim_weight>=ObsBin[x])&(obs_sim_weight<ObsBin[x+1])
            w_sim_mask_obs=w_sim_mask[maskBinObs]
            ObsBin_mean.append(average(obs_sim_weight[maskBinObs],weights=w_sim_mask_obs))
        print  >> fLog, 'total number of objects in bin ',ObsBin[x], ObsBin[x+1],  ': ', len(w_sim_mask[maskBinObs])     
        print >> fLog, '############################################'

        e1_sim_mask_obs=e1_sim_mask[maskBinObs]
        e2_sim_mask_obs=e2_sim_mask[maskBinObs]
        e1_in_sim_mask_obs=e1_in_sim_weight[maskBinObs]
        e2_in_sim_mask_obs=e2_in_sim_weight[maskBinObs]

        #Define the ellipticity by remove self calibration term
        e1_sim_mask_obs_nc=e1_sim_mask[maskBinObs]-e1_corr_sim_weight[maskBinObs]
        e2_sim_mask_obs_nc=e2_sim_mask[maskBinObs]-e2_corr_sim_weight[maskBinObs]

        #snr and R inside the chosen observable bin.
        if(cal==True):
            calibration_snr_sim_mask_obs=snr[maskBinObs]
            calibration_resolution_sim_mask_obs=res[maskBinObs]
        if(file_cal==True):
            cal_m1=m1[maskBinObs]
            cal_m2=m2[maskBinObs]
            cal_c1=c1[maskBinObs]
            cal_c2=c2[maskBinObs]
            #if(mag!=None):
            #    cal_mag=mag[maskBinObs]
            #else:
            #    cal_mag=None

        #print median(obs_sim_weight[maskBinObs]), median(calibration_snr_sim_mask_obs), median(calibration_resolution_sim_mask_obs)

        #Compute weight in the bin

        BinWeight.append(sum(w_sim_mask_obs))

        #Total objects in the bin

        BinCount.append(len(w_sim_mask_obs))

        g1_in_mask_obs=g1_in_mask[maskBinObs]
        g2_in_mask_obs=g2_in_mask[maskBinObs]
    
        g1_out=[]
        g2_out=[]
        g1_out_nc=[]
        g2_out_nc=[]
        g_out_w=[]

        g1_in_used=[]
        g2_in_used=[]
        average_m1_corr=[]
        average_m2_corr=[]

        average_c1_corr=[]
        average_c2_corr=[]

        #print 'num shears..', len(g1_theo)
        for kkk in range(0, len(g1_theo)):
            #print 'input shears:', g1_theo[kkk],g2_theo[kkk], len(g1_in_mask_obs)
            maskShear=(g1_in_mask_obs==g1_theo[kkk])&(g2_in_mask_obs==g2_theo[kkk])
            numMasked=len(e1_sim_mask_obs[maskShear])
            if (numMasked >0):
                #Calculating bin average for calibration quantities
                if(cal==True):
                    #xy_tuple_correction=(calibration_snr_sim_mask_obs[maskShear],calibration_resolution_sim_mask_obs[maskShear],cal_mag[maskShear])
                    xy_tuple_correction=(calibration_snr_sim_mask_obs[maskShear],calibration_resolution_sim_mask_obs[maskShear])
                    if(table_cal==True):
                        xy=array(xy_tuple_correction)
                        m1_sample=[]
                        m2_sample=[]
                        c1_sample=[]
                        c2_sample=[]
                        weight_sample=w_sim_mask_obs[maskShear]
                        for x_quant in xy.T:
                            all_biases = bias_2D_table_query_all(x_quant[0], x_quant[1],m1_table,m2_table,c1_table,c2_table,tree,x_scale,y_scale,x_shift,y_shift)
                            m1_sample.append(all_biases[0])
                            m2_sample.append(all_biases[1])
                        m1_sample=array(m1_sample)
                        m2_sample=array(m2_sample)

                        current_m1=average(m1_sample,weights=weight_sample)
                        current_m2=average(m2_sample,weights=weight_sample)
                        #I NEEDED TO KEEP THIS LINE IN LIKE THIS SINCE SPLINES CANNOT BE CHOSEN TO GIVE 0
                        current_c1=0.0#average(c1_sample,weights=weight_sample)
                        current_c2=0.0#average(c2_sample,weights=weight_sample)
                    elif(bin_cal==True):
                        xy=array(xy_tuple_correction)
                        m1_sample=[]
                        m2_sample=[]
                        c1_sample=[]
                        c2_sample=[]
                        weight_sample=w_sim_mask_obs[maskShear]
                        for x_quant in xy.T:
                            all_biases = bias_2D_bin_query(x_quant[0], x_quant[1],table_data)
                            m1_sample.append(all_biases[0])
                            m2_sample.append(all_biases[1])
                        m1_sample=array(m1_sample)
                        m2_sample=array(m2_sample)

                        current_m1=average(m1_sample,weights=weight_sample)
                        current_m2=average(m2_sample,weights=weight_sample)
                        #I NEEDED TO KEEP THIS LINE IN LIKE THIS SINCE SPLINES CANNOT BE CHOSEN TO GIVE 0
                        current_c1=0.0#average(c1_sample,weights=weight_sample)
                        current_c2=0.0#average(c2_sample,weights=weight_sample)
            
                    elif(spline_cal==True):
                        if(spline_cal2==True):
                            xy=array(xy_tuple_correction)
                            m1_sample=[]
                            m2_sample=[]
                            c1_sample=[]
                            c2_sample=[]
                            weight_sample=w_sim_mask_obs[maskShear]
                            for x_quant in xy.T:
                                m1_sample.append(float(query_spline2_emp(x_quant,m1_splines[0],m1_splines[1])))
                                m2_sample.append(float(query_spline2_emp(x_quant,m2_splines[0],m2_splines[1])))
                                #c1_sample.append(float(query_spline2_emp(x_quant,c1_spline,c1_spline2)))
                                #c2_sample.append(float(query_spline2_emp(x_quant,c2_spline,c2_spline2)))
                            m1_sample=array(m1_sample)
                            m2_sample=array(m2_sample)
                            #c1_sample=array(c1_sample)
                            #c2_sample=array(c2_sample)

                        elif(spline_cal_multi==True):
                            xy=array(xy_tuple_correction)
                            m1_sample=[]
                            m2_sample=[]
                            c1_sample=[]
                            c2_sample=[]
                            weight_sample=w_sim_mask_obs[maskShear]
                            for x_quant in xy.T:
                                (one,two) = x_quant
                                #m1_sample.append(float(query_multi_bin_splines(one,two,m1_splines,spline_thresholds,three,mag_spline1)))
                                m1_sample.append(float(query_multi_bin_splines(one,two,m1_splines,spline_thresholds)))
                                #m2_sample.append(float(query_multi_bin_splines(one,two,m2_splines,spline_thresholds,three,mag_spline2)))
                                m2_sample.append(float(query_multi_bin_splines(one,two,m2_splines,spline_thresholds)))
                                #c1_sample.append(float(query_spline2_emp(x_quant,c1_spline,c1_spline2)))
                                #c2_sample.append(float(query_spline2_emp(x_quant,c2_spline,c2_spline2)))
                            m1_sample=array(m1_sample)
                            m2_sample=array(m2_sample)
                            #c1_sample=array(c1_sample)
                            #c2_sample=array(c2_sample)



                        else:
                            xy=array(xy_tuple_correction)
                            m1_sample=[]
                            m2_sample=[]
                            c1_sample=[]
                            c2_sample=[]
                            weight_sample=w_sim_mask_obs[maskShear]
                            for x_quant in xy.T:
                                m1_sample.append(float(query_spline_emp(x_quant,m1_splines[0])))
                                m2_sample.append(float(query_spline_emp(x_quant,m2_splines[0])))
                                #c1_sample.append(float(query_spline_emp(x_quant,c1_spline)))
                                #c2_sample.append(float(query_spline_emp(x_quant,c2_spline)))
                            m1_sample=array(m1_sample)
                            m2_sample=array(m2_sample)
                            #c1_sample=array(c1_sample)
                            #c2_sample=array(c2_sample)

                                
                        current_m1=average(m1_sample,weights=weight_sample)
                        current_m2=average(m2_sample,weights=weight_sample)
                        #I NEEDED TO KEEP THIS LINE IN LIKE THIS SINCE SPLINES CANNOT BE CHOSEN TO GIVE 0
                        current_c1=0.0#average(c1_sample,weights=weight_sample)
                        current_c2=0.0#average(c2_sample,w            popt1_m, punc1_m, rc1_m, dc1_m,pcov1_m =general_fit(func,g1_in_used,g1_out,p0,g_out_w)eights=weight_sample)

                    elif(file_cal==True):
                        current_m1=average(cal_m1[maskShear],weights=w_sim_mask_obs[maskShear])
                        current_m2=average(cal_m2[maskShear],weights=w_sim_mask_obs[maskShear])
                        current_c1=average(cal_c1[maskShear],weights=w_sim_mask_obs[maskShear])
                        current_c2=average(cal_c2[maskShear],weights=w_sim_mask_obs[maskShear])

                    else:
                        current_m1=average(m_func(xy_tuple_correction,*m1_params),weights=w_sim_mask_obs[maskShear])
                        current_m2=average(m_func(xy_tuple_correction,*m2_params),weights=w_sim_mask_obs[maskShear])
                        current_c1=average(c_func(xy_tuple_correction,*c1_params),weights=w_sim_mask_obs[maskShear])
                        current_c2=average(c_func(xy_tuple_correction,*c2_params),weights=w_sim_mask_obs[maskShear])

                if(Flag=='CH'):
                    if(cal==False):
                        g1_out.append(average(e1_sim_mask_obs[maskShear]-e1_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear]))
                        g_out_w.append(1./(sum(w_sim_mask_obs[maskShear]))**0.5)
                        g1_out_print=average(e1_sim_mask_obs[maskShear]-e1_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                        g1_out_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                        g2_out.append(std(e2_sim_mask_obs[maskShear]-e2_in_sim_mask_obs[maskShear]))
                        g2_out_print=average(e2_sim_mask_obs[maskShear]-e2_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                        g2_out_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                        g1_out_nc.append(average(e1_sim_mask_obs_nc[maskShear]-e1_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear]))
                        g2_out_nc.append(average(e2_sim_mask_obs_nc[maskShear]-e2_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear]))
                        g1_out_nc_print=average(e1_sim_mask_obs_nc[maskShear]-e1_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                        g1_out_nc_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                        g2_out_nc_print=average(e2_sim_mask_obs_nc[maskShear]-e2_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                        g2_out_nc_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                        g1_in_used.append(g1_theo[kkk])
                        g2_in_used.append(g2_theo[kkk])
                        g1_in_print=g1_theo[kkk]
                        g2_in_print=g2_theo[kkk]
                        average_m1_corr.append(0)
                        average_m2_corr.append(0)

                        average_c1_corr.append(0)
                        average_c2_corr.append(0)

                        print >> fOutShear, numpy.median(obs_sim_weight[maskBinObs]), strehl_set,g1_in_print,g2_in_print,g1_out_print,g2_out_print,g1_out_nc_print,g2_out_nc_print,g1_out_err_print,g2_out_err_print,g1_out_nc_err_print,g2_out_nc_err_print,0,0,0,0
                    else:
                        #Compute the average ellipticity (subtracting the intrinsic one..)
                        #...for selfcal
                        value1=average(e1_sim_mask_obs[maskShear]-e1_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                        value2=average(e2_sim_mask_obs[maskShear]-e2_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                        value1_w=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                        #...and removing the self-calibration
                        value1_nc=average(e1_sim_mask_obs_nc[maskShear]-e1_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                        value2_nc=average(e2_sim_mask_obs_nc[maskShear]-e2_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                        #...apply the average c and m correction to the shear value from selfcal
                        value1_c=(value1-current_c1)/(1.+current_m1)
                        value2_c=(value2-current_c2)/(1.+current_m2)
                        #...apply the average c and m correction to the shear value from non-selfcal
                        value1_nc_c=(value1_nc-current_c1)/(1.+current_m1)
                        value2_nc_c=(value2_nc-current_c2)/(1.+current_m2)
                        g1_out.append(value1_c)
                        g_out_w.append(value1_w)
                        g1_out_print=value1_c
                        g1_out_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                        g2_out.append(value2_c)
                        g2_out_print=value2_c
                        g2_out_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                        g1_out_nc.append(value1_nc_c)
                        g2_out_nc.append(value2_nc_c)
                        g1_out_nc_print=value1_nc_c
                        g1_out_nc_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                        g2_out_nc_print=value2_nc_c
                        g2_out_nc_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                        g1_in_used.append(g1_theo[kkk])
                        g2_in_used.append(g2_theo[kkk])
                        g1_in_print=g1_theo[kkk]
                        g2_in_print=g2_theo[kkk]
                        average_m1_corr.append(current_m1)
                        average_m2_corr.append(current_m2)
                        average_c1_corr.append(current_c1)
                        average_c2_corr.append(current_c2)
                        print >> fOutShear, numpy.median(obs_sim_weight[maskBinObs]), strehl_set,g1_in_print,g2_in_print,g1_out_print,g2_out_print,g1_out_nc_print,g2_out_nc_print,g1_out_err_print,g2_out_err_print,g1_out_nc_err_print,g2_out_nc_err_print, current_m1,current_m2, current_c1,current_c2

                else:
                    if(cal==False):
                        g1_out.append(average(e1_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear]))
                        g_out_w.append(1./(sum(w_sim_mask_obs[maskShear]))**0.5)
                        g1_out_print=average(e1_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                        g1_out_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                        g2_out.append(average(e2_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear]))
                        g2_out_print=average(e2_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                        g2_out_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                        g1_out_nc.append(average(e1_sim_mask_obs_nc[maskShear], weights=w_sim_mask_obs[maskShear]))
                        g2_out_nc.append(average(e2_sim_mask_obs_nc[maskShear], weights=w_sim_mask_obs[maskShear]))
                        g1_out_nc_print=average(e1_sim_mask_obs_nc[maskShear], weights=w_sim_mask_obs[maskShear])
                        g1_out_nc_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                        g2_out_nc_print=average(e2_sim_mask_obs_nc[maskShear], weights=w_sim_mask_obs[maskShear])
                        g2_out_nc_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                        g1_in_used.append(g1_theo[kkk])
                        g2_in_used.append(g2_theo[kkk])
                        g1_in_print=g1_theo[kkk]
                        g2_in_print=g2_theo[kkk]
                        average_m1_corr.append(0)
                        average_m2_corr.append(0)
                        average_c1_corr.append(0)
                        average_c2_corr.append(0)
                        print >> fOutShear, numpy.median(obs_sim_weight[maskBinObs]),strehl_set, g1_in_print,g2_in_print,g1_out_print,g2_out_print,g1_out_nc_print,g2_out_nc_print,g1_out_err_print,g2_out_err_print,g1_out_nc_err_print,g2_out_nc_err_print,0,0,0,0
                    else:
                        value1=average(e1_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                        value2=average(e2_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                        value1_w=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                        value1_nc=average(e1_sim_mask_obs_nc[maskShear], weights=w_sim_mask_obs[maskShear])
                        value2_nc=average(e2_sim_mask_obs_nc[maskShear], weights=w_sim_mask_obs[maskShear])
                        value1_c=(value1-current_c1)/(1.+current_m1)
                        value2_c=(value2-current_c2)/(1.+current_m2)
                        value1_nc_c=(value1_nc-current_c1)/(1.+current_m1)
                        value2_nc_c=(value2_nc-current_c2)/(1.+current_m2)
                        #print 'g1',value1, current_c,(1.+current_m),value1_nc_c
                        #print 'g2', value2, current_c,(1.+current_m),value2_nc_c

                        g1_out.append(value1_c)
                        g_out_w.append(value1_w)
                        g1_out_print=value1_c
                        g1_out_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                        g2_out.append(value2_c)
                        g2_out_print=value2_c
                        g2_out_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                        g1_out_nc.append(value1_nc_c)
                        g2_out_nc.append(value2_nc_c)
                        g1_out_nc_print=value1_nc_c
                        g1_out_nc_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                        g2_out_nc_print=value2_nc_c
                        g2_out_nc_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                        g1_in_used.append(g1_theo[kkk])
                        g2_in_used.append(g2_theo[kkk])
                        g1_in_print=g1_theo[kkk]
                        g2_in_print=g2_theo[kkk]
                        average_m1_corr.append(current_m1)
                        average_m2_corr.append(current_m2)
                        average_c1_corr.append(current_c1)
                        average_c2_corr.append(current_c2)

                        print >> fOutShear, numpy.median(obs_sim_weight[maskBinObs]), strehl_set,g1_in_print,g2_in_print,g1_out_print,g2_out_print,g1_out_nc_print,g2_out_nc_print,g1_out_err_print,g2_out_err_print,g1_out_nc_err_print,g2_out_nc_err_print, current_m1,current_m2, current_c1,current_c2

        #print g1_theo[kk], average(e1_sim_mask_obs[maskShear], weights=w_sim_maskobs[maskShear])
        # Check I have at least 3 shear values to do the regression
        numShear=len(g1_out)
        if(numShear<3):
            print 'Cannot do the regression in bin ', ObsBin[x], ObsBin[x+1], ' less than 3 shear values! (', numShear, ')'
            exit()
            continue
    
        else:
            #print 'I am using ', numShear, 'values of the shear'
            g1_in_used=array(g1_in_used)
            g2_in_used=array(g2_in_used)
            g1_out=array(g1_out)
            g2_out=array(g2_out)
            g_out_w=array(g_out_w)
            g1_out_nc=array(g1_out_nc)
            g2_out_nc=array(g2_out_nc)
            average_m1_corr=array(average_m1_corr)
            average_m2_corr=array(average_m2_corr)

            average_c1_corr=array(average_c1_corr)
            average_c2_corr=array(average_c2_corr)

            average_m1_bin.append(average(average_m1_corr))
            average_c1_bin.append(average(average_c1_corr))
            average_m2_bin.append(average(average_m2_corr))
            average_c2_bin.append(average(average_c2_corr))
            popt1_m, punc1_m, rc1_m, dc1_m,pcov1_m =general_fit(func,g1_in_used,g1_out,p0,g_out_w)
            popt2_m, punc2_m, rc2_m, dc2_m,pcov2_m =general_fit(func,g2_in_used,g2_out,p0,g_out_w)

            popt1_m_nc, punc1_m_nc, rc1_m_nc, dc1_m_nc,pcov1_m_nc =general_fit(func,g1_in_used,g1_out_nc,p0,g_out_w)
            popt2_m_nc, punc2_m_nc, rc2_m_nc, dc2_m_nc,pcov2_m_nc =general_fit(func,g2_in_used,g2_out_nc,p0,g_out_w)

            m1_m.append(popt1_m[0])
            m1_m_err.append(punc1_m[0])

            c1_m.append(popt1_m[1])
            c1_m_err.append(punc1_m[1])

            m2_m.append(popt2_m[0])
            m2_m_err.append(punc2_m[0])

            c2_m.append(popt2_m[1])
            c2_m_err.append(punc2_m[1])

            m1_m_nc.append(popt1_m_nc[0])
            m1_m_err_nc.append(punc1_m_nc[0])

            c1_m_nc.append(popt1_m_nc[1])
            c1_m_err_nc.append(punc1_m_nc[1])

            m2_m_nc.append(popt2_m_nc[0])
            m2_m_err_nc.append(punc2_m_nc[0])

            c2_m_nc.append(popt2_m_nc[1])
            c2_m_err_nc.append(punc2_m_nc[1])

            #Performing Bootstrap
            #m1_sample=[]
            #m2_sample=[]
            #c1_sample=[]
            #c2_sample=[]
            #m1_sample_nc=[]
            #m2_sample_nc=[]
            #c1_sample_nc=[]
            #c2_sample_nc=[]
            #for BS_index in range(0,20):
                #BS_g1_in=numpy.empty(numShear)
                #BS_g1_out=numpy.empty(numShear)
                #BS_g1_out_nc=numpy.empty(numShear)
                #BS_g2_in=numpy.empty(numShear)
                #BS_g2_out=numpy.empty(numShear)
                #BS_g2_out_nc=numpy.empty(numShear)
                #BS_g_w=numpy.empty(numShear)
                #Retrieving random shears
                #BSs=numpy.random.randint(0,numShear,numShear)
                #for index in range(0,numShear):
                    #BS_g1_in[index]=g1_in_used[BSs[index]]
                    #BS_g2_in[index]=g2_in_used[BSs[index]]
                    #BS_g1_out[index]=g1_out[BSs[index]]
                    #BS_g2_out[index]=g2_out[BSs[index]]
                    #BS_g1_out_nc[index]=g1_out_nc[BSs[index]]
                    #BS_g2_out_nc[index]=g2_out_nc[BSs[index]]
                    #BS_g_w[index]=g_out_w[BSs[index]]
                #BS_g1_in=array(BS_g1_in)
                #BS_g1_out=array(BS_g1_out)
                #BS_g1_out_nc=array(BS_g1_out_nc)
                #BS_g2_in=array(BS_g2_in)
                #BS_g2_out=array(BS_g2_out)
                #BS_g2_out_nc=array(BS_g2_out_nc)
                #BS_g_w=array(BS_g_w)

                #popt1_m=general_fitBS(func,BS_g1_in,BS_g1_out,p0,BS_g_w)
                #popt2_m=general_fitBS(func,BS_g2_in,BS_g2_out,p0,BS_g_w)
                #popt1_m_nc=general_fitBS(func,BS_g1_in,BS_g1_out_nc,p0,BS_g_w)
                #popt2_m_nc=general_fitBS(func,BS_g2_in,BS_g2_out_nc,p0,BS_g_w)
                #m1_sample.append(popt1_m[0])
                #m2_sample.append(popt2_m[0])
                #c1_sample.append(popt1_m[1])
                #c2_sample.append(popt2_m[1])
                #m1_sample_nc.append(popt1_m_nc[0])
                #m2_sample_nc.append(popt2_m_nc[0])
                #c1_sample_nc.append(popt1_m_nc[1])
                #c2_sample_nc.append(popt2_m_nc[1])
            #m1_sample=array(m1_sample)
            #m2_sample=array(m2_sample)
            #c1_sample=array(c1_sample)
            #c2_sample=array(c2_sample)
            #m1_sample_nc=array(m1_sample_nc)
            #m2_sample_nc=array(m2_sample_nc)
            #c1_sample_nc=array(c1_sample_nc)
            #c2_sample_nc=array(c2_sample_nc)

            m1_m_err_BS.append(1.)
            m2_m_err_BS.append(1.)
            c1_m_err_BS.append(1.)
            c2_m_err_BS.append(1.)
            m1_m_err_BS_nc.append(1.)
            m2_m_err_BS_nc.append(1.)
            c1_m_err_BS_nc.append(1.)
            c2_m_err_BS_nc.append(1.)
                
        
    ObsBin_mean=array(ObsBin_mean)
    m1_m=array(m1_m)
    c1_m=array(c1_m)

    m1_m_err=array(m1_m_err)
    m1_m_err_BS=array(m1_m_err_BS)
    c1_m_err=array(c1_m_err)
    c1_m_err_BS=array(c1_m_err_BS)

    m2_m=array(m2_m)
    c2_m=array(c2_m)

    m2_m_err=array(m2_m_err)
    m2_m_err_BS=array(m2_m_err_BS)
    c2_m_err=array(c2_m_err)

    m1_m_nc=array(m1_m_nc)
    c1_m_nc=array(c1_m_nc)

    m1_m_err_nc=array(m1_m_err_nc)
    m1_m_err_BS_nc=array(m1_m_err_BS_nc)
    c1_m_err_nc=array(c1_m_err_nc)
    c1_m_err_BS_nc=array(c1_m_err_BS_nc)

    m2_m_nc=array(m2_m_nc)
    c2_m_nc=array(c2_m_nc)

    m2_m_err_nc=array(m2_m_err_nc)
    m2_m_err_BS_nc=array(m2_m_err_BS_nc)
    c2_m_err_nc=array(c2_m_err_nc)
    c2_m_err_BS_nc=array(c2_m_err_BS_nc)
    average_m1_bin=array(average_m1_bin)
    average_m2_bin=array(average_m2_bin)

    average_c1_bin=array(average_c1_bin)
    average_c2_bin=array(average_c2_bin)

    print 'I am saving the results here: ', Dir+'/Results/1bin/'+Flag+'_'+str(Psf_set)+'_'+str(subSet)+'_'+fileOut
    fOut=open(Dir+'/Results/1bin/'+Flag+'_'+str(Psf_set)+'_'+str(subSet)+'_'+fileOut,'wb')

    BinWeight=array(BinWeight)
    if(len(m1_m)>0):
        print >> fOut, "#obs_mean" ,"PSFquant", "m1","m1_err","c1","c1_err", "m2","m2_err","c2","c2_err","m1_nc","m1_err_nc","c1_nc","c1_err_nc", "m2_nc","m2_err_nc","c2_nc","c2_err_nc", "wBin/wTot", "psf_e1", "psf_e2", "psf_size", "m1_err_BS", "m2_err_BS", "c1_err_BS", "c2_err_BS", "m1_err_BS_nc", "m2_err_BS_nc", "c1_err_BS_nc", "c2_err_BS_nc", "Nobj", "m1_corr", "m2_corr", "c1_corr", "c2_corr"
        for x in range(0,len(m1_m)):
            print >> fOut, ObsBin_mean[x],strehl_set,m1_m[x],m1_m_err[x],c1_m[x],c1_m_err[x], m2_m[x],m2_m_err[x],c2_m[x],c2_m_err[x],m1_m_nc[x],m1_m_err_nc[x],c1_m_nc[x],c1_m_err_nc[x], m2_m_nc[x],m2_m_err_nc[x],c2_m_nc[x],c2_m_err_nc[x], BinWeight[x]/TotWeight, average(psf_e1_sim_weight),average(psf_e2_sim_weight), average(psf_size_sim_weight), m1_m_err_BS[x], m2_m_err_BS[x], c1_m_err_BS[x], c2_m_err_BS[x], m1_m_err_BS_nc[x], m2_m_err_BS_nc[x], c1_m_err_BS_nc[x], c2_m_err_BS_nc[x], BinCount[x], average_m1_bin[x],average_m2_bin[x], average_c1_bin[x], average_c2_bin[x]


        fac=1.0#BinWeight/TotWeight

        fig=plt.figure(figsize=(10,8))
        fig.add_subplot(2,2,1)

        plt.xlabel(ObsLabel)
        plt.ylabel('m1')
        if((Observable=="resolution")or(Observable=="size_in")or(Observable=="snr")):
            plt.semilogx()
        plt.errorbar(ObsBin_mean,fac*m1_m, yerr=fac*m1_m_err,fmt='o', label='selfcal')
        plt.errorbar(ObsBin_mean,fac*m1_m_nc, yerr=fac*m1_m_err_nc,fmt='o', label='no selfcal')
        plt.axhline(0,linewidth=2,color='black')
        plt.legend(loc=0)
        fig.add_subplot(2,2,2)
        plt.xlabel(ObsLabel)
        plt.ylabel('c1')
        if((Observable=="resolution")or(Observable=="size_in")or(Observable=="snr")):
            plt.semilogx()
        plt.errorbar(ObsBin_mean,fac*c1_m, yerr=fac*c1_m_err,fmt='o', label='selfcal')
        plt.errorbar(ObsBin_mean,fac*c1_m_nc, yerr=fac*c1_m_err_nc,fmt='o', label='no selfcal')
        plt.axhline(0,linewidth=2,color='black')
        #plt.legend()
        fig.add_subplot(2,2,3)
        plt.xlabel(ObsLabel)
        plt.ylabel('m2')
        if((Observable=="resolution")or(Observable=="size_in")or(Observable=="snr")):
            plt.semilogx()
        plt.errorbar(ObsBin_mean,fac*m2_m, yerr=fac*m2_m_err,fmt='o', label='selfcal')
        plt.errorbar(ObsBin_mean,fac*m2_m_nc, yerr=fac*m2_m_err_nc,fmt='o', label='no selfcal')
        plt.axhline(0,linewidth=2,color='black')
        #plt.legend()
        fig.add_subplot(2,2,4)
        plt.xlabel(ObsLabel)
        plt.ylabel('c2')
        if((Observable=="resolution")or(Observable=="size_in")or(Observable=="snr")):
            plt.semilogx()
        plt.errorbar(ObsBin_mean,fac*c2_m, yerr=fac*c2_m_err,fmt='o', label='selfcal')
        plt.errorbar(ObsBin_mean,fac*c2_m_nc, yerr=fac*c2_m_err_nc,fmt='o', label='no selfcal')
        #plt.legend()
        plt.axhline(0,linewidth=2,color='black')
        plt.tight_layout()
        plt.savefig(Dir+'/Results/1bin/'+Flag+'_'+str(Psf_set)+'_'+str(subSet)+'_'+PlotName)

        
    ######################

def DoAll2Bins(Observable,bins, obs1_sim_weight, obs2_sim_weight, mag_in, size_in, mag,  w_sim_mask, e1_sim_mask, e2_sim_mask, e1_in_sim_weight, e2_in_sim_weight, e1_corr_sim_weight,e2_corr_sim_weight, g1_in_mask,g2_in_mask, g1_theo, g2_theo, fileOut, Flag, PlotName, ObsLabel_1, ObsLabel_2,psf_e1_sim_weight,psf_e2_sim_weight,psf_size_sim_weight, Dir,TotWeight,Psf_set,subSet,Nbin,strehl_set,ending_flag,snr=None,res=None,m1=None,m2=None,c1=None,c2=None,m_func=None,m1_params=None,m2_params=None,c_func=None,c1_params=None,c2_params=None,m1_splines=None,m2_splines=None,c1_splines=None,c2_splines=None,spline_thresholds=None, m1_table=None, m2_table=None, c1_table=None, c2_table=None, x_shift=None, y_shift=None, x_scale=None, y_scale=None, tree=None,table_data=None):

    cal=False
    spline_cal=False
    spline_cal2=False
    spline_cal_multi=False
    table_cal=False
    bin_cal=False
    file_cal=False
    if(snr!=None):
        cal=True
        if(m1_splines!=None):
            spline_cal=True
            if(spline_thresholds!=None):
                spline_cal_multi=True
            elif(len(m1_splines)>1):
                spline_cal2=True
        elif(tree!=None):
            table_cal=True
        elif(table_data!=None):
            bin_cal=True
        elif(m1!=None):
            file_cal=True



    if(ending_flag=='nocorr'):
        if(cal==False):
            fOutShear=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear.dat','wb')
        elif(spline_cal==False):
            if(table_cal==True):
                fLog=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_table_cal.log','wb')
                fOutShear=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'shear_table_cal.dat','wb')
            elif(bin_cal==True):
                fLog=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_bin_cal.log','wb')
                fOutShear=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'shear_bin_cal.dat','wb')
            elif(file_cal==True):
                fLog=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_file_cal.log','wb')
                fOutShear=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'shear_file_cal.dat','wb')
            else:
                fLog=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_cal.log','wb')
                fOutShear=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear_cal.dat','wb')
        elif(spline_cal==True):
            if(spline_cal2==True):
                fOutShear=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear_2spline_cal.dat','wb')
            else:
                fOutShear=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear_spline_cal.dat','wb')


    else:
        if(cal==False):
            fOutShear=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear.'+ending_flag+'.dat','wb')
        elif(spline_cal==False):
            if(table_cal==True):
                fLog=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'.'+ending_flag+'_table_cal.log','wb')
                fOutShear=open(Dir+'/Results/1bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear.'+ending_flag+'_table_cal.dat','wb')
            elif(bin_cal==True):
                fLog=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'.'+ending_flag+'_bin_cal.log','wb')
                fOutShear=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear.'+ending_flag+'_bin_cal.dat','wb')
            elif(file_cal==True):
                fLog=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'.'+ending_flag+'_file_cal.log','wb')
                fOutShear=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear.'+ending_flag+'_file_cal.dat','wb')
            else:
                fLog=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'.'+ending_flag+'_cal.log','wb')
                fOutShear=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear.'+ending_flag+'_cal.dat','wb')
        elif(spline_cal==True):
            if(spline_cal2==True):
                fOutShear=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear.'+ending_flag+'_2spline_cal.dat','wb')
            elif(spline_cal_multi==True):
                fOutShear=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear.'+ending_flag+'_multi_spline_cal.dat','wb')
            else:
                fOutShear=open(Dir+'/Results/2bin/'+Flag+'_'+Observable+'_'+str(Psf_set)+'_'+str(subSet)+'_shear.'+ending_flag+'_spline_cal.dat','wb')


    print >> fOutShear, "#Median bin1", "bin1_lower","bin1_upper",  "Median bin2","bin2_lower","bin2_upper" "Median strehl", "g1_in","g2_in","g1_out","g2_out","g1_out_nc","g2_out_nc","g1_out_err","g2_out_err","g1_out_err_nc","g2_out_err_nc"
    ObsBin1_mean=[]
    ObsBin2_mean=[]
    ObsBin1_lower=[]
    ObsBin1_upper=[]
    ObsBin2_lower=[]
    ObsBin2_upper=[]
    avg_mag_in=[]
    avg_size_in=[]
    avg_mag_out=[]
    BinWeight=[]
    BinDim=[]
    p0=[0,0]

    m1_m=[]
    m1_m_err=[]
    m1_m_err_BS=[]
    c1_m=[]
    c1_m_err=[]
    c1_m_err_BS=[]
    m2_m=[]
    m2_m_err=[]
    m2_m_err_BS=[]
    c2_m=[]
    c2_m_err=[]
    c2_m_err_BS=[]

    m1_m_nc=[]
    m1_m_err_nc=[]
    m1_m_err_BS_nc=[]
    c1_m_nc=[]
    c1_m_err_nc=[]
    c1_m_err_BS_nc=[]
    m2_m_nc=[]
    m2_m_err_nc=[]
    m2_m_err_BS_nc=[]
    c2_m_nc=[]
    c2_m_err_nc=[]
    c2_m_err_BS_nc=[]

    x_bin=bins[0]
    y_bin=bins[1]

    
    for x in range(0,len(x_bin[0])-1):
        for j in range(0,len(y_bin[x])-1):
            if((Observable=="mag_nn")or(Observable=="size_nn")):
                maskBinObs=((obs1_sim_weight>=x_bin[0][x])&(obs1_sim_weight<x_bin[0][x+1])&(obs2_sim_weight>y_bin[x][j])&(obs2_sim_weight<y_bin[x][j+1]))
            else:
                maskBinObs=((obs1_sim_weight>=x_bin[0][x])&(obs1_sim_weight<x_bin[0][x+1])&(obs2_sim_weight>=y_bin[x][j])&(obs2_sim_weight<y_bin[x][j+1]))
            #ObsBin1_mean.append(median(obs1_sim_weight[maskBinObs]))
            #if((Observable=="mag_nn")or(Observable=="size_nn")):
                #ObsBin2_mean.append(y_bin[x][j])
            #    ObsBin2_mean.append(median(obs2_sim_weight[maskBinObs]))
            #else:
            #    ObsBin2_mean.append(median(obs2_sim_weight[maskBinObs]))
            ObsBin1_lower.append(x_bin[0][x])
            ObsBin1_upper.append(x_bin[0][x+1])
            ObsBin2_lower.append(y_bin[x][j])
            ObsBin2_upper.append(y_bin[x][j+1])

            w_sim_mask_obs=w_sim_mask[maskBinObs]

            ObsBin1_mean.append(average(obs1_sim_weight[maskBinObs],weights=w_sim_mask_obs))
            if((Observable=="mag_nn")or(Observable=="size_nn")):
                #ObsBin2_mean.append(y_bin[x][j])
                ObsBin2_mean.append(average(obs2_sim_weight[maskBinObs],weights=w_sim_mask_obs))
            else:
                ObsBin2_mean.append(average(obs2_sim_weight[maskBinObs],weights=w_sim_mask_obs))
            
            e1_sim_mask_obs=e1_sim_mask[maskBinObs]
            e2_sim_mask_obs=e2_sim_mask[maskBinObs]
            e1_in_sim_mask_obs=e1_in_sim_weight[maskBinObs]
            e2_in_sim_mask_obs=e2_in_sim_weight[maskBinObs]
            #Define the ellipticity by remove self calibration term
            e1_sim_mask_obs_nc=e1_sim_mask[maskBinObs]-e1_corr_sim_weight[maskBinObs]
            e2_sim_mask_obs_nc=e2_sim_mask[maskBinObs]-e2_corr_sim_weight[maskBinObs]
            cal_mag_in=mag_in[maskBinObs]
            cal_size_in=size_in[maskBinObs]
            cal_mag=mag[maskBinObs]
            if(cal==True):
                calibration_snr_sim_mask_obs=snr[maskBinObs]
                calibration_resolution_sim_mask_obs=res[maskBinObs]

            if(file_cal==True):
                cal_m1=m1[maskBinObs]
                cal_m2=m2[maskBinObs]
                cal_c1=c1[maskBinObs]
                cal_c2=c2[maskBinObs]
            

            avg_mag_in.append(average(cal_mag_in,weights=w_sim_mask_obs))
            avg_size_in.append(average(cal_size_in,weights=w_sim_mask_obs))
            avg_mag_out.append(average(cal_mag,weights=w_sim_mask_obs))
            #Compute weight in the bin

            BinWeight.append(sum(w_sim_mask_obs))
            BinDim.append(len(w_sim_mask_obs))
            g1_in_mask_obs=g1_in_mask[maskBinObs]
            g2_in_mask_obs=g2_in_mask[maskBinObs]
    
            g1_out=[]
            g2_out=[]
            g1_out_nc=[]
            g2_out_nc=[]
            g_out_w=[]
    
            g1_in_used=[]
            g2_in_used=[]
            for kk in range(0, len(g1_theo)):
                maskShear=(g1_in_mask_obs==g1_theo[kk])&(g2_in_mask_obs==g2_theo[kk])
                numMasked=len(e1_sim_mask_obs[maskShear])
                if (numMasked >0):
                    #Calculating bin average for calibration quantities
                    if(cal==True):

                        xy_tuple_correction=(calibration_snr_sim_mask_obs[maskShear],calibration_resolution_sim_mask_obs[maskShear])
                        if(table_cal==True):
                            xy=array(xy_tuple_correction)
                            m1_sample=[]
                            m2_sample=[]
                            c1_sample=[]
                            c2_sample=[]
                            weight_sample=w_sim_mask_obs[maskShear]
                            for x_quant in xy.T:
                                all_biases = bias_2D_table_query_all(x_quant[0], x_quant[1],m1_table,m2_table,c1_table,c2_table,tree,x_scale,y_scale,x_shift,y_shift)
                                m1_sample.append(all_biases[0])
                                m2_sample.append(all_biases[1])
                            m1_sample=array(m1_sample)
                            m2_sample=array(m2_sample)

                            current_m1=average(m1_sample,weights=weight_sample)
                            current_m2=average(m2_sample,weights=weight_sample)
                            #I NEEDED TO KEEP THIS LINE IN LIKE THIS SINCE SPLINES CANNOT BE CHOSEN TO GIVE 0
                            current_c1=0.0#average(c1_sample,weights=weight_sample)
                            current_c2=0.0#average(c2_sample,weights=weight_sample)
                        elif(bin_cal==True):
                            xy=array(xy_tuple_correction)
                            m1_sample=[]
                            m2_sample=[]
                            c1_sample=[]
                            c2_sample=[]
                            weight_sample=w_sim_mask_obs[maskShear]
                            for x_quant in xy.T:
                                all_biases = bias_2D_bin_query(x_quant[0], x_quant[1],table_data)
                                m1_sample.append(all_biases[0])
                                m2_sample.append(all_biases[1])
                            m1_sample=array(m1_sample)
                            m2_sample=array(m2_sample)

                            current_m1=average(m1_sample,weights=weight_sample)
                            current_m2=average(m2_sample,weights=weight_sample)
                            #I NEEDED TO KEEP THIS LINE IN LIKE THIS SINCE SPLINES CANNOT BE CHOSEN TO GIVE 0
                            current_c1=0.0#average(c1_sample,weights=weight_sample)
                            current_c2=0.0#average(c2_sample,weights=weight_sample)

                        elif(spline_cal==True):
                            if(spline_cal2==True):
                                xy=array(xy_tuple_correction)
                                m1_sample=[]
                                m2_sample=[]
                                c1_sample=[]
                                c2_sample=[]
                                weight_sample=w_sim_mask_obs[maskShear]
                                for x_quant in xy.T:
                                    m1_sample.append(float(query_spline2_emp(x_quant,m1_splines[0],m1_splines[1])))
                                    m2_sample.append(float(query_spline2_emp(x_quant,m2_splines[0],m2_splines[1])))
                                    #c1_sample.append(float(query_spline2_emp(x_quant,c1_spline,c1_spline2)))
                                    #c2_sample.append(float(query_spline2_emp(x_quant,c2_spline,c2_spline2)))
                                m1_sample=array(m1_sample)
                                m2_sample=array(m2_sample)
                                #c1_sample=array(c1_sample)
                                #c2_sample=array(c2_sample)

                            elif(spline_cal_multi==True):
                                xy=array(xy_tuple_correction)
                                m1_sample=[]
                                m2_sample=[]
                                c1_sample=[]
                                c2_sample=[]
                                weight_sample=w_sim_mask_obs[maskShear]
                                for x_quant in xy.T:
                                    (one,two) = x_quant
                                    m1_sample.append(float(query_multi_bin_splines(one,two,m1_splines,spline_thresholds)))
                                    m2_sample.append(float(query_multi_bin_splines(one,two,m2_splines,spline_thresholds)))
                                    #c1_sample.append(float(query_spline2_emp(x_quant,c1_spline,c1_spline2)))
                                    #c2_sample.append(float(query_spline2_emp(x_quant,c2_spline,c2_spline2)))
                                m1_sample=array(m1_sample)
                                m2_sample=array(m2_sample)
                                #c1_sample=array(c1_sample)
                                #c2_sample=array(c2_sample)

                            else:
                                xy=array(xy_tuple_correction)
                                m1_sample=[]
                                m2_sample=[]
                                c1_sample=[]
                                c2_sample=[]
                                weight_sample=w_sim_mask_obs[maskShear]
                                for x_quant in xy.T:
                                    m1_sample.append(float(query_spline_emp(x_quant,m1_splines[0])))
                                    m2_sample.append(float(query_spline_emp(x_quant,m2_splines[0])))
                                    #c1_sample.append(float(query_spline_emp(x_quant,c1_spline)))
                                    #c2_sample.append(float(query_spline_emp(x_quant,c2_spline)))
                                m1_sample=array(m1_sample)
                                m2_sample=array(m2_sample)
                                #c1_sample=array(c1_sample)
                                #c2_sample=array(c2_sample)
                                
                            current_m1=average(m1_sample,weights=weight_sample)
                            current_m2=average(m2_sample,weights=weight_sample)
                            #I NEEDED TO KEEP THIS LINE IN LIKE THIS SINCE SPLINES CANNOT BE CHOSEN TO GIVE 0
                            current_c1=0.0#average(c1_sample,weights=weight_sample)
                            current_c2=0.0#average(c2_sample,weights=weight_sample)

                        elif(file_cal==True):
                            current_m1=average(cal_m1[maskShear],weights=w_sim_mask_obs[maskShear])
                            current_m2=average(cal_m2[maskShear],weights=w_sim_mask_obs[maskShear])
                            current_c1=average(cal_c1[maskShear],weights=w_sim_mask_obs[maskShear]) 
                            current_c2=average(cal_c2[maskShear],weights=w_sim_mask_obs[maskShear])

                        else:
                            current_m1=average(m_func(xy_tuple_correction,*m1_params),weights=w_sim_mask_obs[maskShear])
                            current_m2=average(m_func(xy_tuple_correction,*m2_params),weights=w_sim_mask_obs[maskShear])
                            current_c1=average(c_func(xy_tuple_correction,*c1_params),weights=w_sim_mask_obs[maskShear])
                            current_c2=average(c_func(xy_tuple_correction,*c2_params),weights=w_sim_mask_obs[maskShear])


                    if(Flag=='CH'):
                        if(cal==False):
                            g1_out.append(average(e1_sim_mask_obs[maskShear]-e1_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear]))
                            g2_out.append(average(e2_sim_mask_obs[maskShear]-e2_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear]))
                            g_out_w.append(1./(sum(w_sim_mask_obs[maskShear]))**0.5)
                            g1_out_nc.append(average(e1_sim_mask_obs_nc[maskShear]-e1_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear]))
                            g2_out_nc.append(average(e2_sim_mask_obs_nc[maskShear]-e2_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear]))
                            g1_in_used.append(g1_theo[kk])
                            g2_in_used.append(g2_theo[kk])
                            
                            g1_out_print=average(e1_sim_mask_obs[maskShear]-e1_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                            g1_out_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                            g2_out_print=average(e2_sim_mask_obs[maskShear]-e2_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                            g2_out_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                            g1_out_nc_print=average(e1_sim_mask_obs_nc[maskShear]-e1_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                            g1_out_nc_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                            g2_out_nc_print=average(e2_sim_mask_obs_nc[maskShear]-e2_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                            g2_out_nc_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                            g1_in_print=g1_theo[kk]
                            g2_in_print=g2_theo[kk]
                        else:
                            value1=average(e1_sim_mask_obs[maskShear]-e1_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                            value2=average(e2_sim_mask_obs[maskShear]-e2_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                            value1_w=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                            value1_nc=average(e1_sim_mask_obs_nc[maskShear]-e1_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                            value2_nc=average(e2_sim_mask_obs_nc[maskShear]-e2_in_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                            value1_c=(value1-current_c1)/(1.+current_m1)
                            value2_c=(value2-current_c2)/(1.+current_m2)
                            value1_nc_c=(value1_nc-current_c1)/(1.+current_m1)
                            value2_nc_c=(value2_nc-current_c2)/(1.+current_m2)
                            g1_out.append(value1_c)
                            g1_out_print=value1_c
                            g1_out_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                            g_out_w.append(value1_w)
                            g2_out.append(value2_c)
                            g2_out_print=value2_c
                            g2_out_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                            g1_out_nc.append(value1_nc_c)
                            g2_out_nc.append(value2_nc_c)
                            g1_out_nc_print=value1_nc_c
                            g1_out_nc_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                            g2_out_nc_print=value2_nc_c
                            g2_out_nc_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                            g1_in_used.append(g1_theo[kk])
                            g2_in_used.append(g2_theo[kk])
                            g1_in_print=g1_theo[kk]
                            g2_in_print=g2_theo[kk]
                        print >> fOutShear, numpy.median(obs1_sim_weight[maskBinObs]),x_bin[0][x],x_bin[0][x+1], numpy.median(obs2_sim_weight[maskBinObs]),y_bin[x][j],y_bin[x][j+1],strehl_set,g1_in_print,g2_in_print,g1_out_print,g2_out_print,g1_out_nc_print,g2_out_nc_print,g1_out_err_print,g2_out_err_print,g1_out_nc_err_print,g2_out_nc_err_print

                    else:
                        if(cal==False):
                            g1_out.append(average(e1_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear]))
                            g2_out.append(average(e2_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear]))
                            g_out_w.append(1./(sum(w_sim_mask_obs[maskShear]))**0.5)
                            g1_out_nc.append(average(e1_sim_mask_obs_nc[maskShear], weights=w_sim_mask_obs[maskShear]))
                            g2_out_nc.append(average(e2_sim_mask_obs_nc[maskShear], weights=w_sim_mask_obs[maskShear]))
                            g1_in_used.append(g1_theo[kk])
                            g2_in_used.append(g2_theo[kk])
                            g1_out_print=average(e1_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                            g1_out_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                            g2_out_print=average(e2_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                            g2_out_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                            g1_out_nc_print=average(e1_sim_mask_obs_nc[maskShear], weights=w_sim_mask_obs[maskShear])
                            g1_out_nc_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                            g2_out_nc_print=average(e2_sim_mask_obs_nc[maskShear], weights=w_sim_mask_obs[maskShear])
                            g2_out_nc_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                            g1_in_print=g1_theo[kk]
                            g2_in_print=g2_theo[kk]
                        else:
                            value1=average(e1_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                            value2=average(e2_sim_mask_obs[maskShear], weights=w_sim_mask_obs[maskShear])
                            value1_w=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                            value1_nc=average(e1_sim_mask_obs_nc[maskShear], weights=w_sim_mask_obs[maskShear])
                            value2_nc=average(e2_sim_mask_obs_nc[maskShear], weights=w_sim_mask_obs[maskShear])
                            value1_c=(value1-current_c1)/(1.+current_m1)
                            value2_c=(value2-current_c2)/(1.+current_m2)
                            value1_nc_c=(value1_nc-current_c1)/(1.+current_m1)
                            value2_nc_c=(value2_nc-current_c2)/(1.+current_m2)
                            g1_out.append(value1_c)
                            g1_out_print=value1_c
                            g1_out_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                            g_out_w.append(value1_w)
                            g2_out.append(value2_c)
                            g2_out_print=value2_c
                            g2_out_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                            g1_out_nc.append(value1_nc_c)
                            g2_out_nc.append(value2_nc_c)
                            g1_out_nc_print=value1_nc_c
                            g1_out_nc_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                            g2_out_nc_print=value2_nc_c
                            g2_out_nc_err_print=1./(sum(w_sim_mask_obs[maskShear]))**0.5
                            g1_in_used.append(g1_theo[kk])
                            g2_in_used.append(g2_theo[kk])
                            g1_in_print=g1_theo[kk]
                            g2_in_print=g2_theo[kk]
                        print >> fOutShear,numpy.median(obs1_sim_weight[maskBinObs]),x_bin[0][x],x_bin[0][x+1], numpy.median(obs2_sim_weight[maskBinObs]),y_bin[x][j],y_bin[x][j+1], strehl_set,g1_in_print,g2_in_print,g1_out_print,g2_out_print,g1_out_nc_print,g2_out_nc_print,g1_out_err_print,g2_out_err_print,g1_out_nc_err_print,g2_out_nc_err_print


            #print g1_theo[kk], average(e1_sim_mask_obs[maskShear], weights=w_sim_maskobs[maskShear])
            # Check I have at least 3 shear values to do the regression
            numShear=len(g1_out)

            if(numShear<3):
                print 'Cannot do the regression in bin ', x_bin[0][x], x_bin[0][x+1], y_bin[x][j],y_bin[x][j+1], ' less than 3 shear values! (', numShear, ')'
                exit()
            else:
                #print 'I am using ', numShear, 'values of the shear'
                g1_in_used=array(g1_in_used)
                g2_in_used=array(g2_in_used)
                g1_out=array(g1_out)
                g2_out=array(g2_out)
                g1_out_nc=array(g1_out_nc)
                g2_out_nc=array(g2_out_nc)
                g_out_w=array(g_out_w)

                popt1_m, punc1_m, rc1_m, dc1_m,pcov1_m =general_fit(func,g1_in_used,g1_out,p0,g_out_w)
                popt2_m, punc2_m, rc2_m, dc2_m,pcov2_m =general_fit(func,g2_in_used,g2_out,p0,g_out_w)

                popt1_m_nc, punc1_m_nc, rc1_m_nc, dc1_m_nc,pcov1_m_nc =general_fit(func,g1_in_used,g1_out_nc,p0,g_out_w)
                popt2_m_nc, punc2_m_nc, rc2_m_nc, dc2_m_nc,pcov2_m_nc =general_fit(func,g2_in_used,g2_out_nc,p0,g_out_w)

                m1_m.append(popt1_m[0])
                m1_m_err.append(punc1_m[0])

                c1_m.append(popt1_m[1])
                c1_m_err.append(punc1_m[1])

                m2_m.append(popt2_m[0])
                m2_m_err.append(punc2_m[0])

                c2_m.append(popt2_m[1])
                c2_m_err.append(punc2_m[1])

                m1_m_nc.append(popt1_m_nc[0])
                m1_m_err_nc.append(punc1_m_nc[0])

                c1_m_nc.append(popt1_m_nc[1])
                c1_m_err_nc.append(punc1_m_nc[1])

                m2_m_nc.append(popt2_m_nc[0])
                m2_m_err_nc.append(punc2_m_nc[0])

                c2_m_nc.append(popt2_m_nc[1])
                c2_m_err_nc.append(punc2_m_nc[1])



                #Performing Bootstrap
                m1_sample=[]
                m2_sample=[]
                c1_sample=[]
                c2_sample=[]
                m1_sample_nc=[]
                m2_sample_nc=[]
                c1_sample_nc=[]
                c2_sample_nc=[]
                for BS_index in range(0,20):
                    BS_g1_in=empty(numShear)
                    BS_g1_out=empty(numShear)
                    BS_g1_out_nc=empty(numShear)
                    BS_g2_in=empty(numShear)
                    BS_g2_out=empty(numShear)
                    BS_g2_out_nc=empty(numShear)
                    BS_g_out_w=empty(numShear)
                    #Retrieving random shears
                    BSs=numpy.random.randint(0,numShear,numShear)
                    for index in range(0,numShear):
                        BS_g1_in[index]=g1_in_used[BSs[index]]
                        BS_g2_in[index]=g2_in_used[BSs[index]]
                        BS_g1_out[index]=g1_out[BSs[index]]
                        BS_g2_out[index]=g2_out[BSs[index]]
                        BS_g1_out_nc[index]=g1_out_nc[BSs[index]]
                        BS_g2_out_nc[index]=g2_out_nc[BSs[index]]
                        BS_g_out_w[index]=g_out_w[BSs[index]]
                    BS_g1_in=array(BS_g1_in) 
                    BS_g2_in=array(BS_g2_in) 
                    BS_g1_out=array(BS_g1_out) 
                    BS_g2_out=array(BS_g2_out) 
                    BS_g1_out_nc=array(BS_g1_out_nc) 
                    BS_g2_out_nc=array(BS_g2_out_nc) 
                    BS_g_out_w=array(BS_g_out_w)
                        
                    popt1_m=general_fitBS(func,BS_g1_in,BS_g1_out,p0,BS_g_out_w)
                    popt2_m=general_fitBS(func,BS_g2_in,BS_g2_out,p0,BS_g_out_w)
                    popt1_m_nc=general_fitBS(func,BS_g1_in,BS_g1_out_nc,p0)
                    popt2_m_nc=general_fitBS(func,BS_g2_in,BS_g2_out_nc,p0)
                    m1_sample.append(popt1_m[0])
                    m2_sample.append(popt2_m[0])
                    c1_sample.append(popt1_m[1])
                    c2_sample.append(popt2_m[1])
                    m1_sample_nc.append(popt1_m_nc[0])
                    m2_sample_nc.append(popt2_m_nc[0])
                    c1_sample_nc.append(popt1_m_nc[1])
                    c2_sample_nc.append(popt2_m_nc[1])
                m1_sample=array(m1_sample)
                m2_sample=array(m2_sample)
                c1_sample=array(c1_sample)
                c2_sample=array(c2_sample)
                m1_sample_nc=array(m1_sample_nc)
                m2_sample_nc=array(m2_sample_nc)
                c1_sample_nc=array(c1_sample_nc)
                c2_sample_nc=array(c2_sample_nc)
                        
                m1_m_err_BS.append(std(m1_sample))
                m2_m_err_BS.append(std(m2_sample))
                c1_m_err_BS.append(std(c1_sample))
                c2_m_err_BS.append(std(c2_sample))
                m1_m_err_BS_nc.append(std(m1_sample_nc))
                m2_m_err_BS_nc.append(std(m2_sample_nc))
                c1_m_err_BS_nc.append(std(c1_sample_nc))
                c2_m_err_BS_nc.append(std(c2_sample_nc))




        
    ObsBin1_mean=array(ObsBin1_mean)
    ObsBin2_mean=array(ObsBin2_mean)
    ObsBin1_lower=array(ObsBin1_lower)
    ObsBin1_upper=array(ObsBin1_upper)
    ObsBin2_lower=array(ObsBin2_lower)
    ObsBin2_upper=array(ObsBin2_upper)

    m1_m=array(m1_m)
    c1_m=array(c1_m)

    m1_m_err=array(m1_m_err)
    c1_m_err=array(c1_m_err)

    m2_m=array(m2_m)
    c2_m=array(c2_m)

    m2_m_err=array(m2_m_err)
    c2_m_err=array(c2_m_err)

    m1_m_nc=array(m1_m_nc)
    c1_m_nc=array(c1_m_nc)

    m1_m_err_nc=array(m1_m_err_nc)
    c1_m_err_nc=array(c1_m_err_nc)

    m2_m_nc=array(m2_m_nc)
    c2_m_nc=array(c2_m_nc)

    m2_m_err_nc=array(m2_m_err_nc)
    c2_m_err_nc=array(c2_m_err_nc)

    fOut=open(Dir+'/Results/2bin/'+Flag+'_'+str(Psf_set)+'_'+str(subSet)+'_'+fileOut,'wb')

    BinWeight=array(BinWeight)
    avg_mag_in = array(avg_mag_in)
    avg_size_in = array(avg_size_in)
    avg_mag_out = array(avg_mag_out)
    print >> fOut, "#obs1_mean", "obs1_lower","obs1_upper","obs2_mean","obs2_lower","obs2_upper","strehl_set","m1","m1_err","c1","c1_err", "m2","m2_err","c2","c2_err","m1_nc","m1_err_nc","c1_nc","c1_err_nc", "m2_nc","m2_err_nc","c2_nc","c2_err_nc", "wBin/wTot", "psf_e1", "psf_e2", "psf_size", "Nobj", "m1_err_BS", "m2_err_BS", "c1_err_BS", "c2_err_BS", "m1_err_BS_nc", "m2_err_BS_nc", "c1_err_BS_nc", "c2_err_BS_nc", "mag_in_avg","size_in_avg", "mag_out_avg"

    for x in range(0,len(m1_m)):
        print >> fOut, ObsBin1_mean[x],ObsBin1_lower[x],ObsBin1_upper[x], ObsBin2_mean[x],ObsBin2_lower[x],ObsBin2_upper[x],strehl_set,m1_m[x],m1_m_err[x],c1_m[x],c1_m_err[x], m2_m[x],m2_m_err[x],c2_m[x],c2_m_err[x],m1_m_nc[x],m1_m_err_nc[x],c1_m_nc[x],c1_m_err_nc[x], m2_m_nc[x],m2_m_err_nc[x],c2_m_nc[x],c2_m_err_nc[x], BinWeight[x]/TotWeight, average(psf_e1_sim_weight),average(psf_e2_sim_weight), average(psf_size_sim_weight), BinDim[x], m1_m_err_BS[x], m2_m_err_BS[x], c1_m_err_BS[x], c2_m_err_BS[x], m1_m_err_BS_nc[x], m2_m_err_BS_nc[x], c1_m_err_BS_nc[x], c2_m_err_BS_nc[x],avg_mag_in[x], avg_size_in[x], avg_mag_out[x]

#Matches catalogue 2 into catalogue 1 and returns the full kD-tree
#All inputs are simple 1D arrays, number of nearest neighbours is knn
def match_catalogues(x1,y1,x2,y2,knn):

    #Creating tree input formats
    coords1=numpy.dstack([x1,y1])[0]
    coords2=numpy.dstack([x2,y2])[0]

    #Creating tree hull and building the tree
    tree=[]
    Tree=scipy.spatial.cKDTree(coords1,leafsize=128)

    #Querying tree
    for object in coords2:
        tree.append(Tree.query(object,k=knn))

    return tree

#Matches a sextractor and a prior catalogue in a more sophisticated way. 
#First, a naive match of cat2 into cat1 is performed. From the knn nearest 
#neighbours for each match, each with a distance larger than the
#distance threshold is thrown out. This might remove a match completely.
#From the remaning neighbours, the one with the smallest
#difference is the quantity m is chosen. 
#The routine return the indices cat1 and cat2 for the final, remaining matches.
#An offset can be chose which shifts all coordinates in cat2.

def tailored_match_sextractor_grid(x1,y1,m1,x2,y2,m2,knn,threshold,offset):
    #Creating input coordinates for the tree
    x2=x2+offset
    y2=y2+offset
    
    c1=numpy.dstack([x1,y1])[0]
    c2=numpy.dstack([x2,y2])[0]

    #Creating tree hull and building the tree
    tree=[]
    Tree=scipy.spatial.cKDTree(c1,leafsize=128)

    #Querying tree
    for object in c2:
        tree.append(Tree.query(object,k=knn))

    #Creating final output array.
    index_list=[]

    #Going through all matches in the tree
    for i in range(0,len(tree)):
        #Going through each of the nearest neighbours
        mag_sample=[]
        for j in range(0,knn):
            if(tree[i][0][j] < threshold):
                mag_sample.append(abs(m1[tree[i][1][j]]-m2[i]))
        if(len(mag_sample)!=0):
            mag_sample=array(mag_sample)
            element=[tree[i][1][numpy.argmin(mag_sample)],i]
            index_list.append(element)

    index_list=array(index_list)
    return index_list

def new_tailored_match_sextractor_grid(x1,y1,m1,x2,y2,m2,knn,threshold,offset,Tree=None):
    ## Modified by AKJ
    #Creating input coordinates for the tree
    x2=x2+offset
    y2=y2+offset
    
    c1=numpy.dstack([x1,y1])[0]
    c2=numpy.dstack([x2,y2])[0]

    #Creating tree hull and building the tree
    tree=[]
    if Tree is None: ## HACK ALERT
        Tree=scipy.spatial.cKDTree(c1,leafsize=128)

    #Querying tree
    for obj in c2:
        tree.append(Tree.query(obj,k=knn))

    #Creating final output array.
    index_list=[]

    #Going through all matches in the tree
    for i in range(0,len(tree)):
        #Going through each of the nearest neighbours
        mag_sample=[]
        for j in range(0,knn):
            if(tree[i][0][j] < threshold):
                mag_sample.append(abs(m1[tree[i][1][j]]-m2[i]))
        if(len(mag_sample)!=0):
            element=[tree[i][1][numpy.argmin(mag_sample)],i]
        else:
            element=[-1,i] ## could be None
        index_list.append(element)

    index_list=array(index_list)
    return index_list

#This similar to the sextractor matching but is tailored for the ideal
#case that LF is working on a prior grid. We need this function to
#keep the new pipeline consistent and easy to read. 

def tailored_match_prior_grid(x1,y1,x2,y2):

    c1=numpy.dstack([x1,y1])[0]
    c2=numpy.dstack([x2,y2])[0]

    #Creating tree hull and building the tree
    tree=[]
    Tree=scipy.spatial.cKDTree(c1,leafsize=128)

    #Querying tree
    for object in c2:
        tree.append(Tree.query(object,k=1))

    #Creating final output array.
    index_list=[]

    #Going through all matches in the tree
    for i in range(0,len(tree)):
        element=[tree[i][1],i]
        index_list.append(element)
        
    index_list=array(index_list)
    return index_list

def bias_2D_table_query(obs1, obs2, bias_table, table_tree, norm_obs1=1., norm_obs2=1., shift_x = 0., shift_y = 0.):
    obs1 = (obs1-shift_x)/norm_obs1
    obs2 = (obs2-shift_y)/norm_obs2
    index = table_tree.query((obs1,obs2),k=1)

    return bias_table[index[1]]

def bias_2D_table_query_yay(obs1, obs2):
    obs1 = (obs1-shift_x)/norm_obs1
    obs2 = (obs2-shift_y)/norm_obs2
    index = Tree.query((obs1,obs2),k=1)

    return m1[index[1]]

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

    output=array(output)

    return output


def bias_2D_function_query(obs1, obs2,function1,function2,params1,params2):

    value1 = function1((obs1,obs2),*params1)
    value2 = function2((obs1,obs2),*params2)
    output = (value1[0],value2[0],0.0,0.0)
    output=array(output)

    return output



def bias_2D_table_query_all(obs1, obs2, bias_table1,bias_table2,bias_table3,bias_table4, table_tree, norm_obs1=1., norm_obs2=1., shift_x = 0., shift_y = 0.):
    obs1 = (obs1-shift_x)/norm_obs1
    obs2 = (obs2-shift_y)/norm_obs2
    index = table_tree.query((obs1,obs2),k=1)
    #print obs1, obs2, index[1]
    output = (bias_table1[index[1]],bias_table2[index[1]],bias_table3[index[1]],bias_table4[index[1]])
    return array(output)

 




class bias_table:

    def __init__(self,filename):
        data = genfromtxt(filename,comments='#')
        x=data[:,0]
        y=data[:,3]
        self.m1=data[:,7]
        self.c1=data[:,9]
        self.m2=data[:,11]
        self.c2=data[:,13]
        self.shift_x = min(x)
        self.shift_y = min(y)
        x = x - self.shift_x
        y = y - self.shift_y
        self.scale_x = max(x)
        self.scale_y = max(y)
        x=x/self.scale_x
        y=y/self.scale_y
        tree_pos = numpy.dstack([x,y])[0]
        self.Tree=scipy.spatial.cKDTree(tree_pos,leafsize=128)

    def query_m1(self,obs1,obs2):
        obs1 = (obs1-self.shift_x)/self.scale_x
        obs2 = (obs2-self.shift_y)/self.scale_y
        index = self.Tree.query((obs1,obs2),k=1)
        return self.m1[index[1]]


        

    



    



        
