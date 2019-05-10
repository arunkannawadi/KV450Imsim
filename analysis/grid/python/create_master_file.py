import numpy
import pyfits
import sys
import os
import matplotlib
matplotlib.use('Agg')
from MyFunction import *
#import pdb; pdb.set_trace()

Dir = sys.argv[1]
Psf = sys.argv[2]  #Set to -1 if you want all sets into one file
out_base = sys.argv[3]

run_flag = False
if(int(Psf) == -1):
    run_flag = True

if(run_flag):
    add = "_all_sets.fits"
else:
    add = "_set_" +str(Psf) +".fits"

binary_filename = out_base + add


print "Binary FITS files does not exist, creating it."
###########################
list_shear=numpy.genfromtxt(Dir+'/list.t',dtype="|S32")
###########################

#Initialise strings

snr_sim=[]
mode_sim=[]
w_sim=[]
size_sim=[]
nd_sim=[]
nm_sim=[]
fitClass_sim=[]
mag_sim=[]
g1_in=[]
g2_in=[]
e1_sim=[]
e2_sim=[]
strehl_sim=[]
mag_in_sim=[]
size_in_sim=[]
N_in_sim=[]
e1_in_sim=[]
e2_in_sim=[]
psf_size_sim=[]
psf_e1_sim=[]
psf_e2_sim=[]
r_corr_sim=[]
e1_corr_sim=[]
e2_corr_sim=[]
psf_e1_theo=[]
psf_e2_theo=[]
psf_Q11_sim=[]
psf_Q22_sim=[]
psf_Q12_sim=[]
ls_var_sim=[]
rot_sim=[]
ZB_sim=[]
ID_sim=[]

cc=0

#I am starting looping here through the different shear values.

for x in range(0,len(list_shear)): #shear values
    Name_tmp=str(list_shear[x])[0:8]
    if(True): 
        check_string=list_shear[x][10:11]
        if(check_string=='_'):
            Random_string=str(list_shear[x][11:21])
            Psf_type=int(list_shear[x][9:10])
        else:
            Random_string=str(list_shear[x][12:22])
            Psf_type=int(list_shear[x][9:11])
    else:
        Random_string='bbbbbbbb'
        Psf_type=18789
    print Name_tmp, Psf_type, Random_string
    if(Name_tmp[0]=='m'):
        g1_=float(Name_tmp[0:4].replace("m","-0.0"))
    else:
        g1_=float(Name_tmp[0:4].replace("p","0.0"))
    if(Name_tmp[4]=='m'):
        g2_=float(Name_tmp[4:8].replace("m","-0.0"))
    else:
        g2_=float(Name_tmp[4:8].replace("p","0.0"))

    #New main loop over all available sets for each shear value

    for y in range(0,4): #rotation
        cc=cc+1
        if((Psf_type == int(Psf)) or run_flag):
            tbdata_sim = pyfits.getdata(Dir+'/'+Name_tmp+'_'+str(Psf_type)+'_'+Random_string+'/0'+str(y)+'.output.rot.fits.asc.scheme2b_corr.fits')
            #sim=numpy.genfromtxt(Dir+'/'+Name_tmp+'_'+str(Psf_type)+'_'+Random_string+'/0'+str(y)+'.output.rot.fits.asc', comments='#')

            ##MATCHING###
            print 'I am using..', Name_tmp+'_'+str(Psf_type)+'_'+Random_string
            #Checking if LF was working on prior grid or Sextractor catalogue
            grid_check=os.path.isfile(Dir+'/'+str(list_shear[x])+'/sexrot0'+str(y)+'.cat')
            if(grid_check):
                grid_file='/sexrot0'+str(y)+'.cat'
            else:
                #grid_file='/prior.match.asci'
		print 'I cannot find the sextractor catalogue..'
		exit()
            tbdata_grid = pyfits.getdata(Dir+'/'+str(list_shear[x])+grid_file)
            
            #prior_match=numpy.genfromtxt(Dir+'/'+str(list_shear[x])+grid_file, comments='#')
            #prior_match_xpix=prior_match[:,0]
            #prior_match_ypix=prior_match[:,1]
            prior_match_xpix=tbdata_grid.X_IMAGE
            prior_match_ypix=tbdata_grid.Y_IMAGE
            if(grid_check):
                #prior_match_mag=prior_match[:,6]
                prior_match_mag=tbdata_grid.MAG_AUTO
            else:
                prior_match_mag=prior_match[:,4]#??
                
            #Reading simulation prior
            prior_check=os.path.isfile(Dir+'/prior.fits')
            if(prior_check):
                hdu_prior = pyfits.open(Dir+'/prior.fits')
                tbdata_prior=hdu_prior[1].data
                #prior=numpy.genfromtxt(Dir+'/'+str(list_shear[x])+'/prior_new', comments='#')
                print 'Reading: ', Dir+'/prior.fits'
            else:
                print 'I cannot find the prior file!'
                exit()

            #prior_xpix=prior[:,0]
            #prior_ypix=prior[:,1]
            #prior_mag=prior[:,2]
            #prior_size=prior[:,3]
            #prior_fitClass=prior[:,6]

            prior_xpix=tbdata_prior.X_IMAGE#-1627.71
            prior_ypix=tbdata_prior.Y_IMAGE#-646.21
            prior_mag=tbdata_prior.MAG_AUTO
            prior_size=tbdata_prior.RE_GALFIT_HI
	    #Star/gal flag (FLAG==-1: star, FLAG==2: detected in KiDS, FLAG==0: detected in HST only)
            prior_flag=tbdata_prior.FLAG
            prior_q=tbdata_prior.BA_GALFIT_HI
            prior_phi=tbdata_prior.PA_GALFIT_HI#+90.0)*(3.1415/180.)
            prior_N=tbdata_prior.N_GALFIT_HI
            prior_ZB=tbdata_prior.Z_B
            prior_ID=tbdata_prior.SeqNr
            prior_mode=(1.-prior_q)/(1.+prior_q)
            ############################# 

            #Introduce a mask in the prior file (to remove stars, but it can be generalised)
            
            maskPrior=(prior_flag>-99)
            
            #and we redefine new arrays

            prior_xpix_Mask=prior_xpix[maskPrior]
            prior_ypix_Mask=prior_ypix[maskPrior]
            prior_mag_Mask=prior_mag[maskPrior]
            prior_size_Mask=prior_size[maskPrior]
            prior_N_Mask=prior_N[maskPrior]
            prior_ZB_Mask=prior_ZB[maskPrior]
            prior_ID_Mask=prior_ID[maskPrior]
            
            #Note: this is the ellipticity for rotation 00!!!!!!!!
            #prior_e1_00=prior[:,4]
            #prior_e2_00=prior[:,5]
            prior_e1_00=prior_mode*cos(2.*prior_phi)
            prior_e2_00=prior_mode*sin(2.*prior_phi)

            prior_e1_00_Mask=prior_e1_00[maskPrior]
            prior_e2_00_Mask=prior_e2_00[maskPrior]

            #I am rotating here the intrinsic galaxy ellipticities for the 4 runs
            if(y==0):
                prior_e1_Mask=prior_e1_00_Mask
                prior_e2_Mask=prior_e2_00_Mask
            if (y==1):
                prior_e1_Mask=-prior_e2_00_Mask
                prior_e2_Mask=prior_e1_00_Mask
            if(y==2):
                prior_e1_Mask=-prior_e1_00_Mask
                prior_e2_Mask=-prior_e2_00_Mask
            if(y==3):
                prior_e1_Mask=prior_e2_00_Mask
                prior_e2_Mask=-prior_e1_00_Mask

            #Matching the catalogues according to the underlying LF grid
            if(grid_check):
                
                match_map=tailored_match_sextractor_grid(prior_xpix_Mask,prior_ypix_Mask,prior_mag_Mask,prior_match_xpix,prior_match_ypix,prior_match_mag,5,3.0,-0.5)
            else:
                match_map=tailored_match_prior_grid(prior_xpix,prior_ypix,prior_match_xpix,prior_match_ypix)

            e1_LF=tbdata_sim['bias-corrected <e1>']
            e2_LF=tbdata_sim['bias-corrected <e2>']
            SNR_LF=tbdata_sim['model SNratio']
            mode_LF=tbdata_sim['mean-likelihood |e|']
            w_LF=tbdata_sim.weight
            size_LF=tbdata_sim['bias-corrected scalelength /pixels']
            mag_LF=tbdata_sim.catmag
            nd_LF=tbdata_sim['neighbour distance']
            nm_LF=tbdata_sim['neighbour mag']
            strehl_LF=tbdata_sim['PSF-Strehl-ratio']
            PSF_e1=tbdata_sim['PSF-e1']
            PSF_e2=tbdata_sim['PSF-e2']
            PSF_Q11=tbdata_sim['PSF-Q11']
            PSF_Q12=tbdata_sim['PSF-Q12']
            PSF_Q22=tbdata_sim['PSF-Q22']
            fitclass=tbdata_sim.fitclass
            rcorr=tbdata_sim['r correction']
            e1_corr=tbdata_sim['e1 correction']
            e2_corr=tbdata_sim['e2 correction']
            var_LF=tbdata_sim['2D measurement variance']
            print 'looping...'
            for item in match_map:
                #Keep in mind here, the indices for prior related quantities is 0, for Sextractor related ones 1. 
                mag_in_sim.append(prior_mag_Mask[item[0]])
                size_in_sim.append(prior_size_Mask[item[0]])
                e1_in_sim.append(prior_e1_Mask[item[0]])
                e2_in_sim.append(prior_e2_Mask[item[0]])
                N_in_sim.append(prior_N_Mask[item[0]])
                ZB_sim.append(prior_ZB_Mask[item[0]])
                ID_sim.append(prior_ID_Mask[item[0]])
                #Here now all quantities out of the LF catalogue
                e1_sim.append(e1_LF[item[1]])
                e2_sim.append(e2_LF[item[1]])
                snr_sim.append(SNR_LF[item[1]])
                mode_sim.append(mode_LF[item[1]])
                w_sim.append(w_LF[item[1]])
                size_sim.append(size_LF[item[1]])
                mag_sim.append(mag_LF[item[1]])
                nd_sim.append(nd_LF[item[1]])
                nm_sim.append(nm_LF[item[1]])

                strehl_sim.append(strehl_LF[item[1]])
                psf_e1_sim.append(PSF_e1[item[1]])
                psf_e2_sim.append(PSF_e2[item[1]])
                psf_Q11_sim.append(PSF_Q11[item[1]])
                psf_Q22_sim.append(PSF_Q22[item[1]])
                psf_Q12_sim.append(PSF_Q12[item[1]])
   
                fitClass_sim.append(fitclass[item[1]])
                r_corr_sim.append(rcorr[item[1]])
                e1_corr_sim.append(e1_corr[item[1]])
                e2_corr_sim.append(e2_corr[item[1]])
                ls_var_sim.append(var_LF[item[1]])
                rot_sim.append(y)
                
                g1_in.append(g1_)
                g2_in.append(g2_)

snr_sim=array(snr_sim)
mode_sim=array(mode_sim)
w_sim=array(w_sim)
size_sim=array(size_sim)
ZB_sim=array(ZB_sim)
mag_sim=array(mag_sim)
g1_in_tmp=array(g1_in)
g2_in_tmp=array(g2_in)
e1_sim=array(e1_sim)
e2_sim=array(e2_sim)
strehl_sim=array(strehl_sim)
fitClass_sim=array(fitClass_sim)
mag_in_sim=array(mag_in_sim)
size_in_sim=array(size_in_sim)
N_in_sim=array(N_in_sim)
e1_in_sim=array(e1_in_sim)
e2_in_sim=array(e2_in_sim)
psf_e1_sim=array(psf_e1_sim)
psf_e2_sim=array(psf_e2_sim)
r_corr_sim=array(r_corr_sim)
e1_corr_sim=array(e1_corr_sim)
e2_corr_sim=array(e2_corr_sim)
psf_Q11_sim=array(psf_Q11_sim)
psf_Q22_sim=array(psf_Q22_sim)
psf_Q12_sim=array(psf_Q12_sim)
nd_sim=array(nd_sim)
size_raw_sim=size_sim-r_corr_sim
ls_var_sim=array(ls_var_sim)
rot_sim=array(rot_sim)
ID_sim=array(ID_sim)

psf_size_sim=SizeFromMom(psf_Q11_sim,psf_Q22_sim,psf_Q12_sim)

if(len(psf_e1_sim)!=len(e1_sim)):
    print 'mismatch between psf file and catalogue!', len(psf_e1_sim), len(e1_sim)
    exit()

g1_in=g1_in_tmp
g2_in=g2_in_tmp

#CREATE FITS BINARY HERE
print "Writing FITS binary file before cuts"
c1pre=pyfits.Column(name='LFweight',format='D',array=w_sim)
c2pre=pyfits.Column(name='e1',format='D',array=e1_sim)
c3pre=pyfits.Column(name='e2',format='D',array=e2_sim)
c4pre=pyfits.Column(name='mag_out',format='D',array=mag_sim)
c5pre=pyfits.Column(name='g1',format='D',array=g1_in)
c6pre=pyfits.Column(name='g2',format='D',array=g2_in)
c7pre=pyfits.Column(name='strehl',format='D',array=strehl_sim)
c8pre=pyfits.Column(name='snr_model',format='D',array=snr_sim)
c9pre=pyfits.Column(name='size_out',format='D',array=size_sim)
c10pre=pyfits.Column(name='mag_in',format='D',array=mag_in_sim)
c11pre=pyfits.Column(name='size_in',format='D',array=size_in_sim)
c12pre=pyfits.Column(name='N_in',format='D',array=N_in_sim)
c13pre=pyfits.Column(name='e1_in',format='D',array=e1_in_sim)
c14pre=pyfits.Column(name='e2_in',format='D',array=e2_in_sim)
c15pre=pyfits.Column(name='ZB_in',format='D',array=ZB_sim)
c16pre=pyfits.Column(name='Cat_ID',format='D',array=ID_sim)
c17pre=pyfits.Column(name='size_corr',format='D',array=r_corr_sim)
c18pre=pyfits.Column(name='e1_corr',format='D',array=e1_corr_sim)
c19pre=pyfits.Column(name='e2_corr',format='D',array=e2_corr_sim)
c20pre=pyfits.Column(name='nd',format='D',array=nd_sim)
c21pre=pyfits.Column(name='LS variance',format='D',array=ls_var_sim)
c22pre=pyfits.Column(name='fitclass',format='D',array=fitClass_sim)
c23pre=pyfits.Column(name='rotation',format='D',array=rot_sim)
c24pre=pyfits.Column(name='psf_size_in',format='D',array=psf_size_sim)
c25pre=pyfits.Column(name='psf_e1_in',format='D',array=psf_e1_sim)
c26pre=pyfits.Column(name='psf_e2_in',format='D',array=psf_e2_sim)
tbhdu = pyfits.new_table([c1pre, c2pre, c3pre, c4pre, c5pre, c6pre, c7pre, c8pre, c9pre, c10pre, c11pre, c12pre, c13pre, c14pre,c15pre, c16pre, c17pre, c18pre, c19pre, c20pre ,c21pre, c22pre,c23pre,c24pre,c25pre,c26pre])
tbhdu.writeto(binary_filename,clobber=True)
print "FITS binary before matching  created, proceeding with analysis."
