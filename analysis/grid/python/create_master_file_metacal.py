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
list_name=sys.argv[4]
list_shear=numpy.genfromtxt(list_name,dtype="|S64") ## dtype S32 -> S64 by AKJ on Feb 25, 2018
###########################

#Initialise strings

prior_matched=[]
snr_sim=[]
mode_sim=[]
oldw_sim=[]
w_sim=[]
size_sim=[]
fluxradius_sim=[]
fwhmimage_sim=[]
nd_sim=[]
nm_sim=[]
fitClass_sim=[]
mag_sim=[]
btt_sim=[]
star_gal_prob_sim=[]
contamination_radius_sim=[]
metacal_m1_sim=[]
metacal_m2_sim=[]
metacal_c1_sim=[]
metacal_c2_sim=[]
g1_in=[]
g2_in=[]
e1_sim=[]
e2_sim=[]
strehl_sim=[]
mag_in_sim=[]
size_in_sim=[]
N_in_sim=[]
f_in_sim=[]
e1_in_sim=[]
e2_in_sim=[]
x_in_sim=[]
y_in_sim=[]
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
ZB4_sim=[]
ZB9_sim=[]
ID_sim=[]

cc=0

#I am starting looping here through the different shear values.

for x in range(0,len(list_shear)): #shear values
    if(True):
        try:
            Name_tmp, Psf_type, Random_string = list_shear[x].split('_')
            Psf_type = int(Psf_type)
        except:
            continue ## since it is a different directory
        #    check_string=list_shear[x][10:11]
        #    if(check_string=='_'):
        #        Random_string=str(list_shear[x][11:])
        #        Psf_type=int(list_shear[x][9:10])
        #    else:
        #        Random_string=str(list_shear[x][12:])
        #        Psf_type=int(list_shear[x][9:11])
    else:
        Name_tmp=str(list_shear[x])[0:8]
        Random_string='bbbbbbbb'
        Psf_type=18789
    print Name_tmp, Psf_type, Random_string
    if(Name_tmp[0]=='m'):
        g1_=float(Name_tmp[0:4].replace("m","-0.0"))
    elif(Name_tmp[0]=='p'):
        g1_=float(Name_tmp[0:4].replace("p","0.0"))
    else:
        continue
    if(Name_tmp[4]=='m'):
        g2_=float(Name_tmp[4:8].replace("m","-0.0"))
    elif(Name_tmp[4]=='p'):
        g2_=float(Name_tmp[4:8].replace("p","0.0"))
    else:
        continue
    #New main loop over all available sets for each shear value

    for y in range(0,4): #rotation
        cc=cc+1
        if((Psf_type == int(Psf)) or run_flag):
            #tbdata_sim = pyfits.getdata(Dir+'/'+Name_tmp+'_'+str(Psf_type)+'_'+Random_string+'/0'+str(y)+'.output.rot.fits.asc.scheme2b_corr.fits')
            sim=numpy.genfromtxt(Dir+'/'+Name_tmp+'_'+str(Psf_type)+'_'+Random_string+'/0'+str(y)+'.output.rot.fits.asc.scheme2b_corr', comments='#')

            ##MATCHING###
            print 'I am using..', Name_tmp+'_'+str(Psf_type)+'_'+Random_string
            #Checking if LF was working on prior grid or Sextractor catalogue
            grid_check=os.path.isfile(Dir+'/'+str(list_shear[x])+'/sexrot0'+str(y)+'.cat')
            if(grid_check):
                grid_file='/sexrot0'+str(y)+'.cat'
            else:
                grid_file='/sex.cat'
                #grid_file='/prior.match.asci'
		#print 'I cannot find the sextractor catalogue..'
		#exit()
            #tbdata_grid = pyfits.getdata(Dir+'/'+str(list_shear[x])+grid_file)
            catSE = np.loadtxt(Dir+'/'+str(list_shear[x])+grid_file)
            fluxradius_SE = catSE[:,12]
            fwhmimage_SE = catSE[:,14]

            prior_match=numpy.genfromtxt(Dir+'/'+str(list_shear[x])+grid_file, comments='#')
            prior_match_xpix=prior_match[:,0]
            prior_match_ypix=prior_match[:,1]
            #prior_match_xpix=tbdata_grid.X_IMAGE
            #prior_match_ypix=tbdata_grid.Y_IMAGE
            if(grid_check):
                prior_match_mag=prior_match[:,6]
                #prior_match_mag=tbdata_grid.MAG_AUTO
            else:
                prior_match_mag=prior_match[:,6]#??
                
            #Reading simulation prior
            prior_check=os.path.isfile(Dir+'/'+Name_tmp+'_'+str(Psf_type)+'_'+Random_string+'/prior_new')
            print Dir
            
            if(prior_check):
                #hdu_prior = pyfits.open(Dir+'/prior')
                #tbdata_prior=hdu_prior[1].data
                prior=numpy.genfromtxt(Dir+'/'+Name_tmp+'_'+str(Psf_type)+'_'+Random_string+'/prior_new', comments='#')
                print 'Reading: ', Dir+'/prior_new'
            else:
                print 'I cannot find the prior file!', Dir+'/'+Name_tmp+'_'+str(Psf_type)+'_'+Random_string+'/prior_new'
                exit()

            prior_xpix=prior[:,0]
            prior_ypix=prior[:,1]
            prior_mag=prior[:,2]
            prior_size=prior[:,3]
            prior_fitClass=prior[:,6]
            try:
                prior_N=prior[:,7]
                prior_ZB4=prior[:,8]
                try:
                    prior_ZB9=prior[:,9]
                except:
                    prior_ZB9 = -np.ones(len(prior))
                ##AKJ: Carry over a unique ID from prior file, if it exists. Else, number them.
                if prior.shape[1]>10:
                    prior_ID=prior[:,10]
                else:
                    prior_ID=np.arange(0,prior.shape[0],1)
            except:
                print "Didn't find a few columns in the prior. Filling it up with garbage values"
                prior_N = -np.ones_like(prior_mag)
                prior_ZB4 = -4*np.ones_like(prior_mag)
                prior_ZB9 = -9*np.ones_like(prior_mag)
                prior_ID = 100*np.ones_like(prior_mag)

            #prior_xpix=tbdata_prior.X_IMAGE#-1627.71
            #prior_ypix=tbdata_prior.Y_IMAGE#-646.21
            #prior_mag=tbdata_prior.MAG_AUTO
            #prior_size=tbdata_prior.RE_GALFIT_HI
	    #Star/gal flag (FLAG==-1: star, FLAG==2: detected in KiDS, FLAG==0: detected in HST only)
            #prior_flag=tbdata_prior.FLAG
            #prior_q=tbdata_prior.BA_GALFIT_HI
            #prior_phi=tbdata_prior.PA_GALFIT_HI#+90.0)*(3.1415/180.)
            #prior_N=tbdata_prior.N_GALFIT_HI
            #prior_ZB4=tbdata_prior.ZB
            #prior_ID=tbdata_prior.SeqNr
            #prior_mode=(1.-prior_q)/(1.+prior_q)
            ############################# 

            #Introduce a mask in the prior file (to remove stars, but it can be generalised)
            
            #maskPrior=(prior_flag>-99)
            maskPrior=(prior_fitClass>-99)
            #and we redefine new arrays

#AKJ: Where do the fields below come from???
            prior_xpix_Mask=prior_xpix[maskPrior]
            prior_ypix_Mask=prior_ypix[maskPrior]
            prior_mag_Mask=prior_mag[maskPrior]
            prior_size_Mask=prior_size[maskPrior]
            prior_f_Mask=prior_fitClass[maskPrior]
            prior_N_Mask=prior_N[maskPrior]
            prior_ZB4_Mask=prior_ZB4[maskPrior]
            prior_ZB9_Mask=prior_ZB9[maskPrior]
            prior_ID_Mask=prior_ID[maskPrior]
            
            #Note: this is the ellipticity for rotation 00!!!!!!!!
            prior_e1_00=prior[:,4]
            prior_e2_00=prior[:,5]
            #prior_e1_00=prior_mode*cos(2.*prior_phi)
            #prior_e2_00=prior_mode*sin(2.*prior_phi)

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
            #if(grid_check):
                ### AKJ: Offset changed from -0.5 to +0.5 by AKJ for testing
                ### AKJ: Pixel range changed from 3.0 to 5.0 for TSTRGC run
            if os.path.isfile(Dir+'/'+Name_tmp+'_'+str(Psf_type)+'_'+Random_string+'/.tree.p'): ## HACK ALERT
                print "Tree found!"
                import cPickle
                with open(Dir+'/'+Name_tmp+'_'+str(Psf_type)+'_'+Random_string+'/.tree.p','r') as ffff:
                    Tree = cPickle.load(ffff)
                    match_map=new_tailored_match_sextractor_grid(prior_xpix_Mask,prior_ypix_Mask,prior_mag_Mask,prior_match_xpix,prior_match_ypix,prior_match_mag,5,3.0,-0.5,Tree)
            else:
                match_map=new_tailored_match_sextractor_grid(prior_xpix_Mask,prior_ypix_Mask,prior_mag_Mask,prior_match_xpix,prior_match_ypix,prior_match_mag,5,3.0,-0.5)
            #else:
            #match_map=tailored_match_prior_grid(prior_xpix,prior_ypix,prior_match_xpix,prior_match_ypix)

            e1_LF=sim[:,22]#tbdata_sim['bias-corrected <e1>']
            e2_LF=sim[:,23]#tbdata_sim['bias-corrected <e2>']
            SNR_LF=sim[:,10]#tbdata_sim['model SNratio']
            mode_LF=sim[:,21]#tbdata_sim['mean-likelihood |e|'] ## IS THIS WHAT IT SHOULD BE
            oldw_LF=sim[:,59]
            w_LF=sim[:,4]#tbdata_sim.weight
            size_LF=sim[:,19]#tbdata_sim['bias-corrected scalelength /pixels']
            mag_LF=sim[:,26]#tbdata_sim.catmag
            nd_LF=sim[:,25]#tbdata_sim['neighbour distance']
            nm_LF=sim[:,24]#tbdata_sim['neighbour mag']
            strehl_LF=sim[:,14]#tbdata_sim['PSF-Strehl-ratio']
            PSF_e1=sim[:,12]#tbdata_sim['PSF-e1']
            PSF_e2=sim[:,13]#tbdata_sim['PSF-e2']
            PSF_Q11=sim[:,15]#tbdata_sim['PSF-Q11']
            PSF_Q12=sim[:,17]#tbdata_sim['PSF-Q12']
            PSF_Q22=sim[:,16]#tbdata_sim['PSF-Q22']
            fitclass=sim[:,5]#tbdata_sim.fitclass
            rcorr=size_LF-sim[:,6]#tbdata_sim['r correction']
            e1_corr=e1_LF-sim[:,2]#tbdata_sim['e1 correction']
            e2_corr=e2_LF-sim[:,3]#tbdata_sim['e2 correction']
            var_LF=sim[:,20]#tbdata_sim['2D measurement variance']
            bulgefraction_LF=sim[:,7]
            star_gal_prob_LF=sim[:,18]
            contamination_radius_LF=sim[:,11]
            metacal_m1=sim[:,39]
            metacal_m2=sim[:,40]
            metacal_c1=sim[:,41]
            metacal_c2=sim[:,42]
            print 'looping...'
            for item in match_map:
                prior_matched.append( item[0]>-1 )
                ## The matching is bogus if item[0]==-1. It is assigned to keep the following long list of assignments happy.

                #Keep in mind here, the indices for prior related quantities is 0, for Sextractor related ones 1. 
                mag_in_sim.append(prior_mag_Mask[item[0]])
                size_in_sim.append(prior_size_Mask[item[0]])
                e1_in_sim.append(prior_e1_Mask[item[0]])
                e2_in_sim.append(prior_e2_Mask[item[0]])
                N_in_sim.append(prior_N_Mask[item[0]])
                f_in_sim.append(prior_f_Mask[item[0]])
                ZB4_sim.append(prior_ZB4_Mask[item[0]])
                ZB9_sim.append(prior_ZB9_Mask[item[0]])
                ID_sim.append(prior_ID_Mask[item[0]])
                x_in_sim.append(prior_xpix_Mask[item[0]])
                y_in_sim.append(prior_ypix_Mask[item[0]])
                #Here now all quantities out of the LF catalogue
                e1_sim.append(e1_LF[item[1]])
                e2_sim.append(e2_LF[item[1]])
                snr_sim.append(SNR_LF[item[1]])
                mode_sim.append(mode_LF[item[1]])
                oldw_sim.append(oldw_LF[item[1]])
                w_sim.append(w_LF[item[1]])
                size_sim.append(size_LF[item[1]])
                mag_sim.append(mag_LF[item[1]])
                nd_sim.append(nd_LF[item[1]])
                nm_sim.append(nm_LF[item[1]])
                btt_sim.append(bulgefraction_LF[item[1]])
                star_gal_prob_sim.append(star_gal_prob_LF[item[1]])
                contamination_radius_sim.append(contamination_radius_LF[item[1]])

                strehl_sim.append(strehl_LF[item[1]])
                psf_e1_sim.append(PSF_e1[item[1]])
                psf_e2_sim.append(PSF_e2[item[1]])
                psf_Q11_sim.append(PSF_Q11[item[1]])
                psf_Q22_sim.append(PSF_Q22[item[1]])
                psf_Q12_sim.append(PSF_Q12[item[1]])
   
                fluxradius_sim.append(fluxradius_SE[item[1]])
                fwhmimage_sim.append(fwhmimage_SE[item[1]])

                fitClass_sim.append(fitclass[item[1]])
                r_corr_sim.append(rcorr[item[1]])
                e1_corr_sim.append(e1_corr[item[1]])
                e2_corr_sim.append(e2_corr[item[1]])
                ls_var_sim.append(var_LF[item[1]])
                rot_sim.append(y)

                metacal_m1_sim.append(metacal_m1[item[1]])
                metacal_m2_sim.append(metacal_m2[item[1]])
                metacal_c1_sim.append(metacal_c1[item[1]])
                metacal_c2_sim.append(metacal_c2[item[1]])
                
                g1_in.append(g1_)
                g2_in.append(g2_)

prior_matched=array(prior_matched)
snr_sim=array(snr_sim)
mode_sim=array(mode_sim)
oldw_sim=array(oldw_sim)
w_sim=array(w_sim)
size_sim=array(size_sim)
ZB4_sim=array(ZB4_sim)
ZB9_sim=array(ZB9_sim)
mag_sim=array(mag_sim)
g1_in_tmp=array(g1_in)
g2_in_tmp=array(g2_in)
fluxradius_sim=array(fluxradius_sim)
fwhmimage_sim=array(fwhmimage_sim)
e1_sim=array(e1_sim)
e2_sim=array(e2_sim)
strehl_sim=array(strehl_sim)
fitClass_sim=array(fitClass_sim)
mag_in_sim=array(mag_in_sim)
size_in_sim=array(size_in_sim)
N_in_sim=array(N_in_sim)
f_in_sim=array(f_in_sim)
e1_in_sim=array(e1_in_sim)
e2_in_sim=array(e2_in_sim)
x_in_sim=array(x_in_sim)
y_in_sim=array(y_in_sim)
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
btt_sim=array(btt_sim)
star_gal_prob_sim=array(star_gal_prob_sim)
contamination_radius_sim=array(contamination_radius_sim)
metacal_m1_sim=array(metacal_m1_sim)
metacal_m2_sim=array(metacal_m2_sim)
metacal_c1_sim=array(metacal_c1_sim)
metacal_c2_sim=array(metacal_c2_sim)

psf_size_sim=SizeFromMom(psf_Q11_sim,psf_Q22_sim,psf_Q12_sim)

if(len(psf_e1_sim)!=len(e1_sim)):
    print 'mismatch between psf file and catalogue!', len(psf_e1_sim), len(e1_sim)
    exit()

g1_in=g1_in_tmp
g2_in=g2_in_tmp

#CREATE FITS BINARY HERE
print "Writing FITS binary file before cuts"
c0pre=pyfits.Column(name='oldLFweight',format='D',array=oldw_sim)
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
c12bpre=pyfits.Column(name='f_in',format='D',array=f_in_sim)
c13pre=pyfits.Column(name='e1_in',format='D',array=e1_in_sim)
c14pre=pyfits.Column(name='e2_in',format='D',array=e2_in_sim)
c15pre=pyfits.Column(name='ZB4_in',format='D',array=ZB4_sim)
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
c27pre=pyfits.Column(name='X_IMAGE',format='D',array=x_in_sim)
c28pre=pyfits.Column(name='Y_IMAGE',format='D',array=y_in_sim)
c29pre=pyfits.Column(name='prior_matched',format='D',array=prior_matched)
c30pre=pyfits.Column(name='ZB9_in',format='D',array=ZB9_sim)
c31pre=pyfits.Column(name='FLUX_RADIUS',format='D',array=fluxradius_sim)
c32pre=pyfits.Column(name='FWHM_IMAGE',format='D',array=fwhmimage_sim)
c33pre=pyfits.Column(name='bulge_fraction',format='D',array=btt_sim)
c34pre=pyfits.Column(name='star_gal_prob',format='D',array=star_gal_prob_sim)
c35pre=pyfits.Column(name='contamination_radius',format='D',array=contamination_radius_sim)
c36pre=pyfits.Column(name='metacal_m1',format='D',array=metacal_m1_sim)
c37pre=pyfits.Column(name='metacal_m2',format='D',array=metacal_m2_sim)
c38pre=pyfits.Column(name='metacal_c1',format='D',array=metacal_c1_sim)
c39pre=pyfits.Column(name='metacal_c2',format='D',array=metacal_c2_sim)
tbhdu = pyfits.new_table([c0pre,c1pre, c2pre, c3pre, c4pre, c5pre, c6pre, c7pre, c8pre, c9pre, c10pre, c11pre, c12pre, c13pre, c14pre, c15pre, c16pre, c17pre, c18pre, c19pre, c20pre ,c21pre, c22pre,c23pre,c24pre,c25pre,c26pre,c27pre,c28pre,c29pre,c30pre,c31pre,c32pre,c33pre,c34pre,c12bpre,c35pre,c36pre,c37pre,c38pre,c39pre])
tbhdu.writeto(binary_filename,clobber=True)
print "FITS binary before matching  created, proceeding with analysis."
