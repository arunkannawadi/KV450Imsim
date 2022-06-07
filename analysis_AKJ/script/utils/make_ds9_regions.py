from astropy.io import fits
import numpy as np
import sys
import optparse

def make_ds9_regions(cat,cuts=None,filename=None,color='green',label_field=None, shape='elliptical', mode='ang',CS='fk5'):

    if isinstance(cat,str):
        cat_filename = cat
        cat = fits.open(cat_filename)

    if cuts is None:
        dat = cat[1].data
    else:
        dat = cat[1].data[cuts]

    if mode=='ang':
        RADEC_fieldname_opts = [('ALPHA_J2000','DELTA_J2000'),
                                ('RA_THELI', 'DEC_THELI'),
                                ('RA', 'DEC'),
                                ('col1_1','col2_1'), ('col1', 'col2')]
        RAD_fieldname_opts = ['FLUX_RADIUS','FWHM_WORLD_THELI','FWHM_WORLD','col3']

    elif mode=='xy':
        RADEC_fieldname_opts = [('X_IMAGE_THELI','Y_IMAGE_THELI'),
                             ('X_IMAGE','Y_IMAGE'), ('IMAGE_X', 'IMAGE_Y'),
                             ('col1_1','col2_1'), ('col1', 'col2')]

        RAD_fieldname_opts = ['FWHM_IMAGE_THELI', 'FWHM_IMAGE','FLUX_RADIUS']

    POSANG_fieldname ='THETA_WORLD'
    ELLIPTICITY_fieldname = 'ELLIPTICITY'

    found_RA = False
    for RADEC_fieldname in RADEC_fieldname_opts:
        RA_fieldname, DEC_fieldname = RADEC_fieldname
        try:
            RA = dat[RA_fieldname] - 1627.71*0
            DEC = dat[DEC_fieldname] - 646.21*0
            print "Using {0} for RA and {1} for DEC".format(RA_fieldname, DEC_fieldname)
            found_RA = True
            break
        except:
            pass

    if found_RA is False:
        raise RuntimeError("Found no valid RA, DEC fieldnames. Exiting.")

    found_RAD = False
    for RAD_fieldname in RAD_fieldname_opts:
        try:
            RAD = dat[RAD_fieldname]
            RAD_unit = dat.columns[dat.names.index(RAD_fieldname)].unit
            found_RAD = True
            print "Using {0} for the size".format(RAD_fieldname)
            if RAD_fieldname=='size_in':
                RAD = RAD/0.214         ## hack
            break
        except:
            pass

    if found_RAD is False:
        RAD = 5*np.ones_like(RA)
        if mode=='ang':
            RAD_unit = 'arcsec'
        else:
            RAD_unit = 'pix'
        print 'Found no valid size fieldnames. Using 5 {0} as default'.format(RAD_unit)

    if mode=='ang':
        if RAD_unit=='deg' or RAD_unit=='degree' or RAD_unit=='degrees':
            RAD *= 3600. ## Convert to arcsec

    if shape is 'elliptical':
        POSANG = dat[POSANG_fieldname]
#        if 'POSANG' in dat.names:
#            POSANG = dat['POSANG']
#	elif 'POSANG_1' in dat.names:
#	    POSANG = dat['POSANG']
#	elif 'THETA_IMAGE' in dat.names:
#	    POSANG = dat['THETA_IMAGE']
#	elif 'POSANG_THELI' in dat.names:
#	    POSANG = dat['POSANG_THELI']
#	else:
#	    POSANG = dat['PA_GALFIT_HI']
#
#	if 'ELLIPTICITY' in dat.names:
#	    ELLIPTICITY = dat['ELLIPTICITY']
#	elif 'ELLIPTICITY_1' in dat.names:
#	    ELLIPTICITY = dat['ELLIPTICITY_1']
#	elif 'BA_GALFIT_HI' in dat.names:
#	    ELLIPTICITY = 1 - dat['BA_GALFIT_HI']
#	else:
#	    ELLIPTICITY = dat['ELLIPTICITY_THELI']

        ELLIPTICITY = dat[ELLIPTICITY_fieldname]
	Q = (1 - ELLIPTICITY)
        #Q = dat[ELLIPTICITY_fieldname]

        RAD_MAJOR = RAD/(1.+Q)
        RAD_MINOR = RAD*Q/(1.+Q)

    if label_field is not None:
        LABEL = dat[label_field]
    else:
        LABEL = ['']*len(RA)

    if mode=='ang':
        CS = CS
        unit = '"'
    elif mode=='xy':
        CS = 'IMAGE'
        unit = ''

    prefix = '# Region file format: DS9 version 4.1 \n' \
    'global color=green dashlist=8 3 width=1 font="helvetica 1 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n' \
    '{CS} \n'.format(CS=CS)
    prefix = prefix.replace('green', color)

    lines = [prefix]
    for line_no in xrange(len(RA)):
        if shape is 'elliptical':
            line = 'ellipse {ra} {dec} {radius_major}{unit} {radius_minor}{unit} {posang}'.format(\
	            ra=RA[line_no], dec=DEC[line_no], radius_major=RAD_MAJOR[line_no],\
		    radius_minor=RAD_MINOR[line_no], posang=POSANG[line_no], unit=unit)
	else:
            line = 'circle {ra} {dec} {radius}{unit}'.format(\
                   ra=RA[line_no], dec=DEC[line_no],radius=RAD[line_no],unit=unit)

        if label_field is not None:
            line += " #text = {{{0}}}\n".format(LABEL[line_no])
        else:
            line += '\n'

        lines.append( line )

    if filename is None:
        if 'cat_filename' in locals():
            filename = cat_filename.rpartition('.')[0]
        else:
            filename = "ds9_regions"

    with open(filename+'.reg','w') as f:
        f.writelines(lines)


if __name__=='__main__':

    cat_filename = sys.argv[1]

    ## Set up run-time args
    argopts = optparse.OptionParser("usage: %prog [options] arg1")

    argopts.add_option("-m", "--mode",dest="mode",default='xy',type="string",
                        help="Specify mode: 'xy' or 'ang'")
    argopts.add_option("-c", "--color",dest="color",default="green",type="string",
                        help="Specify border color")
    argopts.add_option("-s", "--shape",dest="shape",default="circle",type="string",
                        help="Specify shape of the regions: circle or ellipse")
    argopts.add_option("-C", "--coord",dest="coord",default="fk5",type="string",
                        help="Specify the Coordinate System on the sky (if --mode==ang )")
    argopts.add_option("-l","--label",dest="label",default=None,type="string",
                        help="Specify the column to use as object label")

    (options, args) = argopts.parse_args()

    #check_flux_distribution()
    #get_masked_KiDS_Theli_image()

    ## Make ds9 regions
    from astropy.io import fits
    #cat_filename = '/disks/shear15/KiDS/ImSim/pipeline/archive/m400m000_0_LRIW1MQWKO/sex.cat'
    cuts = None
    #cat = fits.open(cat_filename)
    #dat = cat[1].data
    #RA_Griffith, DEC_Griffith = dat['RA'], dat['DEC']
    #RA_Theli, DEC_Theli = dat['RA_THELI'], dat['DEC_THELI']
    #cuts = (RA_Griffith>=RA_Theli[RA_Theli>0].min())&(RA_Griffith<=RA_Theli.max())&\
    #       (DEC_Griffith>=DEC_Theli[DEC_Theli>0].min())&(DEC_Griffith<=DEC_Theli.max())
    #make_ds9_regions(cat_filename,cuts=cuts,color='blue',mode='xy',shape=None)
    make_ds9_regions(cat_filename, cuts=cuts, color=options.color,mode=options.mode,shape=options.shape,label_field=options.label)

    #readwritetimer()
    #extractDithers()
    #for i in xrange(5):
    #    chopImage('short_expo{0}_noisy.fits'.format(i),output_filename='chopped_short_exp{0}.fits'.format(i))
