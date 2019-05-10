## Author: Arun Kannawadi
## This python routine replaces the measured e1,e2 with true e1,e2 to calculate the selection bias

import os, sys
from astropy.io import fits

randomKey = sys.argv[1]

ARCHDIR = os.path.join('/disks/shear15/KiDS/ImSim/pipeline/archive/',randomKey)
DESTDIR = os.path.join('/disks/shear15/KiDS/ImSim/pipeline/archive/',randomKey+'TrueEllip')

all_items = os.listdir(ARCHDIR)

for item in all_items:
    if 'MasterCat_' in item and '.fits' in item:
        catname = os.path.join(ARCHDIR, item)
        hdulist = fits.open(catname)

        dat = hdulist[1].data

        g1 = dat['g1']
        g2 = dat['g2']

        e1_in = dat['e1_in']
        e2_in = dat['e2_in']

        e_in = e1_in + 1j*e2_in
        g = g1 + 1j*g2

        e_true = (e_in+g)/( 1.+g.conj()*e_in )

        dat['e1'] = e_true.real
        dat['e2'] = e_true.imag

        catname = os.path.join( DESTDIR, item.replace(randomKey,randomKey+'TrueEllip') )

        hdulist.writeto(catname,clobber=True)
        print "Copied ", catname
