from astropy.table import Table, vstack
import sys
import pyfits
import os
import numpy as np

num_PSF=int(sys.argv[1])
DirIni=sys.argv[2]
nameIni=sys.argv[3]
nameOutput=sys.argv[4]

#if not os.path.isfile(DirIni+'/'+nameOutput):
#    print 'Merging master files..'
#    table=[]
##    for kk in range(num_PSF):
#    for kk in [0,1,2,3,4,13,14,15,16,17]:
#        print "HACK ALERT: Running PSF rotated sets only"
#        if(os.path.isfile(DirIni+'/'+nameIni+'_set_'+str(kk)+'.fits')):
#            print 'PSF ', kk
#            table.append(Table.read(DirIni+'/'+nameIni+'_set_'+str(kk)+'.fits'))
#        else:
#            print 'Cannot find master file for PSF ', kk
#        
#    new=vstack(table)
#    new.write(DirIni+'/'+nameOutput)
#else:
#    print 'Merged file already found.'

#Generate tables for each tomographic bin 

#TomoBin=[0.101, 0.301, 0.501, 0.701, 0.901, 10]
TomoBin=[0.101, 0.301, 0.501, 0.701, 0.901, 1.201]
print 'Splitting the merged catalogue into different tomographic bins'
for kk in range(0, len(TomoBin)-1):
    t=pyfits.open(DirIni+'/'+nameOutput)
    tbdata = t[1].data

    ## Old 4-band photometry
    mask=(tbdata['ZB4_in']>=TomoBin[kk])&(tbdata['ZB4_in']<TomoBin[kk+1])
    newtbdata = tbdata[mask]
    hdu = pyfits.BinTableHDU(data=newtbdata)
    print DirIni, kk
    hdu.writeto(DirIni+'/'+'MasterCat_Tomo4Bin_'+str(kk+1)+'.fits', clobber=True)

    ## New 9-bandp photometry
    mask=(tbdata['ZB9_in']>=TomoBin[kk])&(tbdata['ZB9_in']<TomoBin[kk+1])
    newtbdata = tbdata[mask]
    hdu = pyfits.BinTableHDU(data=newtbdata)
    print DirIni, kk
    hdu.writeto(DirIni+'/'+'MasterCat_Tomo9Bin_'+str(kk+1)+'.fits', clobber=True)
