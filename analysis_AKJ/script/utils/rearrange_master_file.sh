#!/bin/bash

nameRun=FC17SB
echo "randomKey: " $nameRun

python_file=/disks/shear15/KiDS/ImSim/analysis/grid/python/rearrange_master_file.py
archdir=/disks/shear15/KiDS/ImSim/pipeline/archive

for psf in 0 1 2 3 4
do
    #python2.7 /disks/shear15/KiDS/ImSim/analysis/grid/python/propagate_chi2nu.py $archdir/$nameRun/MasterCat_$(echo $nameRun)_set_$psf.fits
    python2.7 $python_file $archdir/$nameRun/MasterCat_$(echo $nameRun)_set_$psf.fits
done
