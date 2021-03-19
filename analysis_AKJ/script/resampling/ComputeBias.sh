#!/bin/bash

#####################################################################################
#Resampling the image simulations such that the SNR and R distributions match the 
#one from the data in each tomographic bin. (see Fenech Conti et al. 2017 for details)
#
#Massimo Viola viola@strw.leidenuniv.nl
#
#
#Timing on eemmeer (13 jobs in parallel, 5 tomographic bins)
#
#Resampling stage: 45 minutes
#Calculate bias: 5 minutes
# 
######################################################################################
###############################################################################
name_run=$1 # TSTGr0729trueposCHK
echo "Run ID: " $name_run
###############################################################################

MainDir=/disks/shear15/KiDS/ImSim/analysis_AKJ/resampling/
srcDir=$MainDir/src/

#This is the directory where the master fits file are (post LF run on image simulations)
#DirSim=/disks/shear12/herbonnet/kids/data/vault_03052016_newmagdistr/julian/stars/
DirSim=/disks/shear15/KiDS/ImSim/pipeline/archive/$name_run/

#This is the directory where the KiDS data are (one catalogue per patch)
DirKiDS=/disks/shear14/KIDS/KV450_CATALOGUES_PATCH_V0.5.9/

#Creating new directories
mkdir -p $MainDir/output/$name_run
mkdir -p $MainDir/output/$name_run/PSFset/
mkdir -p $MainDir/output/$name_run/bias/
mkdir -p $MainDir/output/$name_run/bias/average/

DirOut=$MainDir/output/$name_run/PSFset/
DirBias=$MainDir/output/$name_run/bias/
DirBiasFinal=$MainDir/output/$name_run/bias/average/

for patch in G9 G12 G15 G23 GS 
do
    echo 'Doing..' $patch
    for tomo in 0 1 2 3 4
    do 
        python2.7 $srcDir/Resampling_sim_2D_v2.1.py $DirSim $DirKiDS $DirOut $patch $tomo 1 &
    done

    for PSF in 0 1 2 3 4 # 5 6 7 8 9 10 11 12
    do
        python2.7 $srcDir/Resampling_sim_2D_v2.1.py $DirSim $DirKiDS $DirOut $patch $PSF 0 &
    done
    wait
done

wait

#Compute bias for each PSF set and each tomographic bin

for patch in G9 G12 G15 G23 GS 
do
    for PSF in 0 1 2 3 4 5 6 7 8 9 10 11 12
    do
	python2.7 $srcDir/ComputeBias_resampling_Tomoueber_v2.1.py $DirOut $DirBias $patch $PSF &
    done
    wait
done

wait
echo 'Bias for each PSF computed!'

#Compute average bias over all sets

for patch in G9 G12 G15 G23 GS 
do
   python2.7 $srcDir/ComputeAverageBiasSet_v2.1.py $DirBias $DirBiasFinal $patch & 
done

echo 'I am done with resampling!'

