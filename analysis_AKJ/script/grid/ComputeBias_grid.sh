#!/bin/bash

###############################################################################
name_run=$1 # TSTGr0729trueposCHK
echo "Run ID: " $name_run
###############################################################################

MainDir=/disks/shear15/KiDS/ImSim/analysis/grid/
srcDir=$MainDir/python/

#Compute bias as a function of SNR and R (see IFC17)

CatDir=/disks/shear15/KiDS/ImSim/pipeline/archive/$name_run/
strehl=no
rotation=no
Nbin=20 #20
Nbin2=20 #20
cal_type=None
cal_select=snr_res_alt_mv
m_func='empty'
c_func='empty'
TomoFlag=Yes

cd $CatDir ## new, added by AKJ

for PSF in 10 
do
    masterFile=MasterCat_$(echo $name_run)_all_$(echo $PSF)_PSF.fits
    #masterFile=MasterCat_$(echo $name_run)_all_$(echo $PSF)_PSF_clean.fits
    
    #Calculate the surface of doom (using all data)
    if [ ! -f $CatDir/Results/2bin/MV_$(echo $PSF)_100_SignalToNoise_ResolutionAltMV_binning_global.txt ]
    then
        python2.7 $srcDir/Sim_mc_v0.7.py $CatDir $CatDir/$masterFile $strehl $rotation $Nbin $Nbin2 $cal_type $cal_select $m_func $c_func $PSF
    fi
    #Calculate the surface of doom (in tomographic bins)
    if [ "$TomoFlag" == "Yes" ]
    then
	for tomo in  1 2 3 4 5
	do
	    if [ ! -f $CatDir/Results/2bin/MV_Tomo4$(echo $tomo)_100_SignalToNoise_ResolutionAltMV_binning_global.txt ]
	    then
		masterFileTomo=MasterCat_Tomo4Bin_$(echo $tomo).fits
		python2.7 $srcDir/Sim_mc_v0.7.py $CatDir $CatDir/$masterFileTomo $strehl $rotation $Nbin $Nbin2 $cal_type $cal_select $m_func $c_func Tomo4$(echo $tomo) > $CatDir/computebias.txt
	    fi

	    if [ ! -f $CatDir/Results/2bin/MV_Tomo9$(echo $tomo)_100_SignalToNoise_ResolutionAltMV_binning_global.txt ]
	    then
		masterFileTomo=MasterCat_Tomo9Bin_$(echo $tomo).fits
		python2.7 $srcDir/Sim_mc_v0.7.py $CatDir $CatDir/$masterFileTomo $strehl $rotation $Nbin $Nbin2 $cal_type $cal_select $m_func $c_func Tomo9$(echo $tomo) > $CatDir/computebias.txt
	    fi
	done
    fi
done

#Calculate the bias

for PSF in 10
do
    masterFile=MasterCat_$(echo $name_run)_all_$(echo $PSF).fits ## set changed to all by AKJ
    surfaceOfDoom=MV_$(echo $PSF)_100_SignalToNoise_ResolutionAltMV_binning_global.txt
    outputFile=MasterCat_$(echo $name_run)_all_$(echo $PSF).calibrated.fits ## set changed to all by AKJ 

    #Uncomment it if you want to add the calibration to the master file. This is useful
    #if you want to look at residual biases post-calibration
    
    #python2.7 $srcDir/add_cal_to_master.py $CatDir/$masterFile mycalibration MV bin $CatDir/Results/2bin/$surfaceOfDoom $CatDir/$outputFile

    #Add the calibration to the KiDS data
    if [ ! -f $CatDir/Results/Summary_multiplicative.dat ]
    then
        echo 'Add standard calibration to the data with 4 band photo-z'
        python2.7 $srcDir/Apply_multiplicative_patches.py $CatDir/Results/2bin/$surfaceOfDoom $CatDir/Results/ 0 4 > $CatDir/computebias.txt
        echo 'Add standard calibration to the data with 9 band photo-z'
        python2.7 $srcDir/Apply_multiplicative_patches.py $CatDir/Results/2bin/$surfaceOfDoom $CatDir/Results/ 0 9 > $CatDir/computebias.txt
    fi
    
    #Add the calibration to the KiDS data (tomographically)
    if [ "$TomoFlag" == "Yes" ]
    then
	if [ ! -f $CatDir/Results/Summary_multiplicative_tomo4.dat ]
	then
	    echo 'Add tomographic calibration to the data'
	    python2.7 $srcDir/Apply_multiplicative_patches.py $CatDir/Results/2bin/$surfaceOfDoom $CatDir/Results/ 1 4 > $CatDir/computebias4.txt
        python2.7 $srcDir/Money_plot.py $CatDir/Results/ $TomoFlag $CatDir/Results/2bin/$surfaceOfDoom 4 > $CatDir/computebias.txt
	fi
	if [ ! -f $CatDir/Results/Summary_multiplicative_tomo9.dat ]
	then
	    echo 'Add tomographic calibration to the data'
	    python2.7 $srcDir/Apply_multiplicative_patches.py $CatDir/Results/2bin/$surfaceOfDoom $CatDir/Results/ 1 9 > $CatDir/computebias.txt
        python2.7 $srcDir/Money_plot.py $CatDir/Results/ $TomoFlag $CatDir/Results/2bin/$surfaceOfDoom 9 > $CatDir/computebias.txt 
	fi
    fi
    #Make summary plot tomographic bins

    echo $surfaceOfDoom    
    python2.7 $srcDir/Money_plot.py $CatDir/Results/ $TomoFlag $CatDir/Results/2bin/$surfaceOfDoom 9 > $CatDir/computebias.txt 
    python2.7 $srcDir/Money_plot.py $CatDir/Results/ $TomoFlag $CatDir/Results/2bin/$surfaceOfDoom 4 > $CatDir/computebias.txt
#    python2.7 $srcDir/Money_plot.py $CatDir/Results/ "No" $CatDir/Results/2bin/$surfaceOfDoom 9 > $CatDir/computebias.txt 

    #Remove catalogues
#   cd $CatDir/Results
#    rm -f *.cat
    
done
