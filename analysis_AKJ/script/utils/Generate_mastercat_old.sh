#!/bin/bash

################################################################################
name_run=$1 #TSTGr0729trueposCHK
echo "Run ID:" $name_run
################################################################################

MainDir=/disks/shear15/KiDS/ImSim/analysis/grid/
CatDir=/disks/shear15/KiDS/ImSim/pipeline/archive/$name_run/
TmpDir=/disks/shear15/KiDS/ImSim/temp/$name_run/
srcDir=$MainDir/python/
ArchDir=/disks/shear15/KiDS/ImSim/pipeline/archive/

cd $CatDir

#Generate a list of the available runs in the main directory. 

if [ ! -f $ArchDir/list_$(echo $name_run).t ]
then
    echo "Created a new list of directories ... "
    ls -d */ > $ArchDir/$(echo $name_run)/list_$(echo $name_run).t
fi

##Do the weight bias correction (weight recalibration) at a PSF level
for PSF in 28 #0 1 2 3 4
do
    if [ ! "$name_run" == FC17SBmcal ] || [ ! "$name_run" == TSTnewinputmcal ]
    then
        echo "NOT Doing weight recalibration"
        #python2.7 /disks/shear15/KiDS/ImSim/pipeline/utils/recal_LFweights.py $CatDir $CatDir/list_$name_run.t $PSF
    fi

    if [ "$name_run" == FC17SBmcal ] || [ "$name_run" == TSTnewinputmcal ] || [ "$name_run" == FC17gridmcal ] || [ "$name_run" == TSTnewinpGRscramemcal ]
    then
        echo "Doing metacal weight recalibration"
        python2.7 /disks/shear15/KiDS/ImSim/pipeline/utils/recal_LFweights.py $CatDir $CatDir/list_$name_run.t $PSF 1
    fi
done

#Renaming of the unrotated catalogue such that it can be treated in the same
#way as the rotated cataogues in the subsequent analysis. 

while read file
do
    #cp $file/prior $file/prior_new
    if [ "$name_run" == FC17 ] || [ "$name_run" == FC17SB ] || [ "$name_run" == FC17SB77 ] || [ "$name_run" == FC17V23102015 ]
    then
        if [ ! -f $file/prior_new ]
        then
            echo 'expanding prior file'
            python2.7 $srcDir/expand_prior_file.py $file/prior $file/prior_new
        else
            echo 'NOT expanding the prior to prior_new, because it already existed'
        fi
    else
        #mv $file/prior_new $file/prior_old
        #cp $file/prior_test $file/prior_new
        if [ ! -f $file/prior_new ]
        then
            echo 'copying the prior to prior_new'
            cp $file/prior $file/prior_new
            echo 'done'
        else
            echo 'NOT copying the prior to prior_new, because it already existed'
        fi
        #ln -s $file/prior $file/prior_new
    fi

    ## For FC17 alone, or others from the past
    if [ "$name_run" == FC17 ] || [ "$name_run" == FC17V23102015 ]
    then
        ## Create a symoblic link to the one and only SE catalogue file
        ln -s $CatDir/$file/sex.cat $CatDir/$file/sexrot00.cat
        ln -s $CatDir/$file/sex.cat $CatDir/$file/sexrot01.cat
        ln -s $CatDir/$file/sex.cat $CatDir/$file/sexrot02.cat
        ln -s $CatDir/$file/sex.cat $CatDir/$file/sexrot03.cat
    else
        echo "copying the SExtractor catalogues"
        cp $TmpDir/$file/sexrot01.cat $CatDir/$file/
        cp $TmpDir/$file/sexrot02.cat $CatDir/$file/
        cp $TmpDir/$file/sexrot03.cat $CatDir/$file/
        cp $TmpDir/$file/sex.cat $CatDir/$file/sexrot00.cat
    fi

    #if [ ! -f $file/00.output.rot.fits.asc.scheme2b_corr ]
    #then
	#mv $file/output.fits.asc.scheme2b_corr $file/00.output.rot.fits.asc.scheme2b_corr
    #fi
    if [ ! -f $file/00.output.rot.fits.asc ]
    then
        mv $file/output.fits.asc $file/00.output.rot.fits.asc
    fi
done < $ArchDir/$(echo $name_run)/list_$(echo $name_run).t

#Create the master fits files (one per PSF set). These files contain the LF outputs
#as well as the input galaxy properties.

echo "Looping over PSFs now ..."
for PSF in 28 #0 1 2 3 #4
do
    echo MasterCat_$(echo $name_run)_set_$PSF.fits
    #if [ ! -f MasterCat_$(echo $name_run)_set_$PSF.fits ]
    #then
        echo "Executing the Python script"
        echo $srcDir/create_master_file_oldPrior.py $CatDir $PSF  $ArchDir/$(echo $name_run)/MasterCat_$(echo $name_run)  $ArchDir/$(echo $name_run)/list_$(echo $name_run).t
        python2.7 $srcDir/create_master_file_oldPrior.py $CatDir $PSF  $ArchDir/$(echo $name_run)/MasterCat_$(echo $name_run)  $ArchDir/$(echo $name_run)/list_$(echo $name_run).t > $ArchDir/$name_run/mastercat_$PSF.txt
        #python2.7 $srcDir/create_master_file_oldPrior.py $CatDir $PSF  $ArchDir/TSTnewinputDEIMOS/MasterCat_TSTnewinputDEIMOS  $ArchDir/$(echo $name_run)/list_$(echo $name_run).t > $ArchDir/$name_run/mastercat_$PSF.txt
        if [ ! "$name_run" == FC17 ] && [ ! "$name_run" == FC17SB ] && [ ! "$name_run" == FC17gridmcal ]
        then
            python2.7 $srcDir/propagate_chi2nu.py MasterCat_$(echo $name_run)_set_$PSF.fits ## add reduced chi-square from Griffith catalogue
            python2.7 /disks/shear15/KiDS/ImSim/pipeline/utils/assign_ZB.py  $name_run  $PSF ## correct the 9-band ZB in the catalogues from the COSMOS-photoz cat
        fi
    #fi
    python2.7 /disks/shear15/KiDS/ImSim/pipeline/utils/assign_ZB.py  $name_run  $PSF ## correct the 9-band ZB in the catalogues from the COSMOS-photoz cat
done

#Merge the fits files relative to the different PSF sets and generate 
#5 fits files one for each tomographic bin (i.e. we split the galaxies 
#and their LF measurements according to their ZB value from the 
#COSMOS -input- catalogue)

#Check the number of PSF sets available 
numPSF=`ls *set*.fits | wc -l`
echo "Number of PSF sets = " $numPSF

if [ ! -f MasterCat_$(echo $name_run)_all_$(echo $numPSF)_PSF.fits ]
then
    echo 'Merging the catalogues'
    python2.7 $srcDir/merge_tables.py $numPSF $CatDir MasterCat_$(echo $name_run) MasterCat_$(echo $name_run)_all_$(echo $numPSF)_PSF.fits > $ArchDir/$name_run/mastercat.txt
fi

echo 'Done!'
