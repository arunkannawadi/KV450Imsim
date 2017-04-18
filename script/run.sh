#!/bin/bash
# KiDS-450+ Simulations

echo "Starting KiDS-450+ Simulated Field Run. Output logs stored in output.log"

python /disks/shear15/KiDS/ImSim/pipeline/code/pipeline.py --g1=-0.04 --g2=0.00 --psfrange=0 --config=/disks/shear15/KiDS/ImSim/pipeline/config/config.ini --rmtemp=False 2> errorOutput.log > output.log &
