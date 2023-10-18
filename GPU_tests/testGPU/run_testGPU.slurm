#!/bin/tcsh
#SBATCH --partition=GPU
#SBATCH --gres=gpu:2080Ti:1
#SBATCH --time=00:10:00
#SBATCH --ntasks=1

set LogFile=output.log
env > $LogFile
groups >> $LogFile
./testGPU >> $LogFile