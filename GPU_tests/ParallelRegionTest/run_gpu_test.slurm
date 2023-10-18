#!/bin/tcsh
#SBATCH --partition=GPU
#SBATCH --gres=gpu:2080Ti:1
#SBATCH --time=00:05:00
#SBATCH --ntasks=1

set LogFile=gpu_test.log
env > $LogFile
groups >> $LogFile
./parallel_region_test_gpu 256 100 512 >> $LogFile
