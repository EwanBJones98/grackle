#!/bin/tcsh
#SBATCH --partition=GPU
#SBATCH --gres=gpu:2080Ti:1
#SBATCH --time=00:05:00
#SBATCH --ntasks=32

set LogFile=eric_gpu_test.log
env > $LogFile
groups >> $LogFile
./eric-gpu-test >> $LogFile