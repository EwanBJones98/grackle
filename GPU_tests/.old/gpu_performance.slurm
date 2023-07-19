#!/bin/tcsh
#SBATCH --partition=GPU
#SBATCH --gres=gpu:2080Ti:1
#SBATCH --time=00:05:00
#SBATCH --ntasks=1024
#SBATCH -o stdout

set LogFile=Log.$SLURM_JOBID
env > $LogFile
groups >> $LogFile
./deviceQuery >> $LogFile