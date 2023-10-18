#!/bin/tcsh
#SBATCH --job-name=testCPU
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:05:00

setenv OMP_NUM_THREADS $SLURM_JOB_CPUS_PER_NODE

set LogFile=cpu_test.log
./parallel_region_test_cpu 256 1000 8 >> $LogFile
