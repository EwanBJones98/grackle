#!/bin/tcsh
#SBATCH --job-name=testSeri
#SBATCH --time=00:05:00

set LogFile=serial_test.log
./parallel_region_test_serial 1000 20 >> $LogFile