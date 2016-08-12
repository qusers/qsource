#!/bin/bash
#sbatch -n 8 -w machines  -t 6:00:00 -J qtest3 smooth.sh run
module load gcc/5.3.0
sbatch -n 8 -t 24:00:00 -J testdhfr smooth.sh 
