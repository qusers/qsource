#!/bin/bash
#set -x
#trap read debug

export rcran=/home/apps/R
export PATH=$rcran/bin:$PATH

echo "The following is the Q heating time in minutes: "
    grep "Total time of main loop:" eq*.log | awk '{print $6}' | sed 's/(s)//g' \
    | Rscript -e "A <- read.table(pipe('cat /dev/stdin')); (sum(A)/60)"

echo "The following is the Q equilibration time in hours: "
    grep "Total time of main loop:" dc*.log | awk '{print $6}' | sed 's/(s)//g' \
    | Rscript -e "A <- read.table(pipe('cat /dev/stdin')); (sum(A)/60)/60"

echo "The following is the Q heat+equilibration time in hours: "
    grep "Total time of main loop:" *.log | awk '{print $6}' | sed 's/(s)//g' \
    | Rscript -e "A <- read.table(pipe('cat /dev/stdin')); (sum(A)/60)/60"
