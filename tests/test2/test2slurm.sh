#!/bin/bash
# Remember to change here the number tasks (-n 6) to what is available
# and reasonable in the queue you're calling. Perhaps you'll also need
# to issue the --exclusive option.
sbatch -n 16 -A SNIC2014-11-2 -t 02:00:00 -J testsmooth smooth.sh

