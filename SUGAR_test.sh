#!/bin/bash
#SBATCH -o stderr/SG_job.%j.%N.out
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=all
#SBATCH --export=NONE
#SBATCH --time=167:00:00
#SBATCH -A 32cores 
##SBATCH -A Serial
#SBATCH --mem=600000

export OMP_NUM_THREADS=32
ulimit -c unlimited
#time ./SUGAR LGC input/sugar_lc_kitaura.ini input/sugar_lgc_kitaura.param
time ./SUGAR LGC input/sugar_lc.ini input/sugar_lgc.param
