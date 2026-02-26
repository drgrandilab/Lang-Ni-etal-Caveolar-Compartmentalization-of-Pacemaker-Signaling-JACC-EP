#!/bin/bash
#SBATCH --job-name=16_CPU
#SBATCH --output=res_omp.txt
#SBATCH --error=err_omp.txt
##SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# python run_AF_1Hz.py 3


# sh get_all.sh


# sh run_multiple.sh  > sim.log


./SAN_pace
