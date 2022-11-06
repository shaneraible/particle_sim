#!/bin/bash -l

module purge
module load gcc/5.4.0

#SBATCH --time=00:20:00					# adjust the amount of time as necessary
#SBATCH --job-name=serial
#SBATCH --account=ece1570_2022f
#SBATCH --output auto-serial.%j.stdout
#SBATCH --error auto-serial.%j.error
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cluster=smp
#SBATCH --partition=smp

export OMP_NUM_THREADS=24
srun -n 1 -c 24 ./openmp -n 500 -o openmp.txt
