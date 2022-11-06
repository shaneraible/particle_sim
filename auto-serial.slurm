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

rm serial.txt
srun -n 1 -c 1 ./serial -n 500 -no -s serial.txt
srun -n 1 -c 1 ./serial -n 1000 -no -s serial.txt
srun -n 1 -c 1 ./serial -n 2000 -no -s serial.txt
srun -n 1 -c 1 ./serial -n 4000 -no -s serial.txt
srun -n 1 -c 1 ./serial -n 8000 -no -s serial.txt
./autograder -v serial -s serial.txt
