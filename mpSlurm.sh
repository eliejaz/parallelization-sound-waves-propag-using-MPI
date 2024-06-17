#!/bin/bash
# Submission script for NIC5
#SBATCH --job-name=mpi-job
#SBATCH --time=00:01:00 # hh:mm:ss
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1024 # megabytes
module load releases/2021b
module load OpenMPI/4.1.2-GCC-11.2.0

cd $SLURM_SUBMIT_DIR
mpirun -np $SLURM_NTASKS ./a.out param.txt
