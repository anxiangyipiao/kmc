#!/bin/bash
#SBATCH -p cpu
#SBATCH -J KMC++
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -N 1
#SBATCH -n 2


#module load vasp612
#module load vasp6_cconstr

cd $SLURM_SUBMIT_DIR
export I_MPI_ADJUST_REDUCE=3

date
./kmcplus.x  input.yaml > log
date


