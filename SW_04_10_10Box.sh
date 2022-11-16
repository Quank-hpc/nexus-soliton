#!/bin/bash -l
#SBATCH -p batch
#SBATCH -t 03:45:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=4G
#SBATCH -o SW_04_10_10Box.out
 
module load matlab/r2016b
srun matlab_multithread -nosplash -r "SolitonWall_04_10_10_Box(24); exit(0)"
