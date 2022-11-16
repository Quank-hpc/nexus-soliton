#!/bin/bash -l
#SBATCH -p short
#SBATCH -t 00:40:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4G
#SBATCH -o SpinWave.out

module load matlab/r2016b
srun matlab -nosplash -r "SpinDynamics_EigenProblem_PdB_Galerkin; exit(0)" 

