#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p hpxg
SLURM_JOB_NAME=zhouwei_zouyuxiong
export MATLAB_PREFDIR=/tmp/zhouwei_matlab
module load matlab/R2017b 
matlab -nosplash -nodisplay -nodesktop < local_porosity.m