#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p hpxg
SLURM_JOB_NAME=zhouwei_zouyuxiong
export PATH=/home/yxzou/project/voro++/bin:$PATH
export MANPATH=/home/yxzou/project/voro++/man:$MANPATH
export MATLAB_PREFDIR=/tmp/zhouwei_matlab
module load matlab/R2017b 

matlab -nosplash -nodisplay -nodesktop -r "frame=50; Cal_voro_vectors" 

