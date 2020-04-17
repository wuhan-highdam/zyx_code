#!/bin/bash
#SBATCH -N 1
#SBATCH -n 7
#SBATCH -p hpxg
SLURM_JOB_NAME=zhouwei_zouyuxiong
export PATH=/home/yxzou/project/voro++/bin:$PATH
export MANPATH=/home/yxzou/project/voro++/man:$MANPATH
export MATLAB_PREFDIR=/tmp/zhouwei_matlab
module load matlab/R2017b 

mkdir cell_vectors_file
cp -r ./frame_all/*.sample ./cell_vectors_file/
cp -r ./cell_area_file/words4voro.sh ./cell_vectors_file/
matlab -nosplash -nodisplay -nodesktop < sh_change.m
cd ./cell_vectors_file
bash ./words4voro.sh
cd ../
rm -rf ./cell_vectors_file/*.sample

cp -r ./frame_all/*.sample ./cell_area_file/
cd ./cell_area_file
bash ./words4voro.sh
cd ../
rm -rf ./cell_area_file/*.sample

cp -r ./frame_all/*.sample ./surface_node_file/
cd ./surface_node_file
bash ./words4voro.sh
cd ../
rm -rf ./surface_node_file/*.sample
