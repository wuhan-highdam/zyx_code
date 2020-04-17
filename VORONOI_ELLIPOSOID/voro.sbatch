#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p hpxg
SLURM_JOB_NAME=zhouwei_zouyuxiong
export PATH=/home/yxzou/project/voro++/bin:$PATH
export MANPATH=/home/yxzou/project/voro++/man:$MANPATH
export MATLAB_PREFDIR=/tmp/zhouwei_matlab
module load matlab/R2017b 
matlab -nosplash -nodisplay -nodesktop < voro_data.m
cd ./frame_all
cp -r ./*.sample ../surface_node_file
#rm -rf ./*.sample
cd ../
cp -r ./words4voro.sh ./surface_node_file
cd ./surface_node_file
bash words4voro.sh
cp -r ./words4voro.sh ../cell_area_file
cp -r ./*.sample ../cell_area_file
rm -rf ./*.sample
cd ../cell_area_file
bash words4voro.sh
cp -r ./words4voro.sh ../cell_vol_file
cp -r ./*.sample ../cell_vol_file
rm -rf ./*.sample
cd ../cell_vol_file
bash words4voro.sh
rm -rf ./*.sample
cd ../
matlab -nosplash -nodisplay -nodesktop < Cal_voro.m
rm -rf ./cell_area_file/*.vol
rm -rf ./cell_vol_file/*.vol
rm -rf ./surface_node_file/*.vol
rm -rf ./words4voro.sh