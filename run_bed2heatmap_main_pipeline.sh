#!/bin/env bash

#SBATCH --mem=20GB
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=hpc_a10_a
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=nextflow_chip

#source $HOME/.bashrc_rj_test.sh   # use this it works also but not for others

source /lustre/fs4/home/rjohnson/.bashrc_rj_test.sh
# source /ru-auth/local/home/rjohnson/.bashrc_rj_test.sh # or use this, should be the same thing 

conda activate nextflow_three




# so some of the commands should be the same 
# I still want to differentiate between single end, end_seq, spike_in, t7, lambda

nextflow run bed2heatmap_main_pipeline.nf -profile 'bed2heatmap' \
-resume \
--SE \
--spike_in \
--t7 \
--lambda \
--end_seq \
