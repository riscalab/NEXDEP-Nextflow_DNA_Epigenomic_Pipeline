#!/bin/env bash

#SBATCH --mem=20GB
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=lei_ricc
#SBATCH --partition=hpc_l40_a

#source $HOME/.bashrc_rj_test.sh   # use this it works also but not for others

source /lustre/fs4/home/rjohnson/.bashrc_rj_test.sh
# source /ru-auth/local/home/rjohnson/.bashrc_rj_test.sh # or use this, should be the same thing 

conda activate nextflow_three


nextflow run fastq2bam_nextflow_pipeline.nf -profile 'fastq2bam2_pipeline' \
-resume \
--PE \
--BL \
--blacklist_path '/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/blacklist/hg38-blacklist.v2.bed' \
--expr_type 'lei_ricc' \
--paired_end_reads '/lustre/fs8/risc_lab/scratch/lnie/project/h1/ricc/ricc_20260921/fastq/*_{R1,R2}*' \
--experiment_type_field_num 3 \
--condition_type_field_num 0 \
--replicate_type_field_num 4 \
--lane_type_field_num 6 \
--short_reads \
--gloe_seq \
--use_effectiveGenomeSize \
--num_effectiveGenomeSize '2864785220' \
--calc_break_density \
--depth_intersection \
--risc_hotel_bank_account \
--hpc_partition 'hpc_l40_a'