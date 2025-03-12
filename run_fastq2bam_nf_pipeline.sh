#!/bin/env bash

#SBATCH --mem=20GB
#SBATCH --mail-type=FAIL,END
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=nextflow_chip

#source $HOME/.bashrc_rj_test.sh   # use this it works also but not for others

source /lustre/fs4/home/rjohnson/.bashrc_rj_test.sh
# source /ru-auth/local/home/rjohnson/.bashrc_rj_test.sh # or use this, should be the same thing 

conda activate nextflow_three

########## for SE data ###############
# --test: parameter will make the pipeline only take 3 of your fastq files(or fastq pairs in pair end) in your directory of many fastq files. without this the pipeline will run and process all of your fastq files
# --ATAC : if you have atac-seq data, please specify this parameter
# --SE parameter for pair end reads
# when using SE do --single_end_reads and give the path to your single end reads with a glob pattern if you have other files in that directory you dont want ex: path/to/single_end_reads/*file*.fastq
# if you have an adapter sequence use --ada_seq, then specify the sequence with --adapter_seq_str which will take the string sequence
# use --genome and give the path including the file of your genome of choice
# --BL  : this parameter is to tell the pipeline use the process that filters for black list regions
# --blacklist_path : give the path to the blacklist bed file you have and include the file in the path. the defualt used is a path to the hg19 v2 black list. so if using a different species or a different human genome use the correct blacklist and not the default.
# --use_effectiveGenomeSize : this should be called if you want the pipeline to use the path where deeptools bamcoverage will take the effective genome size. this parameter does not take the number see next parameter
# --num_effectiveGenomeSize : if you used the parameter --use_effectiveGenomeSize then you need to use this one also. this one takes the number and you can go to deeptools website to find the correct effective genome size number to use here: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
# --spike_in : lets the pipeline know it should run the workflow for spike ins
# --t7, --lambda : choose one or both and the single end spike in workflow for these will be executed. Have to use with parameter --spike_in
# --yeast : choose the yeast spike in if you have ricc_seq data or data that wants yeast spike in, and the single end spike end will run. must use with parameter --spike_in

# --give_peak_files : this is a parameter for the nasa project. Put all peaks in a directory and give the glob pattern '*.narrowPeak' for the pipeline to find your narrow peak files and use them. give the absolute path followed by the glob pattern. if this is not specified then the pipeline will used the peak files that I choose from the data I have access to.
# --depth_intersection : this new parameter should be for anyone that has alignment bam files and want to check it's depth by seeing how many reads intersect with a given set of already created peak files


########################################

########### for PE data ##############
# I will not put an option to specify adapters and put your own sequence for the Pair End part of this pipeline
# The reason being i specified in fastp that we will look adapters for PE and just trim them. read the parameters used for the fastp tool in the fastp_PE process

# --test: parameter will make the pipeline only take 3 of your fastq files(or fastq pairs in pair end) in your directory of many fastq files. without this the pipeline will run and process all of your fastq files
# --ATAC : if you have atac-seq data, please specify this parameter
# --genome : give the path to the reference genome file you want to use. if you dont specify then the defualt hg19 genome will be used
# --PE : lets the pipeline know you have pair end data
# --paired_end_reads : you need to then also use this parameter to specify the path and the glob pattern to get your forward and reverse reads together in the same input channel. EX: path/to/pair_end_reads/*my_pair_end_file*_{R1,R2}*.fastq 
# --BL : you need to specify BL(black list) if you want the pipeline to do blacklist region filtering.
# --blacklist_path : once you choose --BL, use this parameter to specify the blacklist bed file if you changed the genome from the defualt genome. you dont have to do this if you didnt use the --genome parameter to choose a different genome.
# --use_effectiveGenomeSize : this should be called if you want the pipeline to use the path where deeptools bamcoverage will take the effective genome size. this parameter does not take the number see next parameter
# --num_effectiveGenomeSize : if you used the parameter --use_effectiveGenomeSize then you need to use this one also. this one takes the number as a str and you can go to deeptools website to find the correct effective genome size number to use that here: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
# --spike_in : lets the pipeline know it should run the workflow for spike ins
# --t7, --lambda : choose one or both and the pair end spike in workflow for these will be executed. Have to use with parameter --spike_in
# --yeast : choose the yeast spike in if you have ricc_seq data or data that wants yeast spike in, and the pair end spike end will run. must use with parameter --spike_in

# --give_peak_files : this is a parameter for the nasa project. Put all peaks in a directory and give the glob pattern '*.narrowPeak' for the pipeline to find your narrow peak files and use them. give the absolute path followed by the glob pattern. if this is not specified then the pipeline will used the peak files that I choose from the data I have access to.
# --depth_intersection : this new parameter should be for anyone that has alignment bam files and want to check it's depth by seeing how many reads intersect with a given set of already created peak files
######################################

# NOTE: path to the peak files /lustre/fs4/home/ascortea/store/ascortea/beds
     # then get the different peak directories IMR90, k562, BJ subdirectories hold these peak files

# NEW NOTE: I want to add another process or workflow where i take all the bam_index_tuple_ch that made it to the end of the fastq to bam pipeline and send them to fastqc then multi-qc to get a good html file showing the stats.
#          well neither fastqc nor multiqc takes bam files to be able to do this

#nextflow run fastq2bam_nextflow_pipeline.nf -profile 'nasa_pipeline ' \
#-resume \
#--SE \
#--single_end_reads '/rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs/*_1.fastq.gz' \
#--ada_seq --adapter_seq_str 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA' \
#--BL \
#--spike_in \
#--t7 \
#--lambda




nextflow run fastq2bam_nextflow_pipeline.nf -profile 'fastq2bam2_pipeline' \
-resume \
--PE \
--BL \
--paired_end_reads '/rugpfs/fs0/risc_lab/store/hcanaj/HC_GLOEseq_Novaseq_010925/fastqs_read1_read2/*_{R1,R2}*' \
--use_effectiveGenomeSize \
--num_effectiveGenomeSize '2864785220' \
--calc_break_density \
--spike_in \
--t7 \
--lambda \
--depth_intersection

#--ATAC

 

 # NOTE: If you want to make your own nextflow diagram to see how the pipeline works run this command

# -with-dag fastq2bam_nf_pipeline_flowchart.png{.html,.svg,.pdf}, it seems like this alone overwrites the resume option
# add -preview to render it without having to run the pipeline 
# i like the html best but pdf or png might be a good second