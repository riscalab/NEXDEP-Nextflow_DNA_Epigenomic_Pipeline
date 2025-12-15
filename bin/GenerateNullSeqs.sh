#!/bin/bash



source /ru-auth/local/home/risc_soft/miniconda3/etc/profile.d/conda.sh
conda activate rstudio



bed=$1
if [ ! -f ${bed##*/}.chunked.bed ]
then
	Rscript /lustre/fs4/home/rjohnson/pipelines/hera_pipeline/bin/chunkBedsV2.R $bed 2000 1000
fi
conda activate gkmsvm
Rscript /lustre/fs4/home/rjohnson/pipelines/hera_pipeline/bin/GenerateNullSeqs.R ${bed##*/}.chunked.bed