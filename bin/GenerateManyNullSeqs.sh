#!/bin/bash




for bam in $@
do
     /lustre/fs4/home/rjohnson/pipelines/hera_pipeline/bin/GenerateNullSeqs.sh $bam
done
