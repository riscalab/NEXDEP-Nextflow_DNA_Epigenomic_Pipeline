#!/bin/bash
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -p risc,hpc

#conda env
#source activate rstudio
source /ru-auth/local/home/risc_soft/miniconda3/etc/profile.d/conda.sh

#Main
# Iterate through arguments looking for bams, and for each bam iterate through arguments again looking for beds or peaks to construct command for getting break density within the bed.
#This is a modification of the breakDensityWrapper.sh script which omits the adjusted enrichment, and enrichment values. It will not read normalilze the break counts, but just report the actual break counts themselves.
#breakDensities=()
echo "Making BreakCounts.log\n"
echo "bam	bed	Expected_Density	Observed_Breaks	Read_Norm_Breaks" > densityCalculations.log
for peak in $@
	do 
		echo "iterating"
		if [[ $peak == *.bed ]] || [[ $peak == *Peak ]]
		then
			echo "Working on $peak"
			for bam in $@
				do
					if [[ $bam == *.bam ]]
						then
							echo "Working on $bam"
							echo "expectedDensity"
							conda activate rstudio
							expectedDensity=$(Rscript /lustre/fs4/home/ascortea/Risc_scratch/ascortea/scripts/BreakDensity/getUniformBreakDensity.R $peak hg19)
							echo $expectedDensity
							conda activate fastq2bam
							bash /lustre/fs4/home/ascortea/Risc_scratch/ascortea/scripts/BreakDensity/getBreakDensityInPeaksV3_GLOEseq.sh $bam $peak
							observedBreaks=$(awk '{ total += $4 } END { print total }' ${bam##*/}.${peak##*/}.numBreaksInPeaks.bed)
							echo $observedBreaks
							observedDensity=$(awk '{ total += $4 } END { print total }' ${bam##*/}.${peak##*/}.numBreaksInPeaksNormalized.bed)
							echo $observedDensity
							echo "saving results into densityCalculations.log"
							echo "$bam	$peak	$expectedDensity	$observedBreaks	$observedDensity" >> densityCalculations.log
					fi
				done
		fi
	done
echo "All done!"

# Calculate the uniform break density for each peak file. Hardcoding hg19 for now
# get into correct conda env
#source activate rstudio
#uniformDensities=()
#for peak in $@
#	do
#		if [[ $peak == *.bed ]] || [[ $peak == *Peak ]]
#			then
#				uniformDensities[peak]=$(rscript /lustre/fs4/home/ascortea/Risc_scratch/ascortea/scripts/BreakDensity/getUniformBreakDensity.R $peak hg19)
#		fi
#	done