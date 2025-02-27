#!/bin/bash
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -p risc,hpc

#conda env
#source activate rstudio
source /ru-auth/local/home/risc_soft/miniconda3/etc/profile.d/conda.sh

#Main
# Iterate through arguments looking for bams, and for each bam iterate through arguments again looking for beds or peaks to construct command for getting break density within the bed.
#breakDensities=()
echo "Making densityCalculations.log\n"
echo "bam	bed	Expected_Density	Observed_Density	Enrichment" > densityCalculations.log
for bam in $@
	do 
		echo "iterating"
		if [[ $bam == *.bam ]]
		then
			echo "Working on $bam"
			for peak in $@
				do
					if [[ $peak == *.bed ]] || [[ $peak == *Peak ]]
						then
							echo "Working on $peak"
							echo "expectedDensity"
							conda activate rstudio
							expectedDensity=$(Rscript /lustre/fs4/home/ascortea/Risc_scratch/ascortea/scripts/BreakDensity/getUniformBreakDensity.R $peak hg19)
							echo $expectedDensity
							conda activate fastq2bam
							bash /lustre/fs4/home/ascortea/Risc_scratch/ascortea/scripts/BreakDensity/getBreakDensityInPeaksV3.sh $bam $peak
							observedDensity=$(awk '{ total += $4 } END { print total }' ${bam##*/}.${peak##*/}.numBreaksInPeaksNormalized.bed)
							echo $observedDensity
							enrichment=$(awk -v obs="$observedDensity" -v expt="$expectedDensity" 'BEGIN {print obs/expt}')
							echo $enrichment
							echo "saving results into densityCalculations.log"
							echo "$bam	$peak	$expectedDensity	$observedDensity	$enrichment" >> densityCalculations.log
					fi
				done
		fi
	done

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