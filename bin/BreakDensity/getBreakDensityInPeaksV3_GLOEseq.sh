#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p risc,hpc

#get into correct conda env
source activate fastq2bam

# get arguments to pass into python
bam=$1

# sort bam file in case it isn't already sorted
if [ -f "${bam##*/}.sorted.bam" ]
then
    echo "Sorted bam found, skipping sorting."
else 
	echo "Sorted bam not found. Sorting bam now."
    samtools sort -o ${bam##*/}.sorted.bam -n $bam
fi

# filter the bam for only read 2 (which is the read that contains the break location, 3'OH adapted with i7 adapter -> this is read 2)
samtools view -f 128 -b ${bam##*/}.sorted.bam > ${bam##*/}.sorted.R2.bam

# get bed regions from sorted bam and then crop it to just be intervals of insert locations, skip if already done
if [ -f "${bam##*/}.breaks.bed" ]
then
    echo "Breaks bed file found, skipping."
else 
	echo "Breaks bed file not found, getting break locations now."
    bedtools bamtobed -i ${bam##*/}.sorted.R2.bam > ${bam##*/}.R2.bed
	#touch ${bam##*/}.breaks.bed
	#Extract start postion as break location for + strand ends, and end position as break location for - strand ends, then invert the sign since break occurs on opposite strand of sequencing readout
	awk '{
	if ( $6 == "+" )
		print $1"\t"$2"\t"$2"\tNA\t0\t-"
	else if ( $6 == "-" )
		print $1"\t"$3"\t"$3"\tNA\t0\t+"
	}' ${bam##*/}.R2.bed > ${bam##*/}.breaks.bed
	#awk '{print $1"\t"$2"\t"$2}' ${bam##*/}.PE.bed > ${bam##*/}.breaks.bed
	#awk '{print $1"\t"$6"\t"$6}' ${bam##*/}.PE.bed >> ${bam##*/}.breaks.bed
	bedtools sort -i ${bam##*/}.breaks.bed > ${bam##*/}.breaks.sorted.bed
	#Make separate beds for breaks within + and - strand to look for symmetry or break density variation in a stranded manner
	awk '{
	if ( $6 == "+" )
		print $1"\t"$2"\t"$3"\tNA\t0\t"$6
	}' ${bam##*/}.breaks.sorted.bed > ${bam##*/}.breaks.sorted.plus.bed
	awk '{
	if ( $6 == "-" )
		print $1"\t"$2"\t"$3"\tNA\t0\t"$6
	}' ${bam##*/}.breaks.sorted.bed > ${bam##*/}.breaks.sorted.minus.bed
fi
#bedtools bamtobed -bedpe -i ${bam##*/}.sorted.bam > ${bam##*/}.PE.bed
#awk '{print $1"\t"$2"\t"$2}' ${bam##*/}.PE.bed > ${bam##*/}.breaks.bed
#awk '{print $1"\t"$6"\t"$6}' ${bam##*/}.PE.bed >> ${bam##*/}.breaks.bed

# Find breaks within peaks, any arguments that are .bed or Peak files will be used here.
numBreaks=$(wc -l ${bam##*/}.breaks.bed)
echo "numBreaks"
echo $numBreaks
for peaks in $@
	do 
		if [[ $peaks == *.bed ]] || [[ $peaks == *Peak ]]
		then
			echo $peaks
			bedtools sort -i $peaks > $peaks.sorted
			bedtools intersect -sorted -c -a $peaks.sorted -b ${bam##*/}.breaks.sorted.bed > ${bam##*/}.${peaks##*/}.numBreaksInPeaks.bed
			bedtools intersect -sorted -c -a $peaks.sorted -b ${bam##*/}.breaks.sorted.plus.bed > ${bam##*/}.${peaks##*/}.numBreaksInPeaks.plus.bed
			bedtools intersect -sorted -c -a $peaks.sorted -b ${bam##*/}.breaks.sorted.minus.bed > ${bam##*/}.${peaks##*/}.numBreaksInPeaks.minus.bed
			echo "done counting breaks"
			#awk -v SF=1 '{printf($1"\t"$2"\t"$3"\t"$4*SF"\n")}' ${bam##*/}.${peaks##*/}.numBreaksInPeaks.bed > ${bam##*/}.${peaks##*/}.numBreaksInPeaksCrop.bed
			awk -v SF="$numBreaks" '{printf($1"\t"$2"\t"$3"\t"$NF/SF"\n")}' ${bam##*/}.${peaks##*/}.numBreaksInPeaks.bed > ${bam##*/}.${peaks##*/}.numBreaksInPeaksNormalized.bed
			awk -v SF="$numBreaks" '{printf($1"\t"$2"\t"$3"\t"$NF/SF"\n")}' ${bam##*/}.${peaks##*/}.numBreaksInPeaks.plus.bed > ${bam##*/}.${peaks##*/}.numBreaksInPeaksNormalized.plus.bed
			awk -v SF="$numBreaks" '{printf($1"\t"$2"\t"$3"\t"$NF/SF"\n")}' ${bam##*/}.${peaks##*/}.numBreaksInPeaks.minus.bed > ${bam##*/}.${peaks##*/}.numBreaksInPeaksNormalized.minus.bed
			awk '{print "Number of breaks total " $numBreaks}' > ${bam##*/}.${peaks##*/}OverlapStats.log
			awk '{print "Number of breaks total " $numBreaks}' > ${bam##*/}.${peaks##*/}OverlapStats.plus.log
			awk '{print "Number of breaks total " $numBreaks}' > ${bam##*/}.${peaks##*/}OverlapStats.minus.log
			awk '{ total += $NF } END { print "Number of breaks within peaks " total }' ${bam##*/}.${peaks##*/}.numBreaksInPeaks.bed >> ${bam##*/}.${peaks##*/}OverlapStats.log
			awk '{ total += $NF } END { print "Number of breaks within peaks " total }' ${bam##*/}.${peaks##*/}.numBreaksInPeaks.plus.bed >> ${bam##*/}.${peaks##*/}OverlapStats.plus.log
			awk '{ total += $NF } END { print "Number of breaks within peaks " total }' ${bam##*/}.${peaks##*/}.numBreaksInPeaks.minus.bed >> ${bam##*/}.${peaks##*/}OverlapStats.minus.log
			awk -v SF="$numBreaks" '{ total += $NF } END { print "Number of breaks within peaks normalized by number of breaks " total/SF }' ${bam##*/}.${peaks##*/}.numBreaksInPeaks.bed >> ${bam##*/}.${peaks##*/}OverlapStats.log
			awk -v SF="$numBreaks" '{ total += $NF } END { print "Number of breaks within peaks normalized by number of breaks " total/SF }' ${bam##*/}.${peaks##*/}.numBreaksInPeaks.plus.bed >> ${bam##*/}.${peaks##*/}OverlapStats.plus.log
			awk -v SF="$numBreaks" '{ total += $NF } END { print "Number of breaks within peaks normalized by number of breaks " total/SF }' ${bam##*/}.${peaks##*/}.numBreaksInPeaks.minus.bed >> ${bam##*/}.${peaks##*/}OverlapStats.minus.log
		fi
	done