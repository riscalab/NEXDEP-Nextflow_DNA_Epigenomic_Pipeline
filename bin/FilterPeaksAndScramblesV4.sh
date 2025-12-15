#!/bin/bash


#get into correct conda env
source /ru-auth/local/home/risc_soft/miniconda3/etc/profile.d/conda.sh
conda activate fastq2bam

# get arguments
peak=$1
scramble=$2
blacklist=$3
mappa=$4

#Sort each bed for initial overlap calculations
echo "Sorting"
for bed in $@
do
	bedtools sort -i $bed > ${bed##*/}.st
done

#Calculate Jaccards and save to output file
echo "Calculating initial Jaccard"
peakArea=$(cat $peak | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
scrambleArea=$(cat $scramble | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')

outputStats=${peak##*/}.${scramble##*/}.JaccardCalcs.txt
echo "Initial Jaccard Index before any filtering" > $outputStats
echo "set	filter	intersection	union	jaccard	totalArea	percentWithinFilter" >> $outputStats
#$1/$TA
jaccardPeakBlacklist=$(bedtools jaccard -a ${peak##*/}.st -b ${blacklist##*/}.st | awk -v TA="$peakArea" 'FNR==2 {print $1"\t"$2"\t"$3"\t"TA"\t"$1/TA}')
jaccardPeakMappa=$(bedtools jaccard -a ${peak##*/}.st -b ${mappa##*/}.st | awk -v TA="$peakArea" 'FNR==2 {print $1"\t"$2"\t"$3"\t"TA"\t"$1/TA}')

jaccardScramBlacklist=$(bedtools jaccard -a ${scramble##*/}.st -b ${blacklist##*/}.st | awk -v TA="$scrambleArea" 'FNR==2 {print $1"\t"$2"\t"$3"\t"TA"\t"$1/TA}')
jaccardScramMappa=$(bedtools jaccard -a ${scramble##*/}.st -b ${mappa##*/}.st | awk -v TA="$scrambleArea" 'FNR==2 {print $1"\t"$2"\t"$3"\t"TA"\t"$1/TA}')

echo "${peak##*/}	${blacklist##*/}	$jaccardPeakBlacklist" >> $outputStats
echo "${peak##*/}	${mappa##*/}	$jaccardPeakMappa" >> $outputStats
echo "${scramble##*/}	${blacklist##*/}	$jaccardScramBlacklist" >> $outputStats
echo "${scramble##*/}	${mappa##*/}	$jaccardScramMappa" >> $outputStats

#Run validation plots to check GC distribution and size distribution
echo "Plotting GC and Size distribution"
gen="/lustre/fs4/home/ascortea/store/risc_data/downloaded/hg19/genome/Sequence/WholeGenomeFasta/genome.fa"
/lustre/fs4/home/ascortea/Risc_scratch/ascortea/scripts/GitHub_Scripts/main/ValidateNullSeqs/PlotValidationForNullSeqs.sh ${peak##*/}.st ${scramble##*/}.st /lustre/fs4/home/ascortea/store/risc_data/downloaded/hg19/genome/Sequence/WholeGenomeFasta/genome.fa

#Go through filter process
echo "Filtering"
bedtools subtract -A -a ${peak##*/}.st -b ${blacklist##*/}.st > ${peak##*/}.st.blft
bedtools subtract -A -a ${scramble##*/}.st -b ${blacklist##*/}.st > ${scramble##*/}.st.blft

bedtools intersect -u -f 1 -a ${peak##*/}.st.blft -b ${mappa##*/}.st > ${peak##*/}.st.blft.maft
bedtools intersect -u -f 1 -a ${scramble##*/}.st.blft -b ${mappa##*/}.st > ${scramble##*/}.st.blft.maft

#remove assembly patch regions from scramble
awk -f /lustre/fs4/home/ascortea/Risc_scratch/ascortea/scripts/filterNullSeqs.awk ${peak##*/}.st.blft.maft > ${peak##*/}.st.blft.maft.chft
awk -f /lustre/fs4/home/ascortea/Risc_scratch/ascortea/scripts/filterNullSeqs.awk ${scramble##*/}.st.blft.maft > ${scramble##*/}.st.blft.maft.chft

#Calculate distribution of peak sizes and remove the top and bottom 5% to remove outliers

#Get interval sizes
cat ${peak##*/}.st.blft.maft.chft | awk '{print $1"\t"$2"\t"$3"\t"$3-$2}' > ${peak##*/}.st.blft.maft.chft.sz
cat ${scramble##*/}.st.blft.maft.chft | awk '{print $1"\t"$2"\t" $3"\t"$3-$2}' > ${scramble##*/}.st.blft.maft.chft.sz
#sort rows by interval size
sort -k4 -n ${peak##*/}.st.blft.maft.chft.sz > ${peak##*/}.st.blft.maft.chft.sz.st
sort -k4 -n ${scramble##*/}.st.blft.maft.chft.sz > ${scramble##*/}.st.blft.maft.chft.sz.st
#Calculate how many intervals exist in top and bottom 5%
numPeaks=$(wc -l < ${peak##*/}.st.blft.maft.chft.sz.st)
#numScramble=$(wc -l ${scramble##*/}.st.blft.maft.chft.sz.st)

#fivePercentPeak=$(awk -v num="$numPeaks" "BEGIN {print int(num*0.05)}")
#fivePercentScramble=$(awk -v num="$numScramble" "BEGIN {print int(num*0.05)}")

#Since peaks are ordered by interval size, I can remove the top and bottom 5%
#tail -n +$fivePercentPeak ${peak##*/}.st.blft.maft.chft.sz.st | head -n -$fivePercentPeak > ${peak##*/}.st.blft.maft.chft.sz.st.oo
#tail -n +$fivePercentScramble ${scramble##*/}.st.blft.maft.chft.sz.st | head -n -$fivePercentScramble > ${scramble##*/}.st.blft.maft.chft.sz.st.oo

#Now subset the scramble down to the amount of peaks in the real peak set while avoiding overlapping regions
#numPeaks=$(wc -l < ${peak##*/}.st.blft.maft.chft.sz.st.oo)
>${scramble##*/}.st.blft.maft.chft.sz.st.shuf.subset
subsetCount=0
shuf ${scramble##*/}.st.blft.maft.chft.sz.st > ${scramble##*/}.st.blft.maft.chft.sz.st.shuf
while [ $subsetCount -lt $numPeaks ]
do
	#echo "Looping because $subsetCount is less than $numPeaks"
	#shuf ${scramble##*/}.st.blft.maft.chft.sz.st.oo > ${scramble##*/}.st.blft.maft.chft.sz.st.oo.shuf
	head -n 1 ${scramble##*/}.st.blft.maft.chft.sz.st.shuf > ${scramble##*/}.temp.bed
	sed -i '1d' ${scramble##*/}.st.blft.maft.chft.sz.st.shuf
	bedtools intersect -u -a ${scramble##*/}.st.blft.maft.chft.sz.st.shuf.subset -b ${scramble##*/}.temp.bed > ${scramble##*/}.overlap.bed
	overlapCount=$(wc -l ${scramble##*/}.overlap.bed | awk '{print $1}')
	#echo "There are $overlapCount overlaps choosen"
	if [ $overlapCount -eq 0 ]
	then
		#cat ${scramble##*/}.st.blft.maft.chft.sz.st.oo.shuf.subset temp.bed > subset.bed
		#cat subset.bed > ${scramble##*/}.st.blft.maft.chft.sz.st.oo.shuf.subset
		cat ${scramble##*/}.temp.bed >> ${scramble##*/}.st.blft.maft.chft.sz.st.shuf.subset
	fi	
	subsetCount=$(wc -l ${scramble##*/}.st.blft.maft.chft.sz.st.shuf.subset | awk '{print $1}')
	#bedtools subtract -f 1 -a ${scramble##*/}.st.blft.maft.chft.sz.st.oo -b temp.bed > ${scramble##*/}.st.blft.maft.chft.sz.st.oo2
	#cat ${scramble##*/}.st.blft.maft.chft.sz.st.oo2 > ${scramble##*/}.st.blft.maft.chft.sz.st.oo
done

#${scramble##*/}.st.blft.maft.chft.sz.st.oo.shuf.subset
#resort the peaks and scrambles
bedtools sort -i ${peak##*/}.st.blft.maft.chft.sz.st | awk '{print $1"\t"$2"\t"$3}' > ${peak##*/}.filtered.bed
bedtools sort -i ${scramble##*/}.st.blft.maft.chft.sz.st.shuf.subset | awk '{print $1"\t"$2"\t"$3}' > ${scramble##*/}.filtered.bed

peakArea=$(cat ${peak##*/}.filtered.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
scrambleArea=$(cat ${scramble##*/}.filtered.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')

#recalculate the Jaccards and %blacklist and %mappable
echo "Recalculating Jaccard"
echo "Filtered Jaccard Index" >> $outputStats
echo "set	filter	intersection	union	jaccard	totalArea	percentWithinFilter" >> $outputStats
jaccardPeakBlacklist2=$(bedtools jaccard -a ${peak##*/}.filtered.bed -b ${blacklist##*/}.st | awk -v TA="$peakArea" 'FNR==2 {print $1"\t"$2"\t"$3"\t"TA"\t"$1/TA}')
jaccardPeakMappa2=$(bedtools jaccard -a ${peak##*/}.filtered.bed -b ${mappa##*/}.st | awk -v TA="$peakArea" 'FNR==2 {print $1"\t"$2"\t"$3"\t"TA"\t"$1/TA}')

jaccardScramBlacklist2=$(bedtools jaccard -a ${scramble##*/}.filtered.bed -b ${blacklist##*/}.st | awk -v TA="$scrambleArea" 'FNR==2 {print $1"\t"$2"\t"$3"\t"TA"\t"$1/TA}')
jaccardScramMappa2=$(bedtools jaccard -a ${scramble##*/}.filtered.bed -b ${mappa##*/}.st | awk -v TA="$scrambleArea" 'FNR==2 {print $1"\t"$2"\t"$3"\t"TA"\t"$1/TA}')

echo "${peak##*/}.filtered.bed	${blacklist##*/}	$jaccardPeakBlacklist2" >> $outputStats
echo "${peak##*/}.filtered.bed	${mappa##*/}	$jaccardPeakMappa2" >> $outputStats
echo "${scramble##*/}.filtered.bed	${blacklist##*/}	$jaccardScramBlacklist2" >> $outputStats
echo "${scramble##*/}.filtered.bed	${mappa##*/}	$jaccardScramMappa2" >> $outputStats

#Run validation plots to check GC distribution and size distribution
echo "Plotting GC and Size distribution"
gen="/lustre/fs4/home/ascortea/store/risc_data/downloaded/hg19/genome/Sequence/WholeGenomeFasta/genome.fa"
/lustre/fs4/home/ascortea/Risc_scratch/ascortea/scripts/GitHub_Scripts/main/ValidateNullSeqs/PlotValidationForNullSeqs.sh ${peak##*/}.filtered.bed ${scramble##*/}.filtered.bed /lustre/fs4/home/ascortea/store/risc_data/downloaded/hg19/genome/Sequence/WholeGenomeFasta/genome.fa

echo "Done"