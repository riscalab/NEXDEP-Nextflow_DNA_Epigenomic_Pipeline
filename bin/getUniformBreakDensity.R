#Arguments
args = commandArgs(trailingOnly=TRUE)

bed <- args[1]
genome <- args[2]

# genome options
# values taken from https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
GRCh37=2736124973
GRCh38=2747877777

if(genome == "hg38"){
  #print("Chosen genome is hg38")
  genomeSize <- GRCh38
}else{
  #print("Chosen genome is hg19")
  genomeSize <- GRCh37
}

# Calculating genomic area covered by peaks (area of peak/area of mappable genome)
library(readr)
bedFile <- read_delim(bed,
                      delim = "\t",
                      escape_double = FALSE,
                      col_names = FALSE,
                      trim_ws = TRUE)

widths <- bedFile[3]-bedFile[2]
totalArea <- sum(widths)
genomeCovered <- totalArea/genomeSize

cat(as.numeric(genomeCovered))






