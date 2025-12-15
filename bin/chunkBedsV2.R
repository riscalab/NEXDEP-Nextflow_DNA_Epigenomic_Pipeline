#This script will take an input bed file and for any regions larger than argument $limit will be chopped up into smaller regions set by $chunkSize. The purpose of this is to make it easier to match and filter very large bed regions (LADs or broadPeaks) to a random scramble region.

#packages
library("tidyverse")
options(scipen=999)

# get arguments
args = commandArgs(trailingOnly=TRUE)

bed <- read_table2(args[1],col_names = FALSE)
#bed <- read_table2("E055-H3K36me3.minimal.narrowPeak", col_names = FALSE)
limit <- as.numeric(args[2])
#limit <- 2000
chunkSize <- as.numeric(args[3])
#chunkSize <- 1000
filename <- basename(as.character(args[1]))
#filename <- "LADs_BJ5ta.minimal.bed"

#Main script
print(paste0("Chunking ", filename, " peakset"))
print("Getting sizes")
bed <- bed[,1:3]
bed$X4 <- bed$X3-bed$X2
print("sorting")
bed <- bed[order(bed$X4),]
#remove top and bottom 5% of region sizes to omit outliers
print("Getting 5%")
count5percent <- floor(nrow(bed)*0.05)
print("removing 5%")
bed <- bed[count5percent:(nrow(bed)-count5percent),] 
print("making new bed")
newbed <- data.frame(X1 = "",X2 = "",X3 = "", X4 = "")
newbed$X1 <- as.character(newbed$X1)
newbed$X2 <- as.character(newbed$X2)
newbed$X3 <- as.character(newbed$X3)
newbed$X4 <- as.character(newbed$X4)
print("extracting regions which don't need chunking")
newbed <- rbind(newbed, bed[bed$X4 < limit,])
print("removing regions which don't need chunking")
bed <- bed[!(bed$X4 < limit),1:4]
print("chunking remaining")
if(nrow(bed) > 0){
  print("more than 0 chunks left")
  for (i in 1:nrow(bed)){
    print(i)
    chunks <- floor(as.numeric(bed[[i,4]])/chunkSize)
    for(chunk in 1:chunks){
      print(paste0("Chunks ",chunk," of ",chunks))
      chr <- bed[[i,1]]
      start <- as.numeric(bed[[i,2]])
      end <- as.numeric(bed[[i,3]])
      newStart <- start+((chunk-1)*chunkSize)
      newEnd <- start+((chunk)*chunkSize)
      newChunk <- c(chr,newStart, newEnd, newEnd-newStart)
      #newChunk <- c(chr,start+((chunk-1)*chunkSize), start+((chunk)*chunkSize-1))
      newbed <- rbind(newbed, newChunk)
    }
    runtStart <- start+((chunk)*chunkSize)
    if (end-runtStart > 0){
      runtChunk <- c(chr,runtStart, end, end-runtStart)
      newbed <- rbind(newbed, runtChunk)
    }
    print(paste0("Done with row ", i, " of ", nrow(bed)))
  }
}
#newbed$X4 <- as.numeric(newbed$X3)-as.numeric(newbed$X2)
write_tsv(newbed[-1,1:3], paste0(filename,".chunked.bed"), col_names = FALSE)
print("Done chunking")
