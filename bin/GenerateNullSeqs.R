library(gkmSVM)
library(BSgenome.Hsapiens.UCSC.hg19.masked)

args = commandArgs(trailingOnly=TRUE)
#Should only be 1 argument

name <- basename(as.character(args[1]))
outputName <- paste0("Null_Range_", name)
outputPosFastaName <- paste0("PosFasta_", name)
outputNegFastaName <- paste0("NegFasta_", name)

genNullSeqs(args[1], outputBedFN = outputName, outputPosFastaFN = outputPosFastaName, outputNegFastaFN = outputNegFastaName, nMaxTrials = 1000, xfold = 2)