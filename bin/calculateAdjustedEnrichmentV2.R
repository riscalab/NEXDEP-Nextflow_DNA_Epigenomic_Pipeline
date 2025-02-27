#Takes in TSV of enrichment values and will divide the cells by the PLCs enrichment values to normalize for biases captured in the PLC. Ideally the PLC will be a merge of all respective PLC replicates to have the best dataset.
library(tidyverse)
library(stringr)

#args
args = commandArgs(trailingOnly=TRUE)
#Should only be 1 argument

enrichmentScores <- read_delim(args[1], 
                                                      delim = "\t", escape_double = FALSE, 
                                                      trim_ws = TRUE)

#enrichmentScores <- read_delim("C:/Users/Andrew/Downloads/densityCalculations.log", 
#                                                      delim = "\t", escape_double = FALSE, 
#                                                      trim_ws = TRUE)
#main - This code will divide the non-PLC by the PLC for each peak type. 1st Will have to filter for the different peak types, then within iterate through and do the division.
adjustedEnrichment <- c()
adjustedEnrichmentBams <- data.frame(matrix(ncol=3,nrow=0))
colnames(adjustedEnrichmentBams) <- c("CellsBam", "PLCBam", "peak")
PLCs <- c()
PLCsBams <- c()
Cells <- c()
CellsBams <- c()
for(peak in 1:length(unique(enrichmentScores$bed))){
  for(bam in 1:length(enrichmentScores$bam[enrichmentScores$bed == unique(enrichmentScores$bed)[peak]])){
    if(str_detect(enrichmentScores$bam[enrichmentScores$bed == unique(enrichmentScores$bed)[peak]][bam], "PLC")){
      PLCs <- c(PLCs,enrichmentScores$Enrichment[enrichmentScores$bed == unique(enrichmentScores$bed)[peak]][bam])
      PLCsBams <- c(PLCsBams,enrichmentScores$bam[enrichmentScores$bed == unique(enrichmentScores$bed)[peak]][bam])
      #print(paste0("plcs ",bam," ", PLCsBams))
    }else{
      Cells <- c(Cells,enrichmentScores$Enrichment[enrichmentScores$bed == unique(enrichmentScores$bed)[peak]][bam])
      CellsBams <- c(CellsBams,enrichmentScores$bam[enrichmentScores$bed == unique(enrichmentScores$bed)[peak]][bam])
      #print(paste0("Cells ",bam," ", CellsBams))
    }
  }
}
if(length(PLCs) < length(Cells)){
  PLCs <- rep(PLCs, each = (length(unique(CellsBams))))
  PLCsBams <- rep(PLCsBams, each = (length(unique(CellsBams))))
}
adjustedEnrichment <- c(adjustedEnrichment, (Cells/PLCs))
adjustedEnrichmentBams <- data.frame(CellsBam=CellsBams, PLCBam=PLCsBams, peak = rep(unique(enrichmentScores$bed), each= length(unique(CellsBams))), adjustedEnrichment = adjustedEnrichment)
write.table(adjustedEnrichmentBams, file='adjustedEnrichment.tsv', quote=FALSE, sep='\t', row.names = FALSE)

print(adjustedEnrichmentBams)
adjustedEnrichmentBams$CellsBam <- str_remove(basename(as.character(adjustedEnrichmentBams$CellsBam)),"_001.trim.st.all.blft.qft.rmdup.bam")
adjustedEnrichmentBams$PLCBam <- str_remove(basename(as.character(adjustedEnrichmentBams$PLCBam)),"_001.trim.st.all.blft.qft.rmdup.bam")
adjustedEnrichmentBams$peak <- basename(as.character(adjustedEnrichmentBams$peak))
adjustedEnrichmentBams <- cbind(adjustedEnrichmentBams, paste0(adjustedEnrichmentBams$CellsBam, adjustedEnrichmentBams$PLCBam, adjustedEnrichmentBams$peak))
colnames(adjustedEnrichmentBams) <- c("CellsBam", "PLCBam", "peak", "adjustedEnrichment", "unique")
adjustedEnrichmentBams$unique <- factor(adjustedEnrichmentBams$unique, levels = adjustedEnrichmentBams$unique)

ggplot(adjustedEnrichmentBams)+
  geom_bar(aes(x=CellsBam, y=adjustedEnrichment, fill = peak), position="dodge", stat = "identity")+
  theme_classic()+
  facet_wrap("PLCBam")+
  ggtitle("Adjusted Enrichment of Break Density Within Peaks") +
  ylab("Adjusted Enrichment (Cells/PLC)")+
  xlab("Sample")+
  theme(axis.text.x = element_text(size = 10, vjust = 0.5,), 
        axis.text.y = element_text(size = 20, vjust = 0.5, hjust=1),
        axis.title=element_text(size=20,face="bold"),
        plot.title = element_text(size=22),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))+
  labs(fill="Peaks")+
  scale_y_continuous(breaks = round(seq(min(adjustedEnrichmentBams$adjustedEnrichment)-0.05, max(adjustedEnrichmentBams$adjustedEnrichment)+0.05, by = 0.02),2))+
  coord_cartesian(ylim = c(min(adjustedEnrichmentBams$adjustedEnrichment)-0.05, max(adjustedEnrichmentBams$adjustedEnrichment)+0.05))

ggsave("Adjusted_Enrichment_of_Break_Density_Within_Peaks_Plot.pdf", width = 14, height = 10, device = "pdf")
