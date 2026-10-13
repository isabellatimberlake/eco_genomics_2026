'''
Isabella Timberlake
Ecological Genomics Biol3990
October 2026
Transcriptomics HW 1
'''
## loading everything in ####

# set working directory
setwd("~/eco_genomics_projects/eco_genomics_2026/transcriptomics")

# load packages
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(vsn)  
library("pheatmap")
library("vsn")

# import and convert counts data
countsTable <- read.table("mydata/salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
countsTableRound <- round(countsTable)

# import the sample description table
conds <- read.delim("mydata/ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)

# Look at reads per sample
colSums(countsTableRound)
mean(colSums(countsTableRound))
# 14687388 reads

barplot(colSums(countsTableRound), names.arg=colnames(countsTableRound),cex.names=0.5, las=3,ylim=c(0,21000000))
abline(h=mean(colSums(countsTableRound)), col="blue", lwd=2)

#conds |> filter(treatment %in% c("OWA", "AM"))
#countsTableRoundOWA <- countsTable |> select()


dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ treatment)
# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28)
## suggested by WGCNA on RNAseq FAQ
dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds)

# subset by generation, then only look at the specific results between OWA and AM
dds_F0 <- subset(dds, select = generation == 'F0')



                  
  