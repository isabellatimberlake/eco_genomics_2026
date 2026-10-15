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

# visualizing the samples
barplot(colSums(countsTableRound), names.arg=colnames(countsTableRound),cex.names=0.5, las=3,ylim=c(0,21000000))
abline(h=mean(colSums(countsTableRound)), col="blue", lwd=2)


dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ treatment)
# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28)
## suggested by WGCNA on RNAseq FAQ
dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds)
# 25260

# Generations####
# F4
# subset by generation, then only look at the specific results between OWA and AM
dds_F0 <- subset(dds, select = generation == 'F0')
dim(dds_F0)
# 25260 x 12

# find variable genes
dds_F0 <- DESeq(dds_F0)

## OCEAN WARMING ACIDIFICATION vs AMBIENT
res_OWAvsAM_F0 <- results(dds_F0, name="treatment_OWA_vs_AM", alpha=0.05)
# taking the output and ordering/sorting it by padj value
res_OWAvsAM_F0 <- res_OWAvsAM_F0[order(res_OWAvsAM_F0$padj),]
summary(res_OWAvsAM_F0)



# F4
dds_F4 <- subset(dds, select = generation == 'F4')
dim(dds_F4)
# find variable genes
dds_F4 <- DESeq(dds_F4)
## OCEAN WARMING ACIDIFICATION vs AMBIENT
res_OWAvsAM_F4 <- results(dds_F4, name="treatment_OWA_vs_AM", alpha=0.05)
# taking the output and ordering/sorting it by padj value
res_OWAvsAM_F4 <- res_OWAvsAM_F4[order(res_OWAvsAM_F4$padj),]
summary(res_OWAvsAM_F4)

# F11
dds_F11 <- subset(dds, select = treatment %in% c("AM", "OWA") & generation == 'F11')
dds_F11$treatment <- droplevels(dds_F11$treatment)
dds_F11$treatment <- relevel(dds_F11$treatment, ref = "AM")
table(dds_F11$treatment)

dim(dds_F11)
# 25260 6

# need to also select treatment

# find variable genes
dds_F11 <- DESeq(dds_F11)
## OCEAN WARMING ACIDIFICATION vs AMBIENT
res_OWAvsAM_F11 <- results(dds_F11, name="treatment_OWA_vs_AM", alpha=0.05)
# taking the output and ordering/sorting it by padj value
res_OWAvsAM_F11 <- res_OWAvsAM_F11[order(res_OWAvsAM_F11$padj),]
summary(res_OWAvsAM_F11)


# pulling out DEGs
# F0
res_OWAvsAM_F0 <- res_OWAvsAM_F0[order(res_OWAvsAM_F0$padj),]
res_OWAvsAM_F0 <- res_OWAvsAM_F0[!is.na(res_OWAvsAM_F0$padj),]
degs_OWAvsAM_F0 <- row.names(res_OWAvsAM_F0[res_OWAvsAM_F0$padj < 0.05,])
# F4
res_OWAvsAM_F4 <- res_OWAvsAM_F4[order(res_OWAvsAM_F4$padj),]
res_OWAvsAM_F4 <- res_OWAvsAM_F4[!is.na(res_OWAvsAM_F4$padj),]
degs_OWAvsAM_F4 <- row.names(res_OWAvsAM_F4[res_OWAvsAM_F4$padj < 0.05,])
# F11
res_OWAvsAM_F11 <- res_OWAvsAM_F11[order(res_OWAvsAM_F11$padj),]
res_OWAvsAM_F11 <- res_OWAvsAM_F11[!is.na(res_OWAvsAM_F11$padj),]
degs_OWAvsAM_F11 <- row.names(res_OWAvsAM_F11[res_OWAvsAM_F11$padj < 0.05,])

# making a euler plot to compare DEGs across generations
library(eulerr)

# total DEGs
length(degs_OWAvsAM_F0) # 3918
length(degs_OWAvsAM_F4) # 241
length(degs_OWAvsAM_F11) # 1403

# Intersections
length(intersect(degs_OWAvsAM_F0,degs_OWAvsAM_F4))  # 57
length(intersect(degs_OWAvsAM_F0,degs_OWAvsAM_F11))  # 226
length(intersect(degs_OWAvsAM_F4,degs_OWAvsAM_F11)) # 15    

# Identifying DEGs shared across all groups
ints04 <- intersect(degs_OWAvsAM_F0,degs_OWAvsAM_F4)
length(intersect(degs_OWAvsAM_F11,ints04)) # 2

# number unique to each generation
3918-57-226+2 # 3637 F0
241-57-15+2 # 171 F4
1403-226-15+2 #1164 F11
# Number shared in pairs of treatments
57-2 # 55 FO and F4
226-2 # 224 F0 and F11
15-2 # 13 F4 and F11

# making euler plot
fit1 <- euler(c("F0" = 3637, "F4" = 171, "F11" = 1164, "F0&F4" = 55, "F0&F11" = 224, "F4&F11" = 13, "F0&F4&F11" = 2))
plot(fit1,  lty = 1:3, quantities = TRUE)

# GO analysis setup ####

# F0
# making data frame
F0.df <- as.data.frame(res_OWAvsAM_F0)
F0.df$fullID <- rownames(F0.df)

parts0 <- strsplit(F0.df$fullID, "::")

F0.df$shortID <- sapply(
  parts0,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  F0.df,
  "myresults/OWAvsAM_resultsF0HW.csv",
  row.names = FALSE
)


# F4
# making data frame
F4.df <- as.data.frame(res_OWAvsAM_F4)
F4.df$fullID <- rownames(F4.df)

parts4 <- strsplit(F4.df$fullID, "::")

F4.df$shortID <- sapply(
  parts4,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  F4.df,
  "myresults/OWAvsAM_resultsF4HW.csv",
  row.names = FALSE
)

# F11
# making data frame
F11.df <- as.data.frame(res_OWAvsAM_F11)
F11.df$fullID <- rownames(F11.df)

parts11 <- strsplit(F11.df$fullID, "::")

F11.df$shortID <- sapply(
  parts11,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  F11.df,
  "myresults/OWAvsAM_resultsF11HW.csv",
  row.names = FALSE
)

# GO ANALYSIS ####
library(topGO)
mappingFile <- "mydata/trinotate_annotation_GOblastx_forTopGO.txt"
geneID2GO <- readMappings(
  file = mappingFile
)
# filtering GO terms by size
allGO <- table(unlist(geneID2GO))
keepTerms <- names(allGO)[
  allGO >= 5 &
    allGO <= 500
]

geneID2GO.filtered <- lapply(
  geneID2GO,
  function(x) intersect(x, keepTerms)
)

geneID2GO.filtered <- geneID2GO.filtered[
  lengths(geneID2GO.filtered) > 0
]

## FUNCTION for running Top Go analysis ####
run_topGO_contrast <- function(
    infile,
    outfile,
    ontology = "BP",
    padj.cutoff = 0.05){
  
  deseq <- read.csv(
    infile,
    stringsAsFactors = FALSE
  )
  
  deseq <- subset(
    deseq,
    !is.na(padj)
  )
  
  geneList <- factor(
    as.integer(deseq$padj < padj.cutoff)
  )
  
  names(geneList) <- deseq$shortID
  
  cat("\nGenes tested:",
      length(geneList))
  
  cat("\nSignificant genes:",
      sum(geneList == 1),
      "\n")
  
  GOdata <- new(
    "topGOdata",
    ontology = ontology,
    allGenes = geneList,
    geneSelectionFun = function(x) x == 1,
    annot = annFUN.gene2GO,
    gene2GO = geneID2GO.filtered
  )
  
  resultWeight <- runTest(
    GOdata,
    algorithm = "weight01",
    statistic = "fisher"
  )
  
  GOresults <- GenTable(
    GOdata,
    weightFisher = resultWeight,
    orderBy = "weightFisher",
    topNodes = 100
  )
  write.csv(
    GOresults,
    outfile,
    row.names = FALSE
  )
  
  return(GOresults)
}

# running top go
GO_F0 <- run_topGO_contrast(
  "myresults/OWAvsAM_resultsF0HW.csv",
  "myresults/GO_F0.csv"
)
GO_F4 <- run_topGO_contrast(
  "myresults/OWAvsAM_resultsF4HW.csv",
  "myresults/GO_F4.csv"
)
GO_F11 <- run_topGO_contrast(
  "myresults/OWAvsAM_resultsF11HW.csv",
  "myresults/GO_F11.csv"
)

# F11, F4, F0
# Read TopGO results ###
go <- read.csv(
  "myresults/GO_F0.csv",
  stringsAsFactors = FALSE
)

# Convert p-values to numeric
go$weightFisher <- gsub("^<\\s*", "", go$weightFisher)

go$weightFisher <- as.numeric(go$weightFisher)

# Convert counts to numeric
go$Significant <- as.numeric(go$Significant)

# Create -log10(p)
go$minusLogP <- -log10(go$weightFisher)

# Keep top 10 GO terms
go_top10 <- go %>%
  arrange(weightFisher) %>%
  slice(1:10)

# Order terms for plotting
go_top10$Term <- factor(
  go_top10$Term,
  levels = rev(go_top10$Term)
)

# Bubble plot
ggplot(
  go_top10,
  aes(
    x = minusLogP,
    y = Term
  )
) +
  geom_point(
    aes(
      size = Significant,
      color = minusLogP
    )
  ) +
  scale_color_viridis_c() +
  theme_bw(base_size = 14) +
  labs(
    title = "Top GO Terms F0: OWA vs AM",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Significant\nGenes"
  )

