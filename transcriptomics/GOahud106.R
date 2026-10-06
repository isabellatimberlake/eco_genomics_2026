## Issy Timberlake

# Setting up my data ####
# set working directory
setwd("~/eco_genomics_projects/eco_genomics_2026/transcriptomics")

## Import the libraries that we're likely to need in this session

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(vsn)  

# importing data
countsTable <- read.table("mydata/salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
conds <- read.delim("mydata/ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ treatment)
dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
dds_F0 <- subset(dds, select = generation == 'F0')
dim(dds_F0)
dds_F0 <- DESeq(dds_F0)

## OCEAN WARMING ACIDIFICATION vs AMBIENT
res_OWAvsAM <- results(dds_F0, name="treatment_OWA_vs_AM", alpha=0.05)
# taking the output and ordering/sorting it by padj value
res_OWAvsAM <- res_OWAvsAM[order(res_OWAvsAM$padj),]
head(res_OWAvsAM) 
summary(res_OWAvsAM)
 # making data frame
res_OWAvsAM.df <- as.data.frame(res_OWAvsAM)
res_OWAvsAM.df$fullID <- rownames(res_OWAvsAM.df)

parts <- strsplit(res_OWAvsAM.df$fullID, "::")

res_OWAvsAM.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OWAvsAM.df,
  "myresults/F0_OWAvsAM_results.csv",
  row.names = FALSE
)

## OCEAN WARMING vs AMBIENT
res_OWvsAM <- results(dds_F0, name="treatment_OW_vs_AM", alpha=0.05)
res_OWvsAM <- res_OWvsAM[order(res_OWvsAM$padj),]
head(res_OWvsAM) 
summary(res_OWvsAM)

res_OWvsAM.df <- as.data.frame(res_OWvsAM)
res_OWvsAM.df$fullID <- rownames(res_OWvsAM.df)

parts <- strsplit(res_OWAvsAM.df$fullID, "::")

res_OWvsAM.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OWvsAM.df,
  "myresults/F0_OWvsAM_results.csv",
  row.names = FALSE
)

## OCEAN ACIDIFICATION vs AMBIENT
res_OAvsAM <- results(dds_F0, name = "treatment_OA_vs_AM", alpha = 0.05)
res_OAvsAM <- res_OAvsAM[order(res_OAvsAM$padj)]
head(res_OAvsAM)
summary(res_OAvsAM)

res_OAvsAM.df <- as.data.frame(res_OAvsAM)
res_OAvsAM.df$fullID <- rownames(res_OAvsAM.df)

parts <- strsplit(res_OAvsAM.df$fullID, "::")

res_OAvsAM.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OAvsAM.df,
  "myresults/F0_OAvsAM_results.csv",
  row.names = FALSE
)

# GO ANALYSIS ####
library(topGO)

mappingFile <- "mydata/trinotate_annotation_GOblastx_forTopGO.txt"

geneID2GO <- readMappings(
  file = mappingFile
)

cat("Genes in GO mapping:",
    length(geneID2GO),
    "\n")
# 86453

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

GO_OA <- run_topGO_contrast(
  "myresults/F0_OAvsAM_results.csv",
  "myresults/GO_F0_OAvsAM_BP.csv"
)

GO_OW <- run_topGO_contrast(
  "myresults/F0_OWvsAM_results.csv",
  "myresults/GO_F0_OWvsAM_BP.csv"
)

GO_OWA <- run_topGO_contrast(
  "myresults/F0_OWAvsAM_results.csv",
  "myresults/GO_F0_OWAvsAM_BP.csv"
)

## Visualizing OWA vs AM results

# Read TopGO results
go <- read.csv(
  "myresults/GO_F0_OWAvsAM_BP.csv",
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
    title = "Top GO Terms: OWA vs AM",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Significant\nGenes"
  )

# Read TopGO results
go2 <- read.csv(
  "myresults/GO_F0_OAvsAM_BP.csv",
  stringsAsFactors = FALSE
)

# Convert p-values to numeric
go2$weightFisher <- gsub("^<\\s*", "", go2$weightFisher)

go2$weightFisher <- as.numeric(go2$weightFisher)

# Convert counts to numeric
go2$Significant <- as.numeric(go2$Significant)

# Create -log10(p)
go2$minusLogP <- -log10(go2$weightFisher)

# Keep top 10 GO terms
go2_top10 <- go2 %>%
  arrange(weightFisher) %>%
  slice(1:10)

# Order terms for plotting
go2_top10$Term <- factor(
  go2_top10$Term,
  levels = rev(go2_top10$Term)
)

# Bubble plot
ggplot(
  go2_top10,
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
    title = "Top GO Terms: OA vs AM",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Significant\nGenes"
  )

