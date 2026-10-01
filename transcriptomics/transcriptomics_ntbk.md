# Transcriptomics Notebook

**Course:** Intro to Ecological Genomics - Fall 2026

**Name:** Issy Timberlake

------------------------------------------------------------------------

## 9.15.2026- Setting up lab notebook and learning markdown

-   Setting up transcriptomics notebook

-   Learn how to take notes in markdown

-   Push notes to Github

**Working directory**

`/gpfs1/home/i/t/itimberl/eco_genomics_projects/eco_genomics_2026/transcriptomics`

**Input files**

none

**Output files**

`/gpfs1/home/i/t/itimberl/eco_genomics_projects/eco_genomics_2026/transcriptomics/transcriptomics_ntbk.md`

**Programs/Dependencies**

-   R version 4.5.1

-   R-studio

**Scripts:**

``` {.R .R}
print("Hello World!")
```

**Graphs/Images:**

| Col1 | Col2 | Col3 |
|------|------|------|
| x    | y    | z    |
| x    | y    | z    |
| x    | y    | z    |

: trying a table out

![](markdown-syntax-cheatsheet.webp)

**Notes/Observations:\
**cool graph; good to look back on later

next, I will use this for real data analysis!

------------------------------------------------------------------------

## 9.15.2026- Introducing the study system

**Working directory:**

`/gpfs1/home/i/t/itimberl/eco_genomics_projects/eco_genomics_2026/transcriptomics`

**Input files:**

**Output files:**

**Programs/Dependencies:**

-   R version 4.5.1

-   R-studio

**Notes/Observations:\
**

------------------------------------------------------------------------

## 9.17.2026- Introducing the study system and how the data was prepared through Illumina

**Working directory:**

`/gpfs1/home/i/t/itimberl/eco_genomics_projects/eco_genomics_2026/transcriptomics`

**Input files:**

none

**Output files:**

none

**Programs/Dependencies:**

-   R version 4.5.1

-   R-studio

-   Command Shell

**Notes/Observations:**

-   Questions we could ask/hypotheses...

    -   look for changes in GE in response to stressor

    -   change in gene expression through time/generations

    -   which treatment is most impactful

    -   which changes are most stable

    -   what overlaps occur across treatments

    -   Plasticity vs. Adaptive Evolution

    -   Phenotype Integration and correlation with genes

    -   Additive, synergistic, or antagonistic interactions between treatments?

-   Starting with uploading data through DESeq2

-   Looking at raw data Phred Q score, it is good to see letters (such as I)... shows high score... symbols show low scores

## 9.22.2026- Setting up R environment and loading data

**Working directory:**

`/gpfs1/home/i/t/itimberl/eco_genomics_projects/eco_genomics_2026/transcriptomics`

**Input files:**

counts data

**Output files:**

\~/eco_genomics_projects/eco_genomics_2026/transcriptomics/myscripts/ahud_DESeq_inclass.R

**Programs/Dependencies:**

-   R version 4.5.1

-   R-studio

-   Command Shell

-   DESeq2

**Notes/Observations:**

-   Adding startup to our R directories

    ```         
    setup.sh
    ```

-   Also, add the command "load ecogen-rlibs" under Module when starting an R session

-   Copied data from directory:

    ```         
    /gpfs1/cl/biol3990/Transcriptomics/CountsMatrix
    ```

-   Put data under new mydata folder

-   Created an R script to load data and view it as a matrix using DESeq2

-   Imported libraries and data

    -   needed to get rid of decimal points for DESeq2

    -   reads matrix file containing genes x samples and conditions file listing the treatment and generation of each sample

-   Explored patterns in the data

    -   Saw an average of 14687388 reads per sample

    -   see a variance in the number of count means per sample... for example. AAF11Rep3 has the lowest amount of reads while AAF2Rep3 has the highest... this variance due to sequencing effort and efficiency, show technical artifacts of sequencing depth and output; or could show low quality reads

    -   also see variance in counts per genes... average of 8217.81 reads per gene, but a median of 377 reads per gene

        -   this shows presence of some high outliers

        -   most genes show close to zero (to 1000) counts

-   Start working with DESeq2 and creating a DESeq data object

    -   

        ```{r}
        dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                                      design= ~ generation + treatment)
        ```

-   Filter data... want at least 15 reads in 75% of samples.

    -   Now have 25,260 transcripts

-   Run the DESeq model to test for differential expression using DESeq(dds)

    -   accounts for size differences, dispersions, gene estimates, mean-dispersion

    -   fits a model and applies statistical tests

-   Results of dds:

    ```         
    [1] "Intercept"            "generation_F11_vs_F0" "generation_F2_vs_F0"  "generation_F4_vs_F0" 
    [5] "treatment_OA_vs_AM"   "treatment_OW_vs_AM"   "treatment_OWA_vs_AM" 
    ```

-   Each of those results were tested for

-   Now, visualizing the data

    -   first, transform the data

        -   "normalize" through log transformation

        -   look at meanSD plot

        -   also variance stabilizing transformation

    -   make a heatmap of pairwise differences between samples to look for outliers

        -   see sample (38) x sample (38)

        -   see F2 is a bit of an outlier... looks very light compared to other samples

    -   make a sample tree

        -   look at relationships between samples

        -   again see AH_F2_Rep2 as an outlier... looking back at histogram of average count per sample, see that it has low counts

    -   make a PCA plot

        -   one general, but then split up by generation to make it easier to visualize

        -   Principal Component 1 always describes most of the variance... see most differences along x-axis

    -   ![](myresults/PCA_allGens.png){width="300"}

-   Note on RStudio: dev.off() resets the plot viewing panel

## 9.24.2026- Going over what we have learned

**Working directory**

`/gpfs1/home/i/t/itimberl/eco_genomics_projects/eco_genomics_2026/transcriptomics`

**Input files**

none

**Output files**

`/gpfs1/home/i/t/itimberl/eco_genomics_projects/eco_genomics_2026/transcriptomics/transcriptomics_ntbk.md`

**Programs/Dependencies**

-   R version 4.5.1

-   R-studio

**Notes/Observations**

1)  Where things are

2)  How to move around -\> using bash commands

-   PATHs: "\~" is home directory shortcut... "/users/i/t/itimberl", containing "/eco_genomics_projects/eco_genomics_2026/transcriptomics/mydata" (+ /myresults + myscripts)

-   also: "/gpfs1/cl/biol3990" is class directory, containing "transcriptomics/CountsMatrix.txt"

-   also: "gpfs1/cl/ecogen" containing "/sw/setup.sh"

3)  How to tell the computer what to do (bash, R, etc)

-   cp = copy

-   tab to complete

4)  How to back up and share work -\> Github

Also, played in R to understand basic R functions.

Learned to create sections in an R script

```{r}

# Playing in R ####
# 4 hashtags adds a section! appears at bottom of script (or cntrl+shift+R)

```

Then, looked at DESeq functions that we looked at on Tuesday.

```{r}
# Import the counts matrix
countsTable <- read.table("mydata/salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
# parameters of above read.table() function shows my file has a header, row names are the first column
head(countsTable)
dim(countsTable)
# 67916 genes, 38 samples
```

## 9.29.2026- Continuing analysis of differential expression with DESeq2

-   Loading counts matrix data, filtering out genes with sparse counts, and subsetting the data into different generations
-   Summarize results between different groups in the F0 generation
-   Plot this on scatterplots, MA plots, volcano plots, heatmap
-   Using vsd and heatmap, grouping together genes that are similar
-   Plot Euler plot to look at shared DEGs between different groups (just amounts, not specific genes)

**Working directory**

`/gpfs1/home/i/t/itimberl/eco_genomics_projects/eco_genomics_2026/transcriptomics`

**Input files**

counts data

**Output files**

`/gpfs1/home/i/t/itimberl/eco_genomics_projects/eco_genomics_2026/transcriptomics/transcriptomics_ntbk.md`

\~/eco_genomics_projects/eco_genomics_2026/transcriptomics/myscripts/9-26_10-1_ahud_DESeqpt2.R

**Programs/Dependencies**

-   R version 4.5.1

-   R-studio

-   DESeq2

**Notes/Observations**

-   dim() stands for dimensions

-   Look only at F0 generation:

    ```{r}
    # Subset the DESeqDataSet to the specific level of the "generation" factor
    dds_F0 <- subset(dds, select = generation == 'F0')
    dim(dds_F0)
    # [1] 25260    12
    ```

-   use ?DESeq() to look at the function, its parameters, etc

-   Adjusted p-value vs. regular p-value

    -   adjusts the p-value by weighting it with the number of times that you ran the test... since 5% error on sooooo many samples is many many samples still

-   When running summary on results, receive this output:

    ```         
    out of 25260 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 2343, 9.3%
    LFC < 0 (down)     : 1575, 6.2%
    outliers [1]       : 22, 0.087%
    low counts [2]     : 979, 3.9%
    (mean count < 23)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results
    ```

-   What is an MA plot?

    -   x-axis is mean of counts for each gene

    -   y-axis is log fold change... 0 is the ambient baseline, then each point showing how far off the OW group is

    -   blue genes are significant, gray are not

        ![](myresults/F0MAplot9-29.png){width="561"}

-   What is a volcano plot?

    ![](myresults/9-29volcanoplot_F0OWvsAM.png){width="519"}

    -   y-axis shows the significance factor

    -   x-axis shows the log fold change

    -   also shows down regulation and up regulation through color

    -   seeing lots more significant/great enough log fold changes in upregulated genes than down

-   How to look at and isolate differentially expressed genes in a comparison

    ```{r}
    # For OWA vs AM
    res_OWAvsAM <- results(dds_F0, name="treatment_OWA_vs_AM", alpha=0.05)
    res_OWAvsAM <- res_OWAvsAM[order(res_OWAvsAM$padj),]
    res_OWAvsAM <- res_OWAvsAM[!is.na(res_OWAvsAM$padj),]
    degs_OWAvsAM <- row.names(res_OWAvsAM[res_OWAvsAM$padj < 0.05,])
    ```

-   Heatmap

    ![](myresults/F0heatmap9-29.png){width="505"}

-   Euler plot looks at shared DEGs between groups

    ![](myresults/EulerplotF0.png){width="478"}

-   Upset plot

    -   Displays the same information as Euler plot but in a different context
    -   instead of size of Venn diagram bubble, each bar shows the height/weight of DEGs
    -   Especially useful if there is more than 3 comparisons

    ![](myresults/UpsetplotF0.png){width="497"}

Questions

-   How are all the samples all together in one matrix?

## 10. 1. 2026- Continuing DESeq2 Analysis pt. 3

-   Making scatterplot to compare the LFC response in OWA vs AM and OW vs AM
-   Analyzing the scatterplot
    -   some points could have no significance but large log2fold changes, and can be because there is a lot of variation among biological replicates already (in addition to variation between groups)
    -   See pattern where points significant in both groups follow y=x line, and if they are only significant in one dimension, they will vary across the opposite axis

**Working directory**

`/gpfs1/home/i/t/itimberl/eco_genomics_projects/eco_genomics_2026/transcriptomics`

**Input files**

counts data

**Output files**

`/gpfs1/home/i/t/itimberl/eco_genomics_projects/eco_genomics_2026/transcriptomics/transcriptomics_ntbk.md`

\~/eco_genomics_projects/eco_genomics_2026/transcriptomics/myscripts/9-26_10-1_ahud_DESeqpt2.R

**Programs/Dependencies**

-   R version 4.5.1

-   R-studio

-   DESeq2

**Notes/Observations**

-   What does results output mean?

    -   baseMean shows the mean number of counts

    -   log2FoldChange: shows positive/negative regulation

    -   lfcSE: how much variation in lfc among variants

    -   test statistic: how significant? also +/-

    -   p-value: how significant is that one gene

    -   p-adj: adjusts p-val for correction factor of how big the dataset is

-   Making the scatterplot

    -   First, changing the data

    -   make a new data frame with gene name (from rownames), padj, and LFC for each of the two comparison groups (OWA vs AM and OW vs AM) and then merge them

    -   filter out rows with missing LFC values, mutating a new column for significance in either group or both (use case_when() as a sort of if_else() function to see if padj is \< .05 in either group or both)

    -   Also calculating the correlation between the two groups' LFCs

    -   Sort the data so that it ranges from neither group on bottom to both significance as the top

        ![](myresults/ScatterplotofOWAvsOWvsAM10-1.png){width="455"}

-   annotate adds whatever text

-   important to add the points with significance on the bottom/last, so that these are shown on the top layer of the graph

## 10. .2026-

-   template

**Working directory**

`/gpfs1/home/i/t/itimberl/eco_genomics_projects/eco_genomics_2026/transcriptomics`

**Input files**

none

**Output files**

`/gpfs1/home/i/t/itimberl/eco_genomics_projects/eco_genomics_2026/transcriptomics/transcriptomics_ntbk.md`

**Programs/Dependencies**

-   R version 4.5.1

-   R-studio

**Notes/Observations**
