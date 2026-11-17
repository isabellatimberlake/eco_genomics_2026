'''
Pop Gen HW
looking at genetic diversity between populations ASU and ROS
'''

library(tidyverse)

setwd("/users/i/t/itimberl/eco_genomics_projects/eco_genomics_2026/population_genomics/mydata/GeneticDiversity")
OUTDIR="/users/i/t/itimberl/eco_genomics_projects/eco_genomics_2026/population_genomics/results/PopGenHW"



# read in files
meta <- read_csv("sample_metadata.csv")

sites <- meta %>%
  group_by(Population) %>%
  summarise(Latitude = mean(Latitude))

pi_files <- list.files(path = "HWpi", pattern = "windowed.pi$", full.names = TRUE)
pi <- pi_files %>%
  # name each file by its site, e.g. "pi/BOD.windowed.pi" -> "BOD"
  set_names(str_remove(basename(pi_files), fixed(".windowed.pi"))) %>%
  # read every file and stack them, putting the names in a new column, Site
  map_dfr(read_tsv, .id = "Population", show_col_types = FALSE) %>%
  left_join(sites, by = "Population") %>%
  mutate(Site = fct_reorder(Population, Latitude))
pi_summary <- pi %>%
  group_by(Population, Latitude) %>%
  summarise(n_windows = n(),
            mean_pi   = mean(PI),
            median_pi = median(PI),
            .groups = "drop")
pi_summary

pi_plot <- ggplot(pi, aes(Site, PI, fill = Site)) +
  geom_violin(scale = "width", linewidth = 0.5, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.3) +
  scale_fill_viridis_d(end = 0.9) +
  scale_y_log10() +
  labs(x = "Population", y = "log10(Nucleotide diversity)") +
  theme_classic() +
  theme(legend.position = "none") + 
  coord_flip()
pi_plot

ggsave(paste0(OUTDIR, "pi_by_pop.png"), plot = pi_plot, width = 12, height = 6)

chosenSites <- pi |> filter(Population %in% c("ASU", "ROS"))

#top_scaffold <- chosenSites %>%
 # count(CHROM, sort = TRUE) %>%
  #slice(1) %>%
  #pull(CHROM)

scan_plot <- chosenSites %>%
  #filter(CHROM == top_scaffold) %>%
  ggplot(aes(BIN_START / 1e6, PI, colour = Site)) +
  geom_line(alpha=0.9) +
  #scale_colour_viridis_d(end = 0.9) +
  labs(x = paste( "position (Mb)"), y = "Nucleotide diversity",
       colour = "Site") +
  theme_bw()+
  facet_wrap(~Population)
scan_plot

ggsave(paste0(OUTDIR, "Chrom1scan_plot_wrapped.png"), plot = scan_plot, width = 12, height = 6)

pop1 <- "ASU"   
pop2 <- "ROS"




add_genome_pos <- function(df) {
  scaffolds <- df %>% 
    group_by(CHROM) %>% ##group by CHROM
    summarise(length = max(POS)) %>% ## length = max position within a chromosome
    arrange(match(CHROM, unique(df$CHROM))) %>% ##make sure each is unique, arranged
    mutate(offset = cumsum(as.numeric(length)) - length, ###cumulatively number (helps with labels)
           shade  = factor(row_number() %% 2))    ##value for alternating colour 
  df %>%
    left_join(scaffolds, by = "CHROM") %>% ##join this data by given data with CHROM
    mutate(genome_pos = offset + POS) ###cumalative genome pos 
}

fst_win <- read_tsv("/users/i/t/itimberl/eco_genomics_projects/eco_genomics_2026/population_genomics/mydata/fstHW/ASU_vs_ROS_FST_100kb.windowed.weir.fst")
fst_win <- fst_win %>%
  mutate(POS = (BIN_START + BIN_END) / 2) %>% 
  add_genome_pos() ## this is our new function

# Plotting the Fst outliers
# empirical outlier threshold: the top 1% of windows
fst_cut <- quantile(fst_win$WEIGHTED_FST, probs = 0.99)
fst_cut

fst_win <- mutate(fst_win, outlier = WEIGHTED_FST >= fst_cut)

ggplot(fst_win, aes(x = genome_pos / 1e6, y = WEIGHTED_FST)) +
  geom_point(aes(colour = shade), size = 0.8) +
  geom_point(data = filter(fst_win, outlier), colour = "red", size = 1.2) +
  geom_hline(yintercept = fst_cut, linetype = "dashed") +
  scale_colour_manual(values = c("grey35", "grey70"), guide = "none") +
  labs(x = "Position along the genome (Mb; scaffolds end to end)",
       y = "Weighted FST (10-kb windows)",
       title = paste(pop1, "vs", pop2)) +
  theme_classic()

ggsave(paste0(OUTDIR,"/fst_genome scan_", pop1, "_vs_", pop2, ".png"), width = 10, height = 4)

# looking at XtX

snps <- read_tsv("Lottia.baypass.snps.txt")

snps$MRK<-rownames(snps)
snps$MRK<-as.numeric(snps$MRK)

npop <- 19

xtx <- read.table("Lottia_bp_2_summary_pi_xtx.out", header = TRUE, check.names = FALSE) %>%
  left_join(snps, by = "MRK") %>%
  mutate(lp_high = pchisq(XtXst, npop, lower.tail = FALSE, log.p = TRUE), #
         lp_low = pchisq(XtXst, npop, lower.tail = TRUE, log.p = TRUE),
         direction = if_else(lp_high < lp_low, "high", "low"),
         log10_p = -(pmin(lp_high, lp_low) + log(2)) / log(10),
         pval = 10^(-log10_p),
         qval = p.adjust(pval, method = "BH"))
xtx <- xtx %>%
  add_genome_pos()

sum(xtx$outlier)

threshold <- quantile(xtx$qval, 0.001) # top 0.1% (smallest q)

ggplot(xtx, aes(genome_pos / 1e6, -log10(qval))) +
  geom_point(aes(colour = shade), size = 0.4) +
  geom_hline(yintercept = -log10(threshold), linetype = "dashed", colour = "red") +
  scale_colour_manual(values = c("grey20", "grey80"), guide = "none") +
  labs(x = "Genome pos (Mb)", y = "-log10(q)") +
  theme_classic()
ggsave(paste0(OUTDIR,"/HWbaypass_xtx_scan.png"), width = 10, height = 5)


# PI across genome
pi2 <- chosenSites |> mutate(POS = (BIN_START + BIN_END) / 2) %>% 
  add_genome_pos() 

# taking the top 1% of pi
pi_cut <- quantile(pi2$PI, probs = 0.99)
pi2 <- mutate(pi2, outlier = PI >= pi_cut)

pi_across_genome <- ggplot(pi2, aes(genome_pos / 1e6, PI)) +
      geom_point(aes(colour = shade), size = 0.4) +
      geom_hline(yintercept = pi_cut, linetype = "dashed", colour = "red") +
      scale_colour_manual(values = c("grey20", "grey80"), guide = "none") +
      labs(x = "Genome pos (Mb)", y = "PI") +
      theme_classic()+
    facet_wrap(~Population)
pi_across_genome 
ggsave(paste0(OUTDIR,"pi_across_manhatplot.png"), width = 10, height = 5)
