# density plots for autosomal chromosomes

library(dplyr)
library(ggplot2)
library(here)
library(gridExtra)

multi_counts <- readRDS(here("ckd","multi_aggr_prep","step3_multi_counts.rds"))

# select  barcodes that passed qc
meta <- readRDS(here("ckd","multi_aggr_prep","step3_meta_multi_loy.rds"))

multi_counts <- multi_counts %>%
  dplyr::filter(barcode %in% rownames(meta))


plotA <- multi_counts %>%
  dplyr::filter(seqnames != "chrY") %>%
  ggplot(aes(scaled_atac_frags, alpha=0.5, fill=seqnames)) +
  geom_density(aes(y = after_stat(scaled)), adjust=2) +
  guides(alpha = guide_legend(element_blank())) +
  facet_wrap(~seqnames) +
  xlim(c(0,2)) +
  theme_bw() +
  theme(plot.title = element_text(face="bold"),
        axis.text = element_text(color="black"),
        legend.pos = "none") +
  guides(alpha = "none") +
  guides(fill = guide_legend(title = "Genotype")) +
  ylab("Density") +
  xlab("Scaled RNA Coverage") +
  ggtitle("A)")

plotB <- multi_counts %>%
  dplyr::filter(seqnames != "chrY") %>%
  ggplot(aes(scaled_atac_frags, alpha=0.5, fill=seqnames)) +
  geom_density(aes(y = after_stat(scaled)), adjust=2) +
  guides(alpha = guide_legend(element_blank())) +
  facet_wrap(~seqnames) +
  xlim(c(0,2)) +
  theme_bw() +
  theme(plot.title = element_text(face="bold"),
        axis.text = element_text(color="black"),
        legend.pos = "none") +
  guides(alpha = "none") +
  guides(fill = guide_legend(title = "Genotype")) +
  ylab("Density") +
  xlab("Scaled ATAC Coverage") +
  ggtitle("B)")
  
pdf(here(plots,"sfigure3.pdf"), width = 7.5, height = 10)
 grid.arrange(plotA, plotB, nrow=2)
dev.off()
  
  
  
