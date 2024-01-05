# density plots for m vs. f and different chromosomes

library(dplyr)
library(ggplot2)
library(here)
library(gridExtra)

multi_counts <- readRDS(here("ckd","multi_aggr_prep","step3_multi_counts.rds"))

# male vs female chrY RNA coverage
plotA <- multi_counts %>%
  group_by(seqnames) %>%
  mutate(scaled_rna_counts = scale(corrected_rna_counts, center=0)) %>%
  dplyr::filter(seqnames == "chrY") %>%
  dplyr::mutate(sex = ifelse(malesex == 1,"Male","Female")) %>%
  ggplot(aes(scaled_rna_counts, fill=sex)) +
  geom_density(aes(y = ..scaled..)) +
  facet_wrap(~sex) +
  theme_bw() +
  theme(plot.title = element_text(face="bold"),
        legend.pos = "none",
        axis.text = element_text(color="black")) +
  ggtitle("A)") +
  xlab("chrY RNA coverage") +
  ylab("Density") +
  coord_cartesian(xlim = c(0,5))
  
# male vs female chrY ATAC coverage
plotB <- multi_counts %>%
  group_by(seqnames) %>%
  mutate(scaled_atac_frags = scale(corrected_atac_frags, center=0)) %>%
  dplyr::filter(seqnames == "chrY") %>%
  dplyr::mutate(sex = ifelse(malesex == 1,"Male","Female")) %>%
  ggplot(aes(scaled_atac_frags, fill=sex)) +
  geom_density(aes(y = ..scaled..)) +
  facet_wrap(~sex) +
  theme_bw() +
  theme(plot.title = element_text(face="bold"),
        legend.pos = "none",
        axis.text = element_text(color="black")) +
  ggtitle("B)") +
  xlab("chrY ATAC coverage") +
  ylab("Density") +
  coord_cartesian(xlim = c(0,5))
  
# LOY per library by sequencing depth in PT
# pull meta data with number of counts per cell
meta <- readRDS(here("ckd","multi_aggr_prep","step3_meta_multi_loy.rds"))

multi_counts <- multi_counts %>%
  left_join(data.frame(barcode = rownames(meta), nCount_RNA = meta$nCount_RNA, nCount_ATAC = meta$nCount_ATAC), by = "barcode")

thresholds <- seq(1000, 10000, by = 1000)
loy_rna <- lapply(thresholds, function(threshold) {
  res <- multi_counts %>%
    dplyr::ungroup() %>%
    dplyr::filter(gmm_genotype %in% c("XY","LOY")) %>%
    dplyr::filter(celltype %in% c("PCT","PST")) %>%
    dplyr::filter(nCount_RNA > threshold) %>%
    summarize(genotype = gmm_genotype, library_id = library_id, threshold = threshold) %>%
    table() %>%
    as.data.frame() %>%
    group_by(library_id) %>%
    dplyr::mutate(total_cells = sum(Freq), prop_loy = Freq[genotype == "LOY"] / sum(Freq)) %>%
    distinct(library_id, threshold, total_cells, prop_loy)
  }) %>% bind_rows()

plotC <- loy_rna %>%
  as.data.frame() %>%
  ggplot(aes(threshold, prop_loy, group = library_id)) +
  geom_line(aes(color=library_id)) +
  theme_bw() +
  theme(plot.title = element_text(face="bold"),
        legend.title = element_blank(),
        axis.text = element_text(color="black")) +
  ggtitle("C)") +
  xlab("nCount_RNA Depth Threshold") +
  ylab("Proportion LOY") 

#########
# ATAC depth threshold
thresholds <- seq(10000, 20000, by = 1000)
loy_atac <- lapply(thresholds, function(threshold) {
  res <- multi_counts %>%
    dplyr::ungroup() %>%
    dplyr::filter(gmm_genotype %in% c("XY","LOY")) %>%
    dplyr::filter(celltype %in% c("PCT","PST")) %>%
    dplyr::filter(nCount_ATAC > threshold) %>%
    summarize(genotype = gmm_genotype, library_id = library_id, threshold = threshold) %>%
    table() %>%
    as.data.frame() %>%
    group_by(library_id) %>%
    dplyr::mutate(total_cells = sum(Freq), prop_loy = Freq[genotype == "LOY"] / sum(Freq)) %>%
    distinct(library_id, threshold, total_cells, prop_loy)
  }) %>% bind_rows()

plotD <- loy_atac %>%
  as.data.frame() %>%
  ggplot(aes(threshold, prop_loy, group = library_id)) +
  geom_line(aes(color=library_id)) +
  theme_bw() +
  theme(plot.title = element_text(face="bold"),
        legend.title = element_blank(),
        axis.text = element_text(color="black")) +
  ggtitle("D)") +
  xlab("nCount_ATAC Depth Threshold") +
  ylab("Proportion LOY") 

lay <- rbind(c(1,2),
             c(3,3),
             c(4,4))
gs <- list(plotA, plotB, plotC, plotD)

pdf(here(plots,"sfigure2.pdf"), width = 7.5, height = 10)
 grid.arrange(grobs = gs, layout_matrix = lay)
dev.off()


