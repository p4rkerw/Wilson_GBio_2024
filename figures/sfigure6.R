# density plots for thresholds

library(dplyr)
library(ggplot2)
library(here)
library(gridExtra)

atac_frags <- readRDS(here("ckd","atac_aggr_prep","atac_frags.rds"))
  
# LOY per library by sequencing depth in PT
# pull meta data with number of counts per cell
meta <- readRDS(here("ckd","atac_aggr_prep","step6_meta.rds"))

atac_frags <- atac_frags %>%
  left_join(data.frame(barcode = rownames(meta), nCount_ATAC = meta$nCount_peaks, genotype = meta$akd_genotype), by = "barcode")

#########
# ATAC depth threshold
thresholds <- seq(10000, 20000, by = 1000)
loy_atac <- lapply(thresholds, function(threshold) {
  res <- atac_frags %>%
    dplyr::ungroup() %>%
    dplyr::filter(seqnames == "chrY") %>%
    dplyr::filter(genotype %in% c("XY","LOY")) %>%
    dplyr::filter(celltype %in% c("PCT","PST")) %>%
    dplyr::filter(nCount_ATAC > threshold) %>%
    summarize(genotype = genotype, library_id = library_id, threshold = threshold) %>%
    table() %>%
    as.data.frame() %>%
    group_by(library_id) %>%
    dplyr::mutate(total_cells = sum(Freq), prop_loy = Freq[genotype == "LOY"] / sum(Freq)) %>%
    distinct(library_id, threshold, total_cells, prop_loy)
  }) %>% bind_rows()

plotA <- loy_atac %>%
  dplyr::filter(total_cells > 1000) %>%
  as.data.frame() %>%
  ggplot(aes(threshold, prop_loy, group = library_id)) +
  geom_line(aes(color=library_id)) +
  theme_bw() +
  theme(plot.title = element_text(face="bold"),
        legend.title = element_blank(),
        axis.text = element_text(color="black")) +
  ggtitle("A)") +
  xlab("nCount_ATAC Depth Threshold") +
  ylab("Proportion LOY") 

lay <- rbind(c(1,1))
gs <- list(plotA)

pdf(here(plots,"sfigure6.pdf"), width = 7.5, height = 5)
 grid.arrange(grobs = gs, layout_matrix = lay)
dev.off()
