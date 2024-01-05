# this script will analyze degs for pathway enrichment

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(openxlsx)
library(RColorBrewer)
library(DOSE)
library(here)
library(dplyr)
library(gridExtra)

# create a plot dir
plots <- here("ckd","figures")
dir.create(here(plots))

set.seed(1234)
celltype_sel <- "PT_VCAM1"
gse <- readRDS(here("ckd","rna_aggr_prep","kpmp",paste0("step2_",celltype_sel,"_gsea.rds")))  

terms <- gse@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::mutate(log_pval = -log10(p.adjust)) %>%
  dplyr::mutate(sign = ifelse(NES > 0, 1, 0)) %>%
  group_by(sign) %>%
  top_n(10, log_pval)

# there are identical gene sets that have multiple descriptions
n_distinct(gse@result$core_enrichment)
nrow(gse@result)

# extract first description for each unique gene set
res <- gse@result %>%
  group_by(core_enrichment) %>%
  arrange(core_enrichment) %>%
  filter(row_number ()==1) %>%
  ungroup()
gse@result <- res

plotA <- dotplot(gse, showCategory=terms$Description, split=".sign") + 
  facet_grid(.~.sign) +
  theme_bw() +
  theme(legend.pos = "none", plot.title = element_text(face="bold"), axis.text = element_text(color="black")) +
  ggtitle(paste0("A) ", celltype_sel))


celltype_sel <- "PT_PROM1"
gse <- readRDS(here("ckd","rna_aggr_prep","kpmp",paste0("step2_",celltype_sel,"_gsea.rds")))  

terms <- gse@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::mutate(log_pval = -log10(p.adjust)) %>%
  dplyr::mutate(sign = ifelse(NES > 0, 1, 0)) %>%
  group_by(sign) %>%
  top_n(10, log_pval)

# there are identical gene sets that have multiple descriptions
n_distinct(gse@result$core_enrichment)
nrow(gse@result)

# extract first description for each unique gene set
res <- gse@result %>%
  group_by(core_enrichment) %>%
  arrange(core_enrichment) %>%
  filter(row_number ()==1) %>%
  ungroup()
gse@result <- res

plotB <- dotplot(gse, showCategory=terms$Description, split=".sign") + 
  facet_grid(.~.sign) +
  theme_bw() +
  theme(legend.pos = "none", plot.title = element_text(face="bold"), axis.text = element_text(color="black")) +
  ggtitle(paste0("B) ", celltype_sel))

lay <- rbind(c(1,1),
             c(2,2))
gs <- list(plotA, plotB)

pdf(here(plots,"sfigure9.pdf"), width = 7.5, height = 10)
 grid.arrange(grobs = gs, layout_matrix = lay)
dev.off()
