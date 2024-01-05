# this script will analyze degs for pathway enrichment
# to run locally:
# SCRATCH1=/mnt/g/scratch
# docker run -it --rm \
# --workdir $HOME \
# -v /mnt/g/reference:$HOME/reference \
# -v /mnt/g/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# -v /mnt/g/cellranger_multi_counts:$HOME/cellranger_multi_counts \
# -v /mnt/g/ckd:$HOME/ckd \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/sctools:R4.1.3b R

# 
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(openxlsx)
library(RColorBrewer)
library(DOSE)
library(here)
library(dplyr)
  
file <- here("ckd","atac_aggr_prep","markers","dar.celltype.markers.xlsx")
sheets <- getSheetNames(file)
deg <- lapply(sheets, function(sheet) {
  tryCatch({df <- read.xlsx(file, sheet)
            df$celltype <- sheet
            colnames(df)[1] <- "gene"
            return(df)
            }, error=function(e) return(NULL))
}) %>% bind_rows()
  
sub <- deg %>%
  dplyr::filter(celltype == "PTPROM1", p_val < 0.05)

# we want the log2 fold change 
original_gene_list <- sub$avg_log2FC

# name the vector
names(original_gene_list) <- sub$gene_name

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# perform enrichment
keytypes(org.Hs.eg.db)
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             minGSSize = 5, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

require(DOSE)
pdf(here("ckd","atac_aggr_prep","plots","gsea_ptprom1.pdf"))
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()
  
saveRDS(gse, here("ckd","atac_aggr_prep","step8_gsea_ptprom1.rds"))  

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

pdf(here("ckd","atac_aggr_prep","plots","gsea_ptprom1.pdf"))
plot <- dotplot(gse, showCategory=terms$Description, split=".sign") + 
  facet_grid(.~.sign) +
  theme_bw() +
  theme(legend.pos = "none", plot.title = element_text(face="bold"), axis.text = element_text(color="black")) 
print(plot)  
dev.off()  
