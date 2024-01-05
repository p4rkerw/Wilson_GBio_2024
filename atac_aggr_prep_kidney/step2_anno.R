#!/usr/bin/env Rscript
# this script will preprocess aggregated snATACseq data
# counted and aggregated by cellranger-atac v2.1 without library normalization
#
# to run locally
# counts=/mnt/g/cellranger_atac_counts
# project=/mnt/g/ckd
# SCRATCH1=/mnt/g/scratch
# docker run -it \
# --workdir $HOME \
# -v $project:$HOME/ckd \
# -v $counts:$HOME/cellranger_atac_counts \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/sctools:R4.1.3 R

library(Seurat)
library(Signac)
library(here)
library(ggplot2)
library(dplyr)
library(stringr)
library(data.table)

atac_aggr_prep <- here("ckd","atac_aggr_prep")
atacAggr <- readRDS(here(atac_aggr_prep,"step1f_prep.rds"))
DefaultAssay(atacAggr) <- "peaks"

# finalize the annotations
atacAggr <- RenameIdents(atacAggr,
'0' = 'PCT',
'1' = 'PCT',
'2' = 'DCT1',
'3' = 'PST',
'4' = 'PCT',
'5' = 'TAL2',
'6' = 'ENDO',
'7' = 'PT_VCAM1',
'8' = 'DCT2',
'9' = 'TAL2',
'10' = 'TAL1',
'11' = 'PC',
'12' = 'ICA',
'13' = 'FIB_VSMC_MC',
'14' = 'ATL',
'15' = 'PCT',
'16' = 'PCT',
'17' = 'ICB',
'18' = 'PT_PROM1',
'19' = 'ENDO',
'20' = 'PEC',
'21' = 'PCT',
'22' = 'PCT',
'23' = 'PCT',
'24' = 'TCELL',
'25' = 'PODO',
'26' = 'MONO',
'27' = 'MD',
'28' = 'BCELL'
)

# update levels
levels(atacAggr) <- c("PCT","PST","PT_VCAM1","PT_PROM1","PEC",
    "ATL","TAL1","TAL2","MD","DCT1","DCT2","PC",
    "ICA","ICB","PODO","ENDO","FIB_VSMC_MC",
    "BCELL","TCELL","MONO")
atacAggr@meta.data$celltype <- Idents(atacAggr)
      
# visualize annotations
pdf(here(atac_aggr_prep,"plots","step1l_final.pdf"))
p1 <- DimPlot(atacAggr, label=TRUE, group.by="celltype") + ggtitle("snATAC celltypes")
AugmentPlot(p1)
dev.off()

# save meta data
saveRDS(atacAggr@meta.data, file=here(atac_aggr_prep, "meta.rds"))
 
# save peaks
peaks.df <- as.data.frame((atacAggr[["peaks"]]@ranges))
fwrite(peaks.df, file=here(atac_aggr_prep, "peaks.gr"), sep="\t")

# create a list of cell type barcodes and write to file
dir.create(here("ckd","barcodes"), showWarnings=FALSE)
bc <- atacAggr@meta.data %>%
  dplyr::select(barcode, celltype, library_id) %>%
  dplyr::mutate(barcode = str_split(barcode, pattern="-", simplify=TRUE)[,1]) 
write.table(bc, file=here("ckd","barcodes","atac_barcodes.csv"), sep=",", row.names=FALSE, quote=FALSE)

# save object
saveRDS(atacAggr, here(atac_aggr_prep,"step2_anno.rds"), compress=FALSE)
