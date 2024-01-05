#!/usr/bin/env Rscript
# to run locally:
# SCRATCH1=/mnt/g/scratch
# docker run -it \
# --workdir $HOME \
# -v /mnt/g/ckd:$HOME/ckd \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/sctools:R4.1.0 R

library(Signac) 
library(Seurat) 
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork) 
library(motifmatchr) 
library(here) 
library(chromVAR) 
library(future) 
library(openxlsx) 
set.seed(1234)

multi_aggr_prep <- here("ckd","multi_aggr_prep")
multi <- readRDS(here(multi_aggr_prep,"step3_multi_loy.rds"))
DefaultAssay(multi) <- "ATAC"

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# subset atac object by filtering out alt contigs
# this will result in warnings that can be ignored
# but is required to be compatible with the addmotifs function
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- which(as.character(seqnames(granges(multi))) %in% main.chroms)
multi[["ATAC"]] <- subset(multi[["ATAC"]], features = rownames(multi[["ATAC"]])[keep.peaks])

# add motif information
multi <- AddMotifs(
  object = multi,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

multi <- RegionStats(
  object = multi,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c("-", "-")
)

# compute motif activities using chromvar
multi <- RunChromVAR(
  object = multi,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

saveRDS(multi, here(multi_aggr_prep,"step7_chromVAR.rds"), compress=FALSE)
