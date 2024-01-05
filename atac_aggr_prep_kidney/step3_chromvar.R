#!/usr/bin/env Rscript
# to run locally:
# SCRATCH1=/mnt/g/scratch
# docker run -it \
# --workdir $HOME \
# -v /mnt/g/ckd:$HOME/ckd \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/sctools:R4.1.0
#
# to run interactively on the RIS compute1 cluster:
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/ckd:$HOME/ckd \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=256GB]' -q general-interactive -a 'docker(p4rkerw/sctools:R4.1.0)' /bin/bash

# to run detached:
# git clone https://github.com/p4rkerw/ckd $SCRATCH1/ckd
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/ckd:$HOME/ckd \
# $SCRATCH1:$SCRATCH1"
# bsub -G compute-parkerw -R 'rusage[mem=256GB]' -q general -a 'docker(p4rkerw/sctools:R4.1.0)' \
# -o $SCRATCH1/log_chromvar.out Rscript $SCRATCH1/ckd/atac_aggr_prep/step3_chromvar.R

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

atac_aggr_prep <- here("ckd","atac_aggr_prep")
atacAggr <- readRDS(here(atac_aggr_prep,"step2_anno.rds"))
DefaultAssay(atacAggr) <- "peaks"

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# subset atac object by filtering out alt contigs
# this will result in warnings that can be ignored
# but is required to be compatible with the addmotifs function
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- which(as.character(seqnames(granges(atacAggr))) %in% main.chroms)
atacAggr[["peaks"]] <- subset(atacAggr[["peaks"]], features = rownames(atacAggr[["peaks"]])[keep.peaks])

# add motif information
atacAggr <- AddMotifs(
  object = atacAggr,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

atacAggr <- RegionStats(
  object = atacAggr,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c("-", "-")
)

# compute motif activities using chromvar
atacAggr <- RunChromVAR(
  object = atacAggr,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

saveRDS(atacAggr, here(atac_aggr_prep,"step3_chromVAR.rds"), compress=FALSE)
