# this script will tile windows and count the cellranger-atac bam
# the windows will be hg38 excluding an encode blacklist 
# SCRATCH1=/mnt/g/scratch
# docker run -it \
# --workdir $HOME \
# -v /mnt/g/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# -v /mnt/g/cellranger_rna_counts:$HOME/cellranger_rna_counts \
# -v /mnt/g/reference:$HOME/reference \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# -v /mnt/g/scratch:$HOME/scratch \
# p4rkerw/sctools:R4.1.3 /bin/bash

# # to run RIS interactive
# # # export docker volumes
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# $STORAGE1/reference:$HOME/reference \
# $SCRATCH1:$SCRATCH1"
# # to run interactive
# bsub -Is -G compute-parkerw -R 'rusage[mem=32GB]' -sp 99 -n 4 -q general-interactive -a 'docker(p4rkerw/sctools:R4.1.3)' /bin/bash

library(epiAneufinder)
library(plyranges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
library(utils)
library(GenomicRanges)
library(here)
library(dplyr)
library(tidyr)
library(tibble)

# load blacklist shipped with epianeufinder package
blacklist <- here("reference","hg38-blacklist.v2.bed")
blacklist <- fread(blacklist) %>%
  dplyr::rename(seqnames = V1, start = V2, end = V3) %>%
  makeGRangesFromDataFrame() %>%
  sort()

# tile the genome in 1MB bins
windows <- tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg38), tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
windows <- GenomeInfoDb::keepStandardChromosomes(windows, pruning.mode = "coarse")
windows <- dropSeqlevels(windows, value = "chrM", pruning.mode = 'coarse')
mcols(windows)$wSeq <- as.character(seqnames(windows))
mcols(windows)$wStart <- BiocGenerics::start(windows)
mcols(windows)$wEnd <- BiocGenerics::end(windows)
overlaps <- findOverlaps(windows, blacklist)
idx <- setdiff(1:length(windows), S4Vectors::queryHits(overlaps))
windowsBL <- windows[idx]
names(windowsBL) <- paste0("w",seq_along(windowsBL))
mcols(windowsBL)$name <- names(windowsBL)
windowSplit <- split(windowsBL, as.character(seqnames(windowsBL)))
windowNuc <- lapply(seq_along(windowSplit), function(x){
  message(sprintf("%s of %s", x, length(windowSplit)))
  chrSeq <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38,names(windowSplit)[x])
  grx <- windowSplit[[x]]
  aFreq <- Biostrings::alphabetFrequency(Biostrings::Views(chrSeq, ranges(grx)))
  mcols(grx)$GC <- rowSums(aFreq[, c("G","C")]) / rowSums(aFreq)
  mcols(grx)$AT <- rowSums(aFreq[, c("A","T")]) / rowSums(aFreq)
  return(grx)
}) %>% GRangesList %>% unlist %>% sortSeqlevels %>% sort
windowNuc$N <- 1 - (windowNuc$GC + windowNuc$AT)
windows <- windowNuc

# write windows to file with < 0.1% N content
keep <- which(windows$N < 0.001)
windowSummary <- windows[keep,]
windowSummary %>%
  as.data.frame() %>%
  dplyr::select(seqnames, start, end) %>%
  write.table(file = here("reference","hg38.keep.bed"), sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

############################
# start counting windows from bam files
# library_id <- "Control_1"
library_ids <- c("Control_1","Control_2","Control_3","Control_4","Control_5","Control_6",
                 "CKD_1","CKD_2","CKD_3","CKD_4","CKD_5",
                 "DN_1","DN_2","DN_3","DN_4","DN_5",
                 "SAMN18736215","SAMN18736216","SAMN27505541","SAMN27505542","SAMN27505543","SAMN27505544")
args <- data.frame(library_id = library_ids, version = "version_2.1")

# count the files
lapply(seq(nrow(args)), function(row) {
library_id <- args[row,]$library_id
version <- args[row,]$version
workdir <- here("cellranger_atac_counts", version, library_id, "outs")

# read in library fragments file from cellranger-atac outs folder
fragments <- fread(here(workdir,"fragments.tsv.gz"))
colnames(fragments) <- c('seqnames','start','end','barcode','pcr')
fragments <- plyranges::as_granges(fragments)

# generate count matrix
counts <- epiAneufinder::generateCountMatrix(fragments,
                              windows, 
                              by="barcode", 
                              minFrags = 10000)

# prepare data for GC correction
peaks <- as.data.table(assays(counts)$counts)
colnames(peaks) <- paste0('cell-', colnames(peaks))
rowinfo <- as.data.table(rowRanges(counts))
peaks <- cbind(rowinfo, peaks)

saveRDS(peaks, here(workdir,"cytoband_counts10k.rds"))
})
