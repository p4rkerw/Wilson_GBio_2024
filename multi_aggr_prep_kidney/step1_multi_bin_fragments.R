# this script will tile windows and count the cellranger-multi atac fragments file (including chrY)
# the windows will are hg38 excluding an encode blacklist 
# output is cytoband_counts10k.rds which are 1MB tiles for cells with at least 10k counts in whitelist regions

# SCRATCH1=/mnt/g/scratch
# docker run -it \
# --workdir $HOME \
# -v /mnt/g/cellranger_multi_counts:$HOME/cellranger_multi_counts \
# -v /mnt/g/reference:$HOME/reference \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# -v /mnt/g/scratch:$HOME/scratch \
# p4rkerw/sctools:R4.1.3 /bin/bash

# # to run RIS interactive
# # # export docker volumes
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/cellranger_multi_counts:$HOME/cellranger_multi_counts \
# $STORAGE1/reference:$HOME/reference \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=64GB]' -sp 99 -n 4 -q general-interactive -a 'docker(p4rkerw/sctools:R4.1.3)' /bin/bash

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

# load blacklist
blacklist <- "/opt/epiAneufinder/sample_data/hg38-blacklist.v2.bed"
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
library_ids <- c("1-27Nx","2-15Nx","A2","AJDL105","090922Nx","091422Nx","AIIM164","AIL5160","AJDV174")

# count the files
lapply(library_ids, function(library_id) {
workdir <- here("cellranger_multi_counts", library_id, "outs")

# read in library fragments file from cellranger-atac outs folder
fragments <- fread(here(workdir,"atac_fragments.tsv.gz"))
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
