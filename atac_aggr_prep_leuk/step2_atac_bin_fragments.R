# this script will tile windows and count the cellranger-atac bam
# the windows will be hg38 excluding an encode blacklist 
# SCRATCH1=/mnt/g/scratch
# docker run -it --rm \
# --workdir $HOME \
# -v /mnt/g/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# -v /mnt/c/reference:$HOME/reference \
# -v /mnt/g/ckd:$HOME/ckd \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/c/scratch" \
# -v /mnt/c/scratch:$HOME/scratch \
# p4rkerw/sctools:R4.1.3 R

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
library_ids <- read.csv(here("ckd","cellranger_atac_aggr_rccleuk","outs","aggregation_csv.csv")) %>%
  dplyr::select(library_id) %>%
  unlist()

# count the files
lapply(library_ids, function(library_id) {
workdir <- here("cellranger_atac_counts", "leukocytes", library_id, "outs")

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

# run epianeufinder
input <- here("cellranger_atac_counts", "leukocytes", library_id,"outs","fragments.tsv.gz")
output <- here("cellranger_atac_counts", "leukocytes", library_id,"outs","epiAneufinder_1MB")
dir.create(output, recursive=TRUE)
blacklist <- "/opt/epiAneufinder/sample_data/hg38-blacklist.v2.bed"

epiAneufinder(input=input, #Enter path to your fragments.tsv file or the folder containing bam files
            outdir=output, #Path to the directory where results should be written 
            blacklist=blacklist, #Path to bed file that contains the blacklisted regions of your genome
            windowSize=1e6, 
            genome="BSgenome.Hsapiens.UCSC.hg38", #Substitute with relevant BSgenome
            exclude=c('chrX','chrY','chrM'), 
            reuse.existing=TRUE,
            title_karyo="Karyogram of sample data", 
            ncores=10,
            minFrags=10000)
})
