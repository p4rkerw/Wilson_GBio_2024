# SCRATCH1=/mnt/g/scratch
# docker run -it --rm \
# --workdir $HOME \
# -v /mnt/s:$HOME/data \
# -v /mnt/g/reference:$HOME/reference \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# -v /mnt/g/scratch:$HOME/scratch \
# p4rkerw/sctools:R4.3.2 R

suppressMessages ({
library(Signac)
library(Seurat)
library(hdf5r)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(here)
library(GenomicRanges)
library(dplyr)
})

plan("multicore", workers = 6)
options(future.globals.maxSize = 8000 * 1024^2)

args <- commandArgs(trailingOnly=TRUE)
library_id <- args[1]
datadir <- args[2]
# library_id <- "SAMN20460922"
# datadir <- "data/cellranger_atac_counts/leukocytes/version_2.1"
inputdir <- here(datadir, library_id, "outs")
outputdir <- here(inputdir, "chasm")
dir.create(here(outputdir), recursive=TRUE)

set.seed(1234)
# load the RNA and ATAC data
cellranger_counts <- Read10X_h5(here(inputdir, "filtered_peak_bc_matrix.h5"))
fragpath <- here(inputdir, "fragments.tsv.gz")

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

# create ATAC assay and add it to the object
chrom_assay <- CreateChromatinAssay(
  counts = cellranger_counts,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation,
  min.cells = 10,
  min.features = 200
)

srat <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "cellranger_peaks"
)

srat <- NucleosomeSignal(srat)
srat <- TSSEnrichment(srat, fast = TRUE)

# filter out low quality cells
srat <- subset(
  x = srat,
  subset = nucleosome_signal < 4 &
           TSS.enrichment > 2 
)

# call peaks using MACS2
# macs2 pseudobulk peak calling may bias towards more abundant cell types and miss peaks in less common cell types
# it likely has the same drawbacks as the cellranger-atac / cellranger-arc algorithm
peaks.gr <- CallPeaks(srat, macs2.path = "/usr/local/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks.gr <- keepStandardChromosomes(peaks.gr, pruning.mode = "coarse")
peaks.gr <- subsetByOverlaps(x = peaks.gr, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(srat),
  features = peaks.gr,
  cells = colnames(srat)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
srat[["macs2_peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

# generate a 1mb genome bin matrix
windows <- tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg38), tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
windows <- GenomeInfoDb::keepStandardChromosomes(windows, pruning.mode = "coarse")
windows <- dropSeqlevels(windows, value = "chrM", pruning.mode = 'coarse')

gbin_counts <- FeatureMatrix(
  fragments = Fragments(srat),
  features = windows,
  cells = colnames(srat)
)

# create a new assay using the bin counts and add it to the Seurat object
srat[["gbin_1mb"]] <- CreateChromatinAssay(
  counts = gbin_counts,
  fragments = fragpath,
  annotation = annotation
)

# filter the genome bins by N content and remove blacklist regions
mcols(windows)$wSeq <- as.character(seqnames(windows))
mcols(windows)$wStart <- BiocGenerics::start(windows)
mcols(windows)$wEnd <- BiocGenerics::end(windows)
windows <- subsetByOverlaps(windows, blacklist_hg38_unified, invert = TRUE)
names(windows) <- paste0("w",seq_along(windows))
mcols(windows)$name <- names(windows)
windowSplit <- split(windows, as.character(seqnames(windows)))
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
# keep windows with < 0.1% N content
windows <- windows[which(windows$N < 0.001)]

gbin_counts <- FeatureMatrix(
  fragments = Fragments(srat),
  features = windows,
  cells = colnames(srat)
)

# create a new assay using the bin counts and add it to the Seurat object
srat[["gbin_1mb_filtered"]] <- CreateChromatinAssay(
  counts = gbin_counts,
  fragments = fragpath,
  annotation = annotation
)

# save macs2 peak count matrix
cellxpeak <- srat@assays[["macs2_peaks"]]@data
saveRDS(cellxpeak, here(outputdir, 'cellx_macs2_peak.rds' ))

# save cellranger peak count matrix
cellxpeak <- srat@assays[["cellranger_peaks"]]@data
saveRDS(cellxpeak, here(outputdir, 'cellx_cellranger_peak.rds' ))

# save bin count matrix
cellxbin <- srat@assays[["gbin_1mb"]]@data
saveRDS(cellxbin, here(outputdir, 'cellx_gbin_1mb.rds' ))

# save bin count matrix
cellxbin <- srat@assays[["gbin_1mb_filtered"]]@data
saveRDS(cellxbin, here(outputdir, 'cellx_gbin_1mb_filtered.rds' ))

# save seurat object
saveRDS(srat, here(outputdir, 'srat.rds'), compress=FALSE)
