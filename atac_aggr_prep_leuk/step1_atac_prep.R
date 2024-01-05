#!/usr/bin/env Rscript
# this script will preprocess aggregated snATACseq data
# counted and aggregated by cellranger-atac v2.1 without library normalization
#
# to run locally
# counts=/mnt/g/cellranger_atac_counts
# project=/mnt/g/ckd
# SCRATCH1=/mnt/g/scratch
# docker run -it \
# --memory 100g \
# --cpus 4 \
# --workdir $HOME \
# -v $project:$HOME/ckd \
# -v $counts:$HOME/cellranger_atac_counts \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/sctools:R4.2.3 R

library(Signac) # 
library(Seurat) # 
library(GenomeInfoDb) # 
library(harmony) # 
library(EnsDb.Hsapiens.v86)
library(ggplot2) # 
library(patchwork) # 
library(tibble) # 
library(dplyr) # 
library(here) # 
library(data.table) #
library(stringr)
library(future)

start <- Sys.time()

# create output and plots directory
dir.create(here("ckd","atac_aggr_prep_rccleuk","plots"), recursive=TRUE, showWarnings=FALSE)
atac_aggr_prep <- here("ckd","atac_aggr_prep_rccleuk")

# define cellranger-atac aggregation input dir
aggr_input_dir <- here("ckd","cellranger_atac_aggr_rccleuk","outs")

# define cellranger-atac count dir
# individual counts are subfolders labeled with corresponding library_id
count_input_dir <- here("cellranger_atac_counts","leukocytes")

# load aggregated snATACseq data and create a seurat object
counts <- Read10X_h5(here(aggr_input_dir,"filtered_peak_bc_matrix.h5"))
metadata <- read.csv(here(aggr_input_dir ,"singlecell.csv"), header = TRUE, row.names = 1)
aggcsv <- read.csv(here(aggr_input_dir,"aggregation_csv.csv"), header = TRUE, row.names = 1)

chromassay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = here(aggr_input_dir,"fragments.tsv.gz"),
  min.cells = 10,
  min.features = 200
)
remove(counts)

atacAggr <- CreateSeuratObject(
  counts = chromassay,
  assay = "peaks",
  meta.data = metadata
)
remove(chromassay)
remove(metadata)

# Add sample information to the metadata of the Seurat object
gemgroup <- sapply(strsplit(rownames(atacAggr@meta.data), split="-"), "[[", 2) 
current.gemgroups <- seq(length(rownames(aggcsv))) # no. gemgroups is no. samples
library_ids <- rownames(aggcsv)
sampleID <- plyr::mapvalues(gemgroup, from = current.gemgroups, to = library_ids)
atacAggr <- AddMetaData(atacAggr, metadata=data.frame(library_id=sampleID, row.names=rownames(atacAggr@meta.data)))

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(atacAggr) <- annotations

# calculate nucleosome signal
atacAggr <- NucleosomeSignal(object = atacAggr)

# Add the metadata for mitochondrial fragments from individual snATACseq counts into the aggregated dataset 
# Collect total fragment number data from each of the original CellRangerATAC datasets.
current.gemgroups <- seq(rownames(aggcsv))
metaqc <- 
  lapply(current.gemgroups, 
         function(gemgroup) {
           sampleName <- rownames(aggcsv)[gemgroup]
           file <- here(count_input_dir,sampleName,"outs","singlecell.csv")
           df <- read.csv(file, header=TRUE, row.names=1)
           rownames(df) <- paste(substr(rownames(df), 1, 16), gemgroup, sep = "-") # change gemgroup to reflect sample order
           df <- tibble::rownames_to_column(df, var = "barcode")
           return(df)
         }) %>%
  bind_rows() %>%
  tibble::column_to_rownames(var = "barcode") %>%
  dplyr::select(c("total","mitochondrial"))
atacAggr <- AddMetaData(atacAggr,metaqc)
remove(metaqc)

saveRDS(atacAggr, here(atac_aggr_prep, file = "step1a_prep.rds"), compress=FALSE)

# compute TSS enrichment score per cell
# store TSS matrix to enable plotting with TSSPlot function (fast=FALSE)
atacAggr <- TSSEnrichment(object = atacAggr, fast = TRUE)

# QC metrics and filtering
atacAggr$pct_reads_in_peaks <- atacAggr$peak_region_fragments / atacAggr$passed_filters * 100 # %fragments in peaks
atacAggr$mito_ratio <- atacAggr$mitochondrial / atacAggr$total # %fragments mapping to mitochondrial genome
atacAggr$nucleosome_group <- ifelse(atacAggr$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
atacAggr$high.tss <- ifelse(atacAggr$TSS.enrichment > 2, 'High', 'Low')

# draw QC plots
# Note: pct_read_in_peaks histogram does not look normally distributed and is left-skewed with a bimodal distribution
# this likely represents low quality nuclei
pdf(here(atac_aggr_prep,"plots","step1a_qc.pdf"))
VlnPlot(
  object = atacAggr,
  features = c('pct_reads_in_peaks','peak_region_fragments','mito_ratio','passed_filters'), 
  pt.size = 0,
  ncol = 2, group.by = "library_id") + NoLegend()
hist(atacAggr@meta.data$pct_reads_in_peaks, breaks=100)
hist(atacAggr@meta.data$nFeature_peaks, breaks=100)
hist(atacAggr@meta.data$peak_region_fragments,  breaks=100)
hist(atacAggr@meta.data$mito_ratio,  breaks=100)
hist(atacAggr@meta.data$passed_filters,  breaks=100)
FragmentHistogram(object = atacAggr, group.by = 'nucleosome_group')
# TSSPlot(atacAggr, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.off()

saveRDS(atacAggr, here(atac_aggr_prep, file = "step1b_prep.rds"), compress=FALSE)

# calculate library-specific peak_region_fragments thresholds
prf_thresh <- atacAggr@meta.data %>%
  group_by(library_id) %>%
  summarize(mean = mean(peak_region_fragments),
            median = median(peak_region_fragments),
            sd = sd(peak_region_fragments))

# visualize fragment distribution
pdf(here(atac_aggr_prep,"plots","step1a_qc2.pdf"))
atacAggr@meta.data %>%
  group_by(library_id) %>%
  dplyr::mutate(scaledfrags = scale(peak_region_fragments, center = median(peak_region_fragments))) %>%
  ggplot(aes(scaledfrags)) +
  geom_density(aes(y = ..scaled..)) +
  facet_wrap(~library_id)
dev.off()

bckeep <- atacAggr@meta.data %>%
  rownames_to_column(var = "barcode") %>%
  group_by(library_id) %>%
  dplyr::filter(nucleosome_signal < 4,
                TSS.enrichment > 2,
                pct_reads_in_peaks > 50) %>%
  dplyr::select(barcode) %>%
  unlist()
idy <- rownames(atacAggr@meta.data) %in% bckeep

# subset the object
atacAggr <- atacAggr[,idy]

# perform normalization and dimensional reduction of filtered snATAC object
atacAggr <- RunTFIDF(atacAggr)
atacAggr <- FindTopFeatures(atacAggr, min.cutoff = 'q0')
atacAggr <- RunSVD(atacAggr)

# show that first dimension is related to depth of sequencing
pdf(here(atac_aggr_prep, "plots","step1b_depthCor.pdf"))
DepthCor(atacAggr)
dev.off()

# perform batch correction on lsi embeddings
my_harmony_embeddings <- HarmonyMatrix(
  data_mat  = as.matrix(atacAggr@reductions$lsi@cell.embeddings),
  meta_data = atacAggr@meta.data,
  vars_use  = "library_id",
  do_pca = FALSE
)
rownames(my_harmony_embeddings) <- rownames(atacAggr@reductions$lsi@cell.embeddings)

#store the harmony reduction as a custom dimensional reduction called 'harmony' in the default assay
atacAggr[["harmony"]] <- CreateDimReducObject(embeddings = my_harmony_embeddings, key = "harmony_", assay = DefaultAssay(atacAggr))
atacAggr <- FindNeighbors(object = atacAggr, reduction = "harmony", dims = 2:30)
atacAggr <- FindClusters(object = atacAggr, verbose = TRUE, algorithm = 1) # Louvain algorithm
atacAggr <- RunUMAP(object = atacAggr, reduction = "harmony", dims = 2:30)

# plot the clustering results
pdf(here(atac_aggr_prep,"plots","step1c_cluster.pdf"))
p1 <- DimPlot(atacAggr) + ggtitle("Harmony snATAC clustering")
AugmentPlot(p1)
FeaturePlot(atacAggr, features = "pct_reads_in_peaks")
dev.off()

# read in AMULET designated atac doublet metrics and filter barcodes with FDR < 0.05
library_ids <- rownames(aggcsv)
amulet.df <- lapply(seq(library_ids), function(index) {
  df <- fread(here(count_input_dir,library_ids[index],"amulet", "MultipletProbabilities.txt")) %>%
    as.data.frame()
  df$barcode_update <- paste(substr(df$barcode, 1, 16), index, sep = "-") # change barcode suffix to reflect gemgroup
  df <- dplyr::select(df, barcode_update, "p-value", "q-value") %>%
    dplyr::rename(amulet_pval = "p-value") %>%
    dplyr::rename(amulet_qval = "q-value") %>%
    dplyr::mutate(isdoublet = ifelse(amulet_qval < 0.05, 1, 0))
  return(df)
}) %>% 
bind_rows() %>%
column_to_rownames(var = "barcode_update")

# add amulet annotation to object. 1=doublet, 0=not a doublet
atacAggr <- AddMetaData(atacAggr, amulet.df)
      
pdf(here(atac_aggr_prep,"plots","step1d_amulet.pdf")) 
FeaturePlot(atacAggr, feature="isdoublet", order=TRUE)
dev.off()
      
# remove amulet doublets (0=false, 1=true)
saveRDS(atacAggr, file=here(atac_aggr_prep, "step1c_prep.rds"), compress=FALSE)
atacAggr <- subset(atacAggr, isdoublet == 0)

# create gene activity matrix
# note: long run time for this step
gene.activities <- GeneActivity(atacAggr, process_n=500)
atacAggr[['RNA']] <- CreateAssayObject(counts = gene.activities)
atacAggr <- NormalizeData(
  object = atacAggr,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(atacAggr$nCount_RNA)
)
remove(gene.activities)

# recluster
atacAggr <- FindNeighbors(object = atacAggr, reduction = "harmony", dims = 2:30)
atacAggr <- FindClusters(object = atacAggr, verbose = TRUE, algorithm = 1) # Louvain algorithm
atacAggr <- RunUMAP(object = atacAggr, reduction = "harmony", dims = 2:30)

saveRDS(atacAggr, file=here(atac_aggr_prep, "step1d_prep.rds"), compress=FALSE)

