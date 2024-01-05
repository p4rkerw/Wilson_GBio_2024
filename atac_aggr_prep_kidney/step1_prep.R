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
# p4rkerw/sctools:R4.1.3 /bin/bash

# to run interactively on the RIS compute1 cluster:
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# $STORAGE1/ckd:$HOME/ckd \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=256GB]' -q general-interactive -a 'docker(p4rkerw/sctools:R4.1.3)' /bin/bash

# to run detached:
# git clone https://github.com/p4rkerw/ckd $SCRATCH1/ckd
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# $STORAGE1/ckd:$HOME/ckd \
# $SCRATCH1:$SCRATCH1"
# bsub -G compute-parkerw \
# -R 'rusage[mem=256GB]' \
# -q general \
# -a 'docker(p4rkerw/sctools:R4.1.3)' \
# -o $SCRATCH1/log.atac.step1.out \
# Rscript $SCRATCH1/ckd/atac_aggr_prep/step1_prep.R

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
dir.create(here("ckd","atac_aggr_prep","plots"), recursive=TRUE, showWarnings=FALSE)
atac_aggr_prep <- here("ckd","atac_aggr_prep")

# define cellranger-atac aggregation input dir
aggr_input_dir <- here("ckd","cellranger_atac_aggr","outs")

# define cellranger-atac count dir
# individual counts are subfolders labeled with corresponding library_id
count_input_dir <- here("cellranger_atac_counts","version_2.1")

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
                pct_reads_in_peaks > 30) %>%
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

# perform label transfer from previously annotated snRNAseq object
rnaAggr = readRDS(here("ckd","rna_control_ref","step2_anno.rds"))
Idents(rnaAggr) <- "lowres.celltype"  

transfer.anchors <- FindTransferAnchors(
  reference = rnaAggr,
  query = atacAggr,
  reduction = 'cca',
  normalization.method = 'SCT',
  reference.assay = 'SCT',
  query.assay = 'RNA'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rnaAggr$lowres.celltype,
  weight.reduction = atacAggr[['lsi']],
  dims = 2:30
)

atacAggr <- AddMetaData(object = atacAggr, metadata = predicted.labels)

# visualize label transfer results
pdf(here(atac_aggr_prep,"plots","step1e_label.pdf"))
p1 <- DimPlot(atacAggr, label=TRUE, group.by="predicted.id") + ggtitle("snATAC predicted celltypes")
AugmentPlot(p1)
p2 <- DimPlot(atacAggr, label=TRUE, group.by="seurat_clusters") + ggtitle("snATAC seurat clusters")
AugmentPlot(p2)
dev.off()

saveRDS(atacAggr, file=here(atac_aggr_prep, "step1e_prep.rds"), compress=FALSE)
      
# compute correlation matrix between predicted.id groups to identify unrelated cell types
# eliminate cells with a predicted.id that does not match the most abundant predicted.id within a cluster (ie its unrelated)
DefaultAssay(atacAggr) <- "peaks"
atacAggr <- FindTopFeatures(atacAggr, assay="peaks")
variable_peaks <- VariableFeatures(atacAggr)
av.access <- AverageExpression(atacAggr, features=variable_peaks, group.by="predicted.id", assays="peaks")
cor.access <- as.data.frame(cor(av.access$peaks))
cor.df <- rownames_to_column(cor.access) %>%
  tidyr::gather(celltype, value, -rowname)
remove(av.access)

# plot the correlation results
pdf(here(atac_aggr_prep,"plots","step1f_correlation.pdf"))
p1 <- cor.df %>%
      dplyr::mutate(fill = ifelse(value < 0, 0, 1)) %>%
      ggplot(aes(x = rowname, y = celltype, fill = fill)) + geom_tile()
print(p1)
dev.off()

# loop through each cluster and return barcodes for cells with a predicted.id that does not match the most
# abundant celltype in the cluster && are negatively correlated with that cluster (ie they are unrelated cell types. eg. ENDO-PT)
bc_anno <- lapply(unique(atacAggr@meta.data$seurat_clusters), function(cluster.sel){
    
    meta <- rownames_to_column(atacAggr@meta.data, var="barcode")

    # cells that are predicted to be a different cell type within cluster
    meta <- dplyr::filter(meta, seurat_clusters == cluster.sel) %>%
      add_count(predicted.id, name="num_predicted.id_in_cluster", sort=TRUE) %>%
      mutate(most_abundant_in_cluster = predicted.id[1]) %>%
      mutate(not_most_abundant_cell = ifelse(predicted.id != most_abundant_in_cluster, 1, 0)) 

    correlated_celltypes <- dplyr::filter(cor.df, celltype == unique(meta$most_abundant_in_cluster)) %>%
      dplyr::filter(value > 0) 
    meta$related_celltype <- ifelse(meta$predicted.id %in% correlated_celltypes$rowname, 1, 0)
    meta$filter_barcode <- ifelse(meta$related_celltype == 0 & meta$not_most_abundant_cell == 1, 1, 0)

    # print cell types filtered from cluster
    print(paste0("Filtering for cluster: ", cluster.sel))
    print(table(meta$predicted.id, meta$filter_barcode))

    rownames(meta) <- meta$barcode
    return(meta)
}) %>% bind_rows() 

atacAggr <- AddMetaData(atacAggr, bc_anno)

# visualize predicted doublets
pdf(here(atac_aggr_prep, "plots","step1g_threshold.pdf"))
p1 <- FeaturePlot(atacAggr, features = "filter_barcode", order=TRUE)
AugmentPlot(p1)
dev.off()

# remove predicted doublets / low quality cells
atacAggr <- subset(atacAggr, subset = filter_barcode == 0) # 0=keep, 1=filter

# recluster
atacAggr <- FindNeighbors(object = atacAggr, reduction = "harmony", dims = 2:30)
atacAggr <- FindClusters(object = atacAggr, verbose = TRUE, algorithm = 1) # Louvain algorithm
atacAggr <- RunUMAP(object = atacAggr, reduction = "harmony", dims = 2:30)

# visualize results
pdf(here(atac_aggr_prep,"plots","step1h_filter.pdf"))
p1 <- DimPlot(atacAggr, label=TRUE, group.by="predicted.id") + ggtitle("snATAC predicted celltypes")
AugmentPlot(p1)
FeaturePlot(atacAggr, features = "peak_region_fragments")
hist(atacAggr@meta.data$nFeature_peaks, breaks=100)
VlnPlot(atacAggr, features = "nFeature_peaks", group.by = "predicted.id")
p2 <- DimPlot(atacAggr, label=TRUE, group.by="seurat_clusters") + ggtitle("snATAC seurat clusters")
AugmentPlot(p2)
dev.off()

# visualize marker genes using gene activity
marker.genes <- c("CUBN","HAVCR1","SLC5A1","SLC5A2","VCAM1","PROM1", # PT and PT-VCAM1+ markers
                  "CFH", # PEC
                  "SLC12A1", # TAL NKCC2
                  "CLDN10", #MTAL TAL2
                  "CLDN16", #CTAL TAL1
                  "S100A2", #ATL
                  "SLC12A3","TRPM6", # DCT1 and DCT2 NCC
                  "SCNN1G","TRPV5", # DCT2/CNT ENaC
                  "CALB1", # CNT
                  "AQP2", # PC
                  "ATP6V0D2", # ICA and ICB
                  "SLC4A1","SLC26A7", # ICA
                  "SLC26A4", # ICB
                  "NPHS1","NPHS2", # PODO
                  "PECAM1","FLT1", # ENDO
                  "IGFBP5","IGFBP7", # PTC and AVR
                  "PLVAP", # PTC and AVR https://www.nature.com/articles/s41467-019-12872-5
                  "EHD3", # GEC
                  "SLC6A6","SLC14A1","AQP1", # EA and DVR
                  "NOS1", # MD
                  "ITGA8","PDGFRB","MEIS2","PIEZO2","REN", # MES and JGA
                  "ACTA2","CALD1", # FIB
                  "PROX1","FLT4","PDPN", # Lymphatics
                  "PTPRC","CD3E","MS4A1","CD19","SDC1", # Lymphocytes and plasma cells
                  "FCGR3A","CD14","CSF1R") # Monocyte / Macrophage

DefaultAssay(atacAggr) <- "RNA"
print("Drawing UMAP markers")
pdf(here(atac_aggr_prep,"plots","step1i_markers.pdf")) 
    lapply(marker.genes, function(gene) {
      tryCatch({plot <- FeaturePlot(atacAggr, features=gene, reduction="umap")
                        AugmentPlot(plot) # downsample
                }, warning=function(w) return(NULL)) 
    })
dev.off()
     
pdf(here(atac_aggr_prep,"plots","step1j_dotplot.pdf"), width=10, height=6)
DefaultAssay(atacAggr) <- "RNA"
DotPlot(atacAggr, features=marker.genes) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  
dev.off()

# plots
pdf(here(atac_aggr_prep,"plots","step1k_cluster.pdf"))
p1 <- DimPlot(atacAggr, label=TRUE, group.by="predicted.id") + ggtitle("snATAC predicted celltypes")
AugmentPlot(p1)
p2 <- DimPlot(atacAggr, label=TRUE, group.by="seurat_clusters") + ggtitle("snATAC seurat clusters")
AugmentPlot(p2)
dev.off()
      
DefaultAssay(atacAggr) <- "peaks"

saveRDS(atacAggr, file=here(atac_aggr_prep, "step1f_prep.rds"), compress=FALSE)
     
