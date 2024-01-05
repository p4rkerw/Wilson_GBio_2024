# annotate the aggregated atac obj using bridge integration

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
# -v /mnt/c/reference:$HOME/reference \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/azimuth:1.0 R

BiocManager::install('biovizBase')

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(here)
library(SeuratDisk)

# download 10x multiome bridge from 
# https://www.10xgenomics.com/resources/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0

# the 10x hdf5 file contains both data types.
inputdata.10x <- Read10X_h5(here("reference","pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"))
# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
# Create Seurat object
obj.multi <- CreateSeuratObject(counts = rna_counts)
# Get % of mitochondrial genes
obj.multi[["percent.mt"]] <- PercentageFeatureSet(obj.multi, pattern = "^MT-")

# add the ATAC-seq assay
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
# Get gene annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# Change style to UCSC
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
# File with ATAC per fragment information file
frag.file <- here("reference","pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz")
# Add in ATAC-seq data as ChromatinAssay object
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
# Add the ATAC assay to the multiome object
obj.multi[["ATAC"]] <- chrom_assay
# Filter ATAC data based on QC metrics
obj.multi <- subset(
  x = obj.multi,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)

# load the query assay
query.atac <- readRDS(here("ckd","atac_aggr_prep_rccleuk","step1d_prep.rds"))

# Requantify query ATAC to have same features as multiome ATAC dataset
requant_multiome_ATAC <- FeatureMatrix(
  fragments = Fragments(query.atac),
  features = granges(obj.multi[['ATAC']]),
  cells = Cells(query.atac)
)
# Create assay with requantified ATAC data
ATAC_assay <- CreateChromatinAssay(
  counts = requant_multiome_ATAC,
  fragments = here("ckd","cellranger_atac_aggr_rccleuk","outs","fragments.tsv.gz"),
  annotation = annotations
)
# Create Seurat sbject
obj.atac  <- CreateSeuratObject(counts = ATAC_assay, assay = 'ATAC')
# obj.atac[['peak.orig']] <- query.atac
# obj.atac <- subset(obj.atac, subset = nCount_ATAC < 7e4 & nCount_ATAC > 2000)

# load the reference
obj.rna <- LoadH5Seurat(here("reference","pbmc_multimodal.h5seurat"))

# normalize multiome RNA
DefaultAssay(obj.multi) <- "RNA"
obj.multi <- SCTransform(obj.multi, verbose = FALSE)
# normalize multiome ATAC
DefaultAssay(obj.multi) <- "ATAC"
obj.multi <- RunTFIDF(obj.multi)
obj.multi <- FindTopFeatures(obj.multi, min.cutoff = "q0")
# obj.multi <- FindVariableFeatures(obj.multi)

dims.atac <- 2:50
dims.rna <- 1:50
DefaultAssay(obj.multi) <-  "RNA"
DefaultAssay(obj.rna) <- "SCT"
obj.rna.ext <- PrepareBridgeReference(
  reference = obj.rna,
  bridge = obj.multi,
  reference.reduction = "spca",
  reference.dims = dims.rna,
  normalization.method = "SCT")

# update assay compatibility for new functions in seurat v5
obj.atac[["ATAC"]] <- as(object = obj.atac[["ATAC"]], Class = "Assay5")

bridge.anchor <- FindBridgeTransferAnchors(
  extended.reference = obj.rna.ext,
  query = obj.atac,
  query.assay = "ATAC",
  reduction = "lsiproject",
  dims = dims.atac)

obj.atac <- MapQuery(
  anchorset = bridge.anchor,
  reference = obj.rna.ext,
  query = obj.atac,
  refdata = list(
    l1 = "celltype.l1",
    l2 = "celltype.l2",
    l3 = "celltype.l3"),
  reduction.model = "wnn.umap")

saveRDS(obj.atac@meta.data, here("ckd","atac_aggr_prep_rccleuk","step3_predictions_meta.rds"))

addmeta <- data.frame(predicted.l1 = obj.atac$predicted.l1)
rownames(addmeta) <- rownames(obj.atac@meta.data)
query.atac <- AddMetaData(query.atac, addmeta)

# plot the clustering results
plots <- here("ckd","atac_aggr_prep_rccleuk","plots")

pdf(here(plots,"step4_anno.pdf"))
p1 <- DimPlot(query.atac, group.by = "predicted.l1", label=TRUE) + ggtitle("Label transfer annotation")
AugmentPlot(p1)
dev.off()


