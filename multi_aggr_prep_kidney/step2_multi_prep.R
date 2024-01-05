# this script will preprocess the aggregated multiomes

# to run locally:
# SCRATCH1=/mnt/g/scratch
# docker run -it --rm \
# --workdir $HOME \
# -v /mnt/g/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# -v /mnt/g/cellranger_multi_counts:$HOME/cellranger_multi_counts \
# -v /mnt/g/ckd:$HOME/ckd \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/sctools:R4.1.3 R

# to run interactively on the RIS compute1 cluster:
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/cellranger_multi_counts:$HOME/cellranger_multi_counts \
# $STORAGE1/ckd:$HOME/ckd \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=256GB]' -q general-interactive -a 'docker(p4rkerw/sctools:R4.1.3)' /bin/bash

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
library(stringr)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(here)

set.seed(1234)

# create a folder for qc plots
dir.create(here("ckd","multi_aggr_prep","plots"), recursive = TRUE, showWarnings = FALSE)
plots <- here("ckd","multi_aggr_prep","plots")

# load the RNA and ATAC data
counts <- Read10X_h5(here("ckd","cellranger_arc_aggr","outs","filtered_feature_bc_matrix.h5"))
fragpath <- here("ckd","cellranger_arc_aggr","outs","atac_fragments.tsv.gz")
aggcsv <- fread(here("ckd","cellranger_arc_aggr","outs","aggr.csv")) %>%
  rownames_to_column(var = "gem")

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA data
multi <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
multi[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
rm(counts)

# qc metrics for ATAC assay
DefaultAssay(multi) <- "ATAC"
multi <- NucleosomeSignal(multi)
multi <- TSSEnrichment(multi)

# visualize atac qc metrics
pdf(here(plots, "step2a_qc.pdf"), width=10, height=6)
VlnPlot(
  object = multi,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)
dev.off()

# filter out low quality cells
multi <- subset(
  x = multi,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 20000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 1 &
    TSS.enrichment > 2
)

# normalize rna counts
DefaultAssay(multi) <- "RNA"
multi <- SCTransform(multi)
multi <- RunPCA(multi)

# add library_id annotation to metadata
gem <- str_split(rownames(multi@meta.data), pattern="-", simplify=TRUE)[,2]
multi@meta.data$gem <- gem
meta <- multi@meta.data %>%
  left_join(aggcsv, by = "gem") 
multi@meta.data$library_id <- meta$library_id

# RNA-based doublet identification with the assumption that doublets represent 5% of cells.
FindDoublets <- function(library_id, seurat_aggregate) {
  ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
  rnaAggr <- seurat_aggregate
  seurat_obj <- subset(rnaAggr, idents = library_id)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list_kidney <- paramSweep_v3(seurat_obj, PCs = 1:20, sct = F, num.cores=1)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  pK <- bcmvn_kidney %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    filter(BCmetric == max(BCmetric)) %>%
    select(pK) 
  pK <- as.numeric(as.character(pK[[1]]))
  seurat_doublets <- doubletFinder_v3(seurat_obj, PCs = 1:20, pN = 0.25, pK = pK,
                               nExp = round(0.05*length(seurat_obj@active.ident)), 
                               reuse.pANN = FALSE, sct = F)
 
  # create doublet groupings and visualize results
  DF.class <- names(seurat_doublets@meta.data) %>% str_subset("DF.classifications")
  pANN <- names(seurat_doublets@meta.data) %>% str_subset("pANN")
  
  p1 <- ggplot(bcmvn_kidney, aes(x=pK, y=BCmetric)) +
    geom_bar(stat = "identity") + 
    ggtitle(paste0("pKmax=",pK)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p2 <- DimPlot(seurat_doublets, group.by = DF.class)
  p3 <- FeaturePlot(seurat_doublets, features = pANN)
  
  # visualize doublet metrics
  pdf(here(plots, paste0("step2_doublets.",library_id,".pdf")))
    print(p1)
    print(p2)
    print(p3)
  dev.off()
  
  # create a df of barcodes and doublet designations
  df_doublet_barcodes <- as.data.frame(cbind(rownames(seurat_doublets@meta.data), seurat_doublets@meta.data[[DF.class]]))
  return(df_doublet_barcodes)
}

# take an aggregated snRNA seurat object and a list of library_id to find doublets. return a df of doublet barcodes
# send DimPlot and FeaturePlot of doublets for each library to plots dir
Idents(multi) <- "library_id"
DefaultAssay(multi) <- "RNA"
list.doublet.bc <- lapply(aggcsv$library_id, function(x) {FindDoublets(x, seurat_aggregate = multi)})
df_doublet_id <- list.doublet.bc %>%
  bind_rows() %>%
  dplyr::rename("df_doublet_id" = "V2") %>%
  tibble::column_to_rownames(var = "V1") # this is the barcode column
table(df_doublet_id) # quantify total doublet vs. singlet calls (expect ~6% doublets)
  
# add doublet calls to aggregated snRNA object as doublet_id in meta.data slot
multi <- AddMetaData(multi, df_doublet_id)

# normalize peak counts
DefaultAssay(multi) <- "ATAC"
multi <- FindTopFeatures(multi, min.cutoff = 5)
multi <- RunTFIDF(multi)
multi <- RunSVD(multi)

# build a joint neighbor graph using both assays
multi <- FindMultiModalNeighbors(
  object = multi,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:30, 2:30),
  modality.weight.name = c("SCT.weight","ATAC.weight"), 
  verbose = TRUE
)

# build a joint UMAP visualization
multi <- RunUMAP(
  object = multi,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

# visualize joint UMAP by gem group
# this embedding is not corrected for batch effect and is expected to separate by donor
Idents(multi) <- "gem"
DimPlot(multi, label = TRUE, reduction = "umap")
pdf(here(plots, "step2b_wnn_umap.pdf"))
DimPlot(multi, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()
dev.off()

saveRDS(multi, here("ckd","multi_aggr_prep","multi_withdb.rds"), compress = FALSE)

# split by gem group
obj.list <- SplitObject(multi, split.by = "gem")

# find integration anchors between donors
integration.anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  anchor.features = rownames(multi),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = multi[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)

# visualize integrated umap
pdf(here(plots, "step2c_integrated_umap.pdf"))
DimPlot(integrated, group.by = "gem")
dev.off()

# read in AMULET designated atac doublet metrics 
library_ids <- aggcsv$library_id
amulet.df <- lapply(seq(library_ids), function(index) {
  df <- fread(here("cellranger_multi_counts",library_ids[index],"outs","amulet", "MultipletProbabilities.txt")) %>%
    as.data.frame()
  df$barcode_update <- paste(substr(df$barcode, 1, 16), index, sep = "-") # change barcode suffix to reflect gemgroup
  df <- dplyr::select(df, barcode_update, "p-value", "q-value") %>%
    dplyr::rename(amulet_pval = "p-value") %>%
    dplyr::rename(amulet_qval = "q-value") %>%
    dplyr::mutate(am_doublet_id = ifelse(amulet_pval < 0.05, 1, 0))
  return(df)
}) %>% 
bind_rows() %>%
column_to_rownames(var = "barcode_update")

# add amulet annotation to object. 1=doublet, 0=not a doublet
integrated <- AddMetaData(integrated, amulet.df)

# switch doubletfinder annotation to numeric
integrated@meta.data$df_doublet_id <- ifelse(integrated@meta.data$df_doublet_id == "Doublet",1,0)

# count doublets that were not already filtered and visualize
table(amulet=integrated@meta.data$am_doublet_id, df=integrated@meta.data$df_doublet_id)
pdf(here(plots,"step2d_doublet_umap.pdf"))
p1 <- FeaturePlot(integrated, features="am_doublet_id", reduction="umap", order=TRUE) + ggtitle("Amulet Doublets")
AugmentPlot(p1)
p2 <- FeaturePlot(integrated, features="df_doublet_id", reduction="umap", order=TRUE) + ggtitle("DoubletFinder Doublets")
AugmentPlot(p2)
dev.off()

# remove amulet and doubletfinder doublets
# only retain cells that were annotated as singlets by both tools
integrated <- subset(
  x = integrated,
  subset = am_doublet_id == 0 &
  df_doublet_id == 0)

saveRDS(integrated, here("ckd","multi_aggr_prep","multi_integrate.rds"), compress=FALSE)

# build a joint neighbor graph without the doublets
integrated <- FindNeighbors(integrated, reduction = "integrated_lsi")
integrated <- FindClusters(integrated)

# build a joint UMAP visualization
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)

# plot the updated umap
pdf(here(plots, "step2e_nodb_umap.pdf"))
DimPlot(integrated, label = TRUE, repel = TRUE, reduction = "umap", group.by = "gem") + NoLegend()
dev.off()

# run harmony to batch correct the integrated embeddings 
my_harmony_embeddings <- HarmonyMatrix(
  data_mat  = as.matrix(integrated@reductions$integrated_lsi@cell.embeddings),
  meta_data = integrated@meta.data,
  vars_use  = "gem",
  do_pca = FALSE
)
rownames(my_harmony_embeddings) <- rownames(integrated@reductions$integrated_lsi@cell.embeddings)

#store the harmony reduction as a custom dimensional reduction called 'harmony' in the default assay
integrated[["harmony"]] <- CreateDimReducObject(embeddings = my_harmony_embeddings, key = "harmony_", assay = DefaultAssay(integrated))
integrated <- FindNeighbors(object = integrated, reduction = "harmony", dims = 2:30)
integrated <- FindClusters(object = integrated, verbose = TRUE, algorithm = 1) # Louvain algorithm
integrated <- RunUMAP(object = integrated, reduction = "harmony", dims = 2:30)

# plot the harmony umap
pdf(here(plots, "step2f_harmony_umap.pdf"))
DimPlot(integrated, reduction = "umap")
DimPlot(integrated, reduction = "umap", group.by = "gem")
dev.off()

# read in annotated rna object
rnaAggr <- readRDS(here("ckd","rna_control_ref","step2_anno.rds"))

# label transfer celltype prediction
transfer.anchors <- FindTransferAnchors(
  reference = rnaAggr,
  query = integrated,
  reduction = 'cca',
  normalization.method = 'SCT',
  reference.assay = 'SCT',
  query.assay = 'RNA'
)

# update predicted cell type metadata
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rnaAggr$celltype,
  weight.reduction = integrated[['integrated_lsi']],
  dims = 2:30
)

integrated <- AddMetaData(integrated, predicted.labels)

pdf(here(plots, "step2g_predict_umap.pdf"))
DimPlot(integrated, group.by = "predicted.id", label=TRUE)
DimPlot(integrated, group.by = "seurat_clusters", label=TRUE)
dev.off()

# visualize marker genes using gene expression
marker.genes <- c("CUBN","HAVCR1","SLC5A1","SLC5A2","VCAM1","PROM1", # PT and PT-VCAM1+ markers
                  "CFH", # PEC
                  "SLC12A1", # TAL NKCC2
                  "CLDN10", #MTAL (TAL2)
                  "CLDN16", #CTAL (TAL1)
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
                  "PTPRC","CD3E","MS4A1", # Lymphocytes
                  "FCGR3A","CD14","CSF1R") # Monocyte / Macrophage

DefaultAssay(integrated) <- "SCT"
print("Drawing UMAP markers")
pdf(here(plots,"step2h_markers.pdf")) 
    lapply(marker.genes, function(gene) {
      tryCatch({plot <- FeaturePlot(integrated, features=gene, reduction="umap")
                AugmentPlot(plot) #downsample
                }, warning=function(w) return(NULL)) 
    })
dev.off()
      
pdf(here(plots,"step2i_dotplot.pdf"), width=10, height=6)
DefaultAssay(integrated) <- "RNA"
DotPlot(integrated, features=marker.genes) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  
dev.off()
      
# annotate clusters
integrated <- RenameIdents(integrated,
  '0' = 'PCT',
  '1' = 'DCT1',
  '2' = 'PCT',
  '3' = 'TAL1',
  '4' = 'PCT',
  '5' = 'PT_VCAM1',
  '6' = 'PST',
  '7' = 'TAL2',
  '8' = 'ENDO',
  '9' = 'DCT2',
  '10' = 'PC',
  '11' = 'ICB',
  '12' = 'FIB_VSMC_MC',
  '13' = 'ICA',
  '14' = 'MONO',
  '15' = 'PEC',
  '16' = 'TCELL',
  '17' = 'PCT',
  '18' = 'BCELL',
  '19' = 'PODO',
  '20' = 'ENDO')

# update levels
levels(integrated) <- c("PCT","PST","PT_VCAM1","PEC",
          "TAL1","TAL2","DCT1","DCT2","PC",
          "ICA","ICB","PODO","ENDO","FIB_VSMC_MC",
          "TCELL","BCELL","MONO")


integrated@meta.data$celltype <- Idents(integrated)


pdf(here(plots,"step2j_anno.pdf"))
p1 <- DimPlot(integrated, label=TRUE, group.by="celltype") + ggtitle("snATAC annotated celltypes")
AugmentPlot(p1)
p2 <- DimPlot(integrated, label=TRUE, group.by="predicted.id") + ggtitle("snATAC predicted.id")
AugmentPlot(p2)
dev.off()     
      
# create a list of cell type barcodes and write to file
bc <- integrated@meta.data %>%
  rownames_to_column(var = "barcode") %>%
  dplyr::select(barcode, celltype, library_id) %>%
  dplyr::mutate(barcode = str_split(barcode, pattern="-", simplify=TRUE)[,1]) 
write.table(bc, file=here("ckd","multi_aggr_prep","multi_barcodes.csv"), sep=",", row.names=FALSE, quote=FALSE)
saveRDS(integrated@meta.data, here("ckd","multi_aggr_prep","meta.rds"))

# save object
saveRDS(integrated, here("ckd","multi_aggr_prep","multiome.rds"), compress=FALSE)
