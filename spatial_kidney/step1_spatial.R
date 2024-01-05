# this script will analyze ST libraries from Humphreys lab only
# SCRATCH1=/mnt/g/scratch
# docker run -it --rm \
# --workdir $HOME \
# -v /mnt/g/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# -v /mnt/g/cellranger_multi_counts:$HOME/cellranger_multi_counts \
# -v /mnt/g/spaceranger_visium_counts:$HOME/spaceranger_visium_counts \
# -v /mnt/g/ckd:$HOME/ckd \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/sctools:R4.1.3 R

library(Seurat)
library(tidyr)
library(dplyr)
library(tibble)
library(stringr)
library(data.table)
library(here)
library(ggplot2)
library(harmony)

# note that this older version of Load10X_Spatial only looks for tissue_positions_list.csv in the outs/spatial folder
# so the file must be renamed to tissue_positions.csv for this function to work
library_ids <- c("MGI3535_A1_010322NHK","MGI3535_B1_041921NHK","MGI3535_C1_011319NHK","MGI3535_D1_050619NHK",
"MGI3779_A1_110122_NX","MGI3779_B1_50521_4428AID","MGI3779_C1_110722_463AJKE","MGI3779_D1_070722_446AJGE")
spatial.ls <- lapply(library_ids, function(library_id) {
  
  file <- here("spaceranger_visium_counts",library_id,"outs","spatial","tissue_positions.csv")
  newfile <- here("spaceranger_visium_counts",library_id,"outs","spatial","tissue_positions_list.csv")
  if (file.exists(file)) {
    file.copy(file, newfile)
  }
    
  spatial <- Load10X_Spatial(
  here("spaceranger_visium_counts", library_id,"outs"),
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL)
  
  # update meta data
  spatial$library_id <- library_id

  # preprocess
  spatial <- NormalizeData(spatial)
  spatial <- ScaleData(spatial)
  spatial <- FindVariableFeatures(spatial)
  spatial <- RunPCA(spatial)
})

anchors <- FindIntegrationAnchors(object.list = spatial.ls, reduction = "cca", dims = 1:30)
spatial <- IntegrateData(anchorset = anchors, dims = 1:30)
dir.create(here("ckd","spatial"))

# preprocess reduction
spatial <- FindVariableFeatures(spatial)
spatial <- ScaleData(spatial)
spatial <- RunPCA(spatial, reduction.key = "integrated")
spatial <- RunUMAP(spatial, dims = 1:30)

# visualize prior to batch correction
plots <- here("ckd","spatial","plots")
dir.create(here(plots))
pdf(here(plots, "spatial_integrated_umap.pdf"))
DimPlot(spatial, reduction = "umap")
DimPlot(spatial, reduction = "umap", group.by = "library_id")
FeaturePlot(spatial, features = "nCount_Spatial")
dev.off()

# save obj
saveRDS(spatial, here("ckd","spatial","spatial_step1A.rds"), compress=FALSE)

# run harmony to batch correct the integrated embeddings 
my_harmony_embeddings <- HarmonyMatrix(
  data_mat  = as.matrix(spatial@reductions$pca@cell.embeddings),
  meta_data = spatial@meta.data,
  vars_use  = "library_id",
  do_pca = FALSE
)
rownames(my_harmony_embeddings) <- rownames(spatial@reductions$pca@cell.embeddings)

#store the harmony reduction as a custom dimensional reduction called 'harmony' in the default assay
spatial[["harmony"]] <- CreateDimReducObject(embeddings = my_harmony_embeddings, key = "harmony_", assay = DefaultAssay(spatial))
spatial <- FindNeighbors(object = spatial, reduction = "harmony", dims = 1:30)
spatial <- FindClusters(object = spatial, verbose = TRUE, algorithm = 1) # Louvain algorithm
spatial <- RunUMAP(object = spatial, reduction = "harmony", dims = 1:30)

# plot the harmony umap
plots <- here("ckd","spatial","plots")
dir.create(here(plots))
pdf(here(plots, "spatial_harmony_umap.pdf"))
DimPlot(spatial, reduction = "umap")
DimPlot(spatial, reduction = "umap", group.by = "library_id")
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
      
pdf(here(plots,"spatial_dotplot.pdf"), width=10, height=6)
DotPlot(spatial, features=marker.genes) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  
dev.off()

# annotate
spatial <- RenameIdents(spatial,
                        '0' = 'PT',
                        '1' = 'TAL',
                        '2' = 'PT_INJ',
                        '3' = 'TL',
                        '4' = 'PT',
                        '5' = 'CD',
                        '6' = 'CD',
                        '7' = 'DCT',
                        '8' = 'CD',
                        '9' = 'GLOM',
                        '10' = 'GLOM',
                        '11' = 'FIB_VSMC',
                        '12' = 'PT_INJ',
                        '13' = 'TL')
spatial$celltype <- Idents(spatial)

# viz spot annotations
pdf(here(plots,"spatial_anno.pdf"), width=10, height=6)
DimPlot(spatial, group.by = "celltype", label=TRUE)
dev.off()

levels(spatial) <- c("PT","PT_INJ","TL","TAL","DCT","CD","FIB_VSMC","GLOM")
spatial$celltype <- Idents(spatial)

# save the obj
saveRDS(spatial, here("ckd","spatial","step1b_spatial.rds"), compress=FALSE)
