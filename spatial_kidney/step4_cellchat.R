# cellchat with secreted signaling only

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
# p4rkerw/sctools:R4.1.3d R

library(CellChat)
library(patchwork)
library(here)
library(tibble)
library(tidyr)
library(Seurat)
library(parallel)
library(ggplot2)
options(stringsAsFactors = FALSE)

spatial <- readRDS(here("ckd","spatial","step1b_spatial.rds"))
anno <- data.frame(celltype = spatial$celltype,
                   rownames = rownames(spatial@meta.data),
                   library_id = spatial$library_id)
anno <- anno %>%
  tidyr::separate(rownames, sep = "_", into = c("barcode",NA), remove=FALSE)

# subset spatial obj for single male lib
dir.create(here("ckd","spatial","cellchat"))

library_ids <- c("MGI3535_A1_010322NHK","MGI3535_B1_041921NHK","MGI3535_C1_011319NHK","MGI3535_D1_050619NHK",
"MGI3779_A1_110122_NX","MGI3779_B1_50521_4428AID","MGI3779_C1_110722_463AJKE","MGI3779_D1_070722_446AJGE")
# library_ids <- list.dirs(here("spaceranger_visium_counts"), full.names = FALSE, recursive = FALSE)
cellchat.ls <- mclapply(library_ids, function(library_id) {
  sub <- Load10X_Spatial(
  data.dir = here("spaceranger_visium_counts", library_id, "outs"),
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL)
  
  meta <- sub@meta.data %>% 
  rownames_to_column(var = "barcode") 
  meta$library_id <- library_id
  meta <- meta %>% left_join(anno, by = c("library_id", "barcode"))
  rownames(meta) <- meta$barcode

  sub <- NormalizeData(sub)
  sub <- AddMetaData(sub, meta)
  
  # Prepare input data for CellChat analysis
  Idents(sub) <- "celltype"
  data.input = GetAssayData(sub, slot = "data", assay = "Spatial") # normalized data matrix
  meta = data.frame(labels = Idents(sub), row.names = names(Idents(sub))) # manually create a dataframe consisting of the cell labels
  unique(meta$labels) # check the cell labels

  # load spatial imaging information
  # Spatial locations of spots from full (NOT high/low) resolution images are required
  spatial.locs = GetTissueCoordinates(sub, scale = NULL, cols = c("imagerow", "imagecol")) 

  # Scale factors and spot diameters of the full resolution images 
  scale.factors = jsonlite::fromJSON(txt = here("spaceranger_visium_counts",library_id,"outs","spatial","scalefactors_json.json"))
  scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
                       fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef)

  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                             datatype = "spatial", coordinates = spatial.locs, scale.factors = scale.factors)
  
  CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data

  # use a subset of CellChatDB for cell-cell communication analysis
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
  # use all CellChatDB for cell-cell communication analysis
  # CellChatDB.use <- CellChatDB # simply use the default CellChatDB

  # set the used database in the object
  cellchat@DB <- CellChatDB.use

  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multiprocess", workers = 2) # do parallel
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                                distance.use = TRUE, interaction.length = 200, scale.distance = 0.01)

  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  cellchat <- netAnalysis_computeCentrality(cellchat)
  
  groupSize <- as.numeric(table(cellchat@idents))
  saveRDS(cellchat, here("ckd","spatial","cellchat",paste0(library_id,"_cellchat.rds")), compress=FALSE)
  return(cellchat)
  }, mc.cores=15)

cellchat.ls <- lapply(library_ids, function(library_id) {
  cellchat <- readRDS(here("ckd","spatial","cellchat",paste0(library_id,"_cellchat.rds")))
  cellchat <- netAnalysis_computeCentrality(cellchat)
  return(cellchat)
  })                      

# merge cellchat obj
cellchat <- mergeCellChat(cellchat.ls, add.names = library_ids, cell.prefix=TRUE)

plots <- here("ckd","spatial","plots","cellchat")
dir.create(here(plots))
pdf(here(plots,"netvisual_interaction_v1.pdf"))

  # extract the weights for each sample celltype-celltype interaction
  library(reshape2)
  library(ggplot2)
  df <- lapply(seq_along(cellchat.ls), function(index) {
    df <- cellchat@net[[index]][["weight"]]
    df <- melt(df)
    df$index <- index
    return(df)
    }) %>% bind_rows()

  # compare PT vs PT_VCAM1 interactions weights with other cell types
  res <- dplyr::filter(df, Var1 %in% c("PT","PT_INJ"))

  res %>% ggplot(aes(x=Var2, y=value, fill=Var1)) + geom_boxplot()

dev.off()

saveRDS(cellchat, here("ckd","spatial","step4_cellchat.rds"), compress=FALSE)

