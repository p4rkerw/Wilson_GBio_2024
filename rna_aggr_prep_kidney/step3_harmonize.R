# harmonize kpmp by projecting onto multiome atlas and doing joint deg analysis

# to run locally:
# SCRATCH1=/mnt/g/scratch
# docker run -it --rm \
# --workdir $HOME \
# -v /mnt/g/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# -v /mnt/g/cellranger_multi_counts:$HOME/cellranger_multi_counts \
# -v /mnt/g/kpmp:$HOME/kpmp \
# -v /mnt/g/reference:$HOME/reference \
# -v /mnt/g/scratch:$HOME/scratch \
# -v /mnt/g/ckd:$HOME/ckd \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/sctools:R4.1.3 R

library(Seurat)
library(Signac)
library(data.table)
library(dplyr)
library(here)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)
library(future)
library(openxlsx)

kpmp <- readRDS(here("ckd","rna_aggr_prep","kpmp","kpmp_prep.rds"))

multi <- readRDS(here("ckd","multi_aggr_prep","step3_multi_loy.rds"))
DefaultAssay(multi) <- 'RNA'
DefaultAssay(kpmp) <- 'RNA'

multi[["ATAC"]] <- NULL
multi[["SCT"]] <- NULL

multi <- FindVariableFeatures(multi)
multi <- ScaleData(multi)
multi <- RunPCA(multi)
Idents(multi) <- "celltype"  

# project multiome cell type labels onto kpmp dataset
transfer.anchors <- FindTransferAnchors(
  reference = multi,
  query = kpmp,
  reduction = 'pcaproject',
  normalization.method = 'LogNormalize',
  reference.assay = 'RNA',
  query.assay = 'RNA'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = multi$celltype,
  weight.reduction = kpmp[['pca']],
  dims = 1:30
)

kpmp$old.celltype <- kpmp$celltype
kpmp$celltype <- predicted.labels$predicted.id

# merge rna datasets for joint deg analysis
merged_obj <- merge(multi, kpmp)
Idents(merged_obj) <- "celltype"
levels(merged_obj) <- c("PCT","PST","PT_VCAM1","PEC",
          "TAL1","TAL2","DCT1","DCT2","PC",
          "ICA","ICB","PODO","ENDO","FIB_VSMC_MC",
          "TCELL","BCELL","MONO")

merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)

meta <- merged_obj@meta.data %>%
  dplyr::mutate(genotype = case_when(rkd_genotype == "LOY" ~ "LOY",
                                     rkd_genotype == "XY" ~ "XY",
                                     gmm_genotype == "LOY" ~ "LOY",
                                     gmm_genotype == "XY" ~ "XY"))
merged_obj$genotype <- meta$genotype

# saveRDS(merged_obj, here("ckd","rna_aggr_prep","kpmp","step3_harmonize.rds"), compress=FALSE)

# parallelize
plan("multicore", workers=12)

# joint deg analysis for loy by celltype 
idents <- levels(merged_obj)
mark.ls <- lapply(idents, function(ident) {
  print(ident)
  # script for finding cell-specific deg and deg associated with LOY
  # select cells with loy in cell type
  loy <- merged_obj@meta.data %>%
    dplyr::filter(genotype == "LOY", celltype %in% ident) %>%
    dplyr::distinct(barcode)
  xy <- merged_obj@meta.data %>%
    dplyr::filter(genotype == "XY", celltype %in% ident) %>%
    dplyr::distinct(barcode)

  # find markers associated with loy in the proximal tubule
  markloy <- FindMarkers(merged_obj, ident.1 = loy$barcode, ident.2 = xy$barcode)

}) 

# write to file
markers <- here("ckd","rna_aggr_prep","kpmp","markers")
dir.create(here(markers), recursive=TRUE)
write.xlsx(mark.ls, file = here(markers,"joint_deg_loy_vs_xy.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)

# transform categorical age group into decade numeric variable
meta <- merged_obj@meta.data %>%
  dplyr::mutate(age_decade = case_when(age == "20-29" ~ 3,
                                       age == "30-39" ~ 4,
                                       age == "40-49" ~ 5,
                                       age == "50-59" ~ 6,
                                       age == "60-69" ~ 5,
                                       age == "70-79" ~ 6,
                                       age < 50 && age >= 40 ~ 5,
                                       age < 60 && age >= 50 ~ 6,
                                       age < 70 && age >= 60 ~ 7,
                                       age < 80 && age >= 70 ~ 8,
                                       age < 90 && age >= 80 ~ 9))
merged_obj$age_decade <- meta$age_decade

# deg adjusted for age
mark.ls <- lapply(idents, function(ident) {
  print(ident)
  # script for finding cell-specific deg and deg associated with LOY
  # select cells with loy in cell type
  loy <- merged_obj@meta.data %>%
    dplyr::filter(genotype == "LOY", celltype %in% ident) %>%
    dplyr::distinct(barcode)
  xy <- merged_obj@meta.data %>%
    dplyr::filter(genotype == "XY", celltype %in% ident) %>%
    dplyr::distinct(barcode)

  # find markers associated with loy in the proximal tubule
  markloy <- FindMarkers(merged_obj,
                         ident.1 = loy$barcode,
                         ident.2 = xy$barcode,
                         latent.vars = c("age_decade"),
                         test.use = "LR"
                         )

}) 

# write to file
markers <- here("ckd","rna_aggr_prep","kpmp","markers")
dir.create(here(markers), recursive=TRUE)
write.xlsx(mark.ls, file = here(markers,"joint_deg_loy_vs_xy.age_adjust.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)


                         





