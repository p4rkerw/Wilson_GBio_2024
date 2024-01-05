# analyze cell-specific deg for loy vs. xy in the kpmp dataset

# to run locally:
# SCRATCH1=/mnt/g/scratch
# docker run -it --rm \
# --workdir $HOME \
# -v /mnt/g/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# -v /mnt/g/cellranger_multi_counts:$HOME/cellranger_multi_counts \
# -v /mnt/g/kpmp:$HOME/kpmp \
# -v /mnt/g/diabneph:$HOME/diabneph \
# -v /mnt/g/reference:$HOME/reference \
# -v /mnt/g/scratch:$HOME/scratch \
# -v /mnt/g/ckd:$HOME/ckd \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/sctools:R4.1.3b R

library(Seurat)
library(openxlsx)
library(dplyr)
library(here)
library(future)

plan("multicore", workers=12)

# read in preprocessed obj
kpmp <- readRDS(here("ckd","rna_aggr_prep","kpmp","kpmp_prep.rds"))
Idents(kpmp) <- "celltype"
idents <- levels(kpmp)
levels(idents) <- levels(kpmp$celltype)

mark.ls <- lapply(idents, function(ident) {
  print(ident)
  # script for finding cell-specific deg and deg associated with LOY
  # select cells with loy in cell type
  loy <- kpmp@meta.data %>%
    dplyr::filter(rkd_genotype == "LOY", sex == "Male", celltype %in% ident) %>%
    dplyr::distinct(barcode)
  xy <- kpmp@meta.data %>%
    dplyr::filter(rkd_genotype == "XY", sex == "Male", celltype %in% ident) %>%
    dplyr::distinct(barcode)

  # find markers associated with loy in the proximal tubule
  markloy <- FindMarkers(kpmp, ident.1 = loy$barcode, ident.2 = xy$barcode)

}) 

# write to file
markers <- here("ckd","rna_aggr_prep","kpmp","markers")
dir.create(here(markers), recursive=TRUE)
write.xlsx(mark.ls, file = here(markers,"deg_loy_vs_xy.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)

# transform categorical age group into decade numeric variable
meta <- kpmp@meta.data %>%
  dplyr::mutate(age_decade = case_when(age == "20-29" ~ 3,
                                       age == "30-39" ~ 4,
                                       age == "40-49" ~ 5,
                                       age == "50-59" ~ 6,
                                       age == "60-69" ~ 5,
                                       age == "70-79" ~ 6))

rownames(meta) <- rownames(kpmp@meta.data)
kpmp@meta.data <- meta

# deg adjusted for age
mark.ls <- lapply(idents, function(ident) {
  print(ident)
  # script for finding cell-specific deg and deg associated with LOY
  # select cells with loy in cell type
  loy <- kpmp@meta.data %>%
    dplyr::filter(rkd_genotype == "LOY", sex == "Male", celltype %in% ident) %>%
    dplyr::distinct(barcode)
  xy <- kpmp@meta.data %>%
    dplyr::filter(rkd_genotype == "XY", sex == "Male", celltype %in% ident) %>%
    dplyr::distinct(barcode)

  # find markers associated with loy in the proximal tubule
  markloy <- FindMarkers(kpmp,
                         ident.1 = loy$barcode,
                         ident.2 = xy$barcode,
                         latent.vars = c("age_decade"),
                         test.use = "LR"
                         )

}) 

# write to file
markers <- here("ckd","rna_aggr_prep","kpmp","markers")
dir.create(here(markers), recursive=TRUE)
write.xlsx(mark.ls, file = here(markers,"deg_loy_vs_xy.age_adjust.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)

# deg loy vs XX adjusted for age 
mark.ls <- lapply(idents, function(ident) {
  print(ident)
  # script for finding cell-specific deg and deg associated with LOY
  # select cells with loy in cell type
  loy <- kpmp@meta.data %>%
    dplyr::filter(rkd_genotype == "LOY", sex == "Male", celltype %in% ident) %>%
    dplyr::distinct(barcode)
  xx <- kpmp@meta.data %>%
    dplyr::filter(sex == "Female", celltype %in% ident) %>%
    dplyr::distinct(barcode)

  # find markers associated with loy in the proximal tubule
  markloy <- FindMarkers(kpmp,
                         ident.1 = loy$barcode,
                         ident.2 = xx$barcode,
                         latent.vars = c("age_decade"),
                         test.use = "LR"
                         )

}) 

# write to file
markers <- here("ckd","rna_aggr_prep","kpmp","markers")
dir.create(here(markers), recursive=TRUE)
write.xlsx(mark.ls, file = here(markers,"deg_loy_vs_xx.age_adjust.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)

# deg loy vs XX
mark.ls <- lapply(idents, function(ident) {
  print(ident)
  # script for finding cell-specific deg and deg associated with LOY
  # select cells with loy in cell type
  loy <- kpmp@meta.data %>%
    dplyr::filter(rkd_genotype == "LOY", sex == "Male", celltype %in% ident) %>%
    dplyr::distinct(barcode)
  xx <- kpmp@meta.data %>%
    dplyr::filter(sex == "Female", celltype %in% ident) %>%
    dplyr::distinct(barcode)

  # find markers associated with loy in the proximal tubule
  markloy <- FindMarkers(kpmp,
                         ident.1 = loy$barcode,
                         ident.2 = xx$barcode
                         )

}) 

# write to file
markers <- here("ckd","rna_aggr_prep","kpmp","markers")
dir.create(here(markers), recursive=TRUE)
write.xlsx(mark.ls, file = here(markers,"deg_loy_vs_xx.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)


# deg loy vs no loy (xx or xy)
mark.ls <- lapply(idents, function(ident) {
  print(ident)
  # script for finding cell-specific deg and deg associated with LOY
  # select cells with loy in cell type
  loy <- kpmp@meta.data %>%
    dplyr::filter(rkd_genotype == "LOY", sex == "Male", celltype %in% ident) %>%
    dplyr::distinct(barcode)
  noloy <- kpmp@meta.data %>%
    dplyr::filter(rkd_genotype != "LOY", celltype %in% ident) %>%
    dplyr::distinct(barcode)

  # find markers associated with loy in the proximal tubule
  markloy <- FindMarkers(kpmp,
                         ident.1 = loy$barcode,
                         ident.2 = noloy$barcode
                         )

}) 

# write to file
markers <- here("ckd","rna_aggr_prep","kpmp","markers")
dir.create(here(markers), recursive=TRUE)
write.xlsx(mark.ls, file = here(markers,"deg_loy_vs_noloy.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)
    
