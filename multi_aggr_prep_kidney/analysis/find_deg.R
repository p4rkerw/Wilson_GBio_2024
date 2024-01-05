#!/usr/bin/env Rscript
# to run locally:
# SCRATCH1=/mnt/g/scratch
# docker run -it --rm \
# --workdir $HOME \
# -v /mnt/g/ckd:$HOME/ckd \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/sctools:R4.1.3 R
library(Seurat)
library(dplyr)
library(openxlsx)
library(here)
library(future)

plan("multicore", workers=12)
options(future.globals.maxSize = 1000 * 1024^2)

multi <- readRDS(here("ckd","multi_aggr_prep","step3_multi_loy.rds"))

# define markers directory previously computed files with find_dar.R
markers <- here("ckd","multi_aggr_prep","markers")
dir.create(here(markers))    

# script for finding cell-specific deg and deg associated with LOY
Idents(multi) <- "celltype"
idents <- levels(multi)
# identify cell-specific chromvar activity for loy vs. NO loy
mark.ls <- lapply(idents, function(ident) {
  print(ident)
  ref <- multi@meta.data %>% dplyr::filter(gmm_genotype == "XY", celltype == ident) %>% dplyr::select(barcode) %>% unlist
  test <- multi@meta.data %>% dplyr::filter(gmm_genotype == "LOY", celltype == ident) %>% dplyr::select(barcode) %>% unlist
  
  mark <- tryCatch(FindMarkers(multi,
                              assay="RNA",
                              ident.2 = ref,
                              ident.1 = test,
                              logfc.threshold = 0
                              ), error=function(e) NULL)                                
})  
                        
# write to file                              
write.xlsx(mark.ls, file = here(markers,"deg.celltype.loy_vs_XY.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)

# identify cell-specific deg for loy vs. XY with age adjustment
mark.ls <- lapply(idents, function(ident) {
  print(ident)
  ref <- multi@meta.data %>% dplyr::filter(gmm_genotype == "XY", celltype == ident) %>% dplyr::select(barcode) %>% unlist
  test <- multi@meta.data %>% dplyr::filter(gmm_genotype == "LOY", celltype == ident) %>% dplyr::select(barcode) %>% unlist
  
  mark <- tryCatch(FindMarkers(multi,
                              assay="RNA",
                              ident.2 = ref,
                              ident.1 = test,
                              logfc.threshold = 0,
                              latent.vars = c("age"),
                              test.use = "LR",
                              ), error=function(e) NULL)                                
})  
                        
# write to file                              
write.xlsx(mark.ls, file = here(markers,"deg.celltype.loy_vs_xy.age_adjust.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)  

# identify cell-specific deg for loy vs. XX 
mark.ls <- lapply(idents, function(ident) {
  print(ident)
  ref <- multi@meta.data %>% dplyr::filter(malesex == 0, celltype == ident) %>% dplyr::select(barcode) %>% unlist
  test <- multi@meta.data %>% dplyr::filter(gmm_genotype == "LOY", celltype == ident) %>% dplyr::select(barcode) %>% unlist
  
  mark <- tryCatch(FindMarkers(multi,
                              assay="RNA",
                              ident.2 = ref,
                              ident.1 = test,
                              logfc.threshold = 0
                              ), error=function(e) NULL)                                
})  
                        
# write to file                              
write.xlsx(mark.ls, file = here(markers,"deg.celltype.loy_vs_xx.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)
  
  
# identify cell-specific deg for loy vs. XX with age adjustment
mark.ls <- lapply(idents, function(ident) {
  print(ident)
  ref <- multi@meta.data %>% dplyr::filter(malesex == 0, celltype == ident) %>% dplyr::select(barcode) %>% unlist
  test <- multi@meta.data %>% dplyr::filter(gmm_genotype == "LOY", celltype == ident) %>% dplyr::select(barcode) %>% unlist
  
  mark <- tryCatch(FindMarkers(multi,
                              assay="RNA",
                              ident.2 = ref,
                              ident.1 = test,
                              logfc.threshold = 0,
                              latent.vars = c("age"),
                              test.use = "LR",
                              ), error=function(e) NULL)                                
})  
                        
# write to file                              
write.xlsx(mark.ls, file = here(markers,"deg.celltype.loy_vs_xx.age_adjust.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)

# identify cell-specific deg for XY vs. XX with age adjustment
mark.ls <- lapply(idents, function(ident) {
  print(ident)
  ref <- multi@meta.data %>% dplyr::filter(malesex == 0, celltype == ident) %>% dplyr::select(barcode) %>% unlist
  test <- multi@meta.data %>% dplyr::filter(gmm_genotype == "XY", celltype == ident) %>% dplyr::select(barcode) %>% unlist
  
  mark <- tryCatch(FindMarkers(multi,
                              assay="RNA",
                              ident.2 = ref,
                              ident.1 = test,
                              logfc.threshold = 0,
                              latent.vars = c("age"),
                              test.use = "LR",
                              ), error=function(e) NULL)                                
})  
                        
# write to file                              
write.xlsx(mark.ls, file = here(markers,"deg.celltype.xy_vs_xx.age_adjust.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)

# identify cell-specific deg for loy vs. no loy (xx or xy)
mark.ls <- lapply(idents, function(ident) {
  print(ident)
  ref <- multi@meta.data %>% dplyr::filter(gmm_genotype != "LOY", celltype == ident) %>% dplyr::select(barcode) %>% unlist
  test <- multi@meta.data %>% dplyr::filter(gmm_genotype == "LOY", celltype == ident) %>% dplyr::select(barcode) %>% unlist
  
  mark <- tryCatch(FindMarkers(multi,
                              assay="RNA",
                              ident.2 = ref,
                              ident.1 = test,
                              logfc.threshold = 0
                              ), error=function(e) NULL)                                
})  
                        
# write to file                              
write.xlsx(mark.ls, file = here(markers,"deg.celltype.loy_vs_noloy.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)  

