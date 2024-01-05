# this script will analyze the aggregated snATACseq object to identify cell-specific TF that change in activity for LOY vs no LOY

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

library(Signac) # 
library(Seurat) # 
library(here) # 
library(future) # 
library(openxlsx) # 
library(dplyr) # 
set.seed(1234)

###########################################################################
# annotate chromvar motif id with gene names
AnnoChromvar <- function(df, motifs) {
  motif_lookup <- rownames(df)
  motif_names <- sapply(motif_lookup, function(x) motifs[[x]])
  return(cbind(df, gene = motif_names))
}
##############################################
# read in aggregated snATAC object
atac_aggr_prep <- here("ckd","atac_aggr_prep")
atacAggr <- readRDS(here(atac_aggr_prep,"step6_atac_loy.rds"))
Idents(atacAggr) <- "celltype"
DefaultAssay(atacAggr) <- "chromvar"
                                              
# define markers directory previously computed files with find_dar.R
markers <- here("ckd","atac_aggr_prep","markers")
dir.create(here(markers))                        

# assign idents
idents <- levels(atacAggr)
motif.names <- GetMotifData(object = atacAggr, assay = "peaks", slot = "motif.names")
                        
# identify cell-specific chromvar activity for loy vs. XY
mark.ls <- lapply(idents, function(ident) {
  print(ident)
  ref <- atacAggr@meta.data %>% dplyr::filter(akd_genotype == "XY", celltype == ident) %>% dplyr::select(barcode) %>% unlist
  test <- atacAggr@meta.data %>% dplyr::filter(akd_genotype == "LOY", celltype == ident) %>% dplyr::select(barcode) %>% unlist
  mark <- FindMarkers(atacAggr, assay = "chromvar", ident.2 = ref, ident.1 = test, test.use = 'LR', latent.vars = "peak_region_fragments", logfc.threshold = 0)                        
  mark <- AnnoChromvar(mark, motif.names)
})  
                        
# write to file                              
write.xlsx(mark.ls, file = here(markers,"chromvar.celltype.loy_vs_xy.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)

# add age to metadata                        
age_meta <- read.xlsx(here("ckd","clinical_meta.xlsx")) %>% as.data.frame() %>% dplyr::select(library_id, age)
meta <- atacAggr@meta.data %>% left_join(age_meta, by = "library_id")
rownames(meta) <- rownames(atacAggr@meta.data)
atacAggr@meta.data <- meta
                        
# identify cell-specific chromvar activity for loy vs. XY with age and depth adjustment
mark.ls <- lapply(idents, function(ident) {
  print(ident)
  ref <- atacAggr@meta.data %>% dplyr::filter(akd_genotype == "XY", celltype == ident) %>% dplyr::select(barcode) %>% unlist
  test <- atacAggr@meta.data %>% dplyr::filter(akd_genotype == "LOY", celltype == ident) %>% dplyr::select(barcode) %>% unlist
  mark <- FindMarkers(atacAggr,
                      assay = "chromvar",
                      ident.2 = ref,
                      ident.1 = test,
                      test.use = 'LR',
                      latent.vars = c("age","peak_region_fragments"),
                      logfc.threshold = 0)                        
  mark <- AnnoChromvar(mark, motif.names)
})  
                        
# write to file                              
write.xlsx(mark.ls, file = here(markers,"chromvar.celltype.loy_vs_xy.age_depth_adjust.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)
                                               
# identify cell-specific chromvar activity 
mark.ls <- lapply(idents, function(ident) {
  print(ident)
  mark <- FindMarkers(atacAggr, assay = "chromvar", ident.1 = ident, test.use = 'LR', latent.vars = "peak_region_fragments", logfc.threshold = 0.10, only.pos=TRUE)                        
  mark <- AnnoChromvar(mark, motif.names)
})  
                        
# write to file                              
write.xlsx(mark.ls, file = here(markers,"chromvar.celltype.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)
                       
