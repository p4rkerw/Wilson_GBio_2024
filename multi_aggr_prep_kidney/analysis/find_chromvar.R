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
multi_aggr_prep <- here("ckd","multi_aggr_prep")
multi <- readRDS(here(multi_aggr_prep,"step8_chromVAR.rds"))
Idents(multi) <- "celltype"
DefaultAssay(multi) <- "chromvar"
                                              
# define markers directory previously computed files with find_dar.R
markers <- here("ckd","multi_aggr_prep","markers")
dir.create(here(markers))                        

# assign idents
idents <- levels(multi)
motif.names <- GetMotifData(object = multi, assay = "ATAC", slot = "motif.names")
                        
# identify cell-specific chromvar activity for loy vs. NO loy
mark.ls <- lapply(idents, function(ident) {
  print(ident)
  ref <- multi@meta.data %>% dplyr::filter(gmm_genotype == "XY", celltype == ident) %>% dplyr::select(barcode) %>% unlist
  test <- multi@meta.data %>% dplyr::filter(gmm_genotype == "LOY", celltype == ident) %>% dplyr::select(barcode) %>% unlist
  
  tryCatch({
    mark <- FindMarkers(multi, assay = "chromvar", ident.2 = ref, ident.1 = test, test.use = 'LR', latent.vars = "nCount_ATAC", logfc.threshold = 0)                        
    mark <- AnnoChromvar(mark, motif.names)
    }, error=function(e) return(NULL))
})  
                        
# write to file                              
write.xlsx(mark.ls, file = here(markers,"chromvar.celltype.loy_vs_xy.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)
                        
# identify cell-specific chromvar activity 
mark.ls <- lapply(idents, function(ident) {
  print(ident)
  tryCatch({
    mark <- FindMarkers(multi, assay = "chromvar", ident.1 = ident, test.use = 'LR', latent.vars = "nCount_ATAC", logfc.threshold = 0.10, only.pos=TRUE)                        
    mark <- AnnoChromvar(mark, motif.names)
    }, error=function(e) return(NULL))     
})  
                        
# write to file                              
write.xlsx(mark.ls, file = here(markers,"chromvar.celltype.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)
