# use neighborhoods from step1 to identify neighborhood-specific degs

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

plan("multicore", workers=6)
options(future.globals.maxSize = 1000 * 1024^2)

sobj <- readRDS(here("ckd","spatial","spatial_step1b.rds"))

# define markers directory previously computed files with find_dar.R
markers <- here("ckd","spatial","markers")
dir.create(here(markers))    

# script for finding cell-specific deg and deg associated with LOY
Idents(sobj) <- "celltype"
idents <- levels(sobj)
DefaultAssay(sobj) <- "Spatial"

# identify cell-specific chromvar activity for loy vs. NO loy
mark.ls <- lapply(idents, function(ident) {
  print(ident)
  df <- FindMarkers(sobj, ident.1 = ident)                                
})  
                        
# write to file                              
write.xlsx(mark.ls, file = here(markers,"deg.celltype.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)
