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
library(EnsDb.Hsapiens.v86)
library(stringr)
library(tibble)


atac <- readRDS(here("ckd","atac_aggr_prep","step6_atac_loy.rds"))

plan("multicore", workers=6)
options(future.globals.maxSize = 2000 * 1024^2)

# define markers directory previously computed files with find_dar.R
markers <- here("ckd","atac_aggr_prep","markers")
dir.create(here(markers))    

# script for finding cell-specific deg and deg associated with LOY
Idents(atac) <- "celltype"
idents <- levels(atac)
# identify cell-specific chromvar activity for loy vs. NO loy
mark.ls <- lapply(idents, function(ident) {
  print(ident)
  ref <- atac@meta.data %>% dplyr::filter(akd_genotype == "XY", celltype == ident) %>% dplyr::select(barcode) %>% unlist
  test <- atac@meta.data %>% dplyr::filter(akd_genotype == "LOY", celltype == ident) %>% dplyr::select(barcode) %>% unlist
  
  mark <- tryCatch(FindMarkers(atac,
                              assay="peaks",
                              ident.2 = ref,
                              ident.1 = test,
                              logfc.threshold = 0
                              ), error=function(e) NULL)                                
})

# set up annotation for ClosestFeature 
gene.ranges <- genes(EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(gene.ranges),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(gene.ranges) <- ucsc.levels
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')  
  
mark.anno.ls <- lapply(seq(mark.ls), function(index) {
  print(index)
  dar <- mark.ls[[index]]
  tryCatch({cf <- ClosestFeature(atac, regions=rownames(dar), annotation=gene.ranges) %>%
                  dplyr::rename(peak = "query_region") %>%
                  dplyr::select(gene_name, distance, peak)
  dar <- dar %>% 
         rownames_to_column(var = "peak") %>%
         left_join(cf, by = "peak")
         }, error=function(e) NULL)        
  })
                        
# write to file                              
write.xlsx(mark.anno.ls, file = here(markers,"dar.celltype.loy_vs_xy.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)

# add age
age_meta <- read.xlsx(here("ckd","clinical_meta.xlsx")) %>% as.data.frame() %>% dplyr::select(library_id, age)
meta <- atac@meta.data %>% left_join(age_meta, by = "library_id")
rownames(meta) <- rownames(atac@meta.data)
atac@meta.data <- meta
  
# script for finding cell-specific deg and deg associated with LOY with age and depth adjustment
mark.ls <- lapply(idents, function(ident) {
  print(ident)
  ref <- atac@meta.data %>% dplyr::filter(akd_genotype == "XY", celltype == ident) %>% dplyr::select(barcode) %>% unlist
  test <- atac@meta.data %>% dplyr::filter(akd_genotype == "LOY", celltype == ident) %>% dplyr::select(barcode) %>% unlist
  
  mark <- tryCatch(FindMarkers(atac,
                              assay="peaks",
                              ident.2 = ref,
                              ident.1 = test,
                              logfc.threshold = 0,
                              test.use = "LR",
                              latent.vars = c("nCount_peaks","age"),
                              ), error=function(e) NULL)                                
})

# set up annotation for ClosestFeature 
gene.ranges <- genes(EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(gene.ranges),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(gene.ranges) <- ucsc.levels
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')  
  
mark.anno.ls <- lapply(seq(mark.ls), function(index) {
  print(index)
  dar <- mark.ls[[index]]
  tryCatch({cf <- ClosestFeature(atac, regions=rownames(dar), annotation=gene.ranges) %>%
                  dplyr::rename(peak = "query_region") %>%
                  dplyr::select(gene_name, distance, peak)
  dar <- dar %>% 
         rownames_to_column(var = "peak") %>%
         left_join(cf, by = "peak")
         }, error=function(e) NULL)        
  })
                        
# write to file                              
write.xlsx(mark.anno.ls, file = here(markers,"dar.celltype.loy_vs_xy.age_depth_adjust.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)  
