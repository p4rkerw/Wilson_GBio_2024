# figure qc snATACseq
# supplemental figure qc metrics and lineage markers 

library(Seurat)
library(Signac)
library(here)
library(ggplot2)

# create a plot dir
plots <- here("ckd","figures")
dir.create(here(plots))

atac <- readRDS(here("ckd","atac_aggr_prep","step6_atac_loy.rds"))

plotA <- VlnPlot(atac, features = "nCount_peaks", group.by = "library_id", pt.size = 0) +
 xlab("") +
 ggtitle("A) nCount_peaks") +
 theme_bw() +
 theme(plot.title = element_text(face="bold"), 
       legend.pos = "none", 
       axis.text = element_text(color="black"), 
       axis.text.x = element_text (angle = 90, hjust = 1))

plotB <- VlnPlot(atac, features = "nFeature_peaks", group.by = "library_id", pt.size = 0) +
 xlab("") +
 ggtitle("B) nFeature_peaks") +
 theme_bw() +
 theme(plot.title = element_text(face="bold"), 
       legend.pos = "none", 
       axis.text = element_text(color="black"), 
       axis.text.x = element_text (angle = 90, hjust = 1))
       
## dotplot for lineage markers
marker.genes <- c("CUBN","HAVCR1","SLC5A1","SLC5A2", # PT
                  "VCAM1","PROM1","TNIK","TNFSF10",  # PT-VCAM1+ markers
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
                  "PLVAP", # PTC and AVR https://www.nature.com/articles/s41467-019-12872-5
                  "NOS1", # MD
                  "ITGA8","PDGFRB","MEIS2","PIEZO2","REN", # MES and JGA
                  "ACTA2","CALD1", # FIB
                  "PROX1","FLT4","PDPN", # Lymphatics
                  "PTPRC","CD3E","CD2","MS4A1", # Lymphocytes
                  "FCGR3A","CD14","CSF1R") # Monocyte / Macrophage


DefaultAssay(atac) <- "RNA"
plotC <- DotPlot(atac, features=marker.genes) + 
 ggtitle("C) ATAC Gene Activity Lineage Markers") +
 theme_bw() +
 theme(plot.title = element_text(face="bold"), 
       legend.pos = "none", 
       axis.text = element_text(color="black"), 
       axis.text.x = element_text (angle = 90, hjust = 1)) +
  xlab("") +
  ylab("")
  
  
lay <- rbind(c(1,2),
             c(3,3),
             c(3,3))
gs <- list(plotA, plotB, plotC)

pdf(here(plots,"sfigure5.pdf"), width = 7.5, height = 10)
 grid.arrange(grobs = gs, layout_matrix = lay)
dev.off()
