library(here)
library(ggplot2)

#################
# genes from https://www.mdanderson.org/documents/Labs/Wood-Laboratory/human-dna-repair-genes.html
dna_repair_genes <- c("AAAS","ADA","ADCY6","ADRM1","AK1","AK3","APRT",
                      "ARL6IP1","BCAM","BCAP31","BOLA2","BRF2","BRP44",
                      "CANT1","CCNO","CDA","CETN2","CLP1","CMPK2","COBRA1",
                      "COX17","CSTF3","DAD1","DCTN4","DDB1","DDB2","DFNA5",
                      "DGCR8","DGUOK","DUT","EDF1","EIF1B","EIF2C4","ELL",
                      "ERCC1","ERCC2","ERCC3","ERCC4","ERCC5","ERCC8",
                      "FEN1","GMPR2","GPX4","GTF2A2","GTF2B","GTF2F1",
                      "GTF2H1","GTF2H3","GTF2H5","GTF3C5","GUK1","HCLS1",
                      "HPRT1","IMPDH2","ITPA","LIG1","MPG","MRPL40","NCBP2",
                      "NFX1","NME1","NME3","NME4","NPR2","NT5C","NT5C3","NUDT21","NUDT9","PCNA","PDE4B","PDE6G","PNP","POLA1","POLA2","POLB","POLD1","POLD3","POLD4","POLE4","POLH","POLL","POLR1C","POLR1D","POLR2A","POLR2C","POLR2D","POLR2E","POLR2F","POLR2G","POLR2H","POLR2I","POLR2J","POLR2K","POLR3C","POLR3GL","POM121","PRIM1","RAD51","RAD52","RAE1","RALA","RBX1","RDBP","REV3L","RFC2","RFC3","RFC4","RFC5","RNMT","RPA2","RPA3","RRM2B","SAC3D1","SDCBP","SEC61A1","SF3A3","SMAD5","SNAPC4","SNAPC5","SRSF6","SSRP1","STX3","SUPT4H1","SUPT5H","SURF1","TAF10","TAF12","TAF13","TAF1C","TAF6","TAF9","TARBP2","TCEB3","TH1L","THOC4","TK2","TMED2","TP53","TSG101","TYMS","UMPS","UPF3B","USP11",
                      "VPS28","VPS37B","VPS37D","XPC","ZNF707","ZNRD1","ZWINT")

library(openxlsx)
file <- here("ckd","rna_aggr_prep","kpmp","markers","deg_loy_vs_xy.age_adjust.xlsx")
sheets <- getSheetNames(file)
deg <- lapply(sheets, function(sheet) {
  df <- read.xlsx(file, sheet)
  df$celltype <- sheet
  colnames(df)[1] <- "gene"
  return(df)
  }) %>% bind_rows()

gene_sel <- dna_repair_genes
toplot <- deg %>%
  dplyr::filter(gene %in% gene_sel) %>%
  dplyr::arrange(avg_log2FC) %>%
  dplyr::mutate(star = ifelse(p_val_adj < 0.05, "*","")) %>%
  dplyr::select(avg_log2FC, gene, celltype, star, p_val_adj) %>%
  dplyr::mutate(logFC = cut(avg_log2FC, breaks=c(max(avg_log2FC),0,-1,-2,-5,-10,-Inf)))

# for alphabetical sort
toplot$gene <- factor(toplot$gene, levels=gene_sel)
toplot$celltype <- factor(toplot$celltype, levels = rev(c("PT","PT_VCAM1","PT_PROM1","PT_MT","PEC","LOH","DCT","PC","IC","PODO","ENDO","FIB_VSMC_MC","LEUK")))

library(RColorBrewer)
palette <- rev(brewer.pal(7,"Blues"))
palette[6] <- "#FF7F7F"
plotA <- toplot %>%
  ggplot(aes(celltype, gene, fill=logFC)) +
  geom_tile() + 
  scale_fill_manual(values = palette) +
  xlab("") +
  ylab("") +
  ggtitle("A) DNA Repair") +
  coord_flip() +
  geom_text(aes(label = star)) +
  theme_bw() +
  theme(plot.title = element_text(face="bold"),
        axis.text = element_text(color="black"),
        axis.text.x = element_text(angle=90, vjust=1, hjust=1, size=8)) +
  guides(fill = guide_legend(reverse=TRUE))

############################################
# cell cycle genes
cell_cycle_genes <- c("CDK2","CDK4","CDK6","CDK7","CDKN1A",
                      "CDKN1B","STAG1","CDKN1C","CDKN2A",
                      "CDKN2B","CDKN2C","CDKN2D","ANAPC10",
                      "MAD2L2","STAG2","PTTG2","GADD45G",
                      "DBF4","YWHAQ","CHEK1","CHEK2","CREBBP",
                      "GADD45A","E2F1","E2F2","E2F3","E2F4",
                      "E2F5","EP300","ORC6","ORC3","CDC26",
                      "ABL1","ANAPC13","SMC1B","SFN","GSK3B",
                      "ANAPC2","ANAPC4","HDAC1","HDAC2","MAD2L1",
                      "SMAD2","SMAD3","SMAD4","MCM2","MCM3","MCM4",
                      "MCM5","MCM6","MCM7","MDM2","MYC","GADD45B",
                      "ATM","WEE2","ORC1","ORC2","ORC4","ORC5",
                      "PCNA","FZR1","ANAPC5","ANAPC7","ANAPC11",
                      "PLK1","ATR","PRKDC","RAD21","RB1","RBL1",
                      "RBL2","CCND1","ANAPC1","SKP1","SKP2",
                      "BUB1","BUB1B","TFDP1","TFDP2","TGFB1",
                      "TGFB2","TGFB3","TP53","TTK","SKP1P2",
                      "WEE1","YWHAB","YWHAE","YWHAG","YWHAH","YWHAZ",
                      "ZBTB17","SMC1A","CDC7","CDC45","MAD1L1","CUL1",
                      "CCNB3","CDC14B","CDC14A","CDC23","CDC16","CCNA2",
                      "CCNA1","CCNB1","CCND2","CCND3","CCNE1","CCNH",
                      "PKMYT1","SMC3","CCNB2","CCNE2","BUB3","PTTG1",
                      "ESPL1","CDK1","CDC6","CDC20","CDC25A","CDC25B",
                      "CDC25C","CDC27","RBX1")

gene_sel <- cell_cycle_genes
toplot <- deg %>%
  dplyr::filter(gene %in% gene_sel) %>%
  dplyr::arrange(avg_log2FC) %>%
  dplyr::mutate(star = ifelse(p_val_adj < 0.05, "*","")) %>%
  dplyr::select(avg_log2FC, gene, celltype, star, p_val_adj) %>%
  dplyr::mutate(logFC = cut(avg_log2FC, breaks=c(max(avg_log2FC),0,-1,-2,-5,-10,-Inf)))

# for alphabetical sort
toplot$gene <- factor(toplot$gene, levels=gene_sel)
toplot$celltype <- factor(toplot$celltype, levels = rev(c("PT","PT_VCAM1","PT_PROM1","PT_MT","PEC","LOH","DCT","PC","IC","PODO","ENDO","FIB_VSMC_MC","LEUK")))

library(RColorBrewer)
palette <- rev(brewer.pal(7,"Blues"))
palette[6] <- "#FF7F7F"
plotB <- toplot %>%
  ggplot(aes(celltype, gene, fill=logFC)) +
  geom_tile() + 
  scale_fill_manual(values = palette) +
  xlab("") +
  ylab("") +
  ggtitle("B) Cell Cycle") +
  coord_flip() +
  geom_text(aes(label = star)) +
  theme_bw() +
  theme(plot.title = element_text(face="bold"),
        axis.text = element_text(color="black"),
        axis.text.x = element_text(angle=90, vjust=1, hjust=1, size=8)) +
  guides(fill = guide_legend(reverse=TRUE))

lay <- rbind(c(1,1),
             c(2,2))
gs <- list(plotA, plotB)

pdf(here(plots,"sfigure10.pdf"), width = 7.5, height = 10)
 grid.arrange(grobs = gs, layout_matrix = lay)
dev.off()

