#########################
# figure 2a atac umap
library(Seurat)
library(here)
library(ggplot2)

# create a plot dir
plots <- here("ckd","figures")
dir.create(here(plots))

atac <- readRDS(here("ckd","atac_aggr_prep","step6_atac_loy.rds"))
atac <- RenameIdents(atac,
                     TCELL = "LEUK",
                     BCELL = "LEUK",
                     MONO = "LEUK")
plotA <- DimPlot(atac, label=T, raster=F)

plotA <- plotA + 
  theme_bw() +
  theme(plot.title = element_text(face="bold"), legend.pos = "none", axis.text = element_text(color="black"), panel.grid = element_blank()) +
  ggtitle("A)") +
  annotate(geom="text",
           x = -18,
           y = -16,
           label = paste0("Nuclei = ", nrow(atac@meta.data),"\n","Donors = 22 (12M,10F)"),
           hjust=0)

#########################
# figure 1b 1-dimensional kernel density estimate for chrY for ATAC assay
# with trough identification and shading
library(dplyr)
library(ggplot2)
library(here)
library(gridExtra)

atac_counts <- readRDS(here("ckd","atac_aggr_prep","atac_frags.rds"))

# ATAC kernel density trough
# find the trough of the chrY scaled fragment kernel density estimate to call LOY vs XY genotype
kd <- atac_counts %>%
  group_by(seqnames) %>%
  dplyr::filter(seqnames == "chrY", malesex == 1)
# find the maximum of the tallest peak in the density plot and its
# corresponding x coordinate (this is the median xcoord for XY cells)
density <- density(kd$scaled_atac_frags, n=10000)
global_peak_index <- which.max(density$y)
global_peak_xcoord <- density$x[global_peak_index]
# search up to x coord of the tallest peak for a trough
# note: be careful using the tallest peak without first inspecting the density plot
# In edge cases, the tallest peak may not correspond to the XY genotype
# (ie. LOY is more common than XY genotype). 
x <- density$x
y <- density$y[1:global_peak_index]
# calculate the difference in the y coord for successive steps in density function
diffy <- diff(y)
# find the sign of the change for the difference
signy <- sign(diffy)
# find points where the difference in the sign of change for successive steps in density function
# is negative 2 (ie. the first point was ascending and the second point was descending)
diffpeak <- diff(signy) == -2
# find the first two peak values
local_peak_index <- which(diffpeak)[1]
# set peak vertical line coords
peak1 <- x[local_peak_index]
peak2 <- global_peak_xcoord
# search between the peaks for a trough
trough_index <- which.min(density$y[density$x > peak1 & density$x < peak2]) + local_peak_index
trough <- x[trough_index]
kd <- kd %>%
  dplyr::mutate(genotype = ifelse(scaled_atac_frags < trough, "LOY", "XY"))

# create a shaded density plot using the trough of the chrY atac kernel density estimate
# as a threshold for calling LOY vs XY
p <- kd %>% 
      ggplot(aes(scaled_atac_frags)) + 
      geom_density(aes(y = ..scaled..)) +
      geom_vline(xintercept = trough, linetype = "dotted") + 
      xlab("Scaled coverage") +
      ylab("Density") +
      xlim(c(0,5))
d = ggplot_build(p)$data[[1]]
p = p + geom_area(data = subset(d, x < trough), aes(x=x,y=y), fill = "#00AFBB", alpha = 0.5)
p = p + geom_area(data = subset(d, x > trough), aes(x=x,y=y), fill = "darkgrey", alpha = 0.5)
res <- kd %>%
       summarize(genotype) %>%
       table() %>%
       as.data.frame() %>%
       dplyr::mutate(total = sum(Freq)) %>%
       dplyr::mutate(prop = Freq / total) %>%
       dplyr::filter(genotype == "LOY", seqnames == "chrY")
  
plotB <- p +
      annotate(geom="text",
               x = 3,
               y = 0.5,
               label = paste0("LOY = ", 100*round(res$prop, digits=2),"%")) +
      theme_bw() + 
      theme(plot.title = element_text(face="bold"), axis.text = element_text(color="black"), panel.grid = element_blank()) +
      ggtitle("B)") +
      xlab("Scaled ATAC coverage")

#######################################################
# figure histogram with proportion of each cell type with loy
# calculate aggregate proportion LOY per celltype
library(dplyr)
library(ggplot2)

meta <- readRDS(here("ckd","atac_aggr_prep","step6_meta.rds"))

toplot <- meta %>%
  dplyr::filter(akd_genotype %in% c("XY","LOY")) %>%
  dplyr::distinct(library_id, barcode, celltype, akd_genotype) %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarize(celltype, akd_genotype) %>%
  table() %>%
  as.data.frame() %>%
  group_by(celltype) %>%
  dplyr::mutate(total = sum(Freq)) %>%
  dplyr::mutate(loy_prop = Freq / total) %>%
  as.data.frame() %>%
  dplyr::filter(akd_genotype == "LOY")

# plot aggregate proportion for cell types 
plot_aggregate <- toplot %>% 
  ggplot(aes(x=celltype, y=loy_prop, fill=celltype)) +
  geom_histogram(stat = "identity") + 
  theme(axis.text.x = element_text(angle=90, vjust=1, hjust=1)) +
  xlab("") +
  ylab("Proportion LOY per cell type") +
  NoLegend()

# # plot mean proportion averaged across all male donors
# # group leukocytes together because there are very few of them
toplot <- meta %>%
  dplyr::mutate(celltype = case_when(celltype %in% c("TCELL","BCELL","MONO") ~ "LEUK",
                                     celltype %in% c("PCT","PST") ~ "PT",
                                     celltype %in% c("DCT1","DCT2") ~ "DCT",
                                     celltype %in% c("ATL","TAL1","TAL2","MD") ~ "LOH",
                                     celltype %in% c("ICA","ICB") ~ "IC",
                                     TRUE ~ as.character(celltype))) %>%
  na.omit() %>%
  dplyr::filter(akd_genotype %in% c("XY","LOY")) %>%
  group_by(library_id) %>%
  summarize(celltype, akd_genotype) %>%
  table() %>%
  as.data.frame() %>%
  group_by(library_id, celltype) %>%
  mutate(num_cells = sum(Freq)) %>%
  dplyr::filter(akd_genotype == "LOY") %>%
  mutate(lib_celltype_prop_loy = Freq / num_cells) %>%
  na.omit()
# relevel
toplot$celltype <- factor(toplot$celltype, levels = c("PT","PT_VCAM1","PT_PROM1","PEC","LOH","DCT","PC","IC","PODO","ENDO","FIB_VSMC_MC","LEUK"))

plotC <- toplot %>%
  group_by(celltype) %>%
  ggplot(aes(x=celltype, y=lib_celltype_prop_loy, fill=celltype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(w=0.3), aes(celltype), size=0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") +
  ylab("Proportion LOY") +
  xlab("") +
  ggtitle("C)") + 
  theme(plot.title = element_text(face="bold"), axis.text = element_text(color="black")) 

# # run a kruskal-wallis test for library proportion LOY by celltype
# kruskal.test(lib_celltype_prop_loy ~ celltype, data = toplot)

# pairwise tests
# wilcox.test(toplot[toplot$celltype == "PT",]$lib_celltype_prop_loy, toplot[toplot$celltype == "PT_VCAM1",]$lib_celltype_prop_loy)
# pairwise.t.test(toplot$lib_celltype_prop_loy, toplot$celltype, p.adj = "none")
# pairwise.t.test(toplot$lib_celltype_prop_loy, toplot$celltype, p.adj = "holm")

###############################################
# figure chromvar loy vs xy by cell type
library(here)
library(openxlsx)
library(dplyr)
library(tidyr)
# library(pheatmap)
# library(BuenColors)

xl <- here("ckd","atac_aggr_prep","markers","chromvar.celltype.loy_vs_xy.age_depth_adjust.xlsx")
sheets <- getSheetNames(xl)
chromvar <- lapply(sheets, function(sheet) {
  df <- read.xlsx(xl, sheet = sheet)
  df$celltype <- sheet
  colnames(df)[1] <- "motif"
  return(df)
}) %>% bind_rows

# take motifs that change activity (padj < 0.05) in at least one of the proximal tubule
# cell types and plot them in a heatmap with all of the proximal tubule cell types
# to see if there are any TF that change in the same direction for LOY in PT
celltype_sel <- c("PCT","PST","PT_VCAM1")
toplot <- chromvar %>%
  dplyr::filter(celltype %in% celltype_sel) %>%
  dplyr::group_by(motif) %>%
  dplyr::mutate(sum = sum(p_val_adj < 0.05)) %>%
  dplyr::filter(sum > 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::mutate(star = ifelse(p_val_adj < 0.05, "*","")) %>%
  dplyr::select(avg_log2FC, gene, celltype, star) %>%
  dplyr::mutate(logFC = cut(avg_log2FC, breaks=9, center=0))

# strip the (var.2) suffix from one of the motifs
toplot$gene <- gsub("\\(var.2\\)","",toplot$gene)

# for alphabetical sort
toplot$gene <- as.factor(toplot$gene)
plotD <- toplot %>%
  ggplot(aes(celltype, gene, fill=logFC)) +
  geom_tile() + 
  scale_fill_brewer(palette = "RdBu", direction = -1) +
  xlab("") +
  ylab("") +
  ggtitle("D)") +
  coord_flip() +
  geom_text(aes(label = star)) +
  theme_bw() +
  theme(plot.title = element_text(face="bold"), axis.text = element_text(color="black"),
        axis.text.x = element_text(angle=90, vjust=1, hjust=1)) +
  guides(fill = guide_legend(reverse=TRUE))

###################################################################
# figure forest plot of odds ratios for glmm of loy ~ celltype + (1|library_id)
library(lme4)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggplot2)
library(stringr)
theme_set(theme_sjplot())

df <- atac@meta.data %>%
  dplyr::filter(akd_genotype %in% c("XY","LOY")) %>%
  dplyr::mutate(loy = ifelse(akd_genotype == "LOY",1,0)) %>%
  dplyr::select(library_id, loy, celltype, nCount_peaks)

df <- df %>%
  dplyr::mutate(celltype = case_when(celltype %in% c("TCELL","BCELL","MONO") ~ "LEUK",
                                     celltype %in% c("PCT","PST") ~ "PT",
                                     celltype %in% c("DCT1","DCT2") ~ "DCT",
                                     celltype %in% c("MD","ATL","TAL1","TAL2") ~ "LOH",
                                     celltype %in% c("ICA","ICB") ~ "IC",
                                     TRUE ~ as.character(celltype))) 
# relevel
df$celltype <- factor(df$celltype, levels = c("PT","PT_VCAM1","PT_PROM1","PEC","LOH","DCT","PC","IC","PODO","ENDO","FIB_VSMC_MC","LEUK"))


# glm for effect of celltype on loy
mod1 <- glm(loy ~ celltype, family = binomial(link = "logit"), data = df)
summary(mod1)
# confint(mod1)

# glmm with mixed effect for sample for celltype on loy
mod2 <- lme4::glmer(loy ~ celltype + (1|library_id),
                    family = binomial(link = "logit"),
                    data = df)
summary(mod2)

plotE <- plot_model(mod2, p.adjust=TRUE) 

# swap out yaxis labels
newlabels <- str_replace(plotE$data$term,"celltype","")
newlabels <- factor(newlabels, levels=newlabels)

plotE <- plotE + scale_y_continuous(limits = c(0.01, 3), breaks = c(0.1, 0.5, 1, 3)) +
  ggtitle("E)") +
  geom_hline(yintercept = 1, linetype = "dotted") + 
  theme_bw() +
  theme(plot.title = element_text(face="bold"), axis.text = element_text(color="black")) + 
  ylab("Odds Ratios for LOY relative to PT") + 
  scale_x_discrete(labels=rev(newlabels))

# layout figure
lay <- rbind(c(1,1,2),
             c(1,1,3),
             c(4,4,4),
             c(5,5,NA))
gs <- list(plotA, plotB, plotC, plotD, plotE)

pdf(here(plots,"figure2.pdf"),width = 7.5, height = 10)
 grid.arrange(grobs = gs, layout_matrix = lay, heights = c(0.75,1.25,1,1))
dev.off()

