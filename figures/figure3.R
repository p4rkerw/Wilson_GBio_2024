# kpmp snRNA aggregated analysis loy
#########################
# figure umap
library(Seurat)
library(here)
library(ggplot2)

# create a plot dir
plots <- here("ckd","figures")
dir.create(here(plots))

kpmp <- readRDS(here("ckd","rna_aggr_prep","kpmp","kpmp_prep.rds"))
plotA <- DimPlot(kpmp, label=T, repel=TRUE, raster=F)

plotA <- plotA + 
  theme_bw() +
  theme(plot.title = element_text(face="bold"), legend.pos = "none", axis.text = element_text(color="black"), panel.grid = element_blank()) +
  ggtitle("A)") +
  annotate(geom="text",
           x = 3,
           y = -16,
           label = paste0("Nuclei = ", nrow(kpmp@meta.data),"\n","Donors = 37 (21M,16F)"),
           hjust=0)

#########################
# figure 1b 1-dimensional kernel density estimate for chrY for RNA assay
# with trough identification and shading
library(dplyr)
library(ggplot2)
library(here)
library(gridExtra)

rna_counts <- readRDS(here("ckd","rna_aggr_prep","kpmp","rna_counts.rds"))
addmeta <- kpmp@meta.data %>%
  dplyr::filter(nCount_RNA > 1000) %>%
  dplyr::select(barcode, nCount_RNA)

# kernel density trough
# find the trough of the chrY scaled fragment kernel density estimate to call LOY vs XY genotype
kd <- rna_counts %>%
  left_join(addmeta, by = "barcode") %>%
  dplyr::filter(sex == "Male") %>%
  dplyr::filter(nCount_RNA > 1000)
# find the x coordinate for the peak corresponding to XY (this is the median xcoord)
density <- density(kd$scaled_rna_counts, n=10000)

xy_peak_index <- which(density$x > median(kd$scaled_rna_counts))[1]
xy_peak_xcoord <- density$x[xy_peak_index]

# search up to xy peak for a trough
# note: be careful and first inspect the density plot
x <- density$x
y <- density$y[1:xy_peak_index]
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
peak2 <- xy_peak_xcoord
# search between the peaks for a trough
trough_index <- which.min(density$y[density$x > peak1 & density$x < peak2]) + local_peak_index
trough <- x[trough_index]
kd <- kd %>%
  dplyr::mutate(genotype = ifelse(scaled_rna_counts < trough, "LOY", "XY"))

# create a shaded density plot using the trough of the chrY atac kernel density estimate
# as a threshold for calling LOY vs XY
p <- kd %>% 
      ggplot(aes(scaled_rna_counts)) + 
      geom_density(aes(y = ..scaled..)) +
      geom_vline(xintercept = trough, linetype = "dotted") + 
      xlab("Scaled coverage") +
      ylab("Density") +
      xlim(c(0,5))
d = ggplot_build(p)$data[[1]]
p = p + geom_area(data = subset(d, x < trough), aes(x=x,y=y), fill = "#00AFBB", alpha = 0.5)
p = p + geom_area(data = subset(d, x > trough), aes(x=x,y=y), fill = "darkgrey", alpha = 0.5)
res <- kd %>%
       ungroup() %>%
       summarize(genotype) %>%
       table() %>%
       as.data.frame() %>%
       dplyr::mutate(total = sum(Freq)) %>%
       dplyr::mutate(prop = Freq / total) %>%
       dplyr::filter(genotype == "LOY")
  
plotB <- p +
      annotate(geom="text",
               x = 3,
               y = 0.5,
               label = paste0("LOY = ", 100*round(res$prop, digits=2),"%")) +
      theme_bw() + 
      theme(plot.title = element_text(face="bold"), axis.text = element_text(color="black"), panel.grid = element_blank()) +
      ggtitle("B)") +
      xlab("Scaled RNA coverage")
      
#######################################################
# figure histogram with proportion of each cell type with loy
# calculate aggregate proportion LOY per celltype
library(dplyr)
library(ggplot2)

kpmp <- readRDS(here("ckd","rna_aggr_prep","kpmp","kpmp_prep.rds"))
meta <- kpmp@meta.data

# filter by min nCount_RNA as in the multiome analysis
meta <- meta %>% dplyr::filter(nCount_RNA > 1000)

toplot <- meta %>%
  dplyr::filter(rkd_genotype %in% c("XY","LOY")) %>%
  dplyr::distinct(library_id, barcode, celltype, rkd_genotype) %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarize(celltype, rkd_genotype) %>%
  table() %>%
  as.data.frame() %>%
  group_by(celltype) %>%
  dplyr::mutate(total = sum(Freq)) %>%
  dplyr::mutate(loy_prop = Freq / total) %>%
  as.data.frame() %>%
  dplyr::filter(rkd_genotype == "LOY")

# plot aggregate proportion for cell types with at least 50 cells
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
  dplyr::select(rkd_genotype, library_id, celltype) %>%
  dplyr::filter(rkd_genotype %in% c("XY","LOY")) %>%
  group_by(library_id) %>%
  summarize(celltype, rkd_genotype) %>%
  table() %>%
  as.data.frame() %>%
  group_by(library_id, celltype) %>%
  mutate(num_cells = sum(Freq)) %>%
  dplyr::filter(rkd_genotype == "LOY") %>%
  mutate(lib_celltype_prop_loy = Freq / num_cells) %>%
  na.omit()
# relevel
toplot$celltype <- factor(toplot$celltype, levels = c("PT","PT_VCAM1","PT_PROM1","PT_MT","PEC","LOH","DCT","PC","IC","PODO","ENDO","FIB_VSMC_MC","LEUK"))

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
###############################
library(openxlsx)
library(RColorBrewer)
file <- here("ckd","rna_aggr_prep","kpmp","markers","deg_loy_vs_xy.age_adjust.xlsx")
sheets <- getSheetNames(file)
deg <- lapply(sheets, function(sheet) {
  df <- read.xlsx(file, sheet)
  df$celltype <- sheet
  colnames(df)[1] <- "gene"
  return(df)
}) %>% bind_rows()

# take deg for enriched pathways
# and plot them in a heatmap with selected cell types
# to see if there are deg that change in the same direction for LOY 
dsb_genes <- c("SPIRE1","EPC1","SHLD2","ACTB","SPIDR","SETD2","UBR5")
g1m_genes <- c("WAC","JADE1","RDX","APPL2","ACTB","ARID1B","BCL2","EGFR","CRADD")
gene_sel <- c(dsb_genes, g1m_genes) %>% unique()

toplot <- deg %>%
  dplyr::mutate(gene_group = ifelse(gene %in% dsb_genes, "DSB Repair", NA)) %>%
  dplyr::mutate(gene_group = ifelse(gene %in% g1m_genes, "G1/S Transition", gene_group))

toplot <- toplot %>%
  dplyr::filter(gene %in% gene_sel) %>%
  dplyr::arrange(avg_log2FC) %>%
  dplyr::mutate(star = ifelse(p_val_adj < 0.05, "*","")) %>%
  dplyr::select(avg_log2FC, gene, celltype, star, p_val_adj, gene_group) %>%
  dplyr::mutate(logFC = cut(avg_log2FC, breaks=c(min(avg_log2FC),-5,-2,-0.5,0,0.5,2,5,max(avg_log2FC))))

# for alphabetical sort
toplot$gene <- factor(toplot$gene, levels=gene_sel)
toplot$celltype <- factor(toplot$celltype, levels = rev(levels(kpmp$celltype)))

palette <- rev(brewer.pal(9,"RdBu"))
plotD <- toplot %>%
  filter(celltype %in% c("PT","ENDO","FIB_VSMC_MC","LEUK")) %>%
  ggplot(aes(celltype, gene, fill=logFC)) +
  geom_tile() + 
  scale_fill_manual(values = palette) +
  xlab("") +
  ylab("") +
  ggtitle("D)") +
  coord_flip() +
  geom_text(aes(label = star)) +
  theme_bw() +
  theme(plot.title = element_text(face="bold"), axis.text = element_text(color="black"),
        axis.text.x = element_text(angle=90, vjust=1, hjust=1)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  facet_grid(~gene_group, scales = "free")
##################
library(ggpubr)
# figure loy in ckd vs. non-ckd kpmp samples
meta <- kpmp@meta.data %>%
  dplyr::filter(rkd_genotype %in% c("XY","LOY"), nCount_RNA > 1000) %>%
  dplyr::select(barcode, rkd_genotype, age, specimen_id, celltype) %>%
  dplyr::rename(genotype = rkd_genotype, age_group = age, library_id = specimen_id) 

clinical_meta <- read.csv("G:/kpmp/72ec8f42-750f-4e1b-bda9-7300d4de1365_20221208_OpenAccessClinicalData.csv")
addmeta <- data.frame(library_id = clinical_meta$Participant.ID, Tissue.Type = clinical_meta$Tissue.Type)

meta <- meta %>%
  left_join(addmeta, by = c("library_id")) %>%
  dplyr::mutate(group = case_when(Tissue.Type == "AKI" ~ "Control or AKI",
                                  Tissue.Type == "Healthy Reference" ~ "Control or AKI",
                                  TRUE ~ as.character(Tissue.Type))) 

toplot <- meta %>%
  dplyr::filter(celltype %in% c("PT","PT_VCAM1","PT_MT")) %>%
  dplyr::group_by(library_id) %>%
  dplyr::summarize(group, genotype) %>%
  table() %>%
  as.data.frame() %>%
  group_by(library_id) %>%
  dplyr::mutate(total = sum(Freq)) %>%
  dplyr::mutate(loy_prop = Freq / total) %>%
  as.data.frame() %>%
  dplyr::filter(genotype == "LOY", Freq != 0)

toplot$group <- factor(toplot$group, levels = c("Control or AKI","CKD"))

plotE <- toplot %>%
  group_by(group) %>%
  ggplot(aes(x=group, y=loy_prop, fill=group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(w=0.3), aes(group), size=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Proportion LOY") +
  xlab("") +
  ggtitle("E)") + 
  theme_bw() +
  theme(plot.title = element_text(face="bold"),
        axis.text = element_text(color="black"),
        legend.pos = "none") +
  ggpubr::stat_compare_means(method = 'wilcox.test',
                             label = 'p.signif',
                             show.legend = F,
                             label.y.npc = 0.80)

##############################
# figure forest plot of odds ratios for glmm of loy ~ celltype + (1|library_id)
library(lme4)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggplot2)
library(stringr)
theme_set(theme_sjplot())

meta <- kpmp@meta.data %>%
  dplyr::filter(rkd_genotype %in% c("XY","LOY"), nCount_RNA > 1000) %>%
  dplyr::select(barcode, rkd_genotype, age, specimen_id, celltype) %>%
  dplyr::rename(genotype = rkd_genotype, age_group = age, library_id = specimen_id) 

clinical_meta <- read.csv("G:/kpmp/72ec8f42-750f-4e1b-bda9-7300d4de1365_20221208_OpenAccessClinicalData.csv")
addmeta <- data.frame(library_id = clinical_meta$Participant.ID, Tissue.Type = clinical_meta$Tissue.Type)

meta <- meta %>%
  left_join(addmeta, by = c("library_id")) %>%
  dplyr::mutate(CKD = ifelse(Tissue.Type == "CKD",1,0)) 

df <- meta  %>%
  dplyr::mutate(loy = ifelse(genotype == "LOY",1,0))

# relevel
df$celltype <- factor(df$celltype, levels = c("PT","PT_VCAM1","PT_PROM1","PT_MT","PEC","LOH","DCT","PC","IC","PODO","ENDO","FIB_VSMC_MC","LEUK"))

# glm for effect of celltype on loy
mod1 <- glm(loy ~ celltype, family = binomial(link = "logit"), data = df)
summary(mod1)
# confint(mod1)

# glmm with mixed effect for sample for celltype on loy
mod2 <- lme4::glmer(loy ~ celltype + (1|library_id),
                    family = binomial(link = "logit"),
                    data = df)
summary(mod2)

# glmm with mixed effect for sample for celltype and ckd on loy
mod3 <- lme4::glmer(loy ~ celltype + CKD + (1|library_id),
                    family = binomial(link = "logit"),
                    data = df)
summary(mod3)

plotF <- plot_model(mod3, p.adjust=TRUE)

# swap out yaxis labels
newlabels <- str_replace(plotF$data$term,"celltype","")
newlabels <- factor(newlabels, levels=newlabels)

plotF <- plotF + scale_y_continuous(limits = c(0.01, 5), breaks = c(0.1, 0.5, 1, 5)) +
  ggtitle("F)") +
  geom_hline(yintercept = 1, linetype = "dotted") + 
  theme_bw() +
  theme(plot.title = element_text(face="bold"), axis.text = element_text(color="black")) + 
  ylab("Odds Ratios for LOY relative to PT") + 
  scale_x_discrete(labels=rev(newlabels))

###################
# arrange plots
library(gridExtra)
lay <- rbind(c(1,1,2),
             c(1,1,3),
             c(4,4,4),
             c(5,6,6))
gs <- list(plotA, plotB, plotC, plotD, plotE, plotF)

pdf(here("ckd","figures","figure3.pdf"),width = 7.5, height = 10)
grid.arrange(grobs = gs, layout_matrix = lay, heights=c(0.75,1.25,1,1))
dev.off()


