#########################
# figure multiome umap
library(Seurat)
library(here)
library(ggplot2)

# create a plot dir
plots <- here("ckd","figures")
dir.create(here(plots))

multi <- readRDS(here("ckd","multi_aggr_prep","step3_multi_loy.rds"))
plotA <- DimPlot(multi, label=T)

plotA <- plotA + 
  theme_bw() +
  theme(plot.title = element_text(face="bold"), legend.pos = "none", axis.text = element_text(color="black"), panel.grid = element_blank()) +
  ggtitle("A)") +
  annotate(geom="text",
           x = -18,
           y = -9,
           label = paste0("Nuclei = ", nrow(multi@meta.data),"\n","Donors = 9 (5M,4F)"),
           hjust=0)

#########################
# figure 1-dimensional kernel density estimate for chrY for both RNA and ATAC assays
# with trough identification and shading
library(dplyr)
library(ggplot2)
library(here)
library(gridExtra)

multi_counts <- readRDS(here("ckd","multi_aggr_prep","step3_multi_counts.rds"))

# ATAC kernel density trough
# find the trough of the chrY scaled fragment kernel density estimate to call LOY vs XY genotype
kd <- multi_counts %>%
  dplyr::filter(seqnames == "chrY", malesex == 1) %>%
  ungroup()
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
       dplyr::filter(genotype == "LOY")
  
plotB <- p +
      annotate(geom="text",
               x = 3,
               y = 0.5,
               label = paste0("LOY = ", 100*round(res$prop, digits=2),"%")) +
      theme_bw() + 
      theme(plot.title = element_text(face="bold"), axis.text = element_text(color="black"), panel.grid = element_blank()) +
      ggtitle("B)") +
      xlab("Scaled ATAC coverage")

##########################################
# figure 1-dimensional kernel density estimate
# RNA kernel density trough
# find the trough of the chrY scaled fragment kernel density estimate to call LOY vs XY genotype
kd <- multi_counts %>%
  dplyr::filter(seqnames == "chrY", malesex == 1) %>%
  ungroup()
# find the maximum of the tallest peak in the density plot and its
# corresponding x coordinate (this is the median xcoord for XY cells)
density <- density(kd$scaled_rna_counts, n=10000)
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
       summarize(genotype) %>%
       table() %>%
       as.data.frame() %>%
       dplyr::mutate(total = sum(Freq)) %>%
       dplyr::mutate(prop = Freq / total) %>%
       dplyr::filter(genotype == "LOY")
  
plotC <- p +
      annotate(geom="text",
               x = 3,
               y = 0.5,
               label = paste0("LOY = ", 100*round(res$prop, digits=2),"%")) +
      theme_bw() + 
      theme(plot.title = element_text(face="bold"), axis.text = element_text(color="black"), panel.grid = element_blank()) +
      ggtitle("C)")  +
      xlab("Scaled RNA coverage")

############################
# figure 2-dimensional kernel density estimate for chrY for both RNA and ATAC assays
library(ggside)
library(viridis)
library(dplyr)
library(ggplot2)

multi_counts <- readRDS(here("ckd","multi_aggr_prep","step3_multi_counts.rds"))

# set palette colors for density plot
palette <- viridis(7)
palette[1] <- "white"

toplot <- multi_counts %>%
  na.omit() %>%
  dplyr::filter(scaled_rna_counts < 2, scaled_atac_frags < 2) %>%
  dplyr::filter(seqnames == "chrY", malesex == 1) 
toplot$Genotype <- ifelse(toplot$gmm_genotype == "LOY",1,0)
toplot$Genotype <- as.factor(toplot$Genotype)
levels(toplot$Genotype) <- c("XY","LOY")

plotD <- toplot %>%
  ggplot(aes(scaled_rna_counts, scaled_atac_frags)) +
  geom_density2d_filled(aes(fill = after_stat(level)), bins=7) +
  scale_fill_manual(values = palette) +
  geom_xsidedensity(aes(y = ..scaled.., xfill=Genotype), show.legend = FALSE) +
  geom_ysidedensity(aes(x = ..scaled.., yfill=Genotype)) +
  theme_bw() +
  scale_xfill_manual(values = c("lightgray","#CB181D")) +
  scale_yfill_manual(values = c("lightgray","#CB181D")) +
  labs(fill = "Density") +
  xlab("Scaled RNA coverage") +
  ylab("Scaled ATAC coverage") +
  scale_ysidex_continuous(n.breaks=2) +
  scale_xsidey_continuous(n.breaks=2) +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_fixed() +
  ggtitle("D)") + 
  theme(plot.title = element_text(face="bold"), axis.text = element_text(color="black"))

#######################################################
# figure histogram with proportion of each cell type with loy
# calculate aggregate proportion LOY per celltype
library(dplyr)
library(ggplot2)

multi_counts <- readRDS(here("ckd","multi_aggr_prep","step3_multi_counts.rds"))

toplot <- multi_counts %>%
  dplyr::filter(malesex == 1) %>%
  dplyr::filter(seqnames == "chrY") %>%
  dplyr::distinct(library_id, barcode, celltype, gmm_genotype) %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarize(celltype, gmm_genotype) %>%
  table() %>%
  as.data.frame() %>%
  group_by(celltype) %>%
  dplyr::mutate(total = sum(Freq)) %>%
  dplyr::mutate(loy_prop = Freq / total) %>%
  as.data.frame() %>%
  dplyr::filter(gmm_genotype == "LOY")

# plot aggregate proportion
plot_aggregate <- toplot %>% 
  ggplot(aes(x=celltype, y=loy_prop, fill=celltype)) +
  geom_histogram(stat = "identity") + 
  theme(axis.text.x = element_text(angle=90, vjust=1, hjust=1)) +
  xlab("") +
  ylab("Proportion LOY per cell type") +
  NoLegend()

# # plot mean proportion averaged across all male donors
# # group leukocytes together because there are very few of them
# only consider estimates if there are at least 50 cells for the donor sample
toplot <- multi@meta.data %>%
  dplyr::mutate(celltype = case_when(celltype %in% c("TCELL","BCELL","MONO") ~ "LEUK",
                                     celltype %in% c("PCT","PST") ~ "PT",
                                     celltype %in% c("DCT1","DCT2") ~ "DCT",
                                     celltype %in% c("TAL1","TAL2") ~ "LOH",
                                     celltype %in% c("ICA","ICB") ~ "IC",
                                     TRUE ~ as.character(celltype))) %>%
  na.omit() %>%
  dplyr::filter(gmm_genotype %in% c("XY","LOY")) %>%
  group_by(library_id) %>%
  summarize(celltype, gmm_genotype) %>%
  table() %>%
  as.data.frame() %>%
  group_by(library_id, celltype) %>%
  mutate(num_cells = sum(Freq)) %>%
  dplyr::filter(gmm_genotype == "LOY") %>%
  mutate(lib_celltype_prop_loy = Freq / num_cells) %>%
  na.omit()
# relevel
toplot$celltype <- factor(toplot$celltype, levels = c("PT","PT_VCAM1","PEC","LOH","DCT","PC","IC","PODO","ENDO","FIB_VSMC_MC","LEUK"))

plotE <- toplot %>%
  group_by(celltype) %>%
  ggplot(aes(x=celltype, y=lib_celltype_prop_loy, fill=celltype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(w=0.3), aes(celltype), size=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") +
  ylab("Proportion LOY") +
  xlab("") +
  ggtitle("E)") + 
  theme(plot.title = element_text(face="bold"), axis.text = element_text(color="black"))

# # run a kruskal-wallis test for library proportion LOY by celltype
# kruskal.test(lib_celltype_prop_loy ~ celltype, data = toplot)

# pairwise tests
# wilcox.test(toplot[toplot$celltype == "PT",]$lib_celltype_prop_loy, toplot[toplot$celltype == "PT_VCAM1",]$lib_celltype_prop_loy)
# pairwise.t.test(toplot$lib_celltype_prop_loy, toplot$celltype, p.adj = "none")
# pairwise.t.test(toplot$lib_celltype_prop_loy, toplot$celltype, p.adj = "holm")

# # arrange
# pdf(here(plots, "figure1.pdf"), width=8.5, height=8.5)
# margin = theme(plot.margin = unit(c(2,2,2,2), "cm"))
# grid.arrange(plotA, plotB, plotC, plotD, widths = c(3,2)) 
# grid.arrange(plotD + margin)
# dev.off()
###################################
#######################################################
# # figure 1f gsea for PCT LOY vs PCT XY
library(DOSE)
library(enrichplot)
gse <- readRDS(here("ckd","multi_aggr_prep","step6_gsea.PCT.loy_vs_xy.age_adjust.rds"))

terms <- gse@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::mutate(log_pval = -log10(p.adjust)) %>%
  dplyr::mutate(sign = ifelse(NES > 0, 1, 0)) %>%
  group_by(sign) %>%
  top_n(10, log_pval)

# there are identical gene sets that have multiple descriptions
n_distinct(gse@result$core_enrichment)
nrow(gse@result)
  
# extract first description for each unique gene set
res <- gse@result %>%
  group_by(core_enrichment) %>%
  arrange(core_enrichment) %>%
  filter(row_number ()==1) %>%
  ungroup()
gse@result <- res

plot <- dotplot(gse, showCategory=terms$Description, split=".sign") + 
  facet_grid(.~.sign) +
  theme_bw() +
  theme(legend.pos = "none", plot.title = element_text(face="bold"), axis.text = element_text(color="black")) +
  ggtitle("F)")
print(plot)  

# rename BP to shorten names
levels(plot$data$Description) <- rev(c("growth hormone signaling","apoptotic process","cellular respiration",
                                      "estrogen receptor signaling","symbiotic interaction","signal transduction",
                                      "cytokinesis","endodermal differentiation","heterochromatin organization",
                                      "regulation of gastrulation","histone H3-K9 modification","protein monoubiquitination","mitotic nuclear division",
                                      "endoderm formation","inclusion body assembly","lipoprotein metabolism","ERAD pathway"))
                                    

plotF <- plot 
  
###################################################################
# figure forest plot of odds ratios for glmm of loy ~ celltype + (1|library_id)
library(lme4)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggplot2)
library(stringr)
theme_set(theme_sjplot())

df <- multi@meta.data %>%
  dplyr::filter(malesex == 1, gmm_genotype %in% c("XY","LOY")) %>%
  dplyr::mutate(loy = ifelse(gmm_genotype == "LOY",1,0)) %>%
  dplyr::select(library_id, loy, celltype, nCount_RNA, nCount_ATAC)

df <- df %>%
  dplyr::mutate(celltype = case_when(celltype %in% c("TCELL","BCELL","MONO") ~ "LEUK",
                                     celltype %in% c("PCT","PST") ~ "PT",
                                     celltype %in% c("DCT1","DCT2") ~ "DCT",
                                     celltype %in% c("TAL1","TAL2") ~ "LOH",
                                     celltype %in% c("ICA","ICB") ~ "IC", 
                                     TRUE ~ as.character(celltype))) 
# relevel
df$celltype <- factor(df$celltype, levels = c("PT","PT_VCAM1","PEC","LOH","DCT","PC","IC","PODO","ENDO","FIB_VSMC_MC","LEUK"))


# glm for effect of celltype on loy
mod1 <- glm(loy ~ celltype, family = binomial(link = "logit"), data = df)
summary(mod1)
# confint(mod1)

# glmm with mixed effect for sample for celltype on loy
mod2 <- lme4::glmer(loy ~ celltype + (1|library_id),
                    family = binomial(link = "logit"),
                    data = df)
summary(mod2)

plotG <- plot_model(mod2, p.adjust=TRUE)
saveRDS(mod2, here("ckd","figures","multiomes_mod2.rds"))

# swap out yaxis labels
newlabels <- str_replace(plotG$data$term,"celltype","")
newlabels <- factor(newlabels, levels=newlabels)

plotG <- plotG + scale_y_continuous(limits = c(0.01, 2), breaks = c(0.1, 0.5, 1, 1.5)) + theme_bw() + ggtitle("F)") +
  geom_hline(yintercept = 1, linetype = "dotted") + 
  theme(plot.title = element_text(face="bold"), axis.text = element_text(color="black")) + 
  ylab("Odds Ratios for LOY relative to PT") + 
  scale_x_discrete(labels=rev(newlabels))

# layout figure
lay <- rbind(c(1,1,2),
             c(1,1,3),
             c(4,4,5),
             c(4,4,5),
             c(6,6,6))
gs <- list(plotA, plotB, plotC, plotD, plotE, plotF)

# # arrange
# grob <- grid.arrange(grobs = gs, layout_matrix = lay)

# # add margins
# library(ggpubr)
# ggplot <- as_ggplot(grob) +
#   theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

pdf(here(plots,"figure1.pdf"),width = 7.5, height = 10)
 grid.arrange(grobs = gs, layout_matrix = lay, heights = c(1,1,1,1,1.5))
dev.off()




