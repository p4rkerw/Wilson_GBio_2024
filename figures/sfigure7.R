# LOY analysis in leukocytes from RCC donors by snATAC

#########################
# figure snATAC umap
library(Seurat)
library(here)
library(ggplot2)
library(dplyr)

# create a plot dir
plots <- here("ckd","figures")
dir.create(here(plots))

srat <- readRDS(here("ckd","atac_aggr_prep_rccleuk","step4_anno.rds"))
plotA <- DimPlot(srat, label=T, group.by="predicted.l1")

plotA <- plotA + 
  theme_bw() +
  theme(plot.title = element_text(face="bold"), legend.pos = "none", axis.text = element_text(color="black"), panel.grid = element_blank()) +
  ggtitle("A)") +
  annotate(geom="text",
           x = -19,
           y = -8,
           label = paste0("Nuclei = ", nrow(srat@meta.data),"\n","Donors = 8 (7M,1F)","\n","Libraries = 20"),
           hjust=0)
#########################
# kernel density plot
library(dplyr)
library(ggplot2)
library(here)
library(gridExtra)

atac_counts <- readRDS(here("ckd","atac_aggr_prep_rccleuk","atac_frags.rds"))

# ATAC kernel density trough
# find the trough of the chrY scaled fragment kernel density estimate to call LOY vs XY genotype
kd <- atac_counts %>%
  group_by(seqnames) %>%
  dplyr::filter(seqnames == "chrY", sex == "male")
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

#########################
# proportion celltype by tissue source
clinical_meta <- read.xlsx(here("ckd","supplemental_data_1_clinical_data_qc.xlsx"), sheet = "sample_info") %>%
  dplyr::select(library_id, donor, sex)
meta <- srat@meta.data %>% left_join(clinical_meta, by = "library_id")
rownames(meta) <- rownames(srat@meta.data)

plotC <- meta %>%
  dplyr::select(library_id, predicted.l1, tissue_source, donor) %>%
  dplyr::group_by(library_id, predicted.l1) %>%
  dplyr::add_count(name = "celltype_by_library") %>%
  dplyr::group_by(library_id) %>%
  dplyr::add_count(name = "total_cells_library") %>%
  dplyr::mutate(prop_celltype = celltype_by_library / total_cells_library) %>%
  dplyr::distinct(predicted.l1, prop_celltype, tissue_source, library_id) %>%
  ggplot(aes(predicted.l1, prop_celltype, fill=tissue_source)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=tissue_source)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(face="bold"),
        axis.text = element_text(color="black")) +
  xlab("") +
  ylab("Proportion cell type in sample") +
  guides(fill=guide_legend(title="Source")) + 
  scale_fill_brewer(labels = c("adjacent normal", "pbmc", "tumor")) +
  ggtitle("C)")

#########################
# proportion cell type by LOY
df <- table(library_id = meta$library_id,
            celltype = meta$predicted.l1,
            genotype = meta$akd_genotype,
            source = meta$tissue_source) %>%
  as.data.frame() %>%
  dplyr::group_by(library_id, celltype) %>%
  dplyr::mutate(prop_celltype = Freq / sum(Freq)) %>%
  dplyr::filter(genotype == "LOY") %>%
  mutate_at(vars(prop_celltype), ~replace(., is.nan(.), 0))

plotD <- df  %>%
  ggplot(aes(celltype, prop_celltype, fill=source)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=source)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(face="bold"),
        axis.text = element_text(color="black")) +
  xlab("") +
  ylab("Proportion cell type with LOY") +
  guides(fill=guide_legend(title="Source")) + 
  scale_fill_brewer(labels = c("adjacent normal", "pbmc", "tumor")) +
  ylim(c(0,0.10)) +
  ggtitle("D)")


# layout figure
lay <- rbind(c(1,1,2),
             c(3,3,3),
             c(4,4,4))
gs <- list(plotA, plotB, plotC, plotD)

pdf(here(plots,"sfigure7.pdf"),width = 7.5, height = 10)
 grid.arrange(grobs = gs, layout_matrix = lay)
dev.off()
