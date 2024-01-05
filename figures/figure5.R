# visium spatial analysis
#########################
# figure umap
library(Seurat)
library(here)
library(ggplot2)

plots <- here("ckd","figures")
dir.create(here(plots))

spatial <- readRDS(here("ckd","spatial","step1b_spatial.rds"))
paletteMartin <- c(
"#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
"#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
"#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")

palette <- c("lightgray","#ffff6d", "#009292", "#ff6db6", "#ffb6db", 
"#490092", "#006ddb", "#b66dff", "#6db6ff","#24ff24")

names(palette) <- levels(spatial)
Idents(spatial) <- spatial$celltype

plotA <- DimPlot(spatial, label=F, repel=TRUE)

plotA <- plotA + 
  theme_bw() +
  theme(plot.title = element_text(face="bold"), axis.text = element_text(color="black"), panel.grid = element_blank()) +
  ggtitle("A)") +
  scale_color_manual(values = palette) +
  annotate(geom="text",
           size = 3,
           x = -11,
           y = -9,
           label = paste0("Spots = ", nrow(spatial@meta.data),"\n","Donors = 8"),
           hjust=0) +
  coord_fixed() +
  ylim(c(-10,8))

###########################################
# spatial dimplot with rep annnotations
library(Seurat)

palette <- c("lightgray","#ffff6d", "#009292", "#ff6db6", "#ffb6db", 
"#490092", "#006ddb", "#b66dff")

spatial@images[["slice1.1"]]@coordinates[["tissue"]] <- as.integer(spatial@images[["slice1.1"]]@coordinates[["tissue"]])
spatial@images[["slice1.1"]]@coordinates[["row"]] <- as.integer(spatial@images[["slice1.1"]]@coordinates[["row"]])
spatial@images[["slice1.1"]]@coordinates[["col"]] <- as.integer(spatial@images[["slice1.1"]]@coordinates[["col"]])
spatial@images[["slice1.1"]]@coordinates[["imagerow"]] <- as.integer(spatial@images[["slice1.1"]]@coordinates[["imagerow"]])
spatial@images[["slice1.1"]]@coordinates[["imagecol"]] <- as.integer(spatial@images[["slice1.1"]]@coordinates[["imagecol"]])

names(palette) <- levels(spatial)
Idents(spatial) <- spatial$celltype
plotB <- SpatialDimPlot(spatial, images="slice1.1", cols=palette, stroke=NA, image.alpha=0) +
  theme_bw() +
  theme(plot.title = element_text(face="bold"),
        legend.title = element_blank(),
        legend.pos = "none",
        panel.grid = element_blank()) +
  xlab("") +
  ylab("") +
  ggtitle("B)") +
  coord_fixed()

############################################
# cellchat LR signaling plots
library(CellChat)
library(tidyr)
library(tibble)

library_ids <- unique(spatial$library_id)
cellchat.ls <- lapply(library_ids, function(library_id) {
  cellchat <- readRDS(here("ckd","spatial","cellchat",paste0(library_id,"_cellchat.rds")))
  return(cellchat)
})   

net.aggW <- lapply(seq_along(cellchat.ls), function(index) {
  df <- cellchat.ls[[index]]@net$weight
  df <- as.data.frame(df)
  df$source <- rownames(df)
  df$library_id <- library_ids[index]
  return(df)
}) %>% bind_rows()

net.aggW <- net.aggW %>%
  pivot_longer(cols = !c(source, library_id), names_to = "target", values_to = "weight") %>%
  group_by(source, target) %>%
  na.omit()
net.aggW$source <- factor(net.aggW$source, levels = levels(spatial$celltype))
net.aggW$target <- factor(net.aggW$target, levels = levels(spatial$celltype))

mat <- net.aggW %>%
  summarize(mean_weight = mean(weight), source = source, target = target) %>%
  distinct() %>%
  pivot_wider(names_from = "target", values_from = mean_weight) %>%
  column_to_rownames(var = "source") %>%
  as.matrix()

groupSize <- as.numeric(table(spatial$celltype))

# visualize all intercellular interaction weights
# par(mfrow = c(3,4), xpd=TRUE)
# for (i in 1:nrow(mat)) {
#   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#   mat2[i, ] <- mat[i, ]
#   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
# }

library(viridis)
library(gridGraphics)
# PT_VCAM1 only
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2["PT_INJ", ] <- mat["PT_INJ", ]
# png("G:/scratch/circlize.png")
# netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = "PT_VCAM1")
# dev.off()

# # ## grab the scene as a grid object & save it to P1
# grid.echo()
# plotD <- grid.grab()

library(circlize)
palette <- viridis(8)
circos <- data.frame(from = "PT_INJ", to = colnames(mat2), value = mat2[2,])
groups <- circos$to
names(groups) <- circos$to

pdf(here("ckd","figures","figure5_circlize.pdf"))
circos.clear()
par(mar=c(0.5,0.5,0.5,0.5))
circos.par(gap.after = 5)
chordDiagram(circos,
           grid.col = palette,
           annotationTrack = "grid",
           group=groups,
           col=palette,
           directional = 1,
           direction.type = "arrows",
           preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(circos))))))

# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter,
              CELL_META$ylim[1],
              CELL_META$sector.index,
              facing = "clockwise",
              niceFacing = TRUE,
              cex = 1.2,
              adj = c(0, 0.5))
 }, bg.border = NA) # here set bg.border to N
dev.off()

# available pathways
cellchat.ls[[1]]@netP$pathways

# intercellular pathway signaling for source cell type
net.aggP <- lapply(seq_along(cellchat.ls), function(index) {
  df <- subsetCommunication(cellchat.ls[[index]], slot.name="netP", thresh=1)
  df$library_id <- library_ids[index]
  return(df)
}) %>% bind_rows()
net.aggP$source <- factor(net.aggP$source, levels = levels(spatial$celltype))
net.aggP$target <- factor(net.aggP$target, levels = levels(spatial$celltype))

# mat <- net.aggP %>%
#   group_by(source, target) %>%
#   summarize(median_weight = median(prob), source = source, target = target) %>%
#   distinct() %>%
#   pivot_wider(names_from = "target", values_from = median_weight) %>%
#   column_to_rownames(var = "source") %>%
#   as.matrix()

# groupSize <- as.numeric(table(spatial$celltype))

# par(mfrow = c(3,4), xpd=TRUE)
# for (i in 1:nrow(mat)) {
#   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#   mat2[i, ] <- mat[i, ]
#   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
# }
palette <- c("lightgray","#ffff6d", "#009292", "#ff6db6", "#ffb6db", 
"#490092", "#006ddb", "#b66dff")

pathway="MK"
plotC <- net.aggP %>% 
  filter(pathway_name == pathway) %>%
  filter(source == "PT_INJ") %>%
  ggplot(aes(target, prob, fill=target)) + 
  geom_boxplot() +
  ggtitle("C) Midkine") +
  theme_bw() + 
  theme(plot.title = element_text(face="bold"), 
        axis.text = element_text(color="black"), 
        panel.grid = element_blank(),
        legend.pos = "none",
        axis.text.x = element_text(angle=90, vjust=1, hjust=1)) +
  scale_fill_manual(values = palette) +
  xlab("") +
  ylab("Probability")

pathway="EDN"
plotD <- net.aggP %>% 
  filter(pathway_name == pathway) %>%
  filter(source == "PT_INJ") %>%
  ggplot(aes(target, prob, fill=target)) + 
  geom_boxplot() +
  ggtitle("D) Endothelin") +
  theme_bw() + 
  theme(plot.title = element_text(face="bold"), 
        axis.text = element_text(color="black"), 
        panel.grid = element_blank(),
        legend.pos = "none",
        axis.text.x = element_text(angle=90, vjust=1, hjust=1)) +
  scale_fill_manual(values = palette) +
  xlab("") +
  ylab("Probability")

pathway="EGF"
plotE <- net.aggP %>% 
  filter(pathway_name == pathway) %>%
  filter(source == "PT_INJ") %>%
  ggplot(aes(target, prob, fill=target)) + 
  geom_boxplot() +
  ggtitle("E) EGF") +
  theme_bw() + 
  theme(plot.title = element_text(face="bold"), 
        axis.text = element_text(color="black"), 
        panel.grid = element_blank(),
        legend.pos = "none",
        axis.text.x = element_text(angle=90, vjust=1, hjust=1)) +
  scale_fill_manual(values = palette) +
  xlab("") +
  ylab("Probability")

# pathway="PDGF"
# net.aggP %>% 
#   filter(pathway_name == pathway) %>%
#   filter(source %in% c("FIB_VSMC")) %>%
#   ggplot(aes(source, prob, fill=target)) + geom_boxplot()

#####################################
# figure gsea
library(DOSE)
library(enrichplot)
gse <- readRDS(here("ckd","spatial","step3_gsea.rds"))

terms <- gse@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::mutate(log_pval = -log10(p.adjust)) %>%
  dplyr::mutate(sign = ifelse(NES > 0, 1, 0)) %>%
  dplyr::arrange(desc(log_pval)) %>%
  group_by(sign) %>%
  slice_head(n=5)

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
levels(plot$data$Description) <- rev(c("small molecule catabolism","organic acid metabolism","carboxylic acid metabolism","small molecule metabolism",
                                      "humoral immune response","immune system process","phagocytosis","regulation of immune","transmembrane transport",
                                      "response to stimulus"))
                                    

plotF <- plot 


library(gridExtra)
library(cowplot)
lay <- rbind(c(1,2,2),
             c(3,4,5),
             c(NA,NA,NA))
gs <- list(plotA, plotB, plotC, plotD, plotE, plotF)


row1 <- plotA + plotB
row2 <- cowplot::plot_grid(plotC, plotD, plotE, ncol=3)

pdf(here("ckd","figures","figure5.pdf"), width = 7.5, height = 10)
cowplot::plot_grid(row1, row2, plotF, ncol=1)
dev.off()
