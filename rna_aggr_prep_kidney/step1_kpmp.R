# to run locally:
# SCRATCH1=/mnt/g/scratch
# docker run -it --rm \
# --workdir $HOME \
# -v /mnt/g/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# -v /mnt/g/cellranger_multi_counts:$HOME/cellranger_multi_counts \
# -v /mnt/g/kpmp:$HOME/kpmp \
# -v /mnt/g/diabneph:$HOME/diabneph \
# -v /mnt/g/reference:$HOME/reference \
# -v /mnt/g/scratch:$HOME/scratch \
# -v /mnt/g/ckd:$HOME/ckd \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/sctools:R4.1.3b R

library(SeuratDisk)
library(Seurat)
library(data.table)
library(dplyr)
library(here)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)
library(future)
library(harmony)

plots <- here("ckd","rna_aggr_prep","kpmp","plots")
dir.create(here(plots), recursive=TRUE)

# load the kpmp obj
kpmp <- LoadH5Seurat(here("kpmp","c798e11b-bbde-45dd-bd91-487f27c93f8f_WashU-UCSD_HuBMAP_KPMP-Biopsy_10X-R_12032021.h5Seurat"), assays = list(RNA = "counts"))

# label transfer from reference obj to kpmp obj to harmonize annotations
rna <- readRDS(here("diabneph","analysis","dkd","rna_aggr_prep","step2_anno.rds"))
DefaultAssay(rna) <- "RNA"
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(rna)
rna <- ScaleData(rna, features = all.genes)
rna <- RunPCA(rna, features = VariableFeatures(object = rna))

Idents(rna) <- "celltype"  
transfer.anchors <- FindTransferAnchors(
  reference = rna,
  query = kpmp,
  reduction = 'cca',
  normalization.method = 'LogNormalize',
  reference.assay = 'RNA',
  query.assay = 'RNA'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rna$celltype,
  weight.reduction = kpmp[['pca']],
  dims = 1:30
)

kpmp <- AddMetaData(object = kpmp, metadata = predicted.labels)

pdf(here(plots,"step1_kpmp_label.pdf"))
DimPlot(kpmp, label=TRUE, group.by="subclass.l1") + NoLegend()
DimPlot(kpmp, label=TRUE, group.by="predicted.id") + NoLegend()
FeaturePlot(kpmp, features = "prediction.score.max")
dev.off()

# filter kpmp dataset by prediction scores
kpmp <- subset(kpmp, prediction.score.max > 0.5 & percent.mt < 10 & nCount_RNA < 10000)

# perform batch correction on lsi embeddings
my_harmony_embeddings <- HarmonyMatrix(
  data_mat  = as.matrix(kpmp@reductions$pca@cell.embeddings),
  meta_data = kpmp@meta.data,
  vars_use  = "specimen_id",
  do_pca = FALSE
)
rownames(my_harmony_embeddings) <- rownames(kpmp@reductions$pca@cell.embeddings)

#store the harmony reduction as a custom dimensional reduction called 'harmony' in the default assay
kpmp[["harmony"]] <- CreateDimReducObject(embeddings = my_harmony_embeddings, key = "harmony_", assay = DefaultAssay(kpmp))
kpmp <- FindNeighbors(object = kpmp, reduction = "harmony", dims = 1:30)
kpmp <- FindClusters(object = kpmp, verbose = TRUE, algorithm = 1) # Louvain algorithm
kpmp <- RunUMAP(object = kpmp, reduction = "harmony", dims = 1:26)

pdf(here(plots,"step2_kpmp_harmony.pdf"))
DimPlot(kpmp, label=TRUE, group.by="subclass.l1") + NoLegend()
DimPlot(kpmp, label=TRUE, group.by="subclass.l2") + NoLegend()
DimPlot(kpmp, label=TRUE, group.by="predicted.id") + NoLegend()
FeaturePlot(kpmp, features = "prediction.score.max")
FeaturePlot(kpmp, features = "nCount_RNA")
FeaturePlot(kpmp, features = "percent.mt")
DimPlot(kpmp, label=TRUE, group.by="seurat_clusters") + NoLegend()
VlnPlot(kpmp, features = "percent.mt", group.by = "seurat_clusters", pt.size=0)
VlnPlot(kpmp, features = "nCount_RNA", group.by = "seurat_clusters", pt.size=0)
dev.off()

marker.genes <- c("CUBN","HAVCR1","SLC5A1","SLC5A2","VCAM1","PROM1", # PT and PT-VCAM1+ markers
                  "CFH", # PEC
                  "SLC12A1", # TAL NKCC2
                  "CLDN10", #MTAL TAL2
                  "CLDN16", #CTAL TAL1
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
                  "IGFBP5","IGFBP7", # PTC and AVR
                  "PLVAP", # PTC and AVR https://www.nature.com/articles/s41467-019-12872-5
                  "EHD3", # GEC
                  "SLC6A6","SLC14A1","AQP1", # EA and DVR
                  "NOS1", # MD
                  "ITGA8","PDGFRB","MEIS2","PIEZO2","REN", # MES and JGA
                  "ACTA2","CALD1", # FIB
                  "PROX1","FLT4","PDPN", # Lymphatics
                  "PTPRC","CD3E","MS4A1","CD19","SDC1", # Lymphocytes and plasma cells
                  "FCGR3A","CD14","CSF1R") # Monocyte / Macrophage

pdf(here(plots,"step3_kpmp_dotplot.pdf"), width=10, height=6)
DefaultAssay(kpmp) <- "RNA"
DotPlot(kpmp, features=marker.genes, group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
VlnPlot(kpmp, features = c("VCAM1","PROM1"), group.by="seurat_clusters", pt.size=0)
dev.off()

# annotate the kpmp object
Idents(kpmp) <- "seurat_clusters"
kpmp <- RenameIdents(kpmp,
                     '0' = "PT",
                     '1' = "LOH",
                     '2' = "PT_VCAM1",
                     '3' = "ENDO",
                     '4' = "FIB_VSMC_MC",
                     '5' = "PT",
                     '6' = "PC",
                     '7' = "LOH",
                     '8' = "DCT",
                     '9' = "PT_PROM1",
                     '10' = "DCT",
                     '11' = "IC",
                     '12' = "LEUK",
                     '13' = "ENDO",
                     '14' = "PT_MT",
                     '15' = "FIB_VSMC_MC",
                     '16' = "LEUK",
                     '17' = "PODO",
                     '18' = "ENDO",
                     '19' = "PC",
                     '20' = "IC",
                     '21' = "LEUK",
                     '22' = "PEC",
                     '23' = "FIB_VSMC_MC",
                     '24' = "FIB_VSMC_MC",
                     '25' = "LOH",
                     '26' = "ENDO",
                     '27' = "LEUK"
                    )
                     
pdf(here(plots,"step4_kpmp_anno.pdf"), width=10, height=6)
DimPlot(kpmp, label=TRUE) + NoLegend()
dev.off()

kpmp@meta.data$celltype <- Idents(kpmp)

saveRDS(kpmp, here("ckd","rna_aggr_prep","kpmp","step1_kpmp.rds"), compress=FALSE)

# prepare list of chry features from gtf
gtf <- fread(here("reference","refdata-gex-GRCh38-2020-A","genes","genes.gtf"))
colnames(gtf)[1] <- "chrom"
gtf <- dplyr::filter(gtf, chrom == "chrY", V3 == "transcript")
gtf <- tidyr::separate(gtf, V9, sep = ";", into = c("gene_id","gene_version","transcript_id","transcript_version","gene_type","gene_name","transcript_type","transcript_name","transcript_support","hgnc_id","havana_gene","havana_transcript"))
genes <- gsub(" gene_name ","", gtf$gene_name)
genes <- gsub("\"",replacement = "", genes) %>% sort() %>% unique()
chry_features <- genes

# raw counts are in the "counts" slot
DefaultAssay(kpmp) <- "RNA"
rna_counts <- GetAssayData(kpmp, slot = "counts")

# compute total number of rna counts per cell
# the RNA assay transcript counts will be referred to as "counts" 
total_rna_counts <- colSums(rna_counts)
total_rna_counts <- data.frame(total_rna_counts = total_rna_counts, barcode = names(total_rna_counts))

# subset for chrY counts
chry_rna_counts <- rna_counts[rownames(rna_counts) %in% chry_features,]
chry_rna_counts <- chry_rna_counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "feature") %>%
  pivot_longer(cols = colnames(chry_rna_counts)) %>%
  dplyr::rename(barcode = name)
chry_rna_counts$library_id <- str_split(chry_rna_counts$barcode, pattern = "_", simplify=TRUE)[,1]
chry_rna_counts <- chry_rna_counts %>%
  left_join(total_rna_counts, by = "barcode")

# optional: subset chry features 
# chry_features_keep <- rnacounts %>%
#   group_by(feature, gem) %>%
#   summarize(median_feature_count = median(value), mean_feature_count = mean(value)) %>%
#   dplyr::filter(median_feature_count > 0) %>%
#   dplyr::select(feature) %>%
#   unique()
# rnacounts <- rnacounts %>%
#   dplyr::filter(feature %in% chry_features_keep$feature)

# log normalize the counts after summing all counts on chrY
chry_rna_counts <- chry_rna_counts %>%
  group_by(barcode) %>%
  dplyr::mutate(chrom_rna_counts = sum(value)) %>%
  dplyr::mutate(log_rna_counts = log10(chrom_rna_counts / total_rna_counts + 1)) %>%
  distinct(barcode, chrom_rna_counts, log_rna_counts, library_id, total_rna_counts)

# density plots
pdf(here(plots,"rna_ydensity.pdf"))
chry_rna_counts %>%
  ggplot(aes(log_rna_counts, color=library_id, fill=library_id)) +
  geom_density(aes(y = ..scaled..)) +
  facet_wrap(~library_id) 
dev.off()

# add cell type annotation to visualize loy across cell types
celltype <- data.frame(barcode = rownames(kpmp@meta.data), celltype = kpmp$celltype, specimen_id=kpmp$specimen_id)
chry_rna_counts <- chry_rna_counts %>%
  left_join(celltype, by = "barcode")

chry_rna_counts$seqnames <- "chrY"

# read in metadata
meta <- read.csv(here("kpmp","72ec8f42-750f-4e1b-bda9-7300d4de1365_20221208_OpenAccessClinicalData.csv")) %>%
  dplyr::rename(specimen_id = `Participant.ID`,
                egfr = `Baseline.eGFR..ml.min.1.73m2...Binned.`,
                age = `Age..Years...Binned.`,
                sex = Sex) %>%
  dplyr::select(specimen_id, egfr, age, sex)

meta <- meta %>% dplyr::mutate(age = ifelse(age == "", NA, age))

# join meta data with counts
chry_rna_counts <- chry_rna_counts %>%
  left_join(meta, by = "specimen_id")

# visualize median chrY RNA counts and ATAC frags across celltype
# note that different cell types have variable expression / accessibility levels
# that need to be corrected prior to scaling
# chrY is a special case where we will only consider the male samples to compute a global_median (because female samples are not informative)
pdf(here(plots, "celltype_chry.pdf"))
  toplot <- chry_rna_counts %>%
    group_by(seqnames) %>%
    mutate(global_median = median(log_rna_counts[sex == "Male"])) %>%
    group_by(library_id, celltype) %>%
    summarize(median=median(log_rna_counts), global_median = global_median, sex = sex) %>%
    distinct()
  global_median <- unique(toplot$global_median)
  toplot %>%
    ggplot(aes(celltype, median, fill=as.factor(sex))) + 
    geom_boxplot() +
    geom_hline(yintercept = global_median, color="red")
dev.off()

# compute global median for each chromosome
# for chrY only consider male samples to compute the median
chry_rna_counts <- chry_rna_counts %>%
  group_by(seqnames) %>%
  mutate(global_median_counts = median(log_rna_counts[sex == "Male"]))

# compute median for each celltype
# for chrY only consider male samples to compute the median
chry_rna_counts <- chry_rna_counts %>%
  group_by(seqnames, celltype) %>%
  mutate(celltype_median_counts = median(log_rna_counts[sex == "Male"]))   

# compute celltype correction factor
# and correct the counts
chry_rna_counts <- chry_rna_counts %>%
  group_by(seqnames, celltype) %>%
  mutate(correction_factor_counts = global_median_counts / celltype_median_counts) %>%
  mutate(corrected_rna_counts = log_rna_counts * correction_factor_counts)

# visualize uncorrected and corrected counts and frags for chrY 
pdf(here(plots, "celltype_corrected_chry.pdf"))
  chry_rna_counts %>%
    filter(seqnames == "chrY", sex == "Male") %>%
    ggplot(aes(celltype, log_rna_counts)) + 
    geom_boxplot() 

  chry_rna_counts %>%
    filter(seqnames == "chrY", sex == "Male") %>%
    ggplot(aes(celltype, corrected_rna_counts)) + 
    geom_boxplot() 
dev.off()

# scale corrected rna counts by library_id
# to correct for donor-specific variability
chry_rna_counts <- chry_rna_counts %>%
  group_by(seqnames, library_id) %>%
  dplyr::mutate(scaled_rna_counts = scale(corrected_rna_counts, center=0))

# save the chrY RNA counts
saveRDS(chry_rna_counts, here("ckd","rna_aggr_prep","kpmp","rna_counts.rds"))

chry_rna_counts %>%
  ggplot(aes(scaled_rna_counts, color=sex, fill=sex)) +
  geom_density(aes(y = ..scaled..)) +
  facet_wrap(~sex) +
  xlim(c(-5,5))

chry_rna_counts %>%
  na.omit() %>%
  dplyr::filter(sex == "Male") %>%
  ggplot(aes(scaled_rna_counts, color=age, fill=age)) +
  geom_density(aes(y = ..scaled..)) +
  facet_wrap(~age) +
  xlim(c(0,5))

# select male samples for kd estimate
toplot <- chry_rna_counts %>% 
  dplyr::ungroup() %>%
  dplyr::filter(sex == "Male") 

                                        
# find the maximum of the tallest peak in the density plot and its
# corresponding x coordinate
density <- density(toplot$scaled_rna_counts, n=10000)
# search for the first trough between zero and the maximum x coord
x <- density$x
y <- density$y
# calculate the difference in the y coord for successive steps in density function
diffy <- diff(y)
# find the sign of the change for the difference
signy <- sign(diffy)
# find points where the difference in the sign of change for successive steps in density function
# is negative 2 (ie. the first point was ascending and the second point was descending)
diffpeak <- diff(signy) == 2
# find the first two peak values
local_trough_index <- which(diffpeak)[1]
trough <- x[local_trough_index]
toplot <- toplot %>%
  dplyr::mutate(fill = ifelse(scaled_rna_counts < trough & scaled_rna_counts, "LOY", "Haploid"))
# log-normalize, scale and center on the maximum of the density plot
p2 <- toplot %>% 
  group_by(sex) %>%
  ggplot(aes(scaled_rna_counts, color = sex, fill = sex)) + 
  geom_density(aes(y = ..scaled..)) +
  xlab(paste0("Scaled coverage: ", "chrY")) +
  ylab("Density") +
  facet_wrap( ~ sex) +
  xlim(c(-3,3)) +
  theme(legend.position="none")
# create a shaded density plot
p <-  toplot %>% 
  ggplot(aes(scaled_rna_counts)) + 
  geom_density(aes(y = ..scaled..)) +
  geom_vline(xintercept = trough, linetype = "dotted") + 
  xlab(paste0("Scaled coverage: ", "chrY")) +
  ylab("Density") +
  xlim(c(0,5))
d = ggplot_build(p)$data[[1]]
p = p + geom_area(data = subset(d, x < trough), aes(x=x,y=y), fill = "#00AFBB", alpha = 0.5)
p = p + geom_area(data = subset(d, x > trough), aes(x=x,y=y), fill = "darkgrey", alpha = 0.5)
res <- toplot %>%
  summarize(loy = scaled_rna_counts < trough) %>%
  table() %>%
  as.data.frame() %>%
  dplyr::mutate(total = sum(Freq)) %>%
  dplyr::mutate(loy_prop = Freq / total)
colnames(res)[1] <- "LOY"
res <- res[res$LOY == "TRUE",]

p3 <- p +
  annotate(geom="text",
           x = 1.5,
           y = 0.8,
           label = paste0("LOY: ", round(res$loy_prop, digits=2), "\n", "Cells: ", res$total)) +
  theme_bw()
rm(p)

# summarize results
loy <- toplot %>% 
  na.omit() %>%
  group_by(celltype) %>%
  summarize(loy = scaled_rna_counts < trough) %>%
  table() %>%
  as.data.frame() %>%
  pivot_wider(names_from = loy, values_from = Freq, names_prefix = "LOSS_") %>%
  dplyr::mutate(loy_prop = LOSS_TRUE/(LOSS_TRUE + LOSS_FALSE)) %>%
  dplyr::mutate(total_cells = LOSS_TRUE + LOSS_FALSE) %>%
  dplyr::filter(total_cells > 100)

loy %>%
  ggplot(aes(celltype, loy_prop, color=celltype, fill=celltype, label = total_cells)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=total_cells), vjust=0, nudge_y = 0.001, color = "black") +
  ylab("Proportion LOY in cell type")

pdf(here(plots,"step5_kpmp_loy.pdf"))
DimPlot(kpmp, label=TRUE, group.by="subclass.l1") + NoLegend()
loy <- toplot %>% 
  na.omit() %>%
  group_by(celltype) %>%
  summarize(loy = scaled_rna_counts < trough) %>%
  table() %>%
  as.data.frame() %>%
  pivot_wider(names_from = loy, values_from = Freq, names_prefix = "LOSS_") %>%
  dplyr::mutate(loy_prop = LOSS_TRUE/(LOSS_TRUE + LOSS_FALSE)) %>%
  dplyr::mutate(total_cells = LOSS_TRUE + LOSS_FALSE) %>%
  dplyr::filter(total_cells > 100)

loy %>%
  ggplot(aes(celltype, loy_prop, color=celltype, fill=celltype, label = total_cells)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=total_cells), vjust=0, nudge_y = 0.001, color = "black") +
  ylab("Proportion LOY in cell type")
dev.off()

###########################################
# # plot mean proportion averaged across all male donors
# # group leukocytes together because there are very few of them
# only consider estimates if there are at least 50 cells for the donor sample
loy <- toplot %>%
  dplyr::mutate(loy = scaled_rna_counts < trough) %>%
  na.omit() %>%
  group_by(library_id) %>%
  summarize(celltype, loy) %>%
  table() %>%
  as.data.frame() %>%
  group_by(library_id, celltype) %>%
  mutate(num_cells = sum(Freq)) %>%
  dplyr::filter(loy == "TRUE") %>%
  mutate(lib_celltype_prop_loy = Freq / num_cells) %>%
  dplyr::filter(num_cells > 50) %>%
  na.omit()
# relevel
loy$celltype <- factor(loy$celltype, levels = c("PT","PT_VCAM1","PT_PROM1","PT_MT","PEC","LOH","DCT","PC","IC","PODO","ENDO","FIB_VSMC_MC","LEUK"))

pdf(here(plots,"step6_kpmp_loy_celltype.pdf"))
loy %>%
  group_by(celltype) %>%
  ggplot(aes(x=celltype, y=lib_celltype_prop_loy, fill=celltype)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75), aes(group=celltype)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") +
  ylab("Proportion LOY") +
  xlab("") +
  ggtitle("E)") + 
  theme(plot.title = element_text(face="bold"), axis.text = element_text(color="black"))
dev.off()  

# save kpmp obj with loy anno
toplot <- toplot %>%
  dplyr::mutate(rkd_genotype = ifelse(scaled_rna_counts < trough, "LOY", "XY"))

meta_update <- kpmp@meta.data 
meta_update$barcode <- rownames(kpmp@meta.data)
meta_update <- meta_update %>%
  left_join(meta, by = "specimen_id")

loyanno <- toplot %>% dplyr::select(barcode, rkd_genotype)
meta_update <- meta_update %>%
  left_join(loyanno, by = "barcode")

kpmp@meta.data <- meta_update
rownames(kpmp@meta.data) <- meta_update$barcode

levels(kpmp) <- c("PT","PT_VCAM1","PT_PROM1","PT_MT","PEC","LOH","DCT","PC","IC","PODO","ENDO","FIB_VSMC_MC","LEUK")
kpmp$celltype <- Idents(kpmp)

saveRDS(kpmp, here("ckd","rna_aggr_prep","kpmp","kpmp_prep.rds"), compress=FALSE)
