# add celltype chromosomal RNA and ATAC median correction
# this script will analyze an aggregated multiome and call barcodes for LOY by RNA and ATAC using two finite mixture models
# implemented in mclust
# it will save a new seurat obj with the LOY annotations
# it will also generate a separate file with per barcode chromosome counts

# # to run locally:
# SCRATCH1=/mnt/g/scratch
# docker run -it --rm \
# --workdir $HOME \
# -v /mnt/g/reference:$HOME/reference \
# -v /mnt/g/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# -v /mnt/g/cellranger_multi_counts:$HOME/cellranger_multi_counts \
# -v /mnt/g/ckd:$HOME/ckd \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/sctools:R4.1.3b R
set.seed(1234)

library(Seurat)
library(Signac)
library(here)
library(dplyr)
library(data.table)
library(plyranges)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(mclust)
library(NB.MClust)

multi <- readRDS(here("ckd","multi_aggr_prep","multiome.rds"))

# prepare list of chry features from gtf
gtf <- fread(here("reference","refdata-gex-GRCh38-2020-A","genes","genes.gtf"))
colnames(gtf)[1] <- "chrom"
gtf <- dplyr::filter(gtf, chrom == "chrY", V3 == "transcript")
gtf <- tidyr::separate(gtf, V9, sep = ";", into = c("gene_id","gene_version","transcript_id","transcript_version","gene_type","gene_name","transcript_type","transcript_name","transcript_support","hgnc_id","havana_gene","havana_transcript"))
genes <- gsub(" gene_name ","", gtf$gene_name)
genes <- gsub("\"",replacement = "", genes) %>% sort() %>% unique()
chry_features <- genes

# raw counts are in the "counts" slot
DefaultAssay(multi) <- "RNA"
rna_counts <- GetAssayData(multi, slot = "counts")

# compute total number of rna counts per cell
# the RNA assay transcript counts will be referred to as "counts" and the ATAC assay fragments will be referred to as "frags"
total_rna_counts <- colSums(rna_counts)
total_rna_counts <- data.frame(total_rna_counts = total_rna_counts, barcode = names(total_rna_counts))

# subset for chrY counts
chry_rna_counts <- rna_counts[rownames(rna_counts) %in% chry_features,]
chry_rna_counts <- chry_rna_counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "feature") %>%
  pivot_longer(cols = colnames(chry_rna_counts)) %>%
  dplyr::rename(barcode = name)
chry_rna_counts$gem <- str_split(chry_rna_counts$barcode, pattern = "-", simplify=TRUE)[,2]
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
  distinct(barcode, chrom_rna_counts, log_rna_counts, gem, total_rna_counts)

# density plots
dir.create(here("ckd","multi_aggr_prep","plots","loy"))
plots <- here("ckd","multi_aggr_prep","plots","loy")
pdf(here(plots,"rna_ydensity.pdf"))
chry_rna_counts %>%
  ggplot(aes(log_rna_counts, color=gem, fill=gem)) +
  geom_density(aes(y = ..scaled..)) +
  facet_wrap(~gem) 
dev.off()

# read in cytoband info for hg38
# to process snATAC fragments
cytoband.gr <- fread(here("reference","cytoBand.txt.gz")) %>%
  dplyr::rename(seqnames = V1, start = V2, end = V3, cytoband = V4) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)
cytoband.gr <- keepStandardChromosomes(cytoband.gr, pruning.mode = "coarse")

# create df of args to pass to preprocessing routine
library_ids <- c("1-27Nx","2-15Nx","AJDL105","A2","090922Nx","091422Nx","AIIM164","AIL5160","AJDV174")
malesex <- c(1,0,0,1,1,1,0,0,1)
gem <- seq(1,9)
args <- data.frame(library_id = library_ids, malesex=malesex, gem = gem)

atac_frags <- lapply(seq(nrow(args)), function(row) {
  library_id <- args[row,]$library_id
  malesex <- args[row,]$malesex
  print(library_id)
  df <- readRDS(here("cellranger_multi_counts", library_id,"outs","cytoband_counts10k.rds")) %>%
    dplyr::mutate(bin = paste0(seqnames,"_",start,"_",end)) %>%
    dplyr::select(seqnames, start, end, bin, contains("cell"))
  colnames(df) <- gsub("cell",paste0("cell_",library_id), colnames(df))
  
  # annotate cytobands
  df <- as_granges(df)
  df <- join_overlap_left(df, cytoband.gr)
  df <- as.data.frame(df) 
  barcode_cols <- colnames(df)[grepl("cell",colnames(df))]
  df <- pivot_longer(df, all_of(barcode_cols), names_to = "key", values_to = "counts")
  df <- df %>% 
    dplyr::mutate(cytoband = paste0(seqnames, substr(cytoband,1,1)))
  df$library_id <- library_id
  df$malesex <- malesex
 
  # count bins
  df <- as.data.table(df)
  df <- df[, total_atac_frags:=sum(counts), by=list(key)]
  df <- df[, chrom_atac_frags:=sum(counts), by=list(key, seqnames)]
  df <- df[, cytoband_atac_frags:=sum(counts), by=list(key, cytoband)]
  
  # remove bin designations and keep collapsed barcode counts 
  df <- dplyr::distinct(df, key, seqnames, cytoband, library_id, malesex, total_atac_frags, chrom_atac_frags, cytoband_atac_frags) 
  return(df)
}) %>% dplyr::bind_rows() 

# normalize the counts
atac_frags <- atac_frags %>%
  dplyr::mutate(log_atac_frags = log10(chrom_atac_frags / total_atac_frags + 1)) 
atac_frags$key <- gsub("cell_","",atac_frags$key)

# add library_id to chry rna counts df
chry_rna_counts$gem <- as.numeric(chry_rna_counts$gem)
chry_rna_counts <- left_join(chry_rna_counts, args, by = "gem")
chry_rna_counts <- chry_rna_counts %>%
  dplyr::mutate(bc_prefix = str_split(barcode, pattern="-",simplify=TRUE)[,1]) %>%
  dplyr::mutate(key = paste0(library_id,".",bc_prefix,".1")) %>%
  dplyr::mutate(key = gsub("-",".",key)) %>%
  distinct(key, chrom_rna_counts, log_rna_counts, barcode)

# join the atac fragment counts with the rna transcript counts
multi_counts <- chry_rna_counts %>%
  left_join(atac_frags, by = "key") %>%
  na.omit() %>%
  dplyr::select(key, barcode, log_rna_counts, seqnames, library_id, malesex, log_atac_frags) %>%
  dplyr::distinct()

# density plot for atac frags
pdf(here(plots, "atac_ydensity.pdf"))
multi_counts %>%
  dplyr::filter(seqnames == "chrY") %>%
  ggplot(aes(log_atac_frags, color=library_id, fill=library_id)) +
  geom_density(aes(y = ..scaled..)) +
  facet_wrap(~library_id) +
  xlim(c(-0.5,5))
dev.off()

# scatter and density plots for rna counts and atac frags
pdf(here(plots, "multi_ydensity.pdf"))
multi_counts %>%
  na.omit() %>%
  dplyr::filter(seqnames == "chrY") %>%
  ggplot(aes(log_rna_counts, log_atac_frags)) +
  geom_point() +
  facet_wrap(~library_id) +
  xlim(c(-0.5,5)) +
  ylim(c(-0.5,5))

# simple 2d density plot for rna counts and atac frags in male donors
multi_counts %>%
  dplyr::filter(seqnames == "chrY", malesex == 1) %>%
  ggplot(aes(log_rna_counts, log_atac_frags)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", color="orange", adjust = 2) +
  facet_wrap(~library_id) +
  xlim(c(-0.5,3)) +
  ylim(c(-0.5,3)) +
  ggtitle("Male only")
dev.off()

# add cell type annotation to visualize loy across cell types
celltype <- multi@meta.data %>%
  rownames_to_column(var = "barcode") %>%
  dplyr::select(barcode, celltype)
multi_counts <- multi_counts %>%
  left_join(celltype, by = "barcode")

# visualize median chrY RNA counts and ATAC frags across celltype
# note that different cell types have variable expression / accessibility levels
# that need to be corrected prior to scaling
# chrY is a special case where we will only consider the male samples to compute a global_median (because female samples are not informative)
pdf(here(plots, "celltype_chry.pdf"))
  toplot <- multi_counts %>%
    group_by(seqnames) %>%
    filter(seqnames == "chrY") %>%
    mutate(global_median = median(log_rna_counts[malesex == 1])) %>%
    group_by(library_id, celltype) %>%
    summarize(median=median(log_rna_counts), global_median = global_median, malesex = malesex) %>%
    distinct()
  global_median <- unique(toplot$global_median)
  toplot %>%
    ggplot(aes(celltype, median, fill=as.factor(malesex))) + 
    geom_boxplot() +
    geom_hline(yintercept = as.numeric(global_median), color="red")

  toplot <- multi_counts %>%
    group_by(seqnames) %>%
    filter(seqnames == "chrY") %>%
    mutate(global_median = median(log_atac_frags[malesex == 1])) %>%
    group_by(library_id, celltype) %>%
    summarize(median=median(log_atac_frags), global_median = global_median, malesex = malesex) %>%
    distinct()
  global_median <- unique(toplot$global_median)
  toplot %>%
    ggplot(aes(celltype, median, fill=as.factor(malesex))) + 
    geom_boxplot() +
    geom_hline(yintercept = as.numeric(global_median), color="red")

  toplot <- multi_counts %>%
    group_by(seqnames) %>%
    filter(seqnames == "chr7") %>%
    mutate(global_median = median(log_atac_frags)) %>%
    group_by(library_id, celltype) %>%
    summarize(median=median(log_atac_frags), global_median = global_median) %>%
    distinct()
  global_median <- unique(toplot$global_median)
  toplot %>%
    ggplot(aes(celltype, median)) + 
    geom_boxplot() +
    geom_hline(yintercept = as.numeric(global_median), color="red")
dev.off()

# compute global median for each chromosome
# for chrY only consider male samples to compute the median
# for chrX only consider female samples to compute median
multi_counts <- multi_counts %>%
  group_by(seqnames) %>%
  mutate(global_median_counts = case_when(seqnames == "chrY" ~ median(log_rna_counts[malesex == 1]),
                                          seqnames == "chrX" ~ median(log_rna_counts[malesex == 0]),
                                          TRUE ~ median(log_rna_counts)),
         global_median_frags = case_when(seqnames == "chrY" ~ median(log_atac_frags[malesex == 1]),
                                         seqnames == "chrX" ~ median(log_atac_frags[malesex == 0]),
                                                                          TRUE ~ median(log_atac_frags)))         
# compute median for each celltype
# for chrY only consider male samples to compute the median
multi_counts <- multi_counts %>%
  group_by(seqnames, celltype) %>%
  mutate(celltype_median_counts = case_when(seqnames == "chrY" ~ median(log_rna_counts[malesex == 1]),
                                            seqnames == "chrX" ~ median(log_rna_counts[malesex == 0]),
                                            TRUE ~ median(log_rna_counts)),
         celltype_median_frags = case_when(seqnames == "chrY" ~ median(log_atac_frags[malesex == 1]),
                                           seqnames == "chrX" ~ median(log_atac_frags[malesex == 0]),
                                           TRUE ~ median(log_atac_frags)))   

# compute celltype correction factor
# and correct the counts
multi_counts <- multi_counts %>%
  group_by(seqnames, celltype) %>%
  mutate(correction_factor_counts = global_median_counts / celltype_median_counts,
         correction_factor_frags = global_median_frags / celltype_median_frags) %>%
  mutate(corrected_rna_counts = log_rna_counts * correction_factor_counts,
         corrected_atac_frags = log_atac_frags * correction_factor_frags)

# visualize uncorrected and corrected counts and frags for chrY and chr7
# note that there is a lot less variability in the uncorrected atac frags
pdf(here(plots, "celltype_corrected_chry.pdf"))
  multi_counts %>%
    filter(seqnames == "chrY", malesex == 1) %>%
    ggplot(aes(celltype, log_rna_counts)) + 
    geom_boxplot() 

  multi_counts %>%
    filter(seqnames == "chrY", malesex == 1) %>%
    ggplot(aes(celltype, corrected_rna_counts)) + 
    geom_boxplot() 

  multi_counts %>%
    filter(seqnames == "chrY", malesex == 1) %>%
    ggplot(aes(celltype, log_atac_frags)) + 
    geom_boxplot() 

  multi_counts %>%
    filter(seqnames == "chrY", malesex == 1) %>%
    ggplot(aes(celltype, corrected_atac_frags)) + 
    geom_boxplot() 

  multi_counts %>%
    filter(seqnames == "chr7", malesex == 1) %>%
    ggplot(aes(celltype, log_atac_frags)) + 
    geom_boxplot()

  multi_counts %>%
    filter(seqnames == "chr7", malesex == 1) %>%
    ggplot(aes(celltype, corrected_atac_frags)) + 
    geom_boxplot() 
dev.off()

# scale corrected rna counts and atac frags by chromosome and library_id
# to correct for donor-specific variability
multi_counts <- multi_counts %>%
  group_by(seqnames, library_id) %>%
  dplyr::mutate(scaled_atac_frags = scale(corrected_atac_frags, center=0)) %>%
  dplyr::mutate(scaled_rna_counts = scale(corrected_rna_counts, center=0))

# 2d density plot for rna counts and atac frags in male donors by cell type
pdf(here(plots,"multi_ydensity_celltype.pdf"))
multi_counts %>%
  dplyr::filter(seqnames == "chrY", malesex == 1) %>%
  ggplot(aes(scaled_rna_counts, scaled_atac_frags)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", color="orange", adjust = 2) +
  facet_wrap(~celltype) +
  xlim(c(-5,5)) +
  ylim(c(-5,5))
dev.off()

# use finite mixture models to classify LOY vs XY in male samples
mc <- multi_counts %>%
  dplyr::filter(seqnames == "chrY", malesex == 1)
model.df <- data.frame(counts = mc$scaled_rna_counts, frags = mc$scaled_atac_frags)

# semi-supervised GMM where the cells at the origin are designated as LOY (class=1)
# and cells greater than median are designated XY (class=2)
# remaining cells are designated NA and the GMM assigns them to a cluster and estimates uncertainty
# with 2 components for LOY vs. XY genotypes
model.df <- model.df %>%
  dplyr::mutate(class = ifelse(counts == 0 & frags == 0, 1, NA)) %>%
  dplyr::mutate(class = ifelse(counts > median(counts), 2, class)) %>%
  dplyr::mutate(class = ifelse(frags > median(frags), 2, class))

gmm1 <- MclustSSC(data = model.df[,1:2],
               class = model.df$class,
               G = 2, 
               modelNames = "VII")
saveRDS(gmm1, here("ckd","multi_aggr_prep","step3_gmm.rds"))
summary(gmm1)

# visualize classification for gmm
toplot <- data.frame(class = gmm1$classification,
                     frags = gmm1$data[,2],
                     counts = gmm1$data[,1],
                     label = gmm1$class)
toplot$class <- as.factor(toplot$class)
toplot <- toplot %>%
  dplyr::mutate(label = ifelse(is.na(label), "NA", label))
toplot$label <- as.factor(toplot$label)
levels(toplot$label) <- c("LOY","XY","NA")

# sort for plot order
toplot <- toplot %>%
  dplyr::arrange(desc(label))

pdf(here(plots, "loy_gmm.pdf"))
plot(gmm1, what = "uncertainty", cex=0.1, xlim = c(0,5), ylim = c(0,5), main="uncertainty")
toplot %>%
  ggplot(aes(frags, counts, color=class, fill=class)) +
  geom_point(size = 0.1) +
  facet_wrap(~class, nrow=2) +
  xlim(c(0,5)) +
  ylim(c(0,5)) +
  coord_fixed() +
  ggtitle("Final LOY vs XY classification")

# visualize class labels prior to modeling (3 are cells with NA label prior to model)
toplot %>%
  ggplot(aes(frags, counts, color=label, fill=label)) +
  geom_point(size = 0.1) +
  facet_wrap(~class, nrow=2) +
  xlim(c(0,5)) +
  ylim(c(0,5)) +
  coord_fixed() +
  ggtitle("Labels prior to modeling LOY vs XY classification")

# zoom into the boundary between XY and LOY classes
toplot %>%
  ggplot(aes(frags, counts, color=class, fill=class)) +
  geom_point(size = 0.1) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  coord_fixed() +
  ggtitle("Final LOY vs XY classification")
dev.off()

#####################
# use a NB mixture model 
# similar results but not as clean as gmm
# discretize the continous scaled counts using a scale factor and pseudocount
nbm <- data.frame(round(10000 * model.df[,1:2] + 1))

# run NB model using 2 components
nbm1 <- NB.MClust(nbm, K = 2)
saveRDS(nbm1, here("ckd","multi_aggr_prep","step3_nbm.rds"))

# visualize NB classification
toplot <- data.frame(class = nbm1$cluster,
                     frags=nbm$frags,
                     counts=nbm$counts,
                     pval=nbm1$parameters$posterior[1,])

# assign an alpha based on posterior prob < 0.05 for LOY membership
toplot$alpha <-ifelse(toplot$pval < 0.05, 1, 0.5)

toplot$class <- as.factor(toplot$class)
pdf(here(plots, "loy_nbm.pdf"))
toplot %>%
  ggplot(aes(frags, counts, color=class, fill=class)) +
  geom_point() +
  facet_wrap(~class)

# zoom in to LOY 
toplot %>%
  ggplot(aes(frags, counts, color=class, fill=class, alpha=alpha)) +
  geom_point() +
  xlim(c(0,10000)) +
  ylim(c(0,10000)) +
  facet_wrap(~class)
dev.off()

# add the loy annotation to the seurat obj. Use the pdf files to assign the correct cluster or
# classification for LOY (ie. may not always be class #1)
# male samples are assigned a genotype of LOY or XY
model.df$nbm_genotype <- ifelse(nbm1$cluster == 1, "LOY", "XY")
model.df$gmm_genotype <- ifelse(gmm1$classification == 1, "LOY", "XY")

# compare agreement between the two models
table(nbm=model.df$nbm_genotype, gmm=model.df$gmm_genotype)

# add metadata
meta <- multi@meta.data %>% rownames_to_column(var = "barcode")
addmeta <- data.frame(barcode=mc$barcode, nbm_genotype=model.df$nbm_genotype, gmm_genotype=model.df$gmm_genotype)
meta <- left_join(meta, addmeta, by = "barcode")

# add sex anno
meta <- meta %>%
  left_join(args[,1:2], by = "library_id")

rownames(meta) <- meta$barcode
multi@meta.data <- meta

# save the object with loy annotation from both models
saveRDS(multi@meta.data, here("ckd","multi_aggr_prep","step3_meta_multi_loy.rds"))

# add additional female meta data and age
age_meta <- read.xlsx(here("ckd","clinical_meta.xlsx")) %>% as.data.frame() %>% dplyr::select(library_id, age)
meta <- multi@meta.data %>% left_join(age_meta, by = "library_id")
rownames(meta) <- rownames(multi@meta.data)
multi@meta.data <- meta
saveRDS(multi, here("ckd","multi_aggr_prep","step3_multi_loy.rds"), compress=FALSE)

# join the LOY annotation to the chromosome counts and save in a separate file
# this file only has chrY rna counts, but has atac frag counts for all chroms
addmeta <- data.frame(barcode=rownames(multi@meta.data), nbm_genotype=multi@meta.data$nbm_genotype, gmm_genotype=multi@meta.data$gmm_genotype)
multi_counts <- left_join(multi_counts, addmeta, by="barcode")
saveRDS(multi_counts, here("ckd", "multi_aggr_prep", "step3_multi_counts.rds"), compress=FALSE)

