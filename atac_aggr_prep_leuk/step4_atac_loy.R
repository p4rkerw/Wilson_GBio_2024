# call LOY barcodes for the atac samples and write an annotated seurat obj

# RUN local interactive
# SCRATCH1=/mnt/g/scratch
# docker run -it --rm \
# --workdir $HOME \
# -v /mnt/c/reference:$HOME/reference \
# -v /mnt/g/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# -v /mnt/g/ckd:$HOME/ckd \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/c/scratch" \
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
library(openxlsx)
set.seed=1234

plots <- here("ckd","atac_aggr_prep_rccleuk","plots","loy")
dir.create(here("ckd","atac_aggr_prep_rccleuk","plots","loy"))

# add sex annotation
aggcsv <- read.csv(here("ckd","cellranger_atac_aggr_rccleuk","outs","aggregation_csv.csv")) %>%
  dplyr::select(library_id)
clinical_meta <- read.xlsx(here("ckd","atac_aggr_prep_rccleuk","clinical_meta.xlsx"))
aggcsv <- left_join(aggcsv, clinical_meta, by = "library_id")
aggcsv$gem <- seq(nrow(aggcsv))

# cytoband coords
cytoband.gr <- fread(here("reference","cytoBand.txt.gz")) %>%
  dplyr::rename(seqnames = V1, start = V2, end = V3, cytoband = V4) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)
cytoband.gr <- keepStandardChromosomes(cytoband.gr, pruning.mode = "coarse")

# aggregate the atac chromosome fragment counts
atac_frags <- lapply(seq(nrow(aggcsv)), function(row) {
  library_id <- aggcsv[row,]$library_id
  sex <- aggcsv[row,]$sex
  print(library_id)
  df <- readRDS(here("cellranger_atac_counts", "leukocytes", library_id,"outs","cytoband_counts10k.rds")) %>%
    dplyr::mutate(bin = paste0(seqnames,"_",start,"_",end)) %>%
    dplyr::select(seqnames, start, end, bin, contains("cell"))
  colnames(df) <- gsub("cell",paste0("cell_",library_id), colnames(df))
  
  # annotate chromosome arm cytobands
  df <- as_granges(df)
  df <- join_overlap_left(df, cytoband.gr)
  df <- as.data.frame(df) 
  barcode_cols <- colnames(df)[grepl("cell",colnames(df))]
  df <- pivot_longer(df, all_of(barcode_cols), names_to = "key", values_to = "counts")
  df <- df %>% 
    dplyr::mutate(cytoband = paste0(seqnames, substr(cytoband,1,1)))
  df$library_id <- library_id
  df$sex <- sex
 
  # count bins
  df <- as.data.table(df)
  df <- df[, total_atac_frags:=sum(counts), by=list(key)]
  df <- df[, chrom_atac_frags:=sum(counts), by=list(key, seqnames)]
  df <- df[, cytoband_atac_frags:=sum(counts), by=list(key, cytoband)]
  
  # remove bin designations and keep collapsed barcode counts 
  df <- dplyr::distinct(df, key, seqnames, cytoband, library_id, sex, total_atac_frags, chrom_atac_frags, cytoband_atac_frags) 
  return(df)
}) %>% dplyr::bind_rows()

# normalize the chromosome frags by total frags per cell
atac_frags <- atac_frags %>%
  dplyr::mutate(log_atac_frags = log10(chrom_atac_frags / total_atac_frags + 1)) 

# switch key to barcode with gem suffix
atac_frags$key <- gsub("cell_","",atac_frags$key)
gem <- data.frame(gem=aggcsv$gem, library_id=aggcsv$library_id)
atac_frags <- left_join(atac_frags, gem, by = "library_id")
atac_frags <- atac_frags %>%
  tidyr::separate(key, sep="[.]", into = c(NA,"barcode_prefix",NA))
atac_frags <- atac_frags %>%
  dplyr::mutate(barcode = paste0(barcode_prefix,"-",gem))

# load cell type annotation
meta <- readRDS(here("ckd","atac_aggr_prep_rccleuk","step3_predictions_meta.rds"))

# discard any barcodes that did not meet QC from earlier steps (ie multiplets etc)
# note that the cytoband_counts10k.rds files are generated from raw cellranger-atac count output
atac_frags <- atac_frags[atac_frags$barcode %in% rownames(meta),]

# add cell type annotation
atac_frags <- atac_frags %>%
  left_join(data.frame(barcode=rownames(meta), celltype=meta$predicted.l1), by="barcode")

# do a cell-specific correction to account for variability in median atac coverage across chromosomes
# ie. some cells have more or less coverage independent of chromosomal gains and losses (ie chromatin remodeling etc)
# chrY is a special case where we will only consider the male samples to compute a global_median (because female samples are not informative)
# first visualize uncorrected fragment counts for chrY and chr7
pdf(here(plots, "celltype_chry.pdf"))
  toplot <- atac_frags %>%
    group_by(seqnames) %>%
    filter(seqnames == "chrY") %>%
    mutate(global_median = median(log_atac_frags[sex == "male"])) %>%
    group_by(library_id, celltype) %>%
    summarize(median=median(log_atac_frags), global_median = global_median, sex = sex) %>%
    distinct()
  global_median <- unique(toplot$global_median)
  toplot %>%
    ggplot(aes(celltype, median, fill=as.factor(sex))) + 
    geom_boxplot() +
    geom_hline(yintercept = as.numeric(global_median), color="red") +
    theme_bw() +
    theme(axis.text.x = element_text (angle = 90, hjust = 1))

  toplot <- atac_frags %>%
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
    geom_hline(yintercept = as.numeric(global_median), color="red") +
    theme_bw() +
    theme(axis.text.x = element_text (angle = 90, hjust = 1))

  toplot <- atac_frags %>%
    group_by(seqnames) %>%
    filter(seqnames == "chr10") %>%
    mutate(global_median = median(log_atac_frags)) %>%
    group_by(library_id, celltype) %>%
    summarize(median=median(log_atac_frags), global_median = global_median) %>%
    distinct()
  global_median <- unique(toplot$global_median)
  toplot %>%
    ggplot(aes(celltype, median)) + 
    geom_boxplot() +
    geom_hline(yintercept = as.numeric(global_median), color="red") +
    theme_bw() +
    theme(axis.text.x = element_text (angle = 90, hjust = 1))

dev.off()

# compute global median for each chromosome
# for chrY only consider male samples to compute the median
# for chrX only consider female samples to compute the median
atac_frags <- atac_frags %>%
  group_by(seqnames) %>%
  mutate(global_median_frags = case_when(seqnames == "chrY" ~ median(log_atac_frags[sex == "male"]),
                                         seqnames == "chrX" ~ median(log_atac_frags[sex == "female"]),
                                                                          TRUE ~ median(log_atac_frags)))         
# compute median for each celltype
# for chrY only consider male samples to compute the median
# for chrX only consider female samples to compute the median
atac_frags <- atac_frags %>%
  group_by(seqnames, celltype) %>%
  mutate(celltype_median_frags = case_when(seqnames == "chrY" ~ median(log_atac_frags[sex == "male"]),
                                           seqnames == "chrX" ~ median(log_atac_frags[sex == "female"]),
                                                                          TRUE ~ median(log_atac_frags)))   

# compute celltype correction factor
# and correct the counts
atac_frags <- atac_frags %>%
  group_by(seqnames, celltype) %>%
  mutate(correction_factor_frags = global_median_frags / celltype_median_frags) %>%
  mutate(corrected_atac_frags = log_atac_frags * correction_factor_frags)

# visualize uncorrected and corrected frags for chrY, chrX and chr7
# note that there is a lot less variability in the uncorrected atac frags
pdf(here(plots, "celltype_corrected_chry.pdf"))
  atac_frags %>%
    filter(seqnames == "chrY", sex == "male") %>%
    ggplot(aes(celltype, log_atac_frags)) + 
    geom_boxplot(outlier.shape = NA) +
    ggtitle("Male corrected chrY")

  atac_frags %>%
    filter(seqnames == "chrY", sex == "female") %>%
    ggplot(aes(celltype, corrected_atac_frags)) + 
    geom_boxplot(outlier.shape = NA) +
    ggtitle("Female corrected chrY")

  atac_frags %>%
    filter(seqnames == "chrX") %>%
    ggplot(aes(celltype, log_atac_frags, fill=as.factor(sex))) + 
    geom_boxplot(outlier.shape = NA) 

  atac_frags %>%
    filter(seqnames == "chrX") %>%
    ggplot(aes(celltype, corrected_atac_frags, fill=as.factor(sex))) + 
    geom_boxplot(outlier.shape = NA) 

  atac_frags %>%
    filter(seqnames == "chr7", sex == "male") %>%
    ggplot(aes(celltype, log_atac_frags)) + 
    geom_boxplot(outlier.shape = NA) 

  atac_frags %>%
    filter(seqnames == "chr7", sex == "male") %>%
    ggplot(aes(celltype, corrected_atac_frags)) + 
    geom_boxplot(outlier.shape = NA) 
dev.off()

# scale corrected rna counts and atac frags by chromosome and library_id
# to correct for donor-specific variability
atac_frags <- atac_frags %>%
  group_by(seqnames, library_id) %>%
  dplyr::mutate(scaled_atac_frags = scale(corrected_atac_frags, center=0)) 

# save the scaled atac counts to file
saveRDS(atac_frags, here("ckd","atac_aggr_prep_rccleuk","atac_frags.rds"), compress=FALSE)

# for this analysis only consider entire chromosomes (ie filter out cytoband column)
atac_frags <- atac_frags %>% dplyr::distinct(barcode, seqnames, scaled_atac_frags, log_atac_frags, sex, library_id)

# density plot for atac frags
pdf(here(plots, "atac_ydensity.pdf"))
atac_frags %>%
  dplyr::filter(seqnames == "chrY") %>%
  ggplot(aes(scaled_atac_frags, color=library_id, fill=library_id)) +
  geom_density(aes(y = ..scaled..)) +
  facet_wrap(~library_id) +
  xlim(c(-0.5,5))
atac_frags %>%
  dplyr::filter(seqnames == "chrY") %>%
  ggplot(aes(scaled_atac_frags, color=sex, fill=sex)) +
  geom_density(aes(y = ..scaled..)) +
  facet_wrap(~library_id) +
  xlim(c(-0.5,5))
atac_frags %>%
  dplyr::filter(seqnames == "chrY") %>%
  ggplot(aes(scaled_atac_frags, color=sex, fill=sex)) +
  geom_density(aes(y = ..scaled..)) +
  facet_wrap(~sex) +
  xlim(c(0,5))
dev.off()

############################################################################
# find the trough of the chrY scaled fragment kernel density estimate to call LOY vs XY genotype
kd <- atac_frags %>%
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
      geom_vline(xintercept = peak1, linetype = "dotted") +
      geom_vline(xintercept = peak2, linetype = "dotted") +
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
       dplyr::filter(seqnames == "chrY") %>%
       dplyr::mutate(total = sum(Freq)) %>%
       dplyr::mutate(prop = Freq / total) %>%
       dplyr::filter(genotype == "LOY")
  
p1 <- p +
      annotate(geom="text",
               x = 3,
               y = 0.5,
               label = paste0("LOY: ", round(res$prop, digits=2), "\n", "Cells: ", res$total)) +
      theme_bw()

pdf(here(plots, "step4_akd.pdf"))
print(p1)
dev.off()

# add the kernel density trough LOY estimate to the metadata
meta <- meta %>%
  rownames_to_column(var = "barcode") %>%
  left_join(kd[,c("barcode","genotype")], by = "barcode") %>%
  dplyr::rename(akd_genotype = genotype)

saveRDS(meta, here("ckd","atac_aggr_prep_rccleuk","step4_meta.rds"))

# annotate the object with loy, cell type and tissue source
atacAggr <- readRDS(here("ckd","atac_aggr_prep_rccleuk","step1d_prep.rds"))

# annotate loy genotypes
anno <- meta %>%
  dplyr::select(akd_genotype)
rownames(anno) <- meta$barcode
atacAggr <- AddMetaData(atacAggr, anno)

# add predicted cell types
anno <- meta %>%
  dplyr::select(predicted.l1, predicted.l2, predicted.l3) 
rownames(anno) <- meta$barcode
atacAggr <- AddMetaData(atacAggr, anno)

# add tissue source
clinical_meta <- read.xlsx(here("ckd","supplemental_data_1_clinical_data_qc.xlsx"), sheet = "sample_info") %>%
  dplyr::select(library_id, tissue_source, sex)
meta <- atacAggr@meta.data %>% left_join(tissue_source, by = "library_id")
rownames(meta) <- rownames(atacAggr@meta.data)

atacAggr@meta.data <- meta

saveRDS(atacAggr, here("ckd","atac_aggr_prep_rccleuk","step4_anno.rds"), compress=FALSE)



