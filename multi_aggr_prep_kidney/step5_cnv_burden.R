# this script will compute cnv burden for multiomes using infercnv and epianeufinder results

# SCRATCH1=/mnt/g/scratch
# docker run -it \
# --workdir $HOME \
# -v /mnt/g/ckd:$HOME/ckd \
# -v /mnt/g/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# -v /mnt/g/cellranger_multi_counts:$HOME/cellranger_multi_counts \
# -v /mnt/g/reference:$HOME/reference \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# -v /mnt/g/scratch:$HOME/scratch \
# p4rkerw/sctools:R4.1.3 R

library(data.table)
library(dplyr)
library(tidyr)
library(here)
library(stringr)
library(tidyr)
##################################
# # infercnv cnv burden
# # read gene order from reference dir
# gene_order <- readRDS(here("reference","gene_order.rds"))

# # read gene order from an object
# # obj <- readRDS(here("ckd","multi_aggr_prep","infercnv","A2","run.final.infercnv_obj"))
# # gene_order <- obj@gene_order
# # saveRDS(gene_order, here("reference","gene_order.rds"))

# # make chrom sizes df
# chrom_sizes <- as.data.frame(gene_order) %>%
#   mutate(width = stop - start) %>%
#   mutate(total_size = sum(width)) %>%
#   group_by(chr) %>%
#   mutate(chrom_size = sum(width)) %>%
#   distinct(chr, chrom_size, total_size)

# # get multiome library_ids from aggregation file
# aggcsv <- read.csv(here("ckd","cellranger_arc_aggr","outs","aggr.csv"))
# aggcsv$gem <- seq(nrow(aggcsv))
# aggcsv <- aggcsv %>% dplyr::select(library_id, gem)

# # estimate RNA assay CNV burden per chromosome using infercnv output
# library_ids <- aggcsv$library_id
# rna_cnv_burden <- lapply(library_ids, function(library_id) {
#   print(library_id)
#   df <- fread(here("ckd","multi_aggr_prep","infercnv",library_id,"map_metadata_from_infercnv.txt"))
#   df <- df %>%
#     dplyr::select(V1, contains("proportion_cnv")) %>%
#     pivot_longer(cols = contains("proportion_cnv"), names_to = "chr", values_to = "rna_cnv_chrom_burden") %>%
#     tidyr::separate(chr, sep = "_", into = c(NA, NA, "chr")) %>%
#     left_join(chrom_sizes, by = "chr") %>%
#     mutate(rna_cnv_chrom_bp = rna_cnv_chrom_burden * chrom_size) %>%
#     group_by(V1) %>%
#     mutate(rna_cnv_total_burden = sum(rna_cnv_chrom_bp) / total_size)
#   }) %>% bind_rows()
# colnames(rna_cnv_burden)[1] <- "barcode"
# rna_cnv_burden <- rna_cnv_burden %>%
#   dplyr::select(barcode, chr, rna_cnv_chrom_burden, rna_cnv_total_burden)

# # write cnv results to file
# fwrite(rna_cnv_burden, here("ckd","multi_aggr_prep","step6_rna_cnv_burden.tsv"))

#############################
# get multiome library_ids from aggregation file
aggcsv <- read.csv(here("ckd","cellranger_arc_aggr","outs","aggr.csv"))
aggcsv$gem <- seq(nrow(aggcsv))
aggcsv <- aggcsv %>% dplyr::select(library_id, gem)

# epianeufinder cnv burden
# read in cnv results from epianeufinder
atac_cnv_burden <- lapply(seq_along(aggcsv$library_id), function(gem){
  library_id <- aggcsv$library_id[gem]

  print(library_id)
  df <- read.table(here("cellranger_multi_counts",library_id, "outs","epiAneufinder_1MB",
                        "epiAneufinder_results", "results_table.tsv"))
 
  # update keys to add gem suffix
  keyindex <- grepl("cell", colnames(df))
  keys <- colnames(df)[keyindex]
  bcprefix <- str_split(keys, pattern="[.]", simplify=TRUE)[,2]
  barcodes <- paste0(bcprefix,"-",gem)
  
  # update column names to barcodes with gem suffix
  colnames(df)[keyindex] <- barcodes
  
  # quantify prop cnv per chromosome and total cnv proportion
  # values in data frame are "Normal" = '1', "Loss" = '0', "Gain" = '2'
  # loss and gain are grouped together to quantify proportion
  res <- df %>% 
    select(seq, all_of(barcodes)) %>%
    as.data.frame() %>%
    pivot_longer(cols = all_of(barcodes), names_to = "barcode") %>%
    dplyr::mutate(cnv_logical = ifelse(value != 1, 1, 0)) %>%
    dplyr::group_by(barcode) %>%
    dplyr::mutate(total_bins = n()) %>%
    group_by(seq, barcode) %>%
    dplyr::mutate(chrom_bins = n()) %>%
    dplyr::mutate(atac_cnv_chrom_burden = sum(cnv_logical) / chrom_bins) %>%
    group_by(barcode) %>%
    dplyr:::mutate(atac_cnv_total_burden = sum(cnv_logical) / total_bins) %>%
    ungroup() %>%
    dplyr::select(seq, barcode, atac_cnv_chrom_burden, atac_cnv_total_burden) %>%
    distinct()
    
  return(res)
}) %>% bind_rows()

# write cnv results to file
atac_cnv_burden <- atac_cnv_burden %>% dplyr::rename(chr = seq) %>% dplyr::select(-library_id) %>% distinct()
fwrite(atac_cnv_burden, here("ckd","multi_aggr_prep","step5_atac_cnv_burden.tsv"))

