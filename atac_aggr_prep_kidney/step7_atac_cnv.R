# this script will compare the cnv burden between celltype annotations in the aggregated atac object

library(data.table)
library(dplyr)
library(tidyr)
library(here)
library(stringr)
library(tidyr)

# get library_ids from aggregation file
aggcsv <- read.csv(here("ckd","cellranger_atac_aggr","outs","aggregation_csv.csv"))
aggcsv$gem <- seq(nrow(aggcsv))
aggcsv <- aggcsv %>% dplyr::select(library_id, gem)

# epianeufinder cnv burden
# read in cnv results from epianeufinder
atac_cnv_burden <- lapply(seq_along(aggcsv$library_id), function(gem){
  library_id <- aggcsv$library_id[gem]

  print(library_id)
  df <- read.table(here("cellranger_atac_counts","version_2.1", library_id, "outs","epiAneufinder_1MB",
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
atac_cnv_burden <- atac_cnv_burden %>% dplyr::rename(chr = seq) %>% distinct()
fwrite(atac_cnv_burden, here("ckd","atac_aggr_prep","step7_atac_cnv_burden.tsv"))



