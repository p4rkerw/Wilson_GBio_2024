# this script will run epianeufinder with 1MB windows on the multiome atac fragments files

# SCRATCH1=/mnt/g/scratch
# docker run -it \
# --workdir $HOME \
# -v /mnt/g/ckd:$HOME/project \
# -v /mnt/g/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# -v /mnt/g/cellranger_multi_counts:$HOME/cellranger_multi_counts \
# -v /mnt/g/reference:$HOME/reference \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# -v /mnt/g/scratch:$HOME/scratch \
# p4rkerw/sctools:R4.1.3 /bin/bash

# library_ids=(1-27Nx 2-15Nx AJDL105 A2 090922Nx 091422Nx AIIM164 AIL5160 AJDV174)
# for library_id in ${library_ids[@]}; do
# Rscript $SCRATCH1/ckd/multiomes/step4_multi_epianeufinder.R $library_id
# done

args <- commandArgs(trailingOnly = TRUE)
print(args)
library_id <- args[1]

library(epiAneufinder)
library(here)
library(BSgenome.Hsapiens.UCSC.hg38)

input <- here("cellranger_multi_counts",library_id,"outs","atac_fragments.tsv.gz")
output <- here("cellranger_multi_counts",library_id,"outs","epiAneufinder_1MB")
dir.create(output, recursive=TRUE)
blacklist <- "/opt/epiAneufinder/sample_data/hg38-blacklist.v2.bed"

epiAneufinder(input=input, #Enter path to your fragments.tsv file or the folder containing bam files
              outdir=output, #Path to the directory where results should be written 
              blacklist=blacklist, #Path to bed file that contains the blacklisted regions of your genome
              windowSize=1e6, 
              genome="BSgenome.Hsapiens.UCSC.hg38", #Substitute with relevant BSgenome
              exclude=c('chrX','chrY','chrM'), 
              reuse.existing=TRUE,
              title_karyo="Karyogram of sample data", 
              ncores=10,
              minFrags=10000)
