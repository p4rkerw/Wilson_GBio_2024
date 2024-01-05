# SCRATCH1=/mnt/g/scratch
# docker run -it \
# --workdir $HOME \
# -v /mnt/g/ckd:$HOME/project \
# -v /mnt/g/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# -v /mnt/g/reference:$HOME/reference \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# -v /mnt/g/scratch:$HOME/scratch \
# p4rkerw/sctools:R4.1.3 /bin/bash

# # to run RIS interactive
# # # export docker volumes
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# $STORAGE1/ckd:$HOME/project \
# $STORAGE1/reference:$HOME/reference \
# $SCRATCH1:$SCRATCH1"
# # to run interactive
# bsub -Is -G compute-parkerw -R 'rusage[mem=64GB]' -sp 99 -n 20 -q general-interactive -a 'docker(p4rkerw/sctools:R4.1.3)' /bin/bash

# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# $STORAGE1/ckd:$HOME/project \
# $STORAGE1/reference:$HOME/reference \
# $SCRATCH1:$SCRATCH1"
# bgadd -L 20 /parkerw/epiAneufinder
# # library_ids=(Control_1 Control_2 Control_3 Control_4 Control_5 Control_6 CKD_1 CKD_2 CKD_3 CKD_4 CKD_5 DN_1 DN_2 DN_3 DN_4 DN_5 DN_6 DN_7)
# library_ids=(SAMN18736215 SAMN18736216 SAMN27505541 SAMN27505542 SAMN27505543 SAMN27505544)
# for library_id in ${library_ids[@]}; do
# bsub -G compute-parkerw \
# -g /parkerw/epiAneufinder \
# -R 'rusage[mem=64GB]' \
# -n 10 \
# -q general \
# -a 'docker(p4rkerw/sctools:R4.1.3)' \
# -o $SCRATCH1/log.epianeu.$library_id.out \
# Rscript $SCRATCH1/ckd/atac_aggr_prep/step5_epianeufinder.R $library_id
# done

args <- commandArgs(trailingOnly = TRUE)
print(args)
library_id <- args[1]

library(epiAneufinder)
library(here)
library(BSgenome.Hsapiens.UCSC.hg38)

input <- here("cellranger_atac_counts","version_2.1", library_id,"outs","fragments.tsv.gz")
output <- here("cellranger_atac_counts","version_2.1",library_id,"outs","epiAneufinder_1MB")
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
