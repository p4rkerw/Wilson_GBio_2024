#!/bin/bash
# this script will aggregate snAtacseq output using cellranger-atac 

# # to run locally:
# SCRATCH1=/mnt/c/scratch
# reference=/mnt/c/reference
# docker run -it \
# --workdir $HOME \
# -v $HOME:$HOME \
# -v /mnt/c/scratch:$HOME/counts \
# -v /mnt/g/ckd:$HOME/project \
# -v /mnt/c/reference:$HOME/reference \
# -v /mnt/g/software/ckd/cnv/workflow/renal_carcinoma_pbmc/cellranger:$HOME/aggcsv_dir \
# -v $reference:$HOME/reference \
# -v $SCRATCH1:$SCRATCH1 \
# -e HOME=$HOME \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/cellranger-atac:2.1.0 /bin/bash 

# read in positional args
id=cellranger_atac_aggr_rccleuk
csv=atac_aggr_rccleuk.csv

# create path variables
inputdir=$(pwd)/counts
outputdir=$(pwd)/project
aggcsv_dir=$(pwd)/aggcsv_dir
reference=$(pwd)/reference

# change directory
change_directory () {
  cd $1
}

# set up output dir
mkdir -p $outputdir

# do computation in scratch nvme
change_directory $inputdir

# update file paths in aggregation csv on the fly with sed to match docker mount path
# csv file takes absolute paths 
sed "s|cellranger_atac_counts_dir|${inputdir}|g" $aggcsv_dir/$csv > /tmp/atac_aggr.csv

# aggregate the output without depth normalization
cellranger-atac aggr \
--id=$id \
--csv=/tmp/atac_aggr.csv \
--reference=$reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--nosecondary \
--normalize=none \
--localcores=8 \
--localmem=128
