#!/bin/bash
# this script will aggregate snAtacseq output using cellranger-atac 

# to run LSF detached
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/cellranger_atac_counts/version_2.1:$HOME/counts \
# $STORAGE1/ckd:$HOME/project \
# $SCRATCH1/ckd/:$HOME/github_repository \
# $STORAGE1/reference:$HOME/reference \
# $SCRATCH1:$SCRATCH1"
# bsub -G compute-parkerw \
# -J "cellranger_atac_aggr" \
# -R 'rusage[mem=128GB]' \
# -n20 \
# -q general \
# -o $SCRATCH1/log.cratacagg.out \
# -a 'docker(p4rkerw/cellranger-atac:2.1.0)' \
# bash $SCRATCH1/ckd/atac_aggr_prep/cellranger/cellranger_atac_aggr.sh $id $csv

# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/cellranger_atac_counts/version_2.1:$HOME/counts \
# $STORAGE1/ckd:$HOME/project \
# $SCRATCH1/ckd:$HOME/github_repository \
# $STORAGE1/reference:$HOME/reference \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=128GB]' -q general-interactive -a 'docker(p4rkerw/cellranger-atac:2.1.0)' /bin/bash

# # to run locally:
# SCRATCH1=/mnt/g/scratch
# reference=/mnt/g/reference
# docker run -it \
# --workdir $HOME \
# -v $HOME:$HOME \
# -v /mnt/g/cellranger_atac_counts/version_2.1:$HOME/counts \
# -v /mnt/g/ckd:$HOME/project \
# -v /mnt/g/reference:$HOME/reference \
# -v /mnt/g/github_repository/ckd:$HOME/github_repository \
# -v $reference:$HOME/reference \
# -v $SCRATCH1:$SCRATCH1 \
# -e HOME=$HOME \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/cellranger-atac:2.1.0 /bin/bash 

# read in positional args
id=$1 #eg. cellranger_atac_aggr_24
csv=$2 #eg. atac_aggr_22.csv

# create path variables
inputdir=$(pwd)/counts
outputdir=$(pwd)/project
github_repository=$(pwd)/github_repository
reference=$(pwd)/reference

# change directory
change_directory () {
  cd $1
}

# set up output dir
mkdir -p $outputdir
change_directory $outputdir

# update file paths in aggregation csv on the fly with sed to match docker mount path
# csv file takes absolute paths 
sed "s|cellranger_atac_counts_dir|${inputdir}|g" $github_repository/atac_aggr_prep/cellranger/$csv > /tmp/atac_aggr.csv

# aggregate the output without depth normalization
cellranger-atac aggr \
--id=$id \
--csv=/tmp/atac_aggr.csv \
--reference=$reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--nosecondary \
--normalize=none \
--localcores=16 \
--localmem=128
