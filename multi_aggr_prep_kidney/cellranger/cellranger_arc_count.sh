#! /bin/bash
# this script will count a 10x genomics multiome gex and atac library

# to run LSF detached
# library_id=2-15Nx

# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/reference:$HOME/reference \
# $STORAGE1/cellranger_multi_counts:$HOME/cellranger_multi_counts \
# $STORAGE1/multiFastq:$HOME/multiFastq \
# $SCRATCH1:$SCRATCH1"
# bsub -G compute-parkerw \
# -J "cellranger_multi_count_${library_id}" \
# -R 'rusage[mem=128GB]' \
# -n20 \
# -q general \
# -o $SCRATCH1/log.crmulti.$library_id.out \
# -a 'docker(p4rkerw/cellranger-arc:2.0)' \
# bash $SCRATCH1/ckd/multiomes/cellranger/cellranger_arc_count.sh $library_id

# to run interactive:
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/reference:$HOME/reference \
# $STORAGE1/cellranger_multi_counts:$HOME/cellranger_multi_counts \
# $STORAGE1/multiFastq:$HOME/multiFastq \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=128GB]' -n20 -q general-interactive -a 'docker(p4rkerw/cellranger-arc:2.0)' /bin/bash

# SCRATCH1=/mnt/g/scratch
# docker run -it --rm \
# --workdir $HOME \
# -v /mnt/g/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# -v /mnt/g/cellranger_rna_counts:$HOME/cellranger_rna_counts \
# -v /mnt/g/cellranger_multi_counts:$HOME/cellranger_multi_counts \
# -v /mnt/g/multiFastq:$HOME/multiFastq \
# -v /mnt/g/ckd:$HOME/ckd \
# -v /mnt/g/reference:$HOME/reference \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# -v /mnt/g/scratch:$HOME/scratch \
# p4rkerw/cellranger-arc:2.0 /bin/bash

reference=$(pwd)/reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
libraries=$(pwd)/cellranger_multi_counts/$1.csv

change_directory() {
  cd $1
}

change_directory cellranger_multi_counts
/opt/cellranger-arc-2.0.1/cellranger-arc count --id=$1 \
                       --reference=$reference \
                       --libraries=$libraries \
                       --localcores=20 \
                       --localmem=128
