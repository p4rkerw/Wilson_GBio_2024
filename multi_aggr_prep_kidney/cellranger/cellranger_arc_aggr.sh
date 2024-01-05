#! /bin/bash
# this script will aggregate a 10x genomics multiome gex and atac library

# bgadd -L 10 /parkerw/cellranger

# to run LSF detached
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/reference:$HOME/reference \
# $STORAGE1/cellranger_multi_counts:$HOME/cellranger_multi_counts \
# $SCRATCH1/ckd:$HOME/github_repository \
# $STORAGE1/ckd:$HOME/ckd \
# $SCRATCH1:$SCRATCH1"
# bsub -G compute-parkerw \
# -g /parkerw/cellranger \
# -J "cellranger_arc_aggr" \
# -R 'rusage[mem=128GB]' \
# -n20 \
# -q general \
# -o $SCRATCH1/log.craggr.out \
# -a 'docker(p4rkerw/cellranger-arc:2.0)' \
# bash $SCRATCH1/ckd/multiomes/cellranger/cellranger_arc_aggr.sh

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

# create path variables
inputdir=$(pwd)/cellranger_multi_counts
outputdir=$(pwd)/ckd
githubdir=$(pwd)/github_repository
reference=$(pwd)/reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0

# change directory
change_directory () {
  cd $1
}

# set up output dir
mkdir -p $outputdir
change_directory $outputdir

# update file paths in aggregation csv on the fly with sed to match docker mount path
# csv file takes absolute paths 
sed "s|cellranger_counts_dir|${inputdir}|g" $githubdir/multiomes/cellranger/multi_aggr.csv > /tmp/aggr.csv

/opt/cellranger-arc-2.0.1/cellranger-arc aggr \
                       --id=cellranger_arc_aggr \
                       --reference=$reference \
                       --csv=/tmp/aggr.csv \
                       --localcores=20 \
                       --localmem=128 \
                       --normalize=none
