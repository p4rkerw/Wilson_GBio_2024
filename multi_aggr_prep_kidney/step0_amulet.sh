#!/bin/bash
# this script will detect doublets in individual snATACseq 10X libraries

# to run locally:
# SCRATCH1=/mnt/g/scratch
# docker run \
# --workdir $HOME \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -v /mnt/g/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# -v /mnt/g/cellranger_multi_counts:$HOME/cellranger_multi_counts \
# -v /mnt/g/scratch:$HOME/scratch \
# -e SCRATCH1="/mnt/g/scratch" \
# -it p4rkerw/amulet:1.0 /bin/bash

# to run LSF interactive
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# $STORAGE1/cellranger_multi_counts:$HOME/cellranger_multi_counts \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw \
# -R 'rusage[mem=64GB]' \
# -q general-interactive \
# -a 'docker(p4rkerw/amulet:1.0)' /bin/bash

# to run LSF detached
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# $STORAGE1/cellranger_multi_counts:$HOME/cellranger_multi_counts \
# $SCRATCH1:$SCRATCH1"
# library_ids=(1-27Nx 2-15Nx A2 AJDL105 090922Nx 091422Nx AIIM164 AIL5160 AJDV174)
# for library_id in ${library_ids[*]}; do
# bsub -G compute-parkerw \
# -R 'rusage[mem=64GB]' \
# -q general \
# -o $SCRATCH1/log.amulet.$library_id.out \
# -a 'docker(p4rkerw/amulet:1.0)' $SCRATCH1/ckd/multiomes/step0_amulet.sh $library_id
# done

# update path o/w LSF environ wont see correct python version
export PATH=/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH

# the 1st column (0 position for 0-index) contains the snATAC barcode and is identical to the value in the 2nd column (gex_barcode)
# the value in the 3rd column (is_cell) indicates whether cellranger-arc called this barcode as a cell
# these values are passed as flags to AMULET because the per_barcode_metrics.csv file has a different format than singlecell.csv
# which is output by cellranger-atac 
library_id=$1
outputdir=cellranger_multi_counts/$library_id/outs/amulet
mkdir -p $outputdir
bash /opt/AMULET/AMULET.sh \
--bcidx 0 --cellidx 0 --iscellidx 3 \
cellranger_multi_counts/$library_id/outs/atac_possorted_bam.bam \
cellranger_multi_counts/$library_id/outs/per_barcode_metrics.csv \
/opt/AMULET/human_autosomes.txt \
/opt/AMULET/ENCFF356LFX.bed \
$outputdir \
/opt/AMULET \
--forcesorted
