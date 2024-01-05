#! /bin/bash

# to run LSF detached
# bgadd -L 10 /parkerw/cellranger_atac
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/cellranger_atac_counts/version_2.1:$HOME/counts \
# $STORAGE1/atacFastq:$HOME/fastq \
# $STORAGE1/reference:$HOME/reference \
# $SCRATCH1:$SCRATCH1"
# library_ids=(Control_1 Control_2 Control_3 Control_4 Control_5 CKD_1 CKD_2 CKD_3 CKD_4 CKD_5)
# modality=atac
# for library_id in ${library_ids[*]}; do
# bsub -G compute-parkerw \
# -J "cellranger_atac_count_${library_id}" \
# -R 'rusage[mem=128GB]' \
# -n20 \
# -g /parkerw/cellranger_atac \
# -q general \
# -o $SCRATCH1/log.cratac.$library_id.$modality.out \
# -a 'docker(p4rkerw/cellranger-atac:2.1.0)' \
# bash $SCRATCH1/ckd/atac_aggr_prep/cellranger/cellranger_atac_count.sh $library_id
# done

# to run interactive:
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/reference:$HOME/reference \
# $STORAGE1/ckd:$HOME/project \
# $STORAGE1/atacFastq:$HOME/fastq \
# $STORAGE1/cellranger_atac_counts/version_2.1:$HOME/counts \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=128GB]' -n20 -q general-interactive -a 'docker(p4rkerw/cellranger-atac:2.1.0)' /bin/bash

# TO RUN INTERACTIVE LOCAL
# SCRATCH1=/mnt/g/scratch
# reference=/mnt/g/reference
# fastq=/mnt/g/atacFastq
# counts=/mnt/g/cellranger_atac_counts/version_2.1
# docker run -it --rm \
# --workdir $HOME \
# -v $HOME:$HOME \
# -v $fastq:$HOME/fastq \
# -v $reference:$HOME/reference \
# -v $counts:$HOME/counts \
# -v $SCRATCH1:$SCRATCH1 \
# -e HOME=$HOME \
# -e SCRATCH1="${SCRATCH1}" \
# p4rkerw/cellranger-atac:2.1.0 /bin/bash 

# create path variables
outputdir=$(pwd)/counts
fastq=$(pwd)/fastq
reference=$(pwd)/reference

# assign library_id to positional arg
library_id=$1 # eg. Control_1
# library_id=1153-EO-1

# change directory
change_directory () {
  cd $1
}

# set up output dir
mkdir -p $outputdir
change_directory $outputdir

# count the fastq files
cellranger-atac count \
--id=$library_id \
--reference=$reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--fastqs=$fastq/$library_id \
--localcores=20 \
--localmem=128

