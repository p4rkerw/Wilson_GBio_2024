#!/bin/bash
# this script will detect doublets in individual snATACseq 10X libraries

# to run locally:
# SCRATCH1=/mnt/g/scratch
# docker run \
# --workdir $HOME \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -v /mnt/g/cellranger_atac_counts/version_2.1:$HOME/counts \
# -v /mnt/g/scratch:$HOME/output \
# -e SCRATCH1="/mnt/g/scratch" \
# -it p4rkerw/amulet:1.0 /bin/bash

# to run LSF interactive
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/cellranger_atac_counts/version_2.1:$HOME/counts \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=16GB]' -q general-interactive -a 'docker(p4rkerw/amulet:1.0)' /bin/bash

# to run LSF detached
# bgadd -L 10 /parkerw/amulet
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/cellranger_atac_counts/version_2.1:$HOME/counts \
# $SCRATCH1:$SCRATCH1"
# # library_ids=(Control_1 Control_2 Control_3 Control_4 Control_5 CKD_1 CKD_2 CKD_3 CKD_4 CKD_5)
# library_ids=(Control_6 DN_1 DN_2 DN_3 DN_4 DN_5 DN_6 DN_7 SAMN18736215 SAMN18736216 SAMN27505541 SAMN27505542 SAMN27505543 SAMN27505544)
# modality=atac
# for library_id in ${library_ids[*]}; do
# bsub -G compute-parkerw \
# -R 'rusage[mem=16GB]' \
# -g /parkerw/amulet \
# -q general \
# -o $SCRATCH1/log.amulet.$library_id.$modality.out \
# -a 'docker(p4rkerw/amulet:1.0)' \
# bash $SCRATCH1/ckd/atac_aggr_prep/step0_amulet.sh $library_id
# done

# update path o/w LSF environ wont see correct python version
export PATH=/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH

library_id=$1
outputdir=counts/$library_id/amulet
mkdir -p $outputdir
bash /opt/AMULET/AMULET.sh \
counts/$library_id/outs/possorted_bam.bam \
counts/$library_id/outs/singlecell.csv \
/opt/AMULET/human_autosomes.txt \
/opt/AMULET/ENCFF356LFX.bed \
$outputdir \
/opt/AMULET \
--forcesorted
