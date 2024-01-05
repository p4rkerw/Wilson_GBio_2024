https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181064
scRNA and scATAC in matched tumor and pbmc from donors with RCC

Download raw data files from SRA: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA750487&o=acc_s%3Aa
```
# function for downloading single cell atac from GEO
function atac_fetch_fastq() {
sample=$1
outputdir=$2
accession=$3

prefetch --type all $accession --max-size u
fastq-dump --split-files $accession/$accession.sra --gzip --outdir $accession

# strip last 3 characters from accession and use as lane index
i=$(echo ${accession: -3})

# rename files to conform to bcl2fastq format
mv $accession/${accession}_1.fastq.gz $outputdir/${sample}_S1_L${i}_R1_001.fastq.gz
mv $accession/${accession}_2.fastq.gz $outputdir/${sample}_S1_L${i}_R2_001.fastq.gz
mv $accession/${accession}_3.fastq.gz $outputdir/${sample}_S1_L${i}_R3_001.fastq.gz

# cleanup
rm -rf $accession
}
export -f atac_fetch_fastq

# function for downloading single cell RNA from GEO
function rna_fetch_fastq() {
sample=$1
outputdir=$2
accession=$3

prefetch --type all $accession --max-size u
fastq-dump --split-files $accession/$accession.sra --gzip --outdir $accession

# strip last 3 characters from accession and use as lane index
i=$(echo ${accession: -3})

# rename files to conform to bcl2fastq format
mv $accession/${accession}_1.fastq.gz $outputdir/${sample}_S1_L${i}_R1_001.fastq.gz
mv $accession/${accession}_2.fastq.gz $outputdir/${sample}_S1_L${i}_R2_001.fastq.gz

# cleanup
rm -rf $accession
}
export -f rna_fetch_fastq
```

ATAC Blood immune cell sets (n=8)
```
# pt1002300_PBMC_CD45+_cells_scATACseq
sample=SAMN20460940
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 6 --bar atac_fetch_fastq $sample $outputdir :::: /mnt/g/downloads/$sample.txt

# pt1002310_PBMC_CD45+_cells_scATACseq
sample=SAMN20460945
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 6 --bar atac_fetch_fastq $sample $outputdir :::: /mnt/g/downloads/$sample.txt

# pt2001055_PBMC_CD45+_cells_scATACseq
sample=SAMN20460948
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 6 --bar atac_fetch_fastq $sample $outputdir :::: /mnt/g/downloads/$sample.txt

# pt2001215_PBMC_CD45+_cells_scATACseq
sample=SAMN20460950
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 6 --bar atac_fetch_fastq $sample $outputdir :::: /mnt/g/downloads/$sample.txt

# pt2001077_PBMC_CD45+_cells_scATACseq
sample=SAMN20460943
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 6 --bar atac_fetch_fastq $sample $outputdir :::: /mnt/g/downloads/$sample.txt

# pt2001221_PBMC_CD45+_cells_scATACseq
sample=SAMN20460952
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 6 --bar atac_fetch_fastq $sample $outputdir :::: /mnt/g/downloads/$sample.txt

# pt7001025_PBMC_CD45+_cells_scATACseq
sample=SAMN20460954
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 6 --bar atac_fetch_fastq $sample $outputdir :::: /mnt/g/downloads/$sample.txt

# pt7001031_PBMC_CD45+_cells_scATACseq
sample=SAMN20460923
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 6 --bar atac_fetch_fastq $sample $outputdir :::: /mnt/g/downloads/$sample.txt
```

ATAC Tumor immune cell sets
```
# pt1002300_tumor_CD45+_cells_scATACseq
sample=SAMN20460941
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 6 --bar atac_fetch_fastq $sample $outputdir :::: /mnt/g/downloads/$sample.txt

# pt1002310_tumor_CD45+_cells_scATACseq
sample=SAMN20460946
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 6 --bar atac_fetch_fastq $sample $outputdir :::: /mnt/g/downloads/$sample.txt

# pt2001055_tumor_CD45+_cells_scATACseq
sample=SAMN20460949
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 6 --bar atac_fetch_fastq $sample $outputdir :::: /mnt/g/downloads/$sample.txt

# pt2001077_tumor_CD45+_cells_scATACseq
sample=SAMN20460944
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 6 --bar atac_fetch_fastq $sample $outputdir :::: /mnt/g/downloads/$sample.txt

# pt2001215_tumor_CD45+_cells_scATACseq
sample=SAMN20460951
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 6 --bar atac_fetch_fastq $sample $outputdir :::: /mnt/g/downloads/$sample.txt

# pt2001221_tumor_CD45+_cells_scATACseq
sample=SAMN20460953
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 6 --bar atac_fetch_fastq $sample $outputdir :::: /mnt/g/downloads/$sample.txt

# pt7001025_tumor_CD45+_cells_scATACseq
sample=SAMN20460955
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 6 --bar atac_fetch_fastq $sample $outputdir :::: /mnt/g/downloads/$sample.txt

# pt7001031_tumor_CD45+_cells_scATACseq
sample=SAMN20460924
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 6 --bar atac_fetch_fastq $sample $outputdir :::: /mnt/g/downloads/$sample.txt

```

ATAC Normal adjacent kidney immune cell sets
```
# pt1002300_adjNORM_CD45+_cells_scATACseq
sample=SAMN20460942
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 6 --bar atac_fetch_fastq $sample $outputdir :::: /mnt/g/downloads/$sample.txt

# pt1002310_adjNORM_CD45+_cells_scATACseq
sample=SAMN20460947
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 6 --bar atac_fetch_fastq $sample $outputdir :::: /mnt/g/downloads/$sample.txt

# pt7001025_adjNORM_CD45+_cells_scATACseq
sample=SAMN20460922
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 6 --bar atac_fetch_fastq $sample $outputdir :::: /mnt/g/downloads/$sample.txt

# pt7001031_adjNORM_CD45+_cells_scATACseq
sample=SAMN20460925
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 6 --bar atac_fetch_fastq $sample $outputdir :::: /mnt/g/downloads/$sample.txt

```
