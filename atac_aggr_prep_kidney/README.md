scripts for downloading previously published snATAC-seq kidney datasets

```
# url accessions for downloading accession files

70M control snATAC-seq kidney
GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172008
Accession: SRX10595439
Biosample: SAMN18736215
Sample: GSM5239694
https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRX10595439&o=acc_s%3Aa

52F control snATAC-seq kidney
GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172008
Accession: SRX10595438
Biosample: SAMN18736216
Sample: GSM5239693
https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRX10595438&o=acc_s%3Aa

74M control snATAC-seq kidney
GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200547
Accession: PRJNA825254
Biosample: SAMN27505542
Sample: GSM6037589
https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA825254&o=acc_s%3Aa

65F control snATAC-seq kidney
GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200547
Accession: PRJNA825254
Biosample: SAMN27505541
Sample: GSM6037591
https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA825254&o=acc_s%3Aa

85F control snATAC-seq kidney
GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200547
Accession: PRJNA825254
Biosample: SAMN27505543
Sample: GSM6037588
https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA825254&o=acc_s%3Aa

85F control snATAC-seq kidney
GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200547
Accession: PRJNA825254
Biosample: SAMN27505544
Sample: GSM6037587
https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA825254&o=acc_s%3Aa

```

parallelize the workflow execute from working directory
```

function fetch_fastq() {
sample=$1
outputdir=$2
accession=$3

prefetch --type all $accession
fastq-dump --split-files $accession/$accession.sra --gzip --outdir $accession

# strip last 3 characters from accession and use as lane index
i=$(echo ${accession: -3})

# rename files to conform to bcl2fastq format
mv $accession/${accession}_1.fastq.gz $outputdir/${sample}_S1_L${i}_I1_001.fastq.gz
mv $accession/${accession}_2.fastq.gz $outputdir/${sample}_S1_L${i}_R1_001.fastq.gz
mv $accession/${accession}_3.fastq.gz $outputdir/${sample}_S1_L${i}_R2_001.fastq.gz
mv $accession/${accession}_4.fastq.gz $outputdir/${sample}_S1_L${i}_R3_001.fastq.gz

# cleanup
rm -rf $accession
}
export -f fetch_fastq

export PATH=$PATH:/mnt/g/software/sratoolkit.3.0.0-ubuntu64/bin 

sample=SAMN18736215
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 12 --bar fetch_fastq $sample $outputdir :::: /mnt/g/downloads/ckd-main/ckd-main/cnv/workflow/susztak/$sample.txt

sample=SAMN27505542
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 12 --bar fetch_fastq $sample $outputdir :::: /mnt/g/downloads/ckd-main/ckd-main/cnv/workflow/susztak/$sample.txt 

sample=SAMN18736216
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 12 --bar fetch_fastq $sample $outputdir :::: /mnt/g/downloads/ckd-main/ckd-main/cnv/workflow/susztak/$sample.txt

sample=SAMN27505541
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 12 --bar fetch_fastq $sample $outputdir :::: /mnt/g/downloads/ckd-main/ckd-main/cnv/workflow/susztak/$sample.txt

sample=SAMN27505543
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 12 --bar fetch_fastq $sample $outputdir :::: /mnt/g/downloads/ckd-main/ckd-main/cnv/workflow/susztak/$sample.txt

sample=SAMN27505544
outputdir=/mnt/g/atacFastq/$sample
mkdir -p $outputdir
parallel -j 12 --bar fetch_fastq $sample $outputdir :::: /mnt/g/downloads/ckd-main/ckd-main/cnv/workflow/susztak/$sample.txt

```
