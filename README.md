# **Mosaic loss of Y chromosome is associated with aging and epithelial injury in chronic kidney disease**
__Parker C. Wilson<sup>1</sup>, Amit Verma<sup>1</sup>, Yasuhiro Yoshimura<sup>2</sup>, Yoshiharu Muto<sup>2</sup>, Haikuo Li<sup>2</sup>, Nicole P. Malvin<sup>2</sup>, Eryn
E. Dixon<sup>2</sup>, Benjamin D. Humphreys<sup>2,3</sup>__

<sup>1</sup> Division of Diagnostic Innovation, Department of Pathology and Laboratory Medicine, University of Pennsylvania, Philadelphia, PA, USA </br>
<sup>2</sup> Division of Nephrology, Department of Medicine, Washington University in St. Louis, St. Louis, MO, USA </br>
<sup>3</sup> Department of Developmental Biology, Washington University in St. Louis, St. Louis, MO, USA </br>

If you use any of the code or workflows in this repository please cite our manuscript in Genome Biology [link]()
```

```
The code associated with this publication has been deposited in [Zenodo]()

Single cell multiome and snATAC-seq data generated for this manuscript (multiomes: 6, snATAC-seq: 5 CKD, 2 DKD) and cellranger-arc v2.0 / cellranger-atac v2.1 gene and peak count matrices for can be found in GEO. </br>
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232222

Visium spatial sequencing data generated for this manuscript can be found in a second repository: </br>
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232431

Sequencing data generated for previously-published kidney single cell multiomes (n=3)  can be found at the link below:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220289 <br/>

Sequencing data generated for previously-published kidney snATAC-seq (n=17)  can be found at the links below:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151302 <br/>
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195460 <br/>
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172008 <br/>
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200547 <br/>

[single cell multiome barcodes](https://github.com/p4rkerw/Wilson_GBio_2024/blob/main/multi_aggr_prep_kidney/multi_barcodes.csv)
[snATAC barcodes](https://github.com/p4rkerw/Wilson_GBio_2024/blob/main/atac_aggr_prep_kidney/atac_barcodes.csv) </br>
[snRNA barcodes](https://github.com/p4rkerw/Wilson_GBio_2024/blob/main/rna_aggr_prep/kpmp/rna_barcodes.csv) used for the analysis can be found in this github repository


Welcome to our GitHub repository!  
Here you will find analysis scripts for our manuscript where we use single cell sequencing to detect loss of Y chromosome (LOY) and other mosaic chromosomal alterations (mCA) in chronic kidney disease. Please contact the corresponding author, Dr. Parker Wilson, with questions or comments.  
<br/>
Thanks,  
Parker
<br/><br/>

Visit the Wilson lab website:<br/>
www.parkerwilsonlab.com

Visit the Humphreys lab website:<br/>
www.humphreyslab.com  
<br/>
Check out our interactive datasets with Kidney Interactive mulTiomics (KIT):  
http://humphreyslab.com/SingleCell/
<br/><br/>
Find us on Twitter: 
<br/>
  <a href="https://twitter.com/parkercwilson?ref_src=twsrc%5Etfw" class="twitter-follow-button" data-show-count="false"> @parkercwilson</a>
  <a href="https://twitter.com/HumphreysLab?ref_src=twsrc%5Etfw" class="twitter-follow-button" data-show-count="false"> @HumphreysLab</a>
<br/><br/>
Find us on Docker Hub:  
[p4rkerw@dockerhub](https://hub.docker.com/search?q=p4rkerw&type=image)
<br/>

**Single cell multiome preprocessing and analysis workflow**
1. Align and count each multiome library (multi_aggr_prep_kidney/cellranger/cellranger_arc_count.sh)
Libraries were generated from a nuclear dissociation and aligned to refdata-cellranger-arc-GRCh38-2020-A-2.0.0, which can be downloaded from the 10X genomics website: https://support.10xgenomics.com/.

2. Aggregate the libraries using the multi_aggr_prep_kidney/cellranger/multi_aggr.csv file (multi_aggr_prep_kidney/cellranger/cellranger_arc_aggr.sh)

3. Identify doublet barcodes with AMULET (multi_aggr_prep_kidney/step0_amulet.sh)

4. Bin multiome ATAC fragments into 1Mb bins with epiAneufinder. Exclude barcodes with < 10,000 fragments. Generate cytoband_counts10k.rds file for each library. (multi_aggr_prep_kidney/step1_multi_bin_fragments.R)

5. Run QC and preprocessing routing with Seurat and annotate barcodes (multi_aggr_prep_kidney/step2_multi_prep.R)

6. Count and normalize chrY RNA transcripts. Count and normalize ATAC fragments for all chromosomes. Classify LOY using a gaussian finite mixture model.  (multi_aggr_prep_kidney/step3_multi_loy.R)

7. Find differentially expressed genes for LOY vs XY cells for all cell types with and without age adjustment (multi_aggr_prep_kidney/analysis/find_deg.R)

8. Find differentially accessible chromatin regions for LOY vs XY cells for all cell types with and without age adjustment (multi_aggr_prep_kidney/analysis/find_dar.R)

9. Run epiAneufinder using 1Mb bins with default settings on autosomal chromosomes (multi_aggr_prep_kidney/step4_multi_epianeufinder.R)

10. Estimate autosomal CNV burden using epiAneufinder output (multi_aggr_prep_kidney/step5_cnv_burden.R)

11. Run gene set enrichment analysis on differentially expressed genes that differentiate LOY vs XY cells in the proximal convoluted tubule and other cell types. (multi_aggr_prep_kidney/step6_gsea.R)

12. Run chromVAR to estimate TF motif activities (multi_aggr_prep_kidney/step7_chromvar.R)

13. Find differentially activity of TF motifs with chromVAR for LOY vs XY cells for all cell types (multi_aggr_prep_kidney/analysis/find_chromvar.R)

**snATAC kidney preprocessing and analysis workflow**  
1. Align and count each ATAC library (atac_aggr_prep_kidney/cellranger/cellranger_atac_count.sh)  
Libraries were generated from a nuclear dissociation and aligned to refdata-cellranger-arc-A-2.0.0 which can be downloaded from the 10X genomics website: https://support.10xgenomics.com/. 

2. Aggregate the snATAC libraries using the atac_aggr_prep_kidney/cellranger/atac_aggr_22.csv file (atac_aggr_prep_kidney/cellranger/cellranger_atac_aggr.sh)

3. Identify doublet barcodes with AMULET (atac_aggr_prep_kidney/step0_amulet.sh)

4. Run standard Signac QC on the aggregated snATAC data, perform batch effect correction with Harmony, transfer cell type labels from a previously-published snRNA-seq atlas and visualize cell-specific markers. (atac_aggr_prep_kidney/step1_prep.R)

5. Annotate barcodes (atac_aggr_prep_kidney/step2_anno.R)

6. Run chromVAR to estimate TF motif activities (atac_aggr_prep_kidney/step3_chromvar.R)

8. Bin ATAC fragments into 1Mb bins with epiAneufinder. Exclude barcodes with < 10,000 fragments. Generate cytoband_counts10k.rds file for each library. (atac_aggr_prep_kidney/step4_atac_bin_fragments.R)

9. Run epiAneufinder using 1Mb bins with default settings on autosomal chromosomes (atac_aggr_prep_kidney/step5_multi_epianeufinder.R)

10. Count and normalize ATAC fragments for all chromosomes. Classify LOY using a density threshold model.  (multi_aggr_prep_kidney/step6_atac_loy.R)

12. Estimate autosomal CNV burden using epiAneufinder output (atac_aggr_prep_kidney/step7_cnv_burden.R)

13. Find differentially accessible regions for LOY vs XY cells for all cell types with and without age adjustment (atac_aggr_prep_kidney/analysis/find_dar.R)

14. Find differentially activity of TF motifs with chromVAR for LOY vs XY cells for all cell types (atac_aggr_prep_kidney/analysis/find_chromvar.R)

15. Run gene set enrichment analysis on differentially accessible regions that differentiate LOY vs XY cells in the proximal convoluted tubule and other cell types. (atac_aggr_prep_kidney/step8_gsea.R)

**KPMP sc/snRNA-seq preprocessing and analysis workflow**

1. Download the KPMP dataset in h5seurat format from the KPMP website "c798e11b-bbde-45dd-bd91-487f27c93f8f_WashU-UCSD_HuBMAP_KPMP-Biopsy_10X-R_12032021.h5Seurat"

2. Run QC and preprocessing routing with Seurat, predict cell types using label transfer from a previously-published snRNA-seq atlas, perform batch effect correction with Harmony and annotate barcodes. Count and normalize chrY transcripts and assign LOY using density threshold model (rna_aggr_prep_kidney/step1_kpmp.R)

3. Find differentially expressed genes that differentiate LOY vs XY cells (rna_aggr_prep_kidney/analysis/find_deg.R)

4. Gene set enrichment analysis for LOY vs XY differentially expressed genes (and other comparisons) (rna_aggr_prep_kidney/step2_gsea.R)

5. Project the single cell multiome atlas onto the KPMP atlas to harmonize cell type annotations (rna_aggr_prep_kidney/step3_harmonize.R)


**Visium spatial preprocessing and analysis workflow**


**snATAC leukocyte preprocessing and analysis workflow**  




**Figures**

Scripts for generating figures in the manuscript.

