# Script repository for circRNA analysis

Update on 2022/12/14 by Xianjun Dong, PhD

## Introduction

As an extension of the [BRAINcode](http://www.humanbraincode.org) project, we studied the circRNAs in 197 human samples, including 190 human brain neuronal samples and 7 non-neuronal samples in this study. In detail, dopaminergic neurons (DA) from the midbrain substantia nigra pars compacta of 104 high-quality human brains (DA; HC: n = 59; ILB: n = 27; PD: n = 18), pyramidal neurons from layers V/VI of the middle temporal gyrus of 83 brains (TCPY; HC: n = 40; AD: n = 43) and from the primary motor cortex of three brains (MCPY; HC: n = 3) were laser-captured and pooled for deep total RNA-seq. Human fibroblasts from four (FB) and peripheral blood mononuclear white cells from three individuals (PBMC) were analyzed as non-CNS cell types using the same pipeline. The RNAseq raw data is accessible via GEO accession number GSE218203. 

This is a code repository for circRNA analysis, including:

### 1. circRNA detection and annotation

We ran circExplorer on all 197 BRAINcode RNA-seq samples and 12 RM samples (treated with RNase-R and mock, by Dr. Yunfei Bai), and then merge together into a combined matrix.

#### a. Ran circExplorer on each of the sample

circRNA calling is wrapped as part of our RNA-seq pipeline (see `RNAseq.pipeline.sh` in the [BRAINcode repository](https://github.com/sterding/BRAINcode)). For detail, please refer to the "*#STEP 4.1. mapping & calling circRNA*" in the [`_RNAseq.sh`](https://github.com/sterding/BRAINcode/blob/master/modules/_RNAseq.sh)

- For RM samples: `$HOME/neurogen/pipeline/RNAseq/RNAseq.pipeline.sh $PROJECT_PATH/rawfiles/RNAseq_RNaseR_Mock`
- For BC samples: `$HOME/neurogen/pipeline/RNAseq/RNAseq.pipeline.sh $PROJECT_PATH/rawfiles`

Note that we made two modifications: 

1. add CRD1as into the RefSeq file, e.g. extract CDR1as gtf from UCSC and converted to genePred then manually added to refFlat.txt - 10/23/2017

2. switch to use GENCODE-based refFlat as annotation (see README.txt in `../Annotation/Genes/`) - 8/12/2019


#### b. Merge circRNA called from all samples

We merged CircExplorer output files "circularRNA_known3.txt" to count circRNAs for several samples in one table.  see script `./annotation/main.sh` for details

In the end, it will generate 

- Merge_circexplorer_BC221.annotation.bed14: circRNA annotation (in BED format) for all BRAINcode samples (n=221 before QC).

- Merge_circexplorer_BC221.rawcount.long.txt: circRNA read count in long format for all BRAINcode samples (n=221 before QC).

- Merge_circexplorer_RM.annotation.bed14: circRNA annotation (in BED format) for all RM samples.

- Merge_circexplorer_BC_RM.annotation.bed14: circRNA annotation (in BED format) for all BRAINcode and RM samples.

#### c. circRNA filtering and normalization

We then converted the circRNA read count table from long format to wide format (e.g. rows are circRNAs and columns are samples) and normalized circRNAs expression into RPM (reads per million). 

We also filtered the samples based on BRAINcode sample QC (see Dong et al. *Nature Neuroscience* 2018 doi:10.1038/s41593-018-0223-0) and focused on the 197 samples for downstream analysis. We further filtered circRNAs by requiring it "being expressed"" (i.e. with >=2 unique back-splicing reads in at least 1 sample) and "being enriched" (i.e. at least 20 reads in RNase R treatment AND foldchange(RNase vs Mock) at least 2 in at least one sample). See script `./annotation/main.R` for details. 

In the end, it will generate 

- Merge_circexplorer_BC197.filtered.enriched.annotation.bed14: Annotation (in BED format) of filtered and enriched circRNAs in the BRAINcode samples (n=197 after QC).

- Merge_circexplorer_BC197.filtered.enriched.rawcount.rds: Read count matrix of filtered and enriched circRNAs in the BRAINcode samples (n=197 after QC).

- Merge_circexplorer_BC197.filtered.enriched.normRPM.rds: Normalized expression matrix of filtered and enriched circRNAs in the BRAINcode samples (n=197 after QC).


#### e. Make circRNA tracks for UCSC Genome Browser

See the "*#tracks for UCSC Genome Browser*" code trunk in the script `./annotation/main.sh` to generate trackDb.circRNA.txt file for all filtered and enriched circRNAs in the BRAINcode samples. 

### 2. cell-type specific expression of circRNAs

To call a circRNA "being cell-type specific", we required it with specificity score S >= 0.5 AND mean expression > mean + s.d. of overall expression  in the specific tissue or cell type. See the implementation in the R script `./annotation/getTissueSpecificRNA.R`. 

The specificity score was calculated as cummeRbund (https://github.com/Bioconductor-mirror/cummeRbund/blob/master/R/tools.R), see the library function code in `./annotation/tools.R`. 

### 3. disease/trait GWAS enrichment in circRNAs

First run `./annotation/circRNA.SNP.enrichment.sh` first to generate the input files required for testing the SNP co-localization enrichment analysis, then run `./annotation/circRNA.SNP.enrichment.R` to do the actual test and make plots. 

### 4. disease-associated genes and circRNAs, and pathways

First to make the covariate table by running the `./DE/makeCovariateTable.R`, then to run `./DE/_DE.sh` to call differential expression (DE) analysis for known genes and circRNAs among different covariate groups (e.g. PD vs. HC), and downstream pathway analysis for its DE genes. 

Note that `./DE/_DE.sh` will call several other core scripts, e.g. `./DE/_DE.R` for differential expression analysis based on DEseq2, `./DE/_pathway.R` for pathway/GO enrichment analysis, `_DE/_WGCNA.R` for weighted gene co-expression network analysis (WGCNA) analysis. 

Other scripts for meta analysis and comparative analysis between two DE results can be found in the `./DE/meta.R` and `./DE/comparison.R`. 

## Under development

Some prototype codes for circRNA QTL analysis (incl. eQTL, splicing QTL, circularization QTL) and visualization are at `./qtl/` and `./viz/` folders. 

## System requirements
- R (3.6.1)
- python (2.7.5)
- tophat (2.0.10)
- cufflinks (2.2.1)
- samtools (0.1.19)
- bedtools (2.26.0)
- circExplorer (2.3.0)


## Maintainence and contribution
This code is developed and maintained by Dr Xianjun Dong. Please email to xdong AT bwh.harvard.edu for questions and comments.

## License
Copyright 2022 Xianjun Dong

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

