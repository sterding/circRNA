###########################################
## bash script for running DE analysis using count data and covariance table
###########################################
#!/bin/sh


# begin R script
library("tidyverse")
library(RCurl)
setwd("~/projects/circRNA/results")

## BRAINcode v2
subjectTable_url="https://docs.google.com/spreadsheets/d/e/2PACX-1vQFQ4aQj0sD9oxIqaZ-cEgo7kWcCmNYGBH9emLw8iNu0f6TTjKE5Lte7IBfoMMy57cLjA4pXE0YlPY2/pub?gid=28&single=true&output=tsv"
sampleTable_url="https://docs.google.com/spreadsheets/d/e/2PACX-1vQFQ4aQj0sD9oxIqaZ-cEgo7kWcCmNYGBH9emLw8iNu0f6TTjKE5Lte7IBfoMMy57cLjA4pXE0YlPY2/pub?gid=1049148747&single=true&output=tsv"
subjectTable=read.delim(textConnection(getURL(subjectTable_url)), stringsAsFactors = F) %>% 
  select(CONDITION=DIAGNOSIS, SUBJECT_ID=SOURCE_SUBJECT_ID, AGE, SEX, PMI, PD.pathology.group, MUSS=Modified_Unified_Staging_System, Unified_LB_Stage, AD.pathology.group,	Plaque_density,	CERAD,	Braak_Braak_stage,	NIA_Reagan, MMSE=MMSE_last)
sampleTable=read.delim(textConnection(getURL(sampleTable_url)), stringsAsFactors = F, comment.char = "#") %>%
  filter(BRAINcode2.final.selection==1) %>% 
  select(SAMPLE_ID=SOURCE_SAMPLE_ID, SUBJECT_ID=SOURCE_SUBJECT_ID, CELLTYPE=CELL, BATCH, RIN)
head(subjectTable); head(sampleTable); 
dim(sampleTable) # n=284

## PD pathology group table
left_join(x=sampleTable, y=subjectTable, by = "SUBJECT_ID") %>% filter(CELLTYPE=="SNDA") %>% 
  select(SUBJECT_ID, CONDITION, AGE, SEX, PMI, MUSS, Unified_LB_Stage, PD.pathology.group, SAMPLE_ID, CELLTYPE, BATCH, RIN) %>% 
  mutate(PDpathologygroup=ifelse(PD.pathology.group=="no","hc",ifelse(PD.pathology.group=="late" | PD.pathology.group=="early","pd",PD.pathology.group))) %>%
  mutate(CONDITION2=ifelse(CONDITION=="ILB" | CONDITION=="PD","PD",CONDITION)) %>%
  arrange(CONDITION2, CONDITION) %>% 
  write.table(file="Table.PD.SNDA.pathology.covariates.xls", sep="\t", col.names = TRUE, row.names = F, quote = F)

## AD pathology group table
left_join(x=sampleTable, y=subjectTable, by = "SUBJECT_ID") %>% filter(CELLTYPE=="TCPY") %>% 
  select(SUBJECT_ID, CONDITION, AGE, SEX, PMI, Plaque_density, CERAD, Braak_Braak_stage, NIA_Reagan,  MMSE, AD.pathology.group, SAMPLE_ID, CELLTYPE, BATCH, RIN) %>% 
  arrange(AD.pathology.group, CONDITION) %>% 
  write.table(file="Table.AD.TCPY.pathology.covariates.xls", sep="\t", col.names = TRUE, row.names = F, quote = F)

## PD CSF pathology group table
left_join(x=sampleTable, y=subjectTable, by = "SUBJECT_ID") %>% filter(CELLTYPE=="CSF") %>% 
  select(SUBJECT_ID, CONDITION, AGE, SEX, SAMPLE_ID, CELLTYPE, BATCH) %>% 
  arrange(CONDITION) %>% 
  write.table(file="Table.PD.CSF.pathology.covariates.xls", sep="\t", col.names = TRUE, row.names = F, quote = F)

## Bennett: AD/PD total RNAseq pathology group table
sampleTable_url="https://docs.google.com/spreadsheets/d/e/2PACX-1vQFQ4aQj0sD9oxIqaZ-cEgo7kWcCmNYGBH9emLw8iNu0f6TTjKE5Lte7IBfoMMy57cLjA4pXE0YlPY2/pub?gid=1322938821&single=true&output=tsv"
sampleTable=read.delim(textConnection(getURL(sampleTable_url)), stringsAsFactors = F, comment.char = "#") %>%
  filter(Seq_type=='RNAseq') %>% 
  select(SAMPLE_ID=SampleID, SUBJECT_ID=SubjectID, CELLTYPE=Tissue, CONDITION=Neuropath_Dx, AGE=Age, RIN)
head(sampleTable); 
dim(sampleTable) # n=42
filter(sampleTable, CELLTYPE=="VMB") %>% 
  select(SUBJECT_ID, CONDITION, AGE, SAMPLE_ID, RIN) %>% 
  arrange(CONDITION) %>% 
  write.table(file="Table.Bennett.PD.VMP.pathology.covariates.xls", sep="\t", col.names = TRUE, row.names = F, quote = F)
filter(sampleTable, CELLTYPE=="FCX") %>% 
  select(SUBJECT_ID, CONDITION, AGE, SAMPLE_ID, RIN) %>% 
  arrange(CONDITION) %>% 
  write.table(file="Table.Bennett.AD.FCX.pathology.covariates.xls", sep="\t", col.names = TRUE, row.names = F, quote = F)

cd ~/projects/circRNA/results/
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C CONDITION:PD:HC -s
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C CONDITION:ILB:HC -s
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C CONDITION:PD:ILB -s
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O Mmi -C CONDITION2:PD:HC -s
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C PDpathologygroup:pd:hc -s
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C PD.pathology.group:early:no -s
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C PD.pathology.group:late:early -s
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C PD.pathology.group:late:no -s
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.AD.TCPY.pathology.covariates.xls -o DE_TCPY -O Mmi -C CONDITION:AD:HC -s
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.AD.TCPY.pathology.covariates.xls -o DE_TCPY -O i -C Braak_Braak_stage -s
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O i -C MUSS -s

## PD: collapse to gene
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE2gene_SNDA -O Mmi -C CONDITION:PD:HC -k -s
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE2gene_SNDA -O Mmi -C CONDITION:ILB:HC -k -s
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE2gene_SNDA -O Mmi -C CONDITION2:PD:HC -k -s
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE2gene_SNDA -O i -C MUSS -k -s
## AD: collapse to gene
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.AD.TCPY.pathology.covariates.xls -o DE2gene_TCPY -O Mmi -C CONDITION:AD:HC -k -s
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.AD.TCPY.pathology.covariates.xls -o DE2gene_TCPY -O i -C Braak_Braak_stage -k -s

# for i in DEresult.DE2gene_SNDA.*xls; do awk '$6<=0.05' $i | cut -f1-6; done | cut -f1 | sort -u > DEgenes.list
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE2gene_SNDA -O Mmi -C CONDITION2:PD:HC -k -l DEgenes.list -s
# awk '$2<=0.1 && $4<=0.05 && $6<=0.05{print $1}' summaryTable.VMB.meta.p0.05.xls > summaryTable.VMB.meta.p0.05.list
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE2gene_SNDA -O Mmi -C CONDITION:ILB:HC -k -l summaryTable.VMB.meta.p0.05.list -s

# CSF
Rscript ~/projects/circRNA/src/DE/_DE2gene.R -i ../data/Merge_circexplorer_CSF87.rawcount.rds -c Table.PD.CSF.pathology.covariates.xls -o DE2gene_CSF -O Mmi -C CONDITION2:PD:HC -a ~/projects/circRNA/data/Merge_circexplorer_CSF87.annotation.bed14.rds -k -s

## Bennett's dataset
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_Bennett_VMB.filtered.rawcount.rds -c Table.Bennett.PD.VMP.pathology.covariates.xls -o DE2gene_VMB -O Mmi -C CONDITION:PD:HC -a ~/projects/circRNA/data/Merge_circexplorer_Bennett_VMB.filtered.annotation.bed14.rds -k
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_Bennett_FCX.filtered.rawcount.rds -c Table.Bennett.AD.FCX.pathology.covariates.xls -o DE2gene_FCX -O Mmi -C CONDITION:AD:HC -a ~/projects/circRNA/data/Merge_circexplorer_Bennett_FCX.filtered.annotation.bed14.rds -k

## pathway analysis
cd ~/projects/circRNA/results/DE2gene_SNDA; Rscript ~/projects/circRNA/src/DE/_pathway.R -i DEresult.DE2gene_SNDA.CONDITION2_PD_vs_HC.xls -F loose
cd ~/projects/circRNA/results/DE2gene_TCPY; Rscript ~/projects/circRNA/src/DE/_pathway.R -i DEresult.DE2gene_TCPY.CONDITION_AD_vs_HC.xls -F loose

## WGCNA analysis
cd ~/projects/circRNA/results/DE_SNDA; Rscript ~/projects/circRNA/src/DE/_WGCNA.R -i ~/projects/circRNA/data/Merge_circexplorer_BC197.filtered.enriched.normRPM.rds -c 

## AD vs. PD DE comparison
# see comparison.R


##======================
## gene
##======================
cd ~/neurogen/rnaseq_PD/results/; Rscript ~/projects/circRNA/src/DE/_DE.R -i merged/genes.htseqcount.cufflinks.allSamples.BCv2.uniq.xls -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA.gene -O Mmi -C CONDITION2:PD:HC -t gene
cd ~/neurogen/rnaseq_PD/results/; Rscript ~/projects/circRNA/src/DE/_DE.R -i merged/genes.htseqcount.cufflinks.allSamples.BCv2.uniq.xls -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA.gene -O Mmi -C MUSS -t gene
cd ~/neurogen/rnaseq_PD/results/; Rscript ~/projects/circRNA/src/DE/_DE.R -i merged/genes.htseqcount.cufflinks.allSamples.BCv2.uniq.xls -c Table.AD.TCPY.pathology.covariates.xls -o DE_TCPY.gene -O Mmi -C CONDITION:AD:HC -t gene
cd ~/neurogen/rnaseq_PD/results/; Rscript ~/projects/circRNA/src/DE/_DE.R -i merged/genes.htseqcount.cufflinks.allSamples.BCv2.uniq.xls -c Table.AD.TCPY.pathology.covariates.xls -o DE_TCPY.gene -O Mmi -C Braak_Braak_stage -t gene
cd ~/neurogen/rnaseq_CSF/results; Rscript ~/projects/circRNA/src/DE/_DE.R -i merged/genes.htseqcount.cufflinks.CSF.uniq.xls -c Table.PD.CSF.pathology.covariates.xls -o DE_CSF.gene -O Mmi -C CONDITION:PD:HC -t gene


## pathway analysis
cd ~/neurogen/rnaseq_PD/results/DE_SNDA.gene; Rscript ~/projects/circRNA/src/DE/_pathway.R -i DEresult.DE_SNDA.gene.CONDITION2_PD_vs_HC.xls -F medium
cd ~/neurogen/rnaseq_PD/results/DE_TCPY.gene; Rscript ~/projects/circRNA/src/DE/_pathway.R -i DEresult.DE_TCPY.gene.CONDITION_AD_vs_HC.xls -F medium

## co-expression heatmap
cd ~/neurogen/rnaseq_PD/results/DE_SNDA.gene; Rscript ~/pipeline/modules/_expression2heatmap.R -i genes.htseqcount.cufflinks.allSamples.BCv2.uniq.xls.filtered.vsd_adjusted_log2.xls -l DEresult.DE_SNDA.gene.MUSS.padj0.05.xls