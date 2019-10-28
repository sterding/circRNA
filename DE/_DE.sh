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

#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C CONDITION:PD:HC
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C CONDITION:ILB:HC
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C CONDITION:PD:ILB
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O Mmi -C CONDITION2:PD:HC
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C PDpathologygroup:pd:hc
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C PD.pathology.group:early:no
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C PD.pathology.group:late:early
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C PD.pathology.group:late:no
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.AD.TCPY.pathology.covariates.xls -o DE_TCPY -O Mmi -C CONDITION:AD:HC
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.AD.TCPY.pathology.covariates.xls -o DE_TCPY -O i -C Braak_Braak_stage
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O i -C MUSS

Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE2gene_SNDA -O Mmi -C CONDITION:PD:HC -k
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE2gene_SNDA -O Mmi -C CONDITION:ILB:HC -k

Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE2gene_SNDA -O Mmi -C CONDITION2:PD:HC -k
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.AD.TCPY.pathology.covariates.xls -o DE2gene_TCPY -O Mmi -C CONDITION:AD:HC -k 
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.AD.TCPY.pathology.covariates.xls -o DE2gene_TCPY -O i -C Braak_Braak_stage -k
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE2gene_SNDA -O i -C MUSS -k

# for i in DEresult.DE2gene_SNDA.*xls; do awk '$6<=0.05' $i | cut -f1-6; done | cut -f1 | sort -u > DEgenes.list
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE2gene_SNDA -O Mmi -C CONDITION2:PD:HC -k -l DEgenes.list

# CSF
Rscript ~/projects/circRNA/src/DE/_DE2gene.R -i ../data/Merge_circexplorer_CSF87.rawcount.rds -c Table.PD.CSF.pathology.covariates.xls -o DE2gene_CSF -O Mmi -C CONDITION2:PD:HC -a ~/projects/circRNA/data/Merge_circexplorer_CSF87.annotation.bed14.rds


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
cd ~/neurogen/rnaseq_PD/results/; Rscript ~/projects/circRNA/src/DE/_DE.R -i merged/genes.htseqcount.cufflinks.allSamples.BCv2.uniq.xls -c Table.AD.TCPY.pathology.covariates.xls -o DE_TCPY.gene -O Mmi -C CONDITION:AD:HC -t gene
cd ~/neurogen/rnaseq_PD/results/; Rscript ~/projects/circRNA/src/DE/_DE.R -i merged/genes.htseqcount.cufflinks.allSamples.BCv2.uniq.xls -c Table.AD.TCPY.pathology.covariates.xls -o DE_TCPY.gene -O Mmi -C CONDITION:AD:HC -t gene
cd ~/neurogen/rnaseq_CSF/results; Rscript ~/projects/circRNA/src/DE/_DE.R -i merged/genes.htseqcount.cufflinks.CSF.uniq.xls -c Table.PD.CSF.pathology.covariates.xls -o DE_CSF.gene -O Mmi -C CONDITION:PD:HC -t gene


## pathway analysis
cd ~/neurogen/rnaseq_PD/results/DE_SNDA.gene; Rscript ~/projects/circRNA/src/DE/_pathway.R -i DEresult.DE_SNDA.gene.CONDITION2_PD_vs_HC.xls -F medium
cd ~/neurogen/rnaseq_PD/results/DE_TCPY.gene; Rscript ~/projects/circRNA/src/DE/_pathway.R -i DEresult.DE_TCPY.gene.CONDITION_AD_vs_HC.xls -F medium