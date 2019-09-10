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
dim(sampleTable) # n=197

## PD pathology group table
left_join(x=sampleTable, y=subjectTable, by = "SUBJECT_ID") %>% filter(CELLTYPE=="SNDA") %>% select(SUBJECT_ID, CONDITION, AGE, SEX, PMI, MUSS, Unified_LB_Stage, PD.pathology.group, SAMPLE_ID, CELLTYPE, BATCH, RIN) %>% arrange(PD.pathology.group, CONDITION) %>% write.table(file="Table.PD.SNDA.pathology.covariates.xls", sep="\t", col.names = TRUE, row.names = F, quote = F)

## AD pathology group table
left_join(x=sampleTable, y=subjectTable, by = "SUBJECT_ID") %>% filter(CELLTYPE=="TCPY") %>% select(SUBJECT_ID, CONDITION, AGE, SEX, PMI, Plaque_density, CERAD, Braak_Braak_stage, NIA_Reagan,  MMSE, AD.pathology.group, SAMPLE_ID, CELLTYPE, BATCH, RIN) %>% arrange(AD.pathology.group, CONDITION) %>% write.table(file="Table.AD.TCPY.pathology.covariates.xls", sep="\t", col.names = TRUE, row.names = F, quote = F)

## PD CSF pathology group table
left_join(x=sampleTable, y=subjectTable, by = "SUBJECT_ID") %>% filter(CELLTYPE=="CSF") %>% select(SUBJECT_ID, CONDITION, AGE, SEX, SAMPLE_ID, CELLTYPE, BATCH) %>% arrange(CONDITION) %>% write.table(file="Table.PD.CSF.pathology.covariates.xls", sep="\t", col.names = TRUE, row.names = F, quote = F)

Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C CONDITION:PD:HC
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C CONDITION:ILB:HC
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C CONDITION:PD:ILB
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C PD.pathology.group:early:no
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C PD.pathology.group:late:early
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.AD.TCPY.pathology.covariates.xls -o DE_TCPY -O mi -C CONDITION:AD:HC
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_CSF87.filtered.enriched.rawcount.rds -c Table.PD.CSF.pathology.covariates.xls -o DE_CSF -C CONDITION:PD:HC
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.AD.TCPY.pathology.covariates.xls -o DE_TCPY -O mi -C Braak_Braak_stage
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C MUSS




#### OLD


covariate_table_url="https://docs.google.com/spreadsheets/d/e/2PACX-1vTNGgyiJHx9YItTT3bkI9F4VqYbpbrahX0RPQTiUBhtirwMhUc5bGoTYwbnkLlEkwkTkp1o9TEFQi3o/pub?gid=195725118&single=true&output=tsv" ## covariate_table in BRAINCODE_Sequencing_Log (Freeze to N=142 sample for NN paper)
covariateTable=read.delim(textConnection(getURL(covariate_table_url)), stringsAsFactors = F)
head(covariateTable)
filter(covariateTable, BRAINCODE.final.selection == 1, cellType == 'SNDA', condition %in% c("HC","ILB")) %>% 
  select(SAMPLE_ID=sampleName, SUBJECT_ID=subjectID, CONDITION=condition,	Batch=batch, RIN, PMI, Sex=sex, Age = age) %>%
  mutate(CONDITION=ifelse(CONDITION=="HC","control","treated")) %>%
  write.table(file="HCILB_SNDA.covariates.txt", sep="\t", col.names = TRUE, row.names = F, quote = F)

## HC vs. AD in TCPY
cd ~/projects/circRNA/results
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC106.filtered.enriched.rawcount.rds -c HCILB_SNDA.covariates.txt -o DE_HCvsILB_SNDA



