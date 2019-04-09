###########################################
## bash script for running DE analysis using count data and covariance table
###########################################
#!/bin/sh


# begin R script
library("tidyverse")
library(RCurl)
setwd("~/projects/circRNA/results")
covariate_table_url="https://docs.google.com/spreadsheets/d/e/2PACX-1vTNGgyiJHx9YItTT3bkI9F4VqYbpbrahX0RPQTiUBhtirwMhUc5bGoTYwbnkLlEkwkTkp1o9TEFQi3o/pub?gid=195725118&single=true&output=tsv" ## covariate_table in BRAINCODE_Sequencing_Log (Freeze to N=142 sample for NN paper)
covariateTable=read.delim(textConnection(getURL(covariate_table_url)), stringsAsFactors = F)
head(covariateTable)
filter(covariateTable, BRAINCODE.final.selection == 1, cellType == 'SNDA', condition %in% c("HC","ILB")) %>% 
  select(SAMPLE_ID=sampleName, SUBJECT_ID=subjectID, CONDITION=condition,	Batch=batch, RIN, PMI, Sex=sex, Age = age) %>%
  mutate(CONDITION=ifelse(CONDITION=="HC","control","treated")) %>%
  write.table(file="HCILB_SNDA.covariates.txt", sep="\t", col.names = TRUE, row.names = F, quote = F)
# R end

## HC vs. AD in TCPY
cd ~/projects/circRNA/results
Rscript ~/neurogen/pipeline/RNAseq/modules/_DE.R -i ../data/Merge_circexplorer_BC106.filtered.enriched.rawcount.rds -c HCILB_SNDA.covariates.txt -o DE_HCvsILB_SNDA
