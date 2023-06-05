# R script to make the covariate table
library("tidyverse")
library(RCurl)
setwd("~/projects/circRNA/results")

## BRAINcode v2
# subjectTable_url="https://docs.google.com/spreadsheets/d/e/2PACX-1vQFQ4aQj0sD9oxIqaZ-cEgo7kWcCmNYGBH9emLw8iNu0f6TTjKE5Lte7IBfoMMy57cLjA4pXE0YlPY2/pub?gid=28&single=true&output=tsv"
# sampleTable_url="https://docs.google.com/spreadsheets/d/e/2PACX-1vQFQ4aQj0sD9oxIqaZ-cEgo7kWcCmNYGBH9emLw8iNu0f6TTjKE5Lte7IBfoMMy57cLjA4pXE0YlPY2/pub?gid=1049148747&single=true&output=tsv"
# subjectTable=read.delim(textConnection(getURL(subjectTable_url)), stringsAsFactors = F) %>% 
#   select(CONDITION=DIAGNOSIS, SUBJECT_ID=SOURCE_SUBJECT_ID, AGE, SEX, PMI, PD.pathology.group, MUSS=Modified_Unified_Staging_System, Unified_LB_Stage, AD.pathology.group,	Plaque_density,	CERAD,	Braak_Braak_stage,	NIA_Reagan, MMSE=MMSE_last)
# sampleTable=read.delim(textConnection(getURL(sampleTable_url)), stringsAsFactors = F, comment.char = "#") %>%
#   filter(BRAINcode2.final.selection==1) %>% 
#   select(SAMPLE_ID=SOURCE_SAMPLE_ID, SUBJECT_ID=SOURCE_SUBJECT_ID, CELLTYPE=CELL, BATCH, RIN)

URL="https://docs.google.com/spreadsheets/d/1McWI4zAtC1qIS4wiYMXvToTKrDkFXaiAqhfDPDHyNIA"
library(googlesheets4)
sampleTable = read_sheet(URL, sheet ="RNAseq_statistics", .name_repair = "universal") %>% 
  filter(BRAINcode2.final.selection==1) %>% 
  select(SAMPLE_ID=SOURCE_SAMPLE_ID, SUBJECT_ID=SOURCE_SUBJECT_ID, CELLTYPE=CELL, BATCH, RIN)
subjectTable=read_sheet(URL, sheet ="Subjects_info", .name_repair = "universal") %>% 
  select(CONDITION=DIAGNOSIS, SUBJECT_ID=SOURCE_SUBJECT_ID, AGE, SEX, PMI, PD.pathology.group, MUSS=Modified_Unified_Staging_System, Unified_LB_Stage, AD.pathology.group,	Plaque_density,	CERAD,	Braak_Braak_stage,	NIA_Reagan, MMSE=MMSE_last)
head(subjectTable); head(sampleTable); 
dim(sampleTable) # n=197

## PD pathology group table
## Note: CONDITION includes PD, ILB, and HC; CONDITION2 includes PD (PD+ILB) and HC;
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