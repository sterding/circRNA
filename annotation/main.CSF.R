## Main script to process CSF samples
library(tidyverse)
library(RCurl)

## See ~/projects/circRNA/README.Rmd for how to generate the input files
setwd("~/projects/circRNA/data/")

# load
Merge_circexp_raw <- read.table("Merge_circexplorer_CSF87.rawcount.long.txt", sep="\t",check.names =F, col.names = c("ID","readsCount","sampleID")) %>% 
  group_by(ID, sampleID) %>% summarise(readsCount=sum(readsCount)) %>%
  spread(key = sampleID, value = readsCount, fill = 0) %>%
  column_to_rownames("ID")
head(Merge_circexp_raw); dim(Merge_circexp_raw)

# annotation
annotation<- read.table("Merge_circexplorer_CSF87.annotation.bed14",header=T, check.names = F) # include both ciRNA and circRNA
head(annotation); dim(annotation)
# 42423

# add geneSymbol, geneType
GENCODEv19=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v19.annotation.genes.bed", col.names = c("chr","start","end","geneName","geneID","geneType","strand"), stringsAsFactors = F, header = F)
annotation$geneName=GENCODEv19$geneName[match(annotation$geneID, GENCODEv19$geneID)]
annotation$geneType=GENCODEv19$geneType[match(annotation$geneID, GENCODEv19$geneID)]

# add CDR1as
annotation$geneName[annotation$geneID=='CDR1as'] = 'CDR1as'
annotation$geneType[annotation$geneID=='CDR1as'] = 'antisense'

###########################################
############## normalize to RPM ###########
###########################################
BRAINCODE2_Sequencing_Log.RNAseq_statistics="https://docs.google.com/spreadsheets/d/e/2PACX-1vQFQ4aQj0sD9oxIqaZ-cEgo7kWcCmNYGBH9emLw8iNu0f6TTjKE5Lte7IBfoMMy57cLjA4pXE0YlPY2/pub?gid=1049148747&single=true&output=tsv"
RNAseq_statistics<- read.delim(textConnection(getURL(BRAINCODE2_Sequencing_Log.RNAseq_statistics)), header=T, check.names = T, comment.char = "#", stringsAsFactors = F) 
readsNum_million <- as.numeric(gsub(",","",RNAseq_statistics$total_mapped_reads))/10^6; names(readsNum_million) = RNAseq_statistics$SOURCE_SAMPLE_ID
readsNum_million = readsNum_million[colnames(Merge_circexp_raw)]

Merge_circexp_norm<-sweep(Merge_circexp_raw,2,readsNum_million,"/")

# save
saveRDS(Merge_circexp_raw, file="Merge_circexplorer_CSF87.rawcount.rds")
saveRDS(Merge_circexp_norm, file="Merge_circexplorer_CSF87.normRPM.rds")
saveRDS(annotation, file="Merge_circexplorer_CSF87.annotation.bed14.rds")

#write.table(Merge_circexp_raw, file="Merge_circexplorer_CSF87.rawcount.xls.gz", quote=F, sep="\t", col.names = NA, row.names = TRUE)
#write.table(Merge_circexp_norm, file="Merge_circexplorer_CSF87.normRPM.xls", quote=F, sep="\t", col.names = NA, row.names = TRUE)

## ===============================================
## back-splicing reads vs. number of samples
## similar to Fig. 3E in https://www-sciencedirect-com.ezp-prod1.hul.harvard.edu/science/article/pii/S0092867418316350?via%3Dihub#fig3
## ===============================================
df=Merge_circexp_raw
pdf("../results/Merge_circexplorer_CSF87.nReads_vs_nSample.pdf", width = 6, height = 6)
size=rowSums(df>0)*rowSums(df)
size=(size-min(size))/(max(size)-min(size))
size[size==0]=0.005
plot(jitter(rowSums(df>0)), jitter(rowSums(df)), xaxt='n', cex=sqrt(size)*10+0.01, pch=20, col=rgb(sqrt(size),0,0,sqrt(size)), log='y', xlim=c(0,10), ylim=c(1,400), xlab="Number of samples expressed (out of 87)", ylab="Total number of back-spliced reads", main="Merge_circexplorer_CSF87")
axis(1, at=1:9,labels=1:9)
#axis(2, at=c(1,2,5,10,20,50,100,),labels=x, col.axis="red", las=2)
dev.off()
## only few circRNAs expressed in 90% samples in our case. So, it's applicable to divide into 'high' vs. 'low' like the paper above. 

###########################################
# How many are they also detected in PD SNDA?
###########################################
Merge_circexp_raw = readRDS("Merge_circexplorer_CSF87.rawcount.rds")
dim(Merge_circexp_raw); head(Merge_circexp_raw)
Merge_circexp_raw_PD=select(Merge_circexp_raw, starts_with("PD")); dim(Merge_circexp_raw_PD)
Merge_circexp_raw_filtered_PD=Merge_circexp_raw_PD[rowSums(Merge_circexp_raw_PD)>=2,]; 
dim(Merge_circexp_raw_filtered_PD)
SNDA=readRDS("Merge_circexplorer_BC190.filtered.rawcount.rds")
########## TODO
head(SNDA); dim(SNDA);
Merge_circexp_raw_SNDA=select(SNDA, contains("SNDA")); head(Merge_circexp_raw_SNDA)
Merge_circexp_raw_filtered_SNDA=Merge_circexp_raw_SNDA[rowSums(Merge_circexp_raw_SNDA)>0,]
dim(Merge_circexp_raw_filtered_SNDA)
sum(rownames(Merge_circexp_raw_filtered_PD) %in% rownames(Merge_circexp_raw_filtered_SNDA))
# 148

# collapsed to gene level
table((filter(readRDS("Merge_circexplorer_CSF87.annotation.bed14.rds"),ID %in% rownames(Merge_circexp_raw_filtered_PD)) %>%
  pull(geneName) %>% unique()) %in% 
  (filter(readRDS("Merge_circexplorer_BC197.annotation.bed14.rds"),ID %in% rownames(Merge_circexp_raw_filtered_SNDA)) %>%
                                       pull(geneName) %>% unique()))

# More ciRNAs detected in CSF than in neuron
table(filter(readRDS("Merge_circexplorer_CSF87.annotation.bed14.rds"),ID %in% rownames(Merge_circexp_raw_filtered_PD)) %>% pull(circType))
# circRNA   ciRNA 
# 24    3827 
table(filter(readRDS("Merge_circexplorer_BC197.annotation.bed14.rds"),ID %in% rownames(Merge_circexp_raw_filtered_SNDA)) %>% pull(circType))

###########################################
############ filter cirRNAs     ###########
###########################################

Merge_circexp_raw_filtered <- Merge_circexp_raw[rowSums(Merge_circexp_raw)>=2,]
Merge_circexp_norm_filtered <- Merge_circexp_norm[rownames(Merge_circexp_raw_filtered), ]
annotation_filtered <- annotation[annotation$ID %in% rownames(Merge_circexp_norm_filtered),] 

dim(Merge_circexp_raw_filtered); dim(Merge_circexp_norm_filtered); dim(annotation_filtered)
# 7902

saveRDS(Merge_circexp_raw_filtered, file="Merge_circexplorer_CSF87.filtered.rawcount.rds")
saveRDS(Merge_circexp_norm_filtered, file="Merge_circexplorer_CSF87.filtered.normRPM.rds")
saveRDS(annotation_filtered, file="Merge_circexplorer_CSF87.filtered.annotation.bed14.rds")

# Q: How many are they also detected in TCPY?
Merge_circexp_raw_filtered=readRDS("Merge_circexplorer_CSF87.filtered.rawcount.rds")
Merge_circexp_raw_TCPY=readRDS("Merge_circexplorer_BC190.filtered.rawcount.rds") %>% select(contains("TCPY"))
head(Merge_circexp_raw_TCPY); dim(Merge_circexp_raw_TCPY);
Merge_circexp_raw_filtered_TCPY=Merge_circexp_raw_TCPY[rowSums(Merge_circexp_raw_TCPY)>0,]
dim(Merge_circexp_raw_filtered_TCPY)
sum(rownames(Merge_circexp_raw_filtered) %in% rownames(Merge_circexp_raw_filtered_TCPY))
# 140

# Q: How many are they also detected in brain neurons?
Merge_circexp_raw_filtered=readRDS("Merge_circexplorer_CSF87.filtered.rawcount.rds")
Merge_circexp_raw_brain=readRDS("Merge_circexplorer_BC190.filtered.rawcount.rds") %>% select(contains(c("TCPY","MCPY","SNDA")))
head(Merge_circexp_raw_brain); dim(Merge_circexp_raw_brain);
Merge_circexp_raw_filtered_brain=Merge_circexp_raw_brain[rowSums(Merge_circexp_raw_brain)>0,]
dim(Merge_circexp_raw_filtered_brain)
sum(rownames(Merge_circexp_raw_filtered) %in% rownames(Merge_circexp_raw_filtered_brain))
# 382

########################################################
############ enriched cirRNAs in RNase treatment #######
########################################################

Merge_circexp_raw_long_enriched.RM  = readRDS(file="Merge_circexp_raw_long_filterd.RM.rds")

Merge_circexp_norm_filtered_and_enriched=Merge_circexp_norm_filtered[intersect(rownames(Merge_circexp_raw_filtered), Merge_circexp_raw_long_enriched.RM$ID), ]
Merge_circexp_raw_filtered_and_enriched <- Merge_circexp_raw_filtered[rownames(Merge_circexp_norm_filtered_and_enriched), ]
annotation_filtered_enriched <- annotation[annotation$ID %in% rownames(Merge_circexp_raw_filtered_and_enriched),] 

dim(Merge_circexp_norm_filtered_and_enriched); dim(Merge_circexp_raw_filtered_and_enriched); dim(annotation_filtered_enriched)
# [1] 37

# save
saveRDS(Merge_circexp_raw_filtered_and_enriched, file="Merge_circexplorer_CSF87.filtered.enriched.rawcount.rds")
saveRDS(Merge_circexp_norm_filtered_and_enriched, file="Merge_circexplorer_CSF87.filtered.enriched.normRPM.rds")
saveRDS(annotation_filtered_enriched, file="Merge_circexplorer_CSF87.filtered.enriched.annotation.bed14.rds")

Merge_circexp_norm_enriched=Merge_circexp_norm[intersect(rownames(Merge_circexp_raw), Merge_circexp_raw_long_enriched.RM$ID), ]
Merge_circexp_raw_enriched <- Merge_circexp_raw[rownames(Merge_circexp_norm_enriched), ]
annotation_enriched <- annotation[annotation$ID %in% rownames(Merge_circexp_raw_enriched),] 

dim(Merge_circexp_norm_enriched); dim(Merge_circexp_raw_enriched); dim(annotation_enriched)
table(annotation_enriched$circType)
# circRNA   ciRNA 
# 50       6 



##################################
## pilot (n=10) for grant purpose
##################################
Merge_circexp_raw = readRDS(file="Merge_circexplorer_CSF87.rawcount.rds")
head(Merge_circexp_raw)
plot(colSums(Merge_circexp_raw))
#top10samples = names(head(sort(colSums(Merge_circexp_raw),decreasing = T),10))
#set.seed(123); Merge_circexp_raw_pilot = Merge_circexp_raw[,sample(1:ncol(Merge_circexp_raw),10)]  # random 10 samples
#Merge_circexp_raw_pilot = Merge_circexp_raw[,top10samples]  # first 10 samples with top depth
Merge_circexp_raw_pilot = select(Merge_circexp_raw, contains("HC_"))  # first 10 samples with top depth
head(Merge_circexp_raw_pilot)
Merge_circexp_raw_pilot = Merge_circexp_raw_pilot[rowSums(Merge_circexp_raw_pilot)>=2,]
dim(Merge_circexp_raw_pilot)
# [1] 2496  10
head(sort(rowSums(Merge_circexp_raw_pilot),decreasing = T),20)
# How many are they also detected in SNDA?
sum(rownames(Merge_circexp_raw_pilot) %in% (readRDS("Merge_circexplorer_BC197.annotation.bed14.rds")$ID))
# [1] 275
rowSums(Merge_circexp_raw_pilot['chrX_139865339_139866824',])
# chrX_139865339_139866824 
# 29 

## any of them in PD DE genes?
annotation_filtered = readRDS("Merge_circexplorer_CSF87.filtered.annotation.bed14.rds")
Merge_circexp_raw_filtered_and_enriched = readRDS(file="Merge_circexplorer_CSF87.filtered.enriched.rawcount.rds")
DE=read.table("../results/DE_SNDA/DEresult.DE_SNDA.CONDITION2_PD_vs_HC.xls", header=T, sep="\t", row.names = 1, stringsAsFactors = F)  %>% rownames_to_column("ID")  %>% select(-contains("ind_raw")) # %>% filter(log2FoldChange>1 | log2FoldChange < -1)
head(DE)
# join
inner_join(x=annotation_filtered, y=DE, by="ID") %>% filter(pvalue<=0.05)  
# None

rowSums(Merge_circexp_raw_filtered_and_enriched['chrX_139865339_139866824',])
