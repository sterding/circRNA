## Main script to process Bennett's RNAseq samples
library(tidyverse)
library(RCurl)

## See ~/projects/circRNA/README.Rmd for how to generate the input files
setwd("~/projects/circRNA/data/")

# load
Merge_circexp_raw <- read.table("Merge_circexplorer_Bennett.rawcount.long.txt", sep="\t",check.names =F, col.names = c("ID","readsCount","sampleID")) %>% 
  group_by(ID, sampleID) %>% summarise(readsCount=sum(readsCount)) %>%
  spread(key = sampleID, value = readsCount, fill = 0) %>%
  column_to_rownames("ID")
head(Merge_circexp_raw); dim(Merge_circexp_raw)

# annotation
annotation<- read.table("Merge_circexplorer_Bennett.annotation.bed14",header=T, check.names = F) # include both ciRNA and circRNA
head(annotation); dim(annotation)
# 88698

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
BRAINCODE2_Sequencing_Log.RNAseq_statistics="https://docs.google.com/spreadsheets/d/e/2PACX-1vQFQ4aQj0sD9oxIqaZ-cEgo7kWcCmNYGBH9emLw8iNu0f6TTjKE5Lte7IBfoMMy57cLjA4pXE0YlPY2/pub?gid=1322938821&single=true&output=tsv"
RNAseq_statistics<- read.delim(textConnection(getURL(BRAINCODE2_Sequencing_Log.RNAseq_statistics)), header=T, check.names = T, comment.char = "#", stringsAsFactors = F) 
readsNum_million <- as.numeric(gsub(",","",RNAseq_statistics$total_mapped_reads))/10^6; names(readsNum_million) = RNAseq_statistics$SampleID
readsNum_million = readsNum_million[colnames(Merge_circexp_raw)]

Merge_circexp_norm<-sweep(Merge_circexp_raw,2,readsNum_million,"/")

# save
saveRDS(Merge_circexp_raw, file="Merge_circexplorer_Bennett.rawcount.rds")
saveRDS(Merge_circexp_norm, file="Merge_circexplorer_Bennett.normRPM.rds")
saveRDS(annotation, file="Merge_circexplorer_Bennett.annotation.bed14.rds")

#write.table(Merge_circexp_raw, file="Merge_circexplorer_Bennett.rawcount.xls.gz", quote=F, sep="\t", col.names = NA, row.names = TRUE)
#write.table(Merge_circexp_norm, file="Merge_circexplorer_Bennett.normRPM.xls", quote=F, sep="\t", col.names = NA, row.names = TRUE)

## ===============================================
## back-splicing reads vs. number of samples
## similar to Fig. 3E in https://www-sciencedirect-com.ezp-prod1.hul.harvard.edu/science/article/pii/S0092867418316350?via%3Dihub#fig3
## ===============================================
df=Merge_circexp_raw
pdf("../results/Merge_circexplorer_Bennett.nReads_vs_nSample.pdf", width = 6, height = 6)
plot(jitter(rowMeans(df>0))*100, rowSums(df), xlab="Sample fraction (%) ", ylab="Backspliced reads", pch=21, col='gray', bg='#00000033', log='y')
dev.off()
## only few circRNAs expressed in 90% samples in our case. So, it's applicable to divide into 'high' vs. 'low' like the paper above. 

###########################################
######### filter cirRNAs (PD)    ##########
###########################################

### (0) only selected VMB brain samples  
Bennett_VMB=filter(RNAseq_statistics, Seq_type=="RNAseq", Tissue=='VMB') %>% pull(SampleID); length(Bennett_VMB)
Merge_circexp_raw_filtered <- Merge_circexp_raw[,Bennett_VMB]
Merge_circexp_raw_filtered <- Merge_circexp_raw_filtered[rowSums(Merge_circexp_raw_filtered)>0,]; dim(Merge_circexp_raw_filtered)
saveRDS(Merge_circexp_raw_filtered, file="Merge_circexplorer_Bennett_VMB.rawcount.rds")
saveRDS(annotation[annotation$ID %in% rownames(Merge_circexp_raw_filtered),], file="Merge_circexplorer_Bennett_VMB.annotation.bed14.rds")

### (1) being expressed: with >=2 unique back-splicing reads in all samples (Zheng et al., doi:10.1038/ncomms11215)
### --------------------

Merge_circexp_raw_filtered <- Merge_circexp_raw_filtered[rowSums(Merge_circexp_raw_filtered)>=2,]
Merge_circexp_norm_filtered <- Merge_circexp_norm[rownames(Merge_circexp_raw_filtered), Bennett_VMB]
annotation_filtered <- annotation[annotation$ID %in% rownames(Merge_circexp_norm_filtered),] 

dim(Merge_circexp_raw_filtered); dim(Merge_circexp_norm_filtered); dim(annotation_filtered)
# [1] 37730   23

# save
saveRDS(Merge_circexp_raw_filtered, file="Merge_circexplorer_Bennett_VMB.filtered.rawcount.rds")
saveRDS(Merge_circexp_norm_filtered, file="Merge_circexplorer_Bennett_VMB.filtered.normRPM.rds")
saveRDS(annotation_filtered, file="Merge_circexplorer_Bennett_VMB.filtered.annotation.bed14.rds")
#write.table(annotation_filtered, file="Merge_circexplorer_Bennett_VMB.filtered.annotation.bed14", sep = "\t", quote = F, row.names = F, col.names = F)

table(readRDS("Merge_circexplorer_BC190.filtered.enriched.annotation.bed14.rds")$ID %in% readRDS("Merge_circexplorer_Bennett_VMB.filtered.annotation.bed14.rds")$ID)
# FALSE  TRUE 
# 2207  7975 

###########################################
######### filter cirRNAs (AD)    ##########
###########################################

### (0) only selected FCX brain samples  
Bennett_FCX=filter(RNAseq_statistics, Seq_type=="RNAseq", Tissue=='FCX') %>% pull(SampleID); length(Bennett_FCX)
Merge_circexp_raw_filtered <- Merge_circexp_raw[,Bennett_FCX]
Merge_circexp_raw_filtered <- Merge_circexp_raw_filtered[rowSums(Merge_circexp_raw_filtered)>0,]; dim(Merge_circexp_raw_filtered)
saveRDS(Merge_circexp_raw_filtered, file="Merge_circexplorer_Bennett_FCX.rawcount.rds")
saveRDS(annotation[annotation$ID %in% rownames(Merge_circexp_raw_filtered),], file="Merge_circexplorer_Bennett_FCX.annotation.bed14.rds")

### (1) being expressed: with >=2 unique back-splicing reads in all samples (Zheng et al., doi:10.1038/ncomms11215)
### --------------------

Merge_circexp_raw_filtered <- Merge_circexp_raw_filtered[rowSums(Merge_circexp_raw_filtered)>=2,]
Merge_circexp_norm_filtered <- Merge_circexp_norm[rownames(Merge_circexp_raw_filtered), Bennett_FCX]
annotation_filtered <- annotation[annotation$ID %in% rownames(Merge_circexp_norm_filtered),] 

dim(Merge_circexp_raw_filtered); dim(Merge_circexp_norm_filtered); dim(annotation_filtered)
# [1] 23376   19

# save
saveRDS(Merge_circexp_raw_filtered, file="Merge_circexplorer_Bennett_FCX.filtered.rawcount.rds")
saveRDS(Merge_circexp_norm_filtered, file="Merge_circexplorer_Bennett_FCX.filtered.normRPM.rds")
saveRDS(annotation_filtered, file="Merge_circexplorer_Bennett_FCX.filtered.annotation.bed14.rds")
#write.table(annotation_filtered, file="Merge_circexplorer_Bennett_FCX.filtered.annotation.bed14", sep = "\t", quote = F, row.names = F, col.names = F)

table(readRDS("Merge_circexplorer_BC190.filtered.enriched.annotation.bed14.rds")$ID %in% readRDS("Merge_circexplorer_Bennett_FCX.filtered.annotation.bed14.rds")$ID)
# FALSE  TRUE 
# 3174  7008 
