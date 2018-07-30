# prepare for eQTL using fastQTL
require(RCurl)
library(tidyr)
library(dplyr)

# n=84 HCILB_SNDA samples
sampleID=scan("~/neurogen/rnaseq_PD/results/merged/samplelist.HCILB_SNDA84",character())

# cov
########
covarianceTableURL="https://docs.google.com/spreadsheets/d/1I8nRImE9eJCCuZwpjfrrj-Uwx9bLebnO6o-ph7u6n8s/pub?gid=195725118&single=true&output=tsv"  # for all 140 samples
covs=read.delim(textConnection(getURL(covarianceTableURL)), stringsAsFactors =F)
covs=subset(covs, sampleName %in% sampleID, select = c(subjectID, batch, RIN, sex, age, PMI))
# convert batch to batchxxx
covs$batch=paste0("batch", covs$batch)
write.table(t(covs), "~/projects/circRNA/data/QTL_BC/covs.fastqtl.txt", col.names = F, quote = F, sep="\t", row.names = T)

## expression
########
# genes
df1=read.table("~/neurogen/rnaseq_PD/results/merged/genes.loci.txt", header=T)
df2=read.table("~/neurogen/rnaseq_PD/results/merged/genes.fpkm.cuffnorm.allSamples.uniq.xls", header=T, check.names = F)
df2=df2[,c('tracking_id', sampleID)]
# filter: >0.1 RPKM in at least 10 individuals (32813 out of 57816 remained)
df2=df2[rowSums(df2[,sampleID]>0.1)>=10,]
df = merge(df1,df2,by='tracking_id')[,c(2:4,1,5:88)]
#df$chr=gsub("chr","", df$chr)
df=df[with(df, order(chr, s1, s2)),]
colnames(df)=gsub(".*_(.*)_SNDA_.*","\\1",colnames(df))
colnames(df)[1]=paste0("#",colnames(df)[1]) # vi to add # at the beginning
write.table(df,"~/projects/circRNA/data/QTL_BC/phenotype/genes.expression.bed", row.names=F, quote=F, sep="\t")

# circRNA
########
Merge_circexp_raw_filtered = readRDS("~/projects/circRNA/data/Merge_circexplorer_BC.filtered.rawcount.rds")
Merge_circexp_norm_filtered =readRDS("~/projects/circRNA/data/Merge_circexplorer_BC.filtered.normRPM.rds")
annotation_filtered= readRDS("~/projects/circRNA/data/Merge_circexplorer_BC.filtered.annotation.bed14.rds")

## Filter: expression of >0.001 RPM in at least 10% individuals and ≥1 reads in at least 10% individuals 
# (similar as GTEx: https://www.gtexportal.org/home/documentationPage#staticTextAnalysisMethods)

RPM_threshold = 0.001
raw_threshold = 1
PERCENT_CUTOFF = 0.1

df=Merge_circexp_norm_filtered[rowMeans(Merge_circexp_raw_filtered[,sampleID]) > 0, sampleID]; dim(df) # only those expressed in 84 HCILB_SNDA

df = Merge_circexp_norm_filtered[rowMeans(Merge_circexp_norm_filtered[,sampleID] >= RPM_threshold) >= PERCENT_CUTOFF & rowMeans(Merge_circexp_raw_filtered[,sampleID] >= raw_threshold) >= PERCENT_CUTOFF, sampleID]; dim(df) # n = 374 remained for downstream QTL analysis

annotation_filtered_refilterd4eqtl= annotation_filtered[as.character(annotation_filtered$ID) %in% rownames(df),]
saveRDS(annotation_filtered_refilterd4eqtl, "~/projects/circRNA/data/Merge_circexplorer_BC.filtered.refiltered4eqtl.annotation.bed14.rds")


cairo_pdf("~/projects/circRNA/data/QTL_BC/phenotype/circularRNA.N84.expression.hist.pdf", width = 5, height = 5)
par(mar=c(4,6,4,2)); hist(rowSums(df>0), breaks=100, xlab="Numbers of samples (out of 84 in total)", main='circRNAs expressiong over the cutoff', ylab='# of circRNAs with expression \n ≥0.001 RPM and ≥1 reads in ≥ 10% individuals')
dev.off()

df = df %>% mutate(ID=rownames(.)) %>% 
  separate(ID, c("Chr","Start","End"), remove=F) %>%
  select(Chr, Start, End, ID, dplyr::contains("_")) %>% 
  arrange(Chr, as.numeric(Start), as.numeric(End)) %>%
  setNames(gsub(".*_(.*)_SNDA_.*","\\1",names(.)))
# only circRNA or all circularRNA
colnames(df)[1]=paste0("#",colnames(df)[1]) # vi to add # at the beginning
write.table(df,"~/projects/circRNA/data/QTL_BC/phenotype/circularRNA.expression.bed", row.names=F, quote=F, sep="\t")
write.table(df[df$ID %in% annotation_filtered$ID[annotation_filtered$circType=='circRNA'],],"~/projects/circRNA/data/QTL_BC/phenotype/circRNA.expression.bed", row.names=F, quote=F, sep="\t")

## circularization
df = readRDS("~/projects/circRNA/data/Merge_circexplorer_BC.annotation.bed14.cRatio.rds") # already filtered for 84 HCILB_SNDA samples
df = df[,sampleID]; dim(df)  # 21949    84
# only for the circRNA passed the filter for eQTL above [optional?]
df = df[rownames(df) %in% as.character(annotation_filtered_refilterd4eqtl$ID),]; dim(df)  # 278  84
df = as.data.frame(df) %>% mutate(ID=rownames(df)) %>% 
  separate(ID, c("Chr","Start","End"), remove=F) %>%
  select(Chr, Start, End, ID, dplyr::contains("_")) %>% 
  arrange(Chr, as.numeric(Start), as.numeric(End)) %>%
  setNames(gsub(".*_(.*)_SNDA_.*","\\1",names(.)))
# only circRNA or all circularRNA
colnames(df)[1]=paste0("#",colnames(df)[1]) # vi to add # at the beginning
write.table(df,"~/projects/circRNA/data/QTL_BC/phenotype/circularRNA.PCI.bed", row.names=F, quote=F, sep="\t")
write.table(df[df$ID %in% annotation_filtered$ID[annotation_filtered$circType=='circRNA'],],"~/projects/circRNA/data/QTL_BC/phenotype/circRNA.PCI.bed", row.names=F, quote=F, sep="\t")

## splicing
# see main.sh to get splicing circRNA.PSI.chr*.bed