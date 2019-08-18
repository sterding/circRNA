# prepare for eQTL using fastQTL
require(RCurl)
library(tidyr)
library(dplyr)

# n=84 HCILB_SNDA samples
sampleID=scan("~/neurogen/rnaseq_PD/results/merged/samplelist.HCILB_SNDA84",character())

########
# cov
########
covarianceTableURL="https://docs.google.com/spreadsheets/d/1I8nRImE9eJCCuZwpjfrrj-Uwx9bLebnO6o-ph7u6n8s/pub?gid=195725118&single=true&output=tsv"  # for all 140 samples
covs=read.delim(textConnection(getURL(covarianceTableURL)), stringsAsFactors =F)
covs=subset(covs, sampleName %in% sampleID, select = c(subjectID, batch, RIN, sex, age, PMI))
# convert batch to batchxxx
covs$batch=paste0("batch", covs$batch)
write.table(t(covs), "~/projects/circRNA/data/QTL_BC/covs.fastqtl.txt", col.names = F, quote = F, sep="\t", row.names = T)

########
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

# ########
# # circRNA
# ########
# # only for those high-confident circular RNAs (n=10017)
# Merge_circexp_raw_filtered = readRDS("~/projects/circRNA/data/Merge_circexplorer_BC106.filtered.enriched.rawcount.rds")
# Merge_circexp_norm_filtered =readRDS("~/projects/circRNA/data/Merge_circexplorer_BC106.filtered.enriched.normRPM.rds")
# annotation_filtered= readRDS("~/projects/circRNA/data/Merge_circexplorer_BC106.filtered.enriched.annotation.bed14.rds")
# 
# ## Filter: expression of >0.001 RPM in at least 5% individuals and ≥1 reads in at least 5% individuals 
# # (similar as GTEx: https://www.gtexportal.org/home/documentationPage#staticTextAnalysisMethods)
# 
# RPM_threshold = 0.001
# raw_threshold = 1
# PERCENT_CUTOFF = 0.05
# 
# df = Merge_circexp_norm_filtered[rowMeans(Merge_circexp_norm_filtered[,sampleID] >= RPM_threshold) >= PERCENT_CUTOFF & rowMeans(Merge_circexp_raw_filtered[,sampleID] >= raw_threshold) >= PERCENT_CUTOFF, sampleID]; dim(df) # n = 601 remained for downstream QTL analysis
# 
# annotation_filtered_refilterd4eqtl= annotation_filtered[as.character(annotation_filtered$ID) %in% rownames(df),]
# saveRDS(annotation_filtered_refilterd4eqtl, "~/projects/circRNA/data/Merge_circexplorer_BC.filtered.refiltered4eqtl.annotation.bed14.rds")
# 
# df = df %>% mutate(ID=rownames(.)) %>% 
#   separate(ID, c("Chr","Start","End"), remove=F) %>%
#   select(Chr, Start, End, ID, dplyr::contains("_")) %>% 
#   arrange(Chr, as.numeric(Start), as.numeric(End)) %>%
#   setNames(gsub(".*_(.*)_SNDA_.*","\\1",names(.)))
# # only circRNA or all circularRNA
# colnames(df)[1]=paste0("#",colnames(df)[1]) # add # at the beginning
# write.table(df,"~/projects/circRNA/data/QTL_BC/phenotype/circularRNA.expression.bed", row.names=F, quote=F, sep="\t")
# write.table(df[df$ID %in% annotation_filtered$ID[annotation_filtered$circType=='circRNA'],],"~/projects/circRNA/data/QTL_BC/phenotype/circRNA.expression.bed", row.names=F, quote=F, sep="\t")

########
# circular RNA (from 189128)
########
# from 189128
Merge_circexp_raw = readRDS("~/projects/circRNA/data/Merge_circexplorer_BC106.rawcount.rds")
Merge_circexp_norm= readRDS("~/projects/circRNA/data/Merge_circexplorer_BC106.normRPM.rds")
annotation= readRDS("~/projects/circRNA/data/Merge_circexplorer_BC106.annotation.bed14.rds")
annotation_filtered= readRDS("~/projects/circRNA/data/Merge_circexplorer_BC106.filtered.enriched.annotation.bed14.rds")

## Filter: expression of >0.001 RPM in at least 5% individuals and ≥1 reads in at least 5% individuals 
# (similar as GTEx: https://www.gtexportal.org/home/documentationPage#staticTextAnalysisMethods)

RPM_threshold = 0.001
raw_threshold = 1
PERCENT_CUTOFF = 0.05

df = Merge_circexp_norm[rowMeans(Merge_circexp_norm[,sampleID] >= RPM_threshold) >= PERCENT_CUTOFF & rowMeans(Merge_circexp_raw[,sampleID] >= raw_threshold) >= PERCENT_CUTOFF, sampleID]; dim(df) # n = 1054 remained for downstream QTL analysis

annotation_refilterd4eqtl= annotation[match(rownames(df), as.character(annotation$ID)),]
saveRDS(annotation_refilterd4eqtl, "~/projects/circRNA/data/Merge_circexplorer_BC106.refiltered4eqtl.annotation.bed14.rds")
table(annotation_refilterd4eqtl$circType)
# circRNA   ciRNA 
# 752     302 
## How many of them are also in the filter_enriched set of 10017 circular RNAs?
table(data.frame(circType=annotation_refilterd4eqtl$circType, enriched=annotation_refilterd4eqtl$ID %in% annotation_filtered$ID))
#           enriched
# circType  FALSE TRUE
# circRNA   186  566
# ciRNA     267   35

pdf("~/projects/circRNA/data/QTL_BC/phenotype/circularRNA.N84.expression.hist.pdf", width = 5, height = 3)
d=data.frame(n_samples=rowSums(df>=RPM_threshold), type=annotation_refilterd4eqtl$circType)
library(ggplot2); 
require(scales)
mylog_trans <- function (base = exp(1), from = 0) 
{
  trans <- function(x) log(x, base) - from
  inv <- function(x) base^(x + from)
  trans_new("mylog", trans, inv, log_breaks(base = base), domain = c(base^from, Inf))
}
ggplot(d, aes(x=n_samples)) + geom_histogram(binwidth = 1,color='#fdbb84',fill='#fdbb84') +  geom_histogram(data=subset(d,type=='circRNA'),binwidth = 1,color='#ff0000',fill='#ff0000') + scale_x_continuous(breaks=seq(0,85,5)) + scale_y_continuous(trans = mylog_trans(base=10, from=-1), breaks=c(1,10,100,1000),limits=c(0.1,500)) + geom_text(data=subset(d, n_samples > 40), aes(x=n_samples,y=1,label=rownames(subset(d, n_samples > 40))), angle=90, nudge_y=0.05, hjust=0) + theme_bw() + labs(x="Numbers of samples (out of 84 in total)", title='circular RNAs expressiong over the cutoff', y='# of circular RNAs') 
dev.off()

df = df %>% mutate(ID=rownames(.)) %>% 
  separate(ID, c("Chr","Start","End"), remove=F, convert=T) %>%
  select(Chr, Start, End, ID, dplyr::contains("_")) %>% 
  arrange(Chr, as.numeric(Start), as.numeric(End)) %>%
  setNames(gsub(".*_(.*)_SNDA_.*","\\1",names(.)))

# only circRNA or all circularRNA
colnames(df)[1]=ifelse(grepl("^#", colnames(df)[1]), colnames(df)[1], paste0("#",colnames(df)[1])) # add # at the beginning
write.table(df,"~/projects/circRNA/data/QTL_BC/phenotype/circularRNA.expression.bed", row.names=F, quote=F, sep="\t")
write.table(df[df$ID %in% annotation$ID[annotation$circType=='circRNA'],],"~/projects/circRNA/data/QTL_BC/phenotype/circRNA.expression.bed", row.names=F, quote=F, sep="\t")

## with ciRNAs/intron lariats from the same donor site merged into one group
df2=rbind(df %>% filter(ID %in% annotation$ID[annotation$circType=='circRNA']), df %>% filter(ID %in% annotation$ID[annotation$circType=='ciRNA']) %>% group_by(Chr, Start) %>% mutate(End=max(End)) %>% summarise_at(-c(1,2,4), mean) %>% ungroup() %>% group_by(Chr, End) %>% mutate(Start=min(Start)) %>% summarise_at(-c(1,3), mean) %>% ungroup() %>% unite("ID",c("Chr","Start","End"),remove=F) %>% select(Chr,Start,End,ID,matches(".[0-9]."))) %>% arrange(Chr,Start,End); dim(df2) # n=786
colnames(df2)[1]=ifelse(grepl("^#", colnames(df2)[1]), colnames(df2)[1], paste0("#",colnames(df2)[1])) # add # at the beginning
write.table(df2,"~/projects/circRNA/data/QTL_BC/phenotype/circularRNAmerged.expression.bed", row.names=F, quote=F, sep="\t")

########
## circularization
########
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

########
## splicing
########
# see fastQTL.main.sh to get splicing circRNA.PSI.chr*.bed