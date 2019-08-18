## Main script to work on the output of circExplorer
library('dplyr')
library('reshape2')
library('tidyverse')
library('ggpubr') # install.packages("ggpubr")
setwd("~/projects/circRNA/data/") 

## See ~/projects/circRNA/src/README.Rmd for how to generate the input files: 
# Merge_circexplorer_BC.annotation.bed14
# Merge_circexplorer_BC.rawcount.txt

## see pilot.R in ~/Dropbox/grant/2019R21/pilot.R to make figure for pilot study (for grant purpose)

###########################################
################# load data  ##############
###########################################

## === annotation ===

annotation<- read.table("Merge_circexplorer_BC.annotation.bed14",header=T, check.names = F) # include both ciRNA and circRNA
head(annotation)

## === load raw reads count ===

#removed "#" on first line as well as ful name to just sample name
Merge_circexp_raw <- read.table("Merge_circexplorer_BC.rawcount.txt",sep="\t",check.names =F, header=TRUE, comment.char ='') 
colnames(Merge_circexp_raw)=gsub("#|_candidates.bed_circReads","",colnames(Merge_circexp_raw))

Merge_circexp_raw = mutate(Merge_circexp_raw, concated_column = paste(chrom,start,end, sep = '_'))
row.names(Merge_circexp_raw) = Merge_circexp_raw$concated_column
Merge_circexp_raw=select(Merge_circexp_raw,dplyr::contains('rep'))
dim(Merge_circexp_raw)
#[1] 211772    125

readsNum_filtered<- read.table("~/neurogen/rnaseq_PD/run_output/linescounts.filtered.txt",row.names=1, header=F, check.names = F) 
readsNum_million<-(t(readsNum_filtered)[1,]/10^6)
#only include the 125 (106 HCILB + 19 PD) samples
readsNum_million=readsNum_million[colnames(Merge_circexp_raw)]

###
### below code is to limit to 106 HC+ILB samples, otherwise 125 samples including PD+ILB+HC
###
sample106=scan("~/neurogen/circRNA_seq_Rebeca_HC/BRAINCODE_circexp/BC.n106.samplelist",character())
readsNum_million=readsNum_million[sample106]
Merge_circexp_raw=Merge_circexp_raw[,sample106]
Merge_circexp_raw=Merge_circexp_raw[rowSums(Merge_circexp_raw)>0,]; dim(Merge_circexp_raw)
#[1] 189128    106

# test
Merge_circexp_raw=readRDS(file="~/projects/circRNA/data/Merge_circexplorer_BC106.rawcount.rds")
sample84=scan("~/neurogen/rnaseq_PD/results/merged/samplelist.HCILB_SNDA84",character())
df=read.table("~/projects/circRNA/data/SNPs.on.circRNA.splicingSites.list", header = F, stringsAsFactors = F)
Merge_circexp_raw=Merge_circexp_raw[,sample84]
sum(rowSums(Merge_circexp_raw>0)>4)
Merge_circexp_raw['chr2_40655612_40657444',]
## for Jiajie
df=read.table("~/projects/circRNA/data/QTL_BC/eQTL.nominal.txt.chr15.p1e-5.list", header = F, stringsAsFactors = F)
write.table(Merge_circexp_raw[unique(df$V1),sample84], file="~/tools/SplicePlot2/test_files/chr15.circRNA.rawcount.txt", quote=F, sep="\t", col.names = T)

###########################################
############## normalize to RPM ###########
###########################################

Merge_circexp_norm<-sweep(Merge_circexp_raw,2,readsNum_million,"/")

# save
saveRDS(Merge_circexp_raw, file="Merge_circexplorer_BC106.rawcount.rds")
saveRDS(Merge_circexp_norm, file="Merge_circexplorer_BC106.normRPM.rds")
saveRDS(annotation[annotation$ID %in% rownames(Merge_circexp_norm),], file="Merge_circexplorer_BC106.annotation.bed14.rds")

write.table(Merge_circexp_raw, file="Merge_circexplorer_BC106.rawcount.xls", quote=F, sep="\t", col.names = NA, row.names = TRUE)
write.table(Merge_circexp_norm, file="Merge_circexplorer_BC106.normRPM.xls", quote=F, sep="\t", col.names = NA, row.names = TRUE)

###########################################
############ filter cirRNAs     ###########
###########################################

# What can be a good filter?

### (1) being expressed
### --------------------

# # Definition #1 of being expression: with RPM >= threshold in >=2 samples
# pdf("~/projects/circRNA/results/filter_QC.pdf", width = 6, height = 5)
# boxplot(1/readsNum_million, outcol="NA", ylab="1 RPM cutoff for each sample"); 
# points(jitter(rep(1, length(readsNum_million))), 1/readsNum_million, pch=20, col=rgb(0,0,0,.6)) 
# median(1/readsNum_million) # 0.006118284  --> 0.006
# dev.off()
# threshold<- signif(median(1/readsNum_million), 1)  # 0.006 
# N = ncol(Merge_circexp_norm)
# # with cutoff of 1/(min(readsNum_million)), 4030 remained; while RPM>=0.006, 9161 remained.
# Merge_circexp_norm_filtered <- Merge_circexp_norm[rowSums(Merge_circexp_norm >= threshold) >= 2, ]
# dim(Merge_circexp_norm_filtered)
# # [1] 10279  125

# Final definition of being expression: with >=2 unique back-splicing reads in all samples (Zheng et al., doi:10.1038/ncomms11215)
Merge_circexp_raw_filtered <- Merge_circexp_raw[rowSums(Merge_circexp_raw)>=2, ]
Merge_circexp_norm_filtered <- Merge_circexp_norm[rowSums(Merge_circexp_raw)>=2, ]
annotation_filtered <- annotation[annotation$ID %in% rownames(Merge_circexp_norm_filtered),] 

dim(Merge_circexp_raw_filtered); dim(Merge_circexp_norm_filtered); dim(annotation_filtered)
# [1] 65186   125
# [1] 57130   106


# save
saveRDS(Merge_circexp_raw_filtered, file="Merge_circexplorer_BC106.filtered.rawcount.rds")
saveRDS(Merge_circexp_norm_filtered, file="Merge_circexplorer_BC106.filtered.normRPM.rds")
saveRDS(annotation_filtered, file="Merge_circexplorer_BC106.filtered.annotation.bed14.rds")
write.table(annotation_filtered, file="Merge_circexplorer_BC106.filtered.annotation.bed14", sep = "\t", quote = F, row.names = F, col.names = F)

Merge_circexp_raw=readRDS(file="Merge_circexplorer_BC106.rawcount.rds")
Merge_circexp_norm=readRDS(file="Merge_circexplorer_BC106.normRPM.rds")
Merge_circexp_raw_filtered=readRDS(file="Merge_circexplorer_BC106.filtered.rawcount.rds")
Merge_circexp_norm_filtered=readRDS(file="Merge_circexplorer_BC106.filtered.normRPM.rds")
annotation_filtered=readRDS(file="Merge_circexplorer_BC106.filtered.annotation.bed14.rds")
head(Merge_circexp_raw_filtered)
Merge_circexp_raw_filtered %>% mutate(rowsum_SNDA=rowSums(select(.,contains("_SNDA_")))) %>% filter(rowsum_SNDA>0) %>% dim()

###########################################
############ being enriched     ###########
###########################################

# Definition of being enriched: at least 20 reads in RNase R treatment AND foldchange(RNase vs Mock) at least 2 in at least one sample. 
# Run main.RM.R first
Merge_circexp_raw_long_enriched.RM  = readRDS(file="Merge_circexp_raw_long_filterd.RM.rds")
length(unique(Merge_circexp_raw_long_enriched.RM$name))
# [1] 48016
table(unique(Merge_circexp_raw_long_enriched.RM$name[Merge_circexp_raw_long_enriched.RM$sampleGroup3!='NN']) %in% rownames(Merge_circexp_raw_filtered))
# FALSE  TRUE 
# 28793  9206

Merge_circexp_norm_filtered_and_enriched=Merge_circexp_norm[intersect(rownames(Merge_circexp_raw_filtered), Merge_circexp_raw_long_enriched.RM$name), ]
Merge_circexp_raw_filtered_and_enriched <- Merge_circexp_raw[rownames(Merge_circexp_norm_filtered_and_enriched), ]
annotation_filtered_enriched <- annotation[annotation$ID %in% rownames(Merge_circexp_raw_filtered_and_enriched),] 

dim(Merge_circexp_norm_filtered_and_enriched); dim(Merge_circexp_raw_filtered_and_enriched); dim(annotation_filtered_enriched)
# [1] 10431   125
# [1] 10017   106

# save
saveRDS(Merge_circexp_raw_filtered_and_enriched, file="Merge_circexplorer_BC106.filtered.enriched.rawcount.rds")
saveRDS(Merge_circexp_norm_filtered_and_enriched, file="Merge_circexplorer_BC106.filtered.enriched.normRPM.rds")
saveRDS(annotation_filtered_enriched, file="Merge_circexplorer_BC106.filtered.enriched.annotation.bed14.rds")
annotation_filtered_enriched=readRDS(file="Merge_circexplorer_BC106.filtered.enriched.annotation.bed14.rds")
write.table(annotation_filtered_enriched, file="Merge_circexplorer_BC106.filtered.enriched.annotation.bed14", sep = "\t", quote = F, row.names = F, col.names = F)

#########################################
## Figure 1b: distribution of circRNAs supported by different number of reads
###########################################

# pie chart of circRNA among all circular RNAs
pdf("~/projects/circRNA/results/total_circRNA_raw_reads.BC106.pie.pdf", width=6, height = 2)
par(mfrow=c(1,3))
pie(table(readRDS("Merge_circexplorer_BC106.annotation.bed14.rds")$circType=='circRNA'),labels = NA, border='white', col=c('#fdbb84','red'),main="all")
pie(table(readRDS("Merge_circexplorer_BC106.filtered.annotation.bed14.rds")$circType=='circRNA'),labels = NA, border='white', col=c('#fdbb84','red'),main="filtered")
pie(table(readRDS("Merge_circexplorer_BC106.filtered.enriched.annotation.bed14.rds")$circType=='circRNA'),labels = NA, border='white', col=c('#fdbb84','red'),main="filtered+enriched")
dev.off()

Merge_circexp_raw=readRDS("Merge_circexplorer_BC106.rawcount.rds")
Merge_circexp_raw_filtered_and_enriched=readRDS("Merge_circexplorer_BC106.filtered.enriched.rawcount.rds")
total_circRNA_raw_reads = rowSums(Merge_circexp_raw)
pdf("~/projects/circRNA/results/total_circRNA_raw_reads.BC106.pdf", width = 6, height = 5)
layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(5,1), heights=c(1,1))
ht=hist(total_circRNA_raw_reads[total_circRNA_raw_reads<=200], breaks = 0:200, plot=F)
par(mar=c(4,4,0,1))
plot(log10(ht$counts+.1)+1, lwd=2,xlim=c(0,200), type='h', ylim = c(0.9,6.6),col=c("#000000",rep('#aaaaaa',199)), ylab='Number of circular RNAs', xlab="Number of backspliced reads",yaxt="n",xaxs="i", yaxs="i",bty="n")
y1=floor(range(log10(ht$counts+.1)+1))
pow <- seq(y1[1], y1[2])
ticksat <- as.vector(sapply(pow-1, function(p) log10((2:9)*10^p)))
axis(2, pow, labels = 10^(pow-1), las=1)
axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
# add the remained ones with the #2 filter: enriched
filtered_circRNA_raw_reads = rowSums(Merge_circexp_raw[rownames(Merge_circexp_raw_filtered_and_enriched),])
ht=hist(filtered_circRNA_raw_reads[filtered_circRNA_raw_reads<=200], breaks = 0:200, plot=F)
points(log10(ht$counts+.1)+1, lwd=2, type='h',col='#8B0000')

# add the remained ones with the #3 filter: circRNA only
circRNAonly=as.character(annotation_filtered_enriched$ID[annotation_filtered_enriched$circType=='circRNA'])
filtered_circRNA_raw_reads_circRNAonly = rowSums(Merge_circexp_raw[circRNAonly,])
ht=hist(filtered_circRNA_raw_reads_circRNAonly[filtered_circRNA_raw_reads_circRNAonly<=200], breaks = 0:200, plot=F)
points(log10(ht$counts+.1)+1, lwd=2, type='h',col='#ff0000')

legend("topleft",c(paste0("Distinct circular RNAs (n = ",format(length(total_circRNA_raw_reads),big.mark=","),")"),
                   paste0("- at least 2 reads in overall samples (n = ", format(sum(total_circRNA_raw_reads>1),big.mark=","),")"),
                   paste0("-- enriched in RNase R (n = ", format(length(filtered_circRNA_raw_reads),big.mark=","),")"),
                   paste0("--- circRNA only (n = ", format(length(filtered_circRNA_raw_reads_circRNAonly),big.mark=","),")")),
       col=c("#000000",'#aaaaaa','#8B0000','#ff0000'),
       text.col = c("#000000",'#aaaaaa','#8B0000','#ff0000'),
       bty='n', xpd=T)

## for the 200- region
par(mar=c(4,1,0,1))
n200=total_circRNA_raw_reads[total_circRNA_raw_reads>200]; 
MAX=max(total_circRNA_raw_reads)
ht=hist(n200, breaks = 200:MAX, plot=F)
plot(log10(ht$counts+.1)+1, lwd=2,xlim=range(total_circRNA_raw_reads), type='h', ylim = c(0.9,6.6),col='#aaaaaa', ylab="",xlab="",xaxt="n",yaxt='n',xaxs="i", yaxs="i",bty="n")
y1=floor(range(log10(ht$counts+.1)+1))
pow <- seq(y1[1], y1[2])
ticksat <- as.vector(sapply(pow-1, function(p) log10((2:9)*10^p)))
#axis(2, pow, labels = 10^(pow-1), las=1)
#axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
axis(1, c(0,seq(2000,MAX,10000)-200), labels=c(200,seq(2000,MAX,10000)))

# add the remained ones with the #2 filter
n200=filtered_circRNA_raw_reads[filtered_circRNA_raw_reads>200];
ht=hist(n200, breaks = 200:MAX, plot=F)
points(log10(ht$counts+.1)+1, lwd=2, type='h',col='#8B0000')

# add the remained ones with the #3 filter: circRNA only
ht=hist(filtered_circRNA_raw_reads_circRNAonly[filtered_circRNA_raw_reads_circRNAonly>200], breaks = 200:MAX, plot=F)
points(log10(ht$counts+.1)+1, lwd=2, type='h',col='#ff0000')

dev.off()

## add gene highlight in the figure
annotation_filtered_enriched %>% select(ID, circType, geneName) %>% 
  left_join(rownames_to_column(Merge_circexp_raw_filtered_and_enriched), by = c("ID"="rowname")) %>% 
  mutate(sumVar = rowSums(select(., contains("_")))) %>% select(-contains("_")) %>% 
  filter(circType=='circRNA') %>% arrange(-sumVar) %>% 
  #filter(sumVar<=200, sumVar>150)
  filter(geneName %in% c("SORL1", "DNAJC6","PSEN1","APP", "SNCA","LRRK2","PINK1","RIMS2","ABCA7","CLU","CR1","PICALM","PLD3","TREM2",'PARK7', 'PRKN', 'GBA', 'UCHL1','SRY'))

## ===============================================
## back-splicing reads vs. sample fraction
## similar to Fig. 3E in https://www-sciencedirect-com.ezp-prod1.hul.harvard.edu/science/article/pii/S0092867418316350?via%3Dihub#fig3
## ===============================================
Merge_circexp_raw_filtered_and_enriched=readRDS("Merge_circexplorer_BC106.filtered.enriched.rawcount.rds")
dim(Merge_circexp_raw_filtered_and_enriched)
df=Merge_circexp_raw_filtered_and_enriched
plot(jitter(rowMeans(df>0)), jitter(rowMeans(df)), pch=20, col='#00000033')
## only few circRNAs expressed in 90% samples in our case. So, it's applicable to divide into 'high' vs. 'low' like the paper above. 

## ===============================================
## Number of circRNAs per host gene
## ===============================================
require(scales)
mylog_trans <- function (base = exp(1), from = 0) 
  {
    trans <- function(x) log(x, base) - from
    inv <- function(x) base^(x + from)
    trans_new("mylog", trans, inv, log_breaks(base = base), domain = c(base^from, Inf))
  }
pdf("~/projects/circRNA/results/circRNAs_per_hostgene.histogram.pdf",width = 5, height = 3)
annotation_filtered_enriched_n = annotation_filtered_enriched %>% filter(circType=='circRNA') %>% group_by(geneName) %>% summarise(n=n()) %>% group_by(n) %>% summarise(N=n(), geneNames=paste(geneName, collapse="; "))
ggplot(annotation_filtered_enriched_n, aes(x=n, y=N)) + 
  geom_col() + 
  geom_text(data=subset(annotation_filtered_enriched_n, n > 21), aes(x=n,y=N,label=geneNames), angle=90, nudge_y=0.05, hjust=0) +
  xlab("Number of circRNAs in the host gene") + ylab("Count of host genes") +
  ggtitle("Histogram of number of circRNAs per host gene") +
  scale_y_continuous(trans = mylog_trans(base=10, from=-1), breaks=c(1,10,100,1000),limits=c(0.1,1500)) +
  scale_x_continuous(breaks=pretty_breaks())
dev.off()

## vs. # of exon per host genes
longestTx_annotation=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v19.annotation.longestTx.bed12", header = F, col.names = c('chr','start','end','ID','score','strand','gstart','gend','rgb','nExon','lengths','starts'), check.names = F, stringsAsFactors=F)
longestTx_annotation = separate(longestTx_annotation, ID,c('geneName','ensID',NA,NA),sep='___')
annotation_filtered_enriched$nExon=longestTx_annotation$nExon[match(annotation_filtered_enriched$geneName, longestTx_annotation$geneName)]
head(annotation_filtered_enriched)
# remove NA case caused by inconsistance between RefSeq and GENCODE gene symbol
annotation_filtered_enriched %>% filter(!is.na(nExon), circType=='circRNA') %>% group_by(geneName) %>% summarise(nExon=mean(nExon), nCircRNA=n()) %>% ungroup() %>% group_by(nExon_interval=cut(nExon, breaks = c(0,10,20,30,40,1000))) %>% ggplot(aes(nExon_interval,nCircRNA)) + geom_boxplot() + theme_classic() + ggsave("Merge_circexplorer_BC106.filtered.enriched.nCircRNA.vs.nExon.pdf", width = 3, height = 3)

## ===============================================
## exon and flanking intron length of circRNAs vs. all exons
## ===============================================
annotation_filtered_enriched=readRDS(file="Merge_circexplorer_BC106.filtered.enriched.annotation.bed14.rds") %>% filter(circType=='circRNA') %>% select(1:12)
introns=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/introns.meta.bed", header = F, col.names = c('chr','start','end','ID','score','strand'), check.names = F, stringsAsFactors=F) %>% mutate(score=end-start) %>% unite("chrend",c("chr",'end'), remove =F) %>% unite("chrstart",c("chr",'start'), remove =F)
head(introns)
controls=read.table("Merge_circexplorer_BC.annotation.bed14.matched2", header = F, col.names = c('chrom','start','end','ID','score','strand','thickStart','thickEnd','itemRgb','exonCount','exonSizes','exonOffsets','exonindex','hostgene','matchedCircRNA'), stringsAsFactors = F) %>% filter(matchedCircRNA %in% annotation_filtered_enriched$ID) %>% select(1:12) %>% mutate(itemRgb='0,0,0')
head(controls)
head(annotation_filtered_enriched)

# flanking intron lenth comparsion between called circRNAs and controls
df=bind_rows(mutate(annotation_filtered_enriched,type="circRNA"),mutate(controls,type="control")) %>% 
  unite('left_intron', c("chrom", "start"), remove =F) %>% 
  unite('right_intron', c("chrom", "end"), remove =F) %>% 
  mutate(left_intron_length=introns$score[match(left_intron, introns$chrend)], 
         right_intron_length=introns$score[match(right_intron, introns$chrstart)]) %>% 
  rowwise() %>% mutate(mean_intron_length=mean(c(left_intron_length,right_intron_length), na.rm = T))  
ggplot(df,aes(x=type, y=mean_intron_length, col=type)) + geom_boxplot() + scale_y_log10() + theme_classic() + 
  ggtitle(paste("Wilcox test P value =",signif(wilcox.test(mean_intron_length ~ as.factor(type), df)$p.value,2))) +
  scale_color_manual(guide=FALSE, name="type",values=c("circRNA"="red","control" = "grey"))  
ggsave("../results/Merge_circexplorer_BC106.filtered.enriched.annotation.lengthDistribution.flankingIntron.pdf", width = 2,height = 3)

# single-exon circRNAs lenth vs. random exon controls
exons=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/exons.meta.bed", sep="\t", col.names = c("chrom","start","end","ID","score","strand"), stringsAsFactors = F, header = F) %>% filter(grepl("protein_coding", ID)) %>% select(1:3) %>% sample_n(5000) %>% mutate(type="background")
df=readRDS(file="Merge_circexplorer_BC106.filtered.enriched.annotation.bed14.rds") %>% 
  filter(exonCount==1, circType=='circRNA') %>% select(1:3) %>% mutate(type="circRNA") %>% bind_rows(exons) %>%
  mutate(exon_length=end-start)
ggplot(df, aes(x=type, y=exon_length, col=type)) + geom_boxplot() + scale_y_log10() + theme_classic() + 
  ggtitle(paste("Wilcox test P value =",signif(wilcox.test(exon_length ~ as.factor(type), df)$p.value,2))) +
  scale_color_manual(guide=FALSE, name="type",values=c("circRNA"="red","background" = "grey")) 
ggsave("../results/Merge_circexplorer_BC106.filtered.enriched.annotation.lengthDistribution.single-exon-circRNAs.pdf", width = 2,height = 3)

# region length distribution
readRDS(file="Merge_circexplorer_BC106.filtered.enriched.annotation.bed14.rds") %>% 
  ggplot(aes(x=end-start, fill=factor(circType, levels=c('ciRNA', 'circRNA')))) + 
  geom_histogram(position = 'stack',bins = 60)  + scale_color_brewer(palette="Dark2") + 
  theme_minimal()+theme_classic()+theme(legend.position="top")+ 
  scale_x_log10() + scale_fill_manual(name="circular RNA types",values=c("circRNA"="red","ciRNA" = "orange")) 
ggsave("../results/Merge_circexplorer_BC106.filtered.enriched.annotation.lengthDistribution.pdf", width = 4,height = 3)


# lenth distribution for called circRNAs
filter(df, circRNA %in% BC106$ID, type=="circularized exon") %>% separate(circRNA, c(NA,"start","end"), sep = "_", convert=T) %>% ggplot(aes(x=end-start, fill=as.factor(cirType))) + geom_histogram(position = 'stack',bins = 60)  + scale_color_brewer(palette="Dark2") + theme_minimal()+theme_classic()+theme(legend.position="top")+ scale_x_log10() + scale_fill_manual(values=c("orange", "red")) + ggsave("../results/Merge_circexplorer_BC106.annotation.lengthDistribution.pdf", width = 4,height = 3)


## ===============================================
## host genes function enrichment analysis
## ===============================================
## fold change, meanR, pvalue from T test (see main.RM.R)
Merge_circexp_raw_RM_ttest=readRDS(file="Merge_circexp_raw_RM_ttest.rds")
# get the subset in Merge_circexp_norm_filtered_and_enriched, then convert to gene based format
Merge_circexp_raw_RM_ttest %>% group_by(symbol) %>% top_n(n=1, wt=meanR) %>% top_n(n=-1, wt=pvalue) %>% top_n(n=1, wt=log2FoldChange) %>% ungroup() %>% dplyr::select(symbol, padj=pvalue, log2FoldChange, stat=meanR) %>% distinct() %>% write.table(file="Merge_circexp_raw_RM_ttest.2genes.stat.xls", quote =F, sep="\t", row.names = T)

## then run the pathway analysis
# /data/rnaseq/src/_pathway.R -i Merge_circexp_raw_RM_ttest.2genes.stat.xls --filter medium

## Run topGO analysis for the circRNA-hosting genes
## ------------------------
readRDS("Merge_circexplorer_BC106.filtered.enriched.annotation.bed14.rds") %>% filter(circType=='circRNA') %>% group_by(geneName) %>% summarise(n=n()) %>% dplyr::select(symbol=geneName, stat=n) %>% write.table(file="Merge_circexp_norm_filtered_and_enriched.2genes.stat.xls", quote =F, sep="\t", row.names = T)

# bash
# ~/neurogen/pipeline/RNAseq/modules/_pathway.R -i Merge_circexp_norm_filtered_and_enriched.2genes.stat.xls --filter no

## ORA of circRNA-hosting genes vs. all expressed genes
## ------------------------
genes_fpkm=read.table("~/neurogen/rnaseq_PD/results/merged/genes.fpkm.cufflinks.allSamples.uniq.xls", header=T,row.names = 1, check.names = F)
genes_fpkm=genes_fpkm[,-c(1:7)]; colnames(genes_fpkm)=gsub("FPKM.","",colnames(genes_fpkm))
genes_annotation=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed", header = F, col.names = c('chr','start','end','EnsID','score','strand','symbol','type'), check.names = F, stringsAsFactors=F)
sample106=scan("~/neurogen/circRNA_seq_Rebeca_HC/BRAINCODE_circexp/BC.n106.samplelist",character())
genes_fpkm=genes_fpkm[,sample106]
dim(genes_fpkm)
# fpkm --> tpm
library_size=colSums(genes_fpkm); scaling_factor = 1e-6
genes_tpm=sweep(genes_fpkm, 2, library_size * scaling_factor, "/") 
# filter "expressed genes" with mean tpm > 0.01
genes_expressed=genes_tpm[rowMeans(genes_tpm)>0.01,]
dim(genes_expressed)
# limit to lincRNA and protein-coding
genes_expressed_symbol=subset(genes_annotation, EnsID %in% rownames(genes_expressed) & type %in% c('protein_coding'),select=symbol, drop =T)
length(genes_expressed_symbol)
input_genes=read.delim("Merge_circexp_norm_filtered_and_enriched.2genes.stat.xls", stringsAsFactors=F, row.names = 1, header=T, check.names =F)$symbol

source("~/neurogen/pipeline/RNAseq/bin/lib.R")
gt=ORA(inputGenes=input_genes, allGenes=genes_expressed_symbol, output="Merge_circexp_norm_filtered_and_enriched.2genes.stat")

library(dplyr)
library(ggplot2)
gt %>% arrange(gene_set, -pvalue) %>% group_by(gene_set) %>% top_n(10, -pvalue) %>% ungroup() %>%
  mutate(Term=factor(V1, unique(as.character(V1)))) %>%  
  ggplot(aes(x = Term, y = -log10(pvalue), fill=gene_set)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  theme_bw() + theme_classic() +
  xlab("GO terms") + ylab("-log10(Fisher's test P value)") + 
  ggtitle(paste("Top 10 enriched terms (p <0.01) in"), subtitle="Merge_circexp_norm_filtered_and_enriched.2genes.stat.ORA.xls")
ggsave("Merge_circexp_norm_filtered_and_enriched.2genes.stat.ORA.pdf", width=10, height=8)

## adjusted by length: using GOseq
## only for brain circRNAs
library(goseq) # source("https://bioconductor.org/biocLite.R"); biocLite("goseq")
genes_fpkm=read.table("~/neurogen/rnaseq_PD/results/merged/genes.fpkm.cufflinks.allSamples.uniq.xls", header=T,row.names = 1, check.names = F)
genes_fpkm=genes_fpkm[,-c(1:7)]; colnames(genes_fpkm)=gsub("FPKM.","",colnames(genes_fpkm))
genes_annotation=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed", header = F, col.names = c('chr','start','end','EnsID','score','strand','symbol','type'), check.names = F, stringsAsFactors=F)
samples_neuron=scan("~/neurogen/rnaseq_PD/results/merged/samplelist.HC_Neuron",character()) # 99
genes_fpkm=genes_fpkm[,samples_neuron]
dim(genes_fpkm)
# filter "expressed genes" with at leat 30% samples with rpkm >1 (as https://www.biorxiv.org/content/biorxiv/early/2018/12/19/500991.full.pdf)
genes_expressed=genes_fpkm[rowMeans(genes_fpkm>1)>=0.3,]  
dim(genes_expressed)
# circRNA host genes
Merge_circexplorer_BC_annotation_per_cell=read.table("../results/Merge_circexplorer_BC.annotation_per_cell.xls", sep="\t", header=T, stringsAsFactors = F); head(Merge_circexplorer_BC_annotation_per_cell)
# only circRNAs in SNDA and PY
Merge_circexplorer_BC_annotation_per_cell = filter(Merge_circexplorer_BC_annotation_per_cell, celltype3 != "NN", circType == "circRNA"); dim(Merge_circexplorer_BC_annotation_per_cell) # 6798   16

clip <- pipe("pbcopy", "w")                       
dump(unique(Merge_circexplorer_BC_annotation_per_cell$geneName), file=clip)                               
close(clip)

# limit to lincRNA and protein-coding
genes_expressed_symbol=subset(genes_annotation, EnsID %in% rownames(genes_expressed) & type %in% c('protein_coding', 'lincRNA'),select=symbol, drop =T)
length(genes_expressed_symbol)

source("~/neurogen/pipeline/RNAseq/bin/lib.R")
ORA(inputGenes=unique(Merge_circexplorer_BC_annotation_per_cell$geneName), allGenes=genes_expressed_symbol, output="Merge_circexp_norm_filtered_and_enriched.circRNA.expressedinNeuron")
topGOenrichment(unique(Merge_circexplorer_BC_annotation_per_cell$geneName), allGenes=genes_expressed_symbol, topN=10, pCutoff=0.01, type='all', output="Merge_circexp_norm_filtered_and_enriched.circRNA.expressedinNeuron")

## ===============================================
## expression of linear RNA vs. circRNA
## ===============================================
library(tidyverse)
library(scales)
setwd("~/projects/circRNA/data/") 
nS=read.table("Merge_circexplorer_BC.annotation.bed14.s3s5.gz", col.names=c("ID","type","sample","reads","count"), check.names = F, stringsAsFactors=F)
nU=read.table("Merge_circexplorer_BC.annotation.bed14.u3u5.gz", col.names=c("ID","type","sample","reads","count"), check.names = F, stringsAsFactors=F)
nC=read.table("Merge_circexplorer_BC.annotation.bed14.circReads.txt", header=T, check.names = F, row.names = 1, stringsAsFactors=F)
# nL0=read.table("Merge_circexplorer_BC.annotation.bed14.sum_u3u5s3s5", header=T, check.names = F, row.names = 1, stringsAsFactors=F) # from bedtools coverage (which might be buggy)
sample106=scan("~/neurogen/circRNA_seq_Rebeca_HC/BRAINCODE_circexp/BC.n106.samplelist",character())
filtered_enriched_annotation=readRDS("Merge_circexplorer_BC106.filtered.enriched.annotation.bed14.rds")$ID; length(filtered_enriched_annotation)

nLw=unite(rbind(nS,nU), temp, sample, type, sep = "__") %>% select(ID, temp, reads) %>% spread(key=temp, value=reads, fill=0) %>% column_to_rownames(var = "ID")
nS3=select(nLw, contains("__s3")) %>% rename_all(.funs = funs(sub("__.*", "", .))); dim(nS3)
nS5=select(nLw, contains("__s5")) %>% rename_all(.funs = funs(sub("__.*", "", .))); dim(nS5)
nU3=select(nLw, contains("__u3")) %>% rename_all(.funs = funs(sub("__.*", "", .))); dim(nU3)
nU5=select(nLw, contains("__u5")) %>% rename_all(.funs = funs(sub("__.*", "", .))); dim(nU5)
# nS3=filter(nS,type=='s3') %>% select(ID, sample, reads) %>% spread(key=sample, value=reads, fill = 0) %>% column_to_rownames(var = "ID")
# nS5=filter(nS,type=='s5') %>% select(ID, sample, reads) %>% spread(key=sample, value=reads, fill = 0) %>% column_to_rownames(var = "ID")
# nU3=filter(nU,type=='u3') %>% select(ID, sample, reads) %>% spread(key=sample, value=reads, fill = 0) %>% column_to_rownames(var = "ID")
# nU5=filter(nU,type=='u5') %>% select(ID, sample, reads) %>% spread(key=sample, value=reads, fill = 0) %>% column_to_rownames(var = "ID")
common_Rows=Reduce(intersect, list(filtered_enriched_annotation, nS$ID, nU$ID, rownames(nC))); length(common_Rows)
common_Cols=Reduce(intersect, list(colnames(nC), sample106, nS$sample, nU$sample)); length(common_Cols)
nC=nC[common_Rows, sample106]; nS3=nS3[common_Rows, sample106]; nS5=nS5[common_Rows, sample106]; nU3=nU3[common_Rows, sample106]; nU5=nU5[common_Rows, sample106]
save(nC,nS3,nS5,nU3,nU5, file = "Merge_circexplorer_BC.annotation.bed14.nS3nS5nU3nU5.RData")
load("Merge_circexplorer_BC.annotation.bed14.nS3nS5nU3nU5.RData")
nL=nS3+nS5+nU3+nU5  # not very fair to use the sum
## 3/21/2019: We could change to use max: nL=max(nS3,nS5,nU3,nU5)
# combine
nCL=inner_join(x=rownames_to_column(nC, var = 'ID') %>% gather(key="sample",value=nC, contains("_")),
               y=rownames_to_column(nL, var = 'ID') %>% gather(key="sample",value=nL, contains("_")),
               by=c("ID","sample")) %>% 
  filter(nC>0, nL>0) %>%
  mutate(celltype=gsub(".*_.*_(.*)_.*_.*","\\1",sample))

nCLrange=range(c(nCL$nL,nCL$nC))

# ratio boxplot
head(nCL); dim(nCL)
group_by(nCL, celltype) %>% summarise(rho=cor(nC,nL,method = 'spearman'))
# celltype    rho
# 1 FB       0.0654
# 2 MCPY     0.114 
# 3 PBMC     0.104 
# 4 SNDA     0.113 
# 5 TCPY     0.0895
with(nCL, cor(nC, nL, method = 'spearman')) # 0.1028416

mutate(nCL, rCL=nC/(nC+nL)) %>%
  ggplot(aes(factor(celltype, levels = c("SNDA","TCPY","MCPY","FB","PBMC")),rCL)) + 
  scale_y_log10(breaks = c(0,1e-5,1e-4,1e-3,0.01,.1,.5,1),labels = c(0,1e-5,1e-4,1e-3,0.01,.1,.5,1)) +
  geom_violin(aes(fill = factor(celltype, levels = c("SNDA","TCPY","MCPY","FB","PBMC"))), scale = "width", trim = FALSE) + 
  geom_boxplot(fill='white', width=0.1, outlier.shape = NA) + 
  scale_fill_manual(name="Cell types",
                     values=c("SNDA"="#F22A7B","TCPY" = "#3182bd","MCPY" = "#2659B2","FB" ="#BC9371","PBMC"="#D3CBCB")) +
  labs(x="Cell types", y="Circular-to-linear fraction (ratio = nC/total reads)") +
  theme_bw() + 
  theme(legend.position='none')
ggsave("../results/Merge_circexplorer_BC106.annotation.bed14.nC.vs.nL.boxplot.pdf", width=4.5, height = 4.5)


# nC vs. nL scatter plot
set.seed(1)
ggplot(nCL,aes(x = jitter(nL, amount=0.49), y = jitter(nC, amount=.49))) +
  geom_point(aes(colour = factor(celltype, levels = c("SNDA","TCPY","MCPY","FB","PBMC"))), shape=19, size=.8, alpha=0.8) +
  geom_abline(intercept = 0, slope = 1, linetype=2) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(name="Cell types",
                     values=c("SNDA"="#F22A7B","TCPY" = "#3182bd","MCPY" = "#2659B2","FB" ="#BC9371","PBMC"="#D3CBCB")) +
  coord_fixed() +
  expand_limits(x=nCLrange,y=nCLrange) +
  labs(x="Number of linear reads at junction sites", y="Number of circular reads at junction sites") +
  theme_bw()
ggsave("../results/Merge_circexplorer_BC106.annotation.bed14.nC.vs.nL.pdf", width=6, height = 5, useDingbats=T)
ggsave("../results/Merge_circexplorer_BC106.annotation.bed14.nC.vs.nL.png", width=6, height = 5)

pdf("../results/Merge_circexplorer_BC106.annotation.bed14.nC.vs.nL.sub.pdf", width=6, height = 5)
p=ggplot(nCL,aes(x = nL, y = nC)) +
  geom_point(position=position_jitter(width = 0.1, height = 0.1), color='#aaaaaa', shape=19, size=1, alpha=0.8) + 
  geom_abline(intercept = 0, slope = 1, linetype=2) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed() +
  expand_limits(x=nCLrange,y=nCLrange) + 
  labs(x="Number of linear reads at junction sites", y="Number of circular reads at junction sites") +
  theme_bw()
print(p)
p1=ggplot(filter(nCL,celltype=="SNDA"),aes(x = nL, y = nC)) +
  geom_point(position=position_jitter(width = 0.1, height = 0.1),aes(colour = factor(celltype, levels = c("SNDA"))), shape=19, size=1, alpha=0.8) + 
  geom_abline(intercept = 0, slope = 1, linetype=2) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(name="Cell types",
                     values=c("SNDA" = "#F22A7B")) +
  coord_fixed() +
  expand_limits(x=nCLrange,y=nCLrange) + 
  labs(x="Number of linear reads at junction sites", y="Number of circular reads at junction sites") +
  theme_bw()
print(p1)
p2=ggplot(filter(nCL,celltype %in% c("TCPY","MCPY")),aes(x = nL, y = nC)) +
  geom_point(position=position_jitter(width = 0.1, height = 0.1),aes(colour = factor(celltype, levels = c("TCPY","MCPY"))), shape=19, size=1, alpha=0.8) + 
  geom_abline(intercept = 0, slope = 1, linetype=2) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(name="Cell types",
                     values=c("TCPY" = "#3182bd","MCPY" = "#2659B2")) +
  coord_fixed() +
  expand_limits(x=nCLrange,y=nCLrange) + 
  labs(x="Number of linear reads at junction sites", y="Number of circular reads at junction sites") +
  theme_bw()
print(p2)
p3=ggplot(filter(nCL,celltype %in% c("FB","PBMC")),aes(x = nL, y = nC)) +
  geom_point(position=position_jitter(width = 0.1, height = 0.1),aes(colour = factor(celltype, levels = c("FB","PBMC"))), shape=19, size=1, alpha=0.8) + 
  geom_abline(intercept = 0, slope = 1, linetype=2) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(name="Cell types",
                     values=c("FB" ="#BC9371","PBMC"="#D3CBCB")) +
  coord_fixed() +
  expand_limits(x=nCLrange,y=nCLrange) +
  labs(x="Number of linear reads at junction sites", y="Number of circular reads at junction sites") +
  theme_bw()
print(p3)
dev.off()

## Q: which circRNAs have significantly more circular reads than linear reads?
# use Wilcoxon Signed-Rank Test
length(grep("_SNDA_",names(nS3)))
readsNum_filtered<- read.table("~/neurogen/rnaseq_PD/run_output/linescounts.filtered.txt",row.names=1, header=F, check.names = F) 
readsNum_million<-(t(readsNum_filtered)[1,]/10^6)
# convert to RPM
SNDAsamples = names(nS3)[grep("_SNDA_",names(nS3))]
normalized_factor = readsNum_million[SNDAsamples]
rpmL<-sweep(nL[,SNDAsamples],2,normalized_factor,"/"); dim(rpmL)
rpmC<-sweep(nC[,SNDAsamples],2,normalized_factor,"/"); dim(rpmC)
# filter: at least 10% samples with RPM >0 in both circular and linear reads
expressed_rows = rowMeans(rpmL>0)>.01 & rowMeans(rpmC>0)>.01; sum(expressed_rows)
do.call(rbind, apply(rpmC[expressed_rows,] - rpmL[expressed_rows, ], 1, wilcox.test, alternative = "g", exact = F, correct=F)) %>% data.frame() %>% rownames_to_column() %>% filter(p.value<0.05)


# ## into 2D density (e.g. https://www.r-graph-gallery.com/2d-density-plot-with-ggplot2/)
# ggplot(nCL, aes(x=nL, y=nC) ) +
#   stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
#   scale_fill_distiller(palette="Spectral", direction=-1) +
#   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),expand = c(0, 0),
#                 labels = trans_format("log10", math_format(10^.x))) +
#   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),expand = c(0, 0),
#                 labels = trans_format("log10", math_format(10^.x))) +
#   coord_fixed() +
#   theme(legend.position='none')
