## Main script to work on the output of circExplorer
library('dplyr')
library('reshape2')
setwd("~/projects/circRNA/data/") 

## See ~/projects/circRNA/src/README.Rmd for how to generate the input files: 
# Merge_circexplorer_BC.annotation.bed14
# Merge_circexplorer_BC.rawcount.txt

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

###########################################
############## normalize to RPM ###########
###########################################

Merge_circexp_norm<-sweep(Merge_circexp_raw,2,readsNum_million,"/")

# save
saveRDS(Merge_circexp_raw, file="Merge_circexplorer_BC106.rawcount.rds")
saveRDS(Merge_circexp_norm, file="Merge_circexplorer_BC106.normRPM.rds")
saveRDS(annotation[annotation$ID %in% rownames(Merge_circexp_norm),], file="Merge_circexplorer_BC106.annotation.bed14.rds")

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

###########################################
############ being riched     ###########
###########################################

# Definition of being enriched: at least 20 reads in RNase R treatment AND foldchange(RNase vs Mock) at least 2 in at least one sample. 
# Run main.RM.R first
Merge_circexp_raw_long_enriched.RM  = readRDS(file="Merge_circexp_raw_long_filterd.RM.rds")
length(unique(Merge_circexp_raw_long_enriched.RM$name))
# [1] 48016

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

###########################################
## Figure 1b: distribution of circRNAs supported by different number of reads
###########################################

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


# ################################################### 
# ## For grant: subset of 20 samples as pilot study
# ###################################################
# Merge_circexp_raw0=Merge_circexp_raw; readsNum_million0=readsNum_million;
# set.seed(42); sample20=sample(sample106, 20); sample20
# readsNum_million=readsNum_million[sample20]
# Merge_circexp_raw=Merge_circexp_raw[,sample20]
# Merge_circexp_raw=Merge_circexp_raw[rowSums(Merge_circexp_raw)>0,]; dim(Merge_circexp_raw)
# #[1] 59571    20
# Merge_circexp_norm<-sweep(Merge_circexp_raw,2,readsNum_million,"/")
# Merge_circexp_raw_filtered <- Merge_circexp_raw[rowSums(Merge_circexp_raw)>=2, ]
# Merge_circexp_norm_filtered <- Merge_circexp_norm[rowSums(Merge_circexp_raw)>=2, ]
# 
# total_circRNA_raw_reads = rowSums(Merge_circexp_raw)
# 
# pdf("~/projects/circRNA/results/total_circRNA_raw_reads.BC20.pdf", width = 6, height = 5)
# layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(5,1), heights=c(1,1))
# ht=hist(total_circRNA_raw_reads[total_circRNA_raw_reads<=200], breaks = 0:200, plot=F)
# par(mar=c(4,4,0,1))
# plot(log10(ht$counts+.1)+1, lwd=2,xlim=c(0,200), type='h', ylim = c(0.9,6.6),col=c("#000000",rep('#ff0000',199)), ylab='Number of circular RNAs', xlab="Number of backspliced reads",yaxt="n",xaxs="i", yaxs="i",bty="n")
# y1=floor(range(log10(ht$counts+.1)+1))
# pow <- seq(y1[1], y1[2])
# ticksat <- as.vector(sapply(pow-1, function(p) log10((2:9)*10^p)))
# axis(2, pow, labels = 10^(pow-1), las=1)
# axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
# # add the remained ones with the #2 filter
# filtered_circRNA_raw_reads = rowSums(Merge_circexp_raw[rownames(Merge_circexp_raw_filtered),])
# ht=hist(filtered_circRNA_raw_reads[filtered_circRNA_raw_reads<=200], breaks = 0:200, plot=F)
# points(log10(ht$counts+.1)+1, lwd=2, type='h',col='#ff0000')
# legend("topleft",c(paste0("Distinct circular RNAs (n = ",format(length(total_circRNA_raw_reads),big.mark=","),")"),
#                    paste0("+ at least 2 reads in overall samples (n = ", format(sum(total_circRNA_raw_reads>1),big.mark=","),")")),
#        col=c("#000000",'#ff0000'),
#        text.col = c("#000000",'#ff0000'),
#        bty='n', xpd=T)
# ## for the 200- region
# par(mar=c(4,1,0,1))
# total_circRNA_raw_reads = rowSums(Merge_circexp_raw)
# MAX=max(total_circRNA_raw_reads)
# ht=hist(total_circRNA_raw_reads[total_circRNA_raw_reads>200], breaks = 200:MAX, plot=F)
# plot(log10(ht$counts+.1)+1, lwd=2,xlim=range(total_circRNA_raw_reads), type='h', ylim = c(0.9,6.5),col='#ff0000', ylab="",xlab="",xaxt="n",yaxt='n',xaxs="i", yaxs="i",bty="n")
# y1=floor(range(log10(ht$counts+.1)+1))
# pow <- seq(y1[1], y1[2])
# ticksat <- as.vector(sapply(pow-1, function(p) log10((2:9)*10^p)))
# #axis(2, pow, labels = 10^(pow-1), las=1)
# #axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
# axis(1, c(0,seq(2000,MAX,5000)-200), labels=c(200,seq(2000,MAX,5000)))
# # add the remained ones with the #2 filter
# filtered_circRNA_raw_reads = rowSums(Merge_circexp_raw[rownames(Merge_circexp_norm_filtered),])
# ht=hist(filtered_circRNA_raw_reads[filtered_circRNA_raw_reads>200], breaks = 200:MAX, plot=F)
# points(log10(ht$counts+.1)+1, lwd=2,xlim=range(total_circRNA_raw_reads), type='h',col='#ff0000')
# dev.off()

## ===============================================
## Number of circRNAs per host gene
## ===============================================
require(scales)
mylog_trans <- 
  function (base = exp(1), from = 0) 
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

## ===============================================
## host genes function enrichment analysis
## ===============================================
## fold change, meanR, pvalue from T test (see main.RM.R)
Merge_circexp_raw_RM_ttest=readRDS(file="Merge_circexp_raw_RM_ttest.rds")
# get the subset in Merge_circexp_norm_filtered_and_enriched, then convert to gene based format
Merge_circexp_raw_RM_ttest %>% group_by(symbol) %>% top_n(n=1, wt=meanR) %>% top_n(n=-1, wt=pvalue) %>% top_n(n=1, wt=log2FoldChange) %>% ungroup() %>% dplyr::select(symbol, padj=pvalue, log2FoldChange, stat=meanR) %>% distinct() %>% write.table(file="Merge_circexp_raw_RM_ttest.2genes.stat.xls", quote =F, sep="\t", row.names = T)

## then run the pathway analysis
# /data/rnaseq/src/_pathway.R -i Merge_circexp_raw_RM_ttest.2genes.stat.xls --filter medium

readRDS("Merge_circexplorer_BC106.filtered.enriched.annotation.bed14.rds") %>% filter(circType=='circRNA') %>% group_by(geneName) %>% summarise(n=n()) %>% dplyr::select(symbol=geneName, n=n) %>% write.table(file="Merge_circexp_norm_filtered_and_enriched.2genes.stat.xls", quote =F, sep="\t", row.names = T)

## then run the pathway analysis
# /data/rnaseq/src/_pathway.R -i Merge_circexp_norm_filtered_and_enriched.2genes.stat.xls --filter no

## ===============================================
## get synaptic genes and test its enrichment
## ===============================================
## (1) 1461 genes in Synapse Proteomics Datasets in G2C (http://synsysdb.genes2cognition.org/db/GeneList/L00000069)
G2C_PSP=read.delim("SYNSYSDBdb_L00000069_BAYES-COLLINS-HUMAN-PSD-FULL.txt", header = T, stringsAsFactors = F)$Symbol %>% unique()

## (2) genes with GO term including key words "synapse" or "synaptic";
library(biomaRt) #source("https://bioconductor.org/biocLite.R"); biocLite("biomaRt")
library(dplyr)
library(RCurl)
ensembl_hs=useEnsembl(biomart="ensembl",GRCh=37,dataset='hsapiens_gene_ensembl')
#listAttributes(ensembl_hs)
listAttributes(ensembl_hs, page='feature_page')
# get GO for all genes (if having GO annotation): N = 20475
genesGO_all=getBM(attributes=c('ensembl_gene_id','hgnc_symbol', 'go_id','name_1006','namespace_1003'),
                  filters =c('biotype','transcript_biotype','with_go'), 
                  values =list("protein_coding","protein_coding",T),
                  mart = ensembl_hs) %>% 
  rename(go_name=name_1006, go_domain=namespace_1003) %>% distinct
write.table(genesGO_all, file='genesGO_all.tab',sep = "\t", row.names = F)
genes_with_synapse_GO = genesGO_all %>% subset(grepl("synaptic|synapse", go_name), select=hgnc_symbol, drop=T) %>% unique()

# join "synapse" column
annotation$synapse=ifelse(annotation$geneName %in% c(G2C_PSP, genes_with_synapse_GO), "yes", "no")

total_circRNA_raw_reads_filtered = rowSums(Merge_circexp_raw[rowSums(Merge_circexp_raw)>=10, ])
total_circRNA_raw_reads_filtered_annotation=cbind(annotation[match(names(total_circRNA_raw_reads_filtered), annotation$ID),c('ID','geneName','synapse')], total_circRNA_raw_reads_filtered)

# all genes expressed in the 20 samples (RPKM>0.1)
fpkm=read.table("~/neurogen/rnaseq_PD/results/merged/genes.fpkm.cuffnorm.allSamples.uniq.xls", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
expressed_genes = gsub("\\..*","",rownames(fpkm)[rowMeans(fpkm[,sample20])>1])
expressed_genes_name = unique(genesGO_all$hgnc_symbol[genesGO_all$ensembl_gene_id %in% expressed_genes])
expressed_genes_ID = unique(genesGO_all$ensembl_gene_id[genesGO_all$ensembl_gene_id %in% expressed_genes])
write(expressed_genes_ID, file="pilot20_expressed_hostgenes.txt")
# write gene name to file for GSEA analysis
circRNA_hostgene = unique(as.character(total_circRNA_raw_reads_filtered_annotation$geneName))
write(circRNA_hostgene[circRNA_hostgene %in% expressed_genes_name], file="pilot20_circRNA_hostgenes.txt")

mytable=table(total_circRNA_raw_reads_filtered_annotation$synapse)
lbls=paste(names(mytable), "\n", mytable, sep="")
pie(mytable, labels = lbls, main="Pie chart of circRNAs of host genes with synaptic function")

mytable=table(unique(total_circRNA_raw_reads_filtered_annotation[,2:3])$synapse)
mytable

allgenes = unique(genesGO_all$hgnc_symbol)
synapse=ifelse(allgenes %in% c(G2C_PSP, genes_with_synapse_GO), "yes", "no")
circRNA=ifelse(allgenes %in% circRNA_hostgene, "yes", "no")
fisher.test(table(synapse, circRNA))


hist(total_circRNA_raw_reads_filtered_annotation$total_circRNA_raw_reads_filtered[total_circRNA_raw_reads_filtered_annotation$total_circRNA_raw_reads_filtered<=200], col="Red")
hist(total_circRNA_raw_reads_filtered_annotation$total_circRNA_raw_reads_filtered[total_circRNA_raw_reads_filtered_annotation$synapse=='yes'], breaks=seq(1,15000,1000), col="Blue", add=TRUE)
