## Main script to work on the output of circExplorer
library('dplyr')
library('reshape2')
setwd("~/projects/circRNA/data/") 

## See ~/projects/circRNA/README.Rmd for how to generate the input files

###########################################
################# load data  ##############
###########################################

## === annotation ===

annotation<- read.table("Merge_circexplorer_MSBB.annotation.bed14",header=T, check.names = F) # include both ciRNA and circRNA
head(annotation)

## === covariates ===
ID_DNA_RNA <- read.table("~/neurogen/ROSMAP/MSBB/Samples_hasGenoExpr.ID.DNA.RNA.tab", header = T, stringsAsFactors = F)

RNAseq_cov<- read.csv("~/neurogen/ROSMAP/MSBB/Metadata/MSBB_RNAseq_covariates.csv",header=T, check.names = F) 
readsNum_million<-RNAseq_cov$TotalReads[grepl("accepted_hits", RNAseq_cov$fileName)]/10^6;
names(readsNum_million)=RNAseq_cov$sampleIdentifier[grepl("accepted_hits", RNAseq_cov$fileName)]
#only include the 220 samples with both DNA and RNA
readsNum_million=setNames(readsNum_million[ID_DNA_RNA$RNA_SampleName], ID_DNA_RNA$individualIdentifier); length(readsNum_million)

## === load raw reads count ===

#removed "#" on first line as well as ful name to just sample name
Merge_circexp_raw <- read.table("Merge_circexplorer_MSBB.rawcount.txt",sep="\t",check.names =F, header=TRUE, comment.char ='') 
colnames(Merge_circexp_raw)=gsub("#|_candidates.bed_circReads","",colnames(Merge_circexp_raw)); dim(Merge_circexp_raw)

Merge_circexp_raw = mutate(Merge_circexp_raw, concated_column = paste(chrom,start,end, sep = '_'))
row.names(Merge_circexp_raw) = Merge_circexp_raw$concated_column
# change to use individualIdentifier as names
Merge_circexp_raw=Merge_circexp_raw[,ID_DNA_RNA$RNA_SampleName] %>% setNames(ID_DNA_RNA$individualIdentifier)  
dim(Merge_circexp_raw)
#[1] 206377    220


###########################################
############## normalize to RPM ###########
###########################################

Merge_circexp_norm<-sweep(Merge_circexp_raw,2,readsNum_million,"/")

###########################################
############ filter cirRNAs     ###########
###########################################

# Definition of being expression: with >=2 unique back-splicing reads in all samples (Zheng et al., doi:10.1038/ncomms11215)
Merge_circexp_raw_filtered <- Merge_circexp_raw[rowSums(Merge_circexp_raw)>=2, ]
Merge_circexp_norm_filtered <- Merge_circexp_norm[rowSums(Merge_circexp_raw)>=2, ]
annotation_filtered <- annotation[annotation$ID %in% rownames(Merge_circexp_norm_filtered),] 

dim(Merge_circexp_raw_filtered); dim(Merge_circexp_norm_filtered); dim(annotation_filtered)
# [1] 65119   125

# save
saveRDS(Merge_circexp_raw_filtered, file="Merge_circexplorer_MSBB.filtered.rawcount.rds")
saveRDS(Merge_circexp_norm_filtered, file="Merge_circexplorer_MSBB.filtered.normRPM.rds")
saveRDS(annotation_filtered, file="Merge_circexplorer_MSBB.filtered.annotation.bed14.rds")


## Filter: expression of >CUTOFF RPM in at least 10 individuals and ≥6 reads in at least 10 individuals 
# (like GTEx: https://www.gtexportal.org/home/documentationPage#staticTextAnalysisMethods)

Merge_circexp_norm_filtered = Merge_circexp_norm[rowSums(Merge_circexp_norm >= RPM_threshold) >= 10 & rowSums(Merge_circexp_raw >= raw_threshold) >= 10, ]


pdf("~/projects/circRNA/results/MSBB.total_circRNA_raw_reads.pdf", width = 6, height = 5)
layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(5,1), heights=c(1,1))
total_circRNA_raw_reads = rowSums(Merge_circexp_raw)
ht=hist(total_circRNA_raw_reads[total_circRNA_raw_reads<=200], breaks = 0:200, plot=F)
par(mar=c(4,4,0,1))
plot(log10(ht$counts+.1)+1, lwd=2,xlim=c(0,200), type='h', ylim = c(0.9,6.5),col="#aaaaaa", ylab='Number of circular RNAs', xlab="Number of backspliced reads",yaxt="n",xaxs="i", yaxs="i",bty="n")
y1=floor(range(log10(ht$counts+.1)+1))
pow <- seq(y1[1], y1[2])
ticksat <- as.vector(sapply(pow-1, function(p) log10((2:9)*10^p)))
axis(2, pow, labels = 10^(pow-1), las=1)
axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
# add the remained ones with the #2 filter
filtered_circRNA_raw_reads = rowSums(Merge_circexp_raw[rownames(Merge_circexp_norm_filtered),])
ht=hist(filtered_circRNA_raw_reads[filtered_circRNA_raw_reads<=200], breaks = 0:200, plot=F)
points(log10(ht$counts+.1)+1, lwd=2, type='h',col='#ff0000')
legend("topleft",c(paste(format(length(total_circRNA_raw_reads),big.mark=","),"distinct circular RNAs"),
                   paste(format(nrow(Merge_circexp_norm_filtered),big.mark=","),"with >=", signif(RPM_threshold,2), "RPM and >=", raw_threshold, "reads in >=10 samples")),
       col=c('#aaaaaa','#ff0000'),
       text.col = c('#aaaaaa','#ff0000'),
       bty='n', xpd=TRUE)

## for the 200- region
par(mar=c(4,1,0,1))
total_circRNA_raw_reads = rowSums(Merge_circexp_raw)
MAX=max(total_circRNA_raw_reads)
ht=hist(total_circRNA_raw_reads[total_circRNA_raw_reads>200], breaks = 200:MAX, plot=F)
plot(log10(ht$counts+.1)+1, lwd=2,xlim=range(total_circRNA_raw_reads), type='h', ylim = c(0.9,6.5),col='#aaaaaa', ylab="",xlab="",xaxt="n",yaxt='n',xaxs="i", yaxs="i",bty="n")
pow <- seq(y1[1], y1[2])
ticksat <- as.vector(sapply(pow-1, function(p) log10((2:9)*10^p)))
#axis(2, pow, labels = 10^(pow-1), las=1)
#axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
axis(1, c(0,seq(2000,MAX,5000)-200), labels=c(200,seq(2000,MAX,5000)))
# add the remained ones with the #2 filter
filtered_circRNA_raw_reads = rowSums(Merge_circexp_raw[rownames(Merge_circexp_norm_filtered),])
ht=hist(filtered_circRNA_raw_reads[filtered_circRNA_raw_reads>200], breaks = 200:MAX, plot=F)
points(log10(ht$counts+.1)+1, lwd=2,xlim=range(total_circRNA_raw_reads), type='h',col='#ff0000')

dev.off()