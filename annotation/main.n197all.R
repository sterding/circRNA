## Main script to work on the output of circExplorer
library('dplyr')
library('reshape2')
library('tidyverse')
library('ggpubr') # install.packages("ggpubr")
setwd("~/projects/circRNA/data/") 

## See ~/projects/circRNA/src/main.sh for how to generate the input files: 
# Merge_circexplorer_BC221.annotation.bed14
# Merge_circexplorer_BC221.rawcount.long.txt

## see pilot.R in ~/Dropbox/grant/2019R21/pilot.R to make figure for pilot study (for grant purpose)

###########################################
################# load data  ##############
###########################################

## === annotation ===

annotation<- read.table("Merge_circexplorer_BC221.annotation.bed14",header=T, check.names = F) # include both ciRNA and circRNA
head(annotation); dim(annotation)

# add geneSymbol, geneType
GENCODEv19=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v19.annotation.genes.bed", col.names = c("chr","start","end","geneName","geneID","geneType","strand"), stringsAsFactors = F, header = F)
annotation$geneName=GENCODEv19$geneName[match(annotation$geneID, GENCODEv19$geneID)]
annotation$geneType=GENCODEv19$geneType[match(annotation$geneID, GENCODEv19$geneID)]

# add CDR1as
annotation$geneName[annotation$geneID=='CDR1as'] = 'CDR1as'
annotation$geneType[annotation$geneID=='CDR1as'] = 'antisense'

## === load raw reads count ===

Merge_circexp_raw <- read.table("Merge_circexplorer_BC221.rawcount.long.txt",sep="\t",check.names =F, col.names = c("ID","readsCount","sampleID")) %>% 
  group_by(ID, sampleID) %>% summarise(readsCount=sum(readsCount)) %>%
  spread(key = sampleID, value = readsCount, fill = 0) %>%
  column_to_rownames("ID")
head(Merge_circexp_raw); dim(Merge_circexp_raw)

###########################################
############## normalize to RPM ###########
###########################################
BRAINCODE2_Sequencing_Log.RNAseq_statistics="https://docs.google.com/spreadsheets/d/e/2PACX-1vQFQ4aQj0sD9oxIqaZ-cEgo7kWcCmNYGBH9emLw8iNu0f6TTjKE5Lte7IBfoMMy57cLjA4pXE0YlPY2/pub?gid=1049148747&single=true&output=tsv"
require(RCurl)
RNAseq_statistics<- read.delim(textConnection(getURL(BRAINCODE2_Sequencing_Log.RNAseq_statistics)), header=T, check.names = T, comment.char = "#", stringsAsFactors = F) 
readsNum_million <- as.numeric(gsub(",","",RNAseq_statistics$total_mapped_reads))/10^6; names(readsNum_million) = RNAseq_statistics$SOURCE_SAMPLE_ID
readsNum_million = readsNum_million[colnames(Merge_circexp_raw)]

Merge_circexp_norm<-sweep(Merge_circexp_raw,2,readsNum_million,"/")

# save
saveRDS(Merge_circexp_raw, file="Merge_circexplorer_BC221.rawcount.rds")
saveRDS(Merge_circexp_norm, file="Merge_circexplorer_BC221.normRPM.rds")
saveRDS(annotation, file="Merge_circexplorer_BC221.annotation.bed14.rds")

write.table(Merge_circexp_raw, file="Merge_circexplorer_BC221.rawcount.xls.gz", quote=F, sep="\t", col.names = NA, row.names = TRUE)
write.table(Merge_circexp_norm, file="Merge_circexplorer_BC221.normRPM.xls", quote=F, sep="\t", col.names = NA, row.names = TRUE)

## SN vs. SNDA for UWA616 
# Rscript ~/pipeline/modules/_pairwise_compare.R HC_UWA616_SNDA_2_rep1 HC_UWA616_SN_6_rep1.amplified Merge_circexplorer_BC221.normRPM.rds

###########################################
############ filter cirRNAs     ###########
###########################################

### (0) only selected 197 samples 
BC197=filter(RNAseq_statistics, BRAINcode2.final.selection==1, CELL!='CSF') %>% pull(SOURCE_SAMPLE_ID)
Merge_circexp_raw_filtered <- Merge_circexp_raw[,BC197]
Merge_circexp_raw_filtered <- Merge_circexp_raw_filtered[rowSums(Merge_circexp_raw_filtered)>0,]; dim(Merge_circexp_raw_filtered)
saveRDS(Merge_circexp_raw_filtered, file="Merge_circexplorer_BC197.rawcount.rds")
saveRDS(annotation[annotation$ID %in% rownames(Merge_circexp_raw_filtered),], file="Merge_circexplorer_BC197.annotation.bed14.rds")

### (1) being expressed: with >=2 unique back-splicing reads in all samples (Zheng et al., doi:10.1038/ncomms11215)
### --------------------

Merge_circexp_raw_filtered <- Merge_circexp_raw_filtered[rowSums(Merge_circexp_raw_filtered)>=2,]
Merge_circexp_norm_filtered <- Merge_circexp_norm[rownames(Merge_circexp_raw_filtered), BC197]
annotation_filtered <- annotation[annotation$ID %in% rownames(Merge_circexp_norm_filtered),] 

dim(Merge_circexp_raw_filtered); dim(Merge_circexp_norm_filtered); dim(annotation_filtered)
# [1] 124771   197

# save
saveRDS(Merge_circexp_raw_filtered, file="Merge_circexplorer_BC197.filtered.rawcount.rds")
saveRDS(Merge_circexp_norm_filtered, file="Merge_circexplorer_BC197.filtered.normRPM.rds")
saveRDS(annotation_filtered, file="Merge_circexplorer_BC197.filtered.annotation.bed14.rds")
#write.table(annotation_filtered, file="Merge_circexplorer_BC197.filtered.annotation.bed14", sep = "\t", quote = F, row.names = F, col.names = F)

###########################################
############ being enriched     ###########
###########################################

# Definition of being enriched: at least 20 reads in RNase R treatment AND foldchange(RNase vs Mock) at least 2 in at least one sample. 
# Run main.RM.R first
Merge_circexp_raw_long_enriched.RM  = readRDS(file="Merge_circexp_raw_long_filterd.RM.rds")
length(unique(Merge_circexp_raw_long_enriched.RM$ID))
# [1] 50900
table(unique(Merge_circexp_raw_long_enriched.RM$ID) %in% rownames(Merge_circexp_raw_filtered))
# FALSE  TRUE 
# 37499  13401

Merge_circexp_norm_filtered_and_enriched=Merge_circexp_norm_filtered[intersect(rownames(Merge_circexp_raw_filtered), Merge_circexp_raw_long_enriched.RM$ID), ]
Merge_circexp_raw_filtered_and_enriched <- Merge_circexp_raw_filtered[rownames(Merge_circexp_norm_filtered_and_enriched), ]
annotation_filtered_enriched <- annotation[annotation$ID %in% rownames(Merge_circexp_raw_filtered_and_enriched),] 

dim(Merge_circexp_norm_filtered_and_enriched); dim(Merge_circexp_raw_filtered_and_enriched); dim(annotation_filtered_enriched)
# [1] 13401   197
# [1] 10431   125
# [1] 10017   106

# save
saveRDS(Merge_circexp_raw_filtered_and_enriched, file="Merge_circexplorer_BC197.filtered.enriched.rawcount.rds")
saveRDS(Merge_circexp_norm_filtered_and_enriched, file="Merge_circexplorer_BC197.filtered.enriched.normRPM.rds")
saveRDS(annotation_filtered_enriched, file="Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds")
annotation_filtered_enriched=readRDS(file="Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds")
write.table(annotation_filtered_enriched, file="Merge_circexplorer_BC197.filtered.enriched.annotation.bed14", sep = "\t", quote = F, row.names = F, col.names = F)

#########################################
## Figure 1b: distribution of circRNAs supported by different number of reads
###########################################

# pie chart of circRNA among all circular RNAs
pdf("~/projects/circRNA/results/total_circRNA_raw_reads.BC197.pie.pdf", width=6, height = 2)
par(mfrow=c(1,3))
pie(table(readRDS("Merge_circexplorer_BC197.annotation.bed14.rds")$circType=='circRNA'),labels = NA, border='white', col=c('#fdbb84','red'),main="all")
pie(table(readRDS("Merge_circexplorer_BC197.filtered.annotation.bed14.rds")$circType=='circRNA'),labels = NA, border='white', col=c('#fdbb84','red'),main="filtered")
pie(table(readRDS("Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds")$circType=='circRNA'),labels = NA, border='white', col=c('#fdbb84','red'),main="filtered+enriched")
dev.off()

table(readRDS("Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds")[,c('geneType','circType')])
# geneType                   circRNA ciRNA
# 3prime_overlapping_ncrna       2     0
# antisense                     36     2
# lincRNA                       38     1
# polymorphic_pseudogene         1     0
# processed_transcript          20     3
# protein_coding             12833   405
# pseudogene                    58     0
# sense_intronic                 1     0
# sense_overlapping              1     0

Merge_circexp_raw=readRDS("Merge_circexplorer_BC197.rawcount.rds")
Merge_circexp_raw_filtered_and_enriched=readRDS("Merge_circexplorer_BC197.filtered.enriched.rawcount.rds")
total_circRNA_raw_reads = rowSums(Merge_circexp_raw)
pdf("~/projects/circRNA/results/total_circRNA_raw_reads.BC197.pdf", width = 6, height = 5)
BREAK=300
layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(5,1), heights=c(1,1))
ht=hist(total_circRNA_raw_reads[total_circRNA_raw_reads<=BREAK], breaks = 0:BREAK, plot=F)
par(mar=c(4,4,0,1))
plot(log10(ht$counts+.1)+1, lwd=2,xlim=c(0,BREAK), type='h', ylim = c(0.9,6.6),col=c("#000000",rep('#aaaaaa',199)), ylab='Number of circular RNAs', xlab="Number of backspliced reads",yaxt="n",xaxs="i", yaxs="i",bty="n")
y1=floor(range(log10(ht$counts+.1)+1))
pow <- seq(y1[1], y1[2])
ticksat <- as.vector(sapply(pow-1, function(p) log10((2:9)*10^p)))
axis(2, pow, labels = 10^(pow-1), las=1)
axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
# add the remained ones with the #2 filter: enriched
filtered_circRNA_raw_reads = rowSums(Merge_circexp_raw[rownames(Merge_circexp_raw_filtered_and_enriched),])
ht=hist(filtered_circRNA_raw_reads[filtered_circRNA_raw_reads<=BREAK], breaks = 0:BREAK, plot=F)
points(log10(ht$counts+.1)+1, lwd=2, type='h',col='#8B0000')

# add the remained ones with the #3 filter: circRNA only
circRNAonly=as.character(annotation_filtered_enriched$ID[annotation_filtered_enriched$circType=='circRNA'])
filtered_circRNA_raw_reads_circRNAonly = rowSums(Merge_circexp_raw[circRNAonly,])
ht=hist(filtered_circRNA_raw_reads_circRNAonly[filtered_circRNA_raw_reads_circRNAonly<=BREAK], breaks = 0:BREAK, plot=F)
points(log10(ht$counts+.1)+1, lwd=2, type='h',col='#ff0000')

legend("topleft",c(paste0("Distinct circular RNAs (n = ",format(length(total_circRNA_raw_reads),big.mark=","),")"),
                   paste0("- at least 2 reads in overall samples (n = ", format(sum(total_circRNA_raw_reads>1),big.mark=","),")"),
                   paste0("-- enriched in RNase R (n = ", format(length(filtered_circRNA_raw_reads),big.mark=","),")"),
                   paste0("--- circRNA only (n = ", format(length(filtered_circRNA_raw_reads_circRNAonly),big.mark=","),")")),
       col=c("#000000",'#aaaaaa','#8B0000','#ff0000'),
       text.col = c("#000000",'#aaaaaa','#8B0000','#ff0000'),
       bty='n', xpd=T)

## for the BREAK- region
par(mar=c(4,1,0,1))
n200=total_circRNA_raw_reads[total_circRNA_raw_reads>BREAK]; 
MAX=max(total_circRNA_raw_reads)
ht=hist(n200, breaks = BREAK:MAX, plot=F)
plot(log10(ht$counts+.1)+1, lwd=2,xlim=range(total_circRNA_raw_reads), type='h', ylim = c(0.9,6.6),col='#aaaaaa', ylab="",xlab="",xaxt="n",yaxt='n',xaxs="i", yaxs="i",bty="n")
y1=floor(range(log10(ht$counts+.1)+1))
pow <- seq(y1[1], y1[2])
ticksat <- as.vector(sapply(pow-1, function(p) log10((2:9)*10^p)))
#axis(2, pow, labels = 10^(pow-1), las=1)
#axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
axis(1, c(0,seq(10000,MAX,20000)-BREAK), labels=c(BREAK,seq(10000,MAX,20000)))

# add the remained ones with the #2 filter
n200=filtered_circRNA_raw_reads[filtered_circRNA_raw_reads>BREAK];
ht=hist(n200, breaks = BREAK:MAX, plot=F)
points(log10(ht$counts+.1)+1, lwd=2, type='h',col='#8B0000')

# add the remained ones with the #3 filter: circRNA only
ht=hist(filtered_circRNA_raw_reads_circRNAonly[filtered_circRNA_raw_reads_circRNAonly>BREAK], breaks = BREAK:MAX, plot=F)
points(log10(ht$counts+.1)+1, lwd=2, type='h',col='#ff0000')

dev.off()

## add gene highlight in the figure
annotation_filtered_enriched %>% select(ID, circType, geneName) %>% 
  left_join(rownames_to_column(Merge_circexp_raw_filtered_and_enriched), by = c("ID"="rowname")) %>% 
  mutate(sumVar = rowSums(select(., contains("_")))) %>% select(-contains("_")) %>% 
  filter(circType=='circRNA') %>% arrange(-sumVar) %>% 
  filter(sumVar>10) %>%
  filter(geneName %in% c("SLC8A1", 'SNAP25', 'SV2A', 'RAB3A', 'SYNGR1','NRXN3', 'SYNPR','STXBP1', 'ERC2','CASK','PPFIA1','PPFIA4','UNC13B', "DNAJC6","RIMS1","RIMS2"))
filter(geneName %in% c("SORL1","PSEN1","APP", "ABCA7","CLU","CR1","PICALM","PLD3","TREM2",'UCHL1','SRY',"SNCA","LRRK2","PINK1","ABCA7","CLU","CR1","PLD3","TREM2",'PARK7', 'PRKN', 'GBA'))

# ID circType geneName sumVar
# 1     chr2_40655612_40657444  circRNA   SLC8A1   9885
# 2     chr2_40655612_40657441  circRNA   SLC8A1   1385
# 3     chr6_73005639_73043538  circRNA    RIMS1   1011
# 4     chr6_73016960_73043538  circRNA    RIMS1    977
# 5    chr11_85707868_85714494  circRNA   PICALM    951
# 6     chr1_65830317_65831879  circRNA   DNAJC6    710
# 7    chr11_85707868_85742653  circRNA   PICALM    637
# 8     chr2_40366540_40405633  circRNA   SLC8A1    226
# 9     chr6_72960032_72961071  circRNA    RIMS1    201
# 10  chr8_105080739_105161076  circRNA    RIMS2    196
# 11    chr6_72960032_73043538  circRNA    RIMS1    159
# 12   chr11_85701292_85742653  circRNA   PICALM    152
# 13   chr11_85685750_85742653  circRNA   PICALM    142
# 14   chr14_73614502_73640415  circRNA    PSEN1    125
# 15  chr8_105105698_105161076  circRNA    RIMS2    123
# 16   chr11_85685750_85695016  circRNA   PICALM    109
# 17   chr14_73614502_73614814  circRNA    PSEN1    109
# 18    chr6_73001636_73043538  circRNA    RIMS1    102
# 19  chr8_105053553_105161076  circRNA    RIMS2     93
# 20    chr1_65830317_65860715  circRNA   DNAJC6     84
# 21   chr21_27326903_27354790  circRNA      APP     70
# 22   chr11_85707868_85712201  circRNA   PICALM     55
# 23   chr14_73614502_73614802  circRNA    PSEN1     55
# 24   chr11_85718584_85742653  circRNA   PICALM     52
# 25    chr8_27462440_27464041  circRNA      CLU     33
# 26   chr14_73614502_73659572  circRNA    PSEN1     29
# 27   chr11_85692171_85742653  circRNA   PICALM     27
# 28   chr11_85722072_85742653  circRNA   PICALM     27
# 29  chr8_105025669_105161076  circRNA    RIMS2     25
# 30   chr14_73614502_73664837  circRNA    PSEN1     23
# 31   chr11_85707868_85718626  circRNA   PICALM     22
# 32   chr21_27347382_27372497  circRNA      APP     22
# 33    chr6_73001636_73023375  circRNA    RIMS1     21
# 34    chr6_72952016_73043538  circRNA    RIMS1     20
# 35   chr11_85685750_85714494  circRNA   PICALM     18
# 36    chr2_40387904_40405633  circRNA   SLC8A1     17
# 37  chr1_155207131_155209868  circRNA      GBA     16
# 38   chr11_85692787_85722179  circRNA   PICALM     14
# 39  chr8_104943490_105161076  circRNA    RIMS2     14
# 40    chr8_27462440_27462852  circRNA      CLU     14
# 41   chr21_27326903_27372497  circRNA      APP     13
# 42 chr11_121348826_121367758  circRNA    SORL1     12
# 43    chr1_65830317_65845207  circRNA   DNAJC6     12
# 44    chr1_65830317_65871816  circRNA   DNAJC6     12
# 45  chr8_104922361_105026843  circRNA    RIMS2     12
# 46  chr8_105025669_105080840  circRNA    RIMS2     12
# 47   chr11_85707868_85726006  circRNA   PICALM     11
# 48   chr12_40696590_40697936  circRNA    LRRK2     11

## ===============================================
## back-splicing reads vs. sample fraction
## similar to Fig. 3E in https://www-sciencedirect-com.ezp-prod1.hul.harvard.edu/science/article/pii/S0092867418316350?via%3Dihub#fig3
## ===============================================
Merge_circexp_raw_filtered_and_enriched=readRDS("Merge_circexplorer_BC197.filtered.enriched.rawcount.rds")
dim(Merge_circexp_raw_filtered_and_enriched)
df=Merge_circexp_raw_filtered_and_enriched
pdf("../results/Merge_circexplorer_BC197.filtered.enriched.nReads_vs_nSample.pdf", width = 3, height = 3)
plot(jitter(rowMeans(df>0)), jitter(rowMeans(df)), pch=20, col='#00000033')
dev.off()
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
annotation_filtered_enriched=readRDS("Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds")
annotation_filtered_enriched_n = annotation_filtered_enriched %>% filter(circType=='circRNA') %>% group_by(geneName) %>% summarise(n=n()) %>% group_by(n) %>% summarise(N=n(), geneNames=paste(geneName, collapse="; "))
ggplot(annotation_filtered_enriched_n, aes(x=n, y=N)) + 
  geom_col() + 
  geom_text(data=subset(annotation_filtered_enriched_n, n > 17), aes(x=n,y=N,label=geneNames), size=2,angle=90, nudge_y=0.05, hjust=0) +
  xlab("Number of circRNAs in the host gene") + ylab("Count of host genes") +
  ggtitle("Histogram of number of circRNAs per host gene") +
  scale_y_continuous(trans = mylog_trans(base=10, from=-1), breaks=c(1,10,100,1000),limits=c(0.1,1500)) +
  scale_x_continuous(breaks=pretty_breaks())
ggsave("../results/Merge_circexplorer_BC197.filtered.enriched.circRNAs_per_hostgene.hist.pdf",width = 5, height = 3)

## vs. # of exon per host genes
longestTx_annotation=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v19.annotation.longestTx.bed12", header = F, col.names = c('chr','start','end','ID','score','strand','gstart','gend','rgb','nExon','lengths','starts'), check.names = F, stringsAsFactors=F)
longestTx_annotation=separate(longestTx_annotation, ID,c('geneName','geneID',NA,NA),sep='___')
annotation_filtered_enriched$nExon=longestTx_annotation$nExon[match(annotation_filtered_enriched$geneID, longestTx_annotation$geneID)]
annotation_filtered_enriched$nExon[annotation_filtered_enriched$geneID=='CDR1as']=1
head(annotation_filtered_enriched)
# remove NA case caused by inconsistance between RefSeq and GENCODE gene symbol
annotation_filtered_enriched %>% filter(!is.na(nExon), circType=='circRNA') %>% 
  group_by(geneID) %>% summarise(nExon=mean(nExon), nCircRNA=n()) %>% ungroup() %>% 
  group_by(nExon_interval=cut(nExon, breaks = c(0,10,20,30,1000))) %>% 
  ggplot(aes(nExon_interval,nCircRNA)) + geom_boxplot() + theme_classic() 
ggsave("../results/Merge_circexplorer_BC197.filtered.enriched.nCircRNA_vs_nExon.pdf", width = 3, height = 3)

## ===============================================
## exon and flanking intron length of circRNAs vs. all exons
## ===============================================
annotation_filtered_enriched=readRDS(file="Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds") %>% filter(circType=='circRNA') %>% select(1:12)
introns=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/introns.meta.bed", header = F, col.names = c('chr','start','end','ID','score','strand'), check.names = F, stringsAsFactors=F) %>% mutate(score=end-start) %>% unite("chrend",c("chr",'end'), remove =F) %>% unite("chrstart",c("chr",'start'), remove =F)
head(introns)
controls=read.table("Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.matched2", header = F, col.names = c('chrom','start','end','ID','score','strand','thickStart','thickEnd','itemRgb','exonCount','exonSizes','exonOffsets','exonindex','hostgene','matchedCircRNA'), stringsAsFactors = F) %>% filter(matchedCircRNA %in% annotation_filtered_enriched$ID) %>% select(1:12) %>% mutate(itemRgb='0,0,0')
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
ggsave("../results/Merge_circexplorer_BC197.filtered.enriched.length.flankingIntron.pdf", width = 2,height = 3)

## number of repetitive elements comparsion between called circRNAs and controls
# bash
# slopBed -g ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/hg19.chrom.size -i Merge_circexplorer_BC197.filtered.enriched.annotation.bed14 -b 1 | bedtools intersect -a ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/introns.meta.bed -b - -wo | awk '$23==1' | cut -f1-6,10 | bedtools intersect -a - -b ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Variation/rmsk.hg19.bed -wo | cut -f7 | sort | uniq -c | sed -e 's/^[ \t]*//' > Merge_circexplorer_BC197.filtered.enriched.nRepeats_in_flankingIntron.txt
# awk '{OFS="\t"; $4=$15; print $0,$15}' Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.matched2 | slopBed -g ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/hg19.chrom.size -b 1 | bedtools intersect -a ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/introns.meta.bed -b - -wo | awk '$23==1' | cut -f1-6,10 | bedtools intersect -a - -b ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Variation/rmsk.hg19.bed -wo | cut -f7 | sort | uniq -c  | sed -e 's/^[ \t]*//' > Merge_circexplorer_BC197.filtered.enriched.matched2.nRepeats_in_flankingIntron.txt
df=bind_rows(mutate(read.table('Merge_circexplorer_BC197.filtered.enriched.nRepeats_in_flankingIntron.txt', header = F, col.names = c('nRepeats','circRNAID')),type="circRNA"),
             mutate(read.table('Merge_circexplorer_BC197.filtered.enriched.matched2.nRepeats_in_flankingIntron.txt', header = F, col.names = c('nRepeats','circRNAID')),type="control"))
ggplot(df,aes(x=type, y=nRepeats, col=type)) + geom_boxplot() + scale_y_log10() + theme_classic() + 
  ggtitle(paste("Wilcox test P value =",signif(wilcox.test(nRepeats ~ as.factor(type), df)$p.value,2))) +
  scale_color_manual(guide=FALSE, name="type",values=c("circRNA"="red","control" = "grey"))
ggsave("../results/Merge_circexplorer_BC197.filtered.enriched.nRepeats_in_flankingIntron.pdf", width = 2,height = 3)

# single-exon circRNAs lenth vs. random exon controls
exons=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/exons.meta.bed", sep="\t", col.names = c("chrom","start","end","ID","score","strand"), stringsAsFactors = F, header = F) %>% filter(grepl("protein_coding", ID)) %>% select(1:3) %>% sample_n(5000) %>% mutate(type="background")
df=readRDS(file="Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds") %>% 
  filter(exonCount==1, circType=='circRNA') %>% select(1:3) %>% mutate(type="circRNA") %>% bind_rows(exons) %>%
  mutate(exon_length=end-start)
ggplot(df, aes(x=type, y=exon_length, col=type)) + geom_boxplot() + scale_y_log10() + theme_classic() + 
  ggtitle(paste("Wilcox test P value =",signif(wilcox.test(exon_length ~ as.factor(type), df)$p.value,2))) +
  scale_color_manual(guide=FALSE, name="type",values=c("circRNA"="red","background" = "grey")) 
ggsave("../results/Merge_circexplorer_BC197.filtered.enriched.length.single-exon-circRNAs.pdf", width = 2,height = 3)

# region length distribution
readRDS(file="Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds") %>% 
  ggplot(aes(x=end-start, fill=factor(circType, levels=c('ciRNA', 'circRNA')))) + 
  geom_histogram(position = 'stack',bins = 60)  + scale_color_brewer(palette="Dark2") + 
  theme_minimal()+theme_classic()+theme(legend.position="top")+ 
  scale_x_log10() + scale_fill_manual(name="circular RNA types",values=c("circRNA"="red","ciRNA" = "orange")) 
ggsave("../results/Merge_circexplorer_BC197.filtered.enriched.lengthDistribution.pdf", width = 4,height = 3)


# lenth distribution for called circRNAs
readRDS(file="Merge_circexplorer_BC197.annotation.bed14.rds") %>% 
  ggplot(aes(x=end-start, fill=factor(circType, levels=c('ciRNA', 'circRNA')))) + 
  geom_histogram(position = 'stack',bins = 60)  + scale_color_brewer(palette="Dark2") + 
  theme_minimal()+theme_classic()+theme(legend.position="top")+ 
  scale_x_log10() + scale_fill_manual(name="circular RNA types",values=c("circRNA"="red","ciRNA" = "orange"))
ggsave("../results/Merge_circexplorer_BC197.unfiltered.lengthDistribution.pdf", width = 4,height = 3)


## ===============================================
## host genes function enrichment analysis
## ===============================================

# ### tryout 1:  --- ABANDONED
# ## For each host gene, pick the circRNA with strongest meanR (then smallest pvalue, bigger fold change) from the R-vs-M T-test (see main.RM.R)
# ## ------------------------
# readRDS(file="Merge_circexp_raw_RM_ttest.rds") %>% group_by(symbol) %>% 
#   top_n(n=1, wt=meanR) %>% top_n(n=-1, wt=pvalue) %>% top_n(n=1, wt=log2FoldChange) %>% 
#   ungroup() %>% dplyr::select(symbol, padj=pvalue, log2FoldChange, stat=meanR) %>% distinct() %>% 
#   write.table(file="Merge_circexp_raw_RM_ttest.2genes.stat.xls", quote =F, sep="\t", row.names = T)
# ## then run the pathway analysis
# # /data/rnaseq/src/_pathway.R -i Merge_circexp_raw_RM_ttest.2genes.stat.xls --filter medium
# # see the output Merge_circexp_raw_RM_ttest.2genes.stat.xls.GOenrichment.pdf

# ### tryout 2:  --- ABANDONED
# ## Run topGO analysis for the circRNA-hosting genes
# ## ------------------------
# readRDS("Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds") %>% filter(circType=='circRNA') %>% 
#   group_by(geneName) %>% summarise(n=n()) %>% dplyr::select(symbol=geneName, stat=n) %>% 
#   write.table(file="Merge_circexplorer_BC197.filtered.enriched.annotation.2genes.stat.xls", quote =F, sep="\t", row.names = T)
# # bash
# # ~/neurogen/pipeline/RNAseq/modules/_pathway.R -i Merge_circexplorer_BC197.filtered.enriched.annotation.2genes.stat.xls --filter no
# 
# ### tryout 3:  --- ABANDONED
# ## ORA of circRNA-hosting genes vs. all expressed genes (similar to You et al. Nature Neuroscience, 2015)
# ## ------------------------
# genes_fpkm=read.table("~/neurogen/rnaseq_PD/results/merged/genes.fpkm.cufflinks.allSamples.BCv2.uniq.xls", header=T,row.names = 1, check.names = F)
# genes_annotation=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed", header = F, col.names = c('chr','start','end','EnsID','score','strand','symbol','type'), check.names = F, stringsAsFactors=F)
# genes_fpkm=genes_fpkm[,BC197]
# dim(genes_fpkm)
# # fpkm --> tpm
# library_size=colSums(genes_fpkm); scaling_factor = 1e-6
# genes_tpm=sweep(genes_fpkm, 2, library_size * scaling_factor, "/") 
# # filter "expressed genes" with mean tpm > 0.01
# genes_expressed=genes_tpm[rowMeans(genes_tpm)>0.01,]
# dim(genes_expressed)
# # limit to lincRNA and protein-coding
# genes_expressed_symbol=subset(genes_annotation, EnsID %in% rownames(genes_expressed) & type %in% c('protein_coding'),select=symbol, drop =T)
# length(genes_expressed_symbol)
# input_genes=read.delim("Merge_circexplorer_BC197.filtered.enriched.annotation.2genes.stat.xls", stringsAsFactors=F, row.names = 1, header=T, check.names =F)$symbol
# 
# source("~/neurogen/pipeline/RNAseq/bin/lib.R")
# gt=ORA(inputGenes=input_genes, allGenes=genes_expressed_symbol, output="Merge_circexplorer_BC197.filtered.enriched.annotation.2genes.stat")
# 
# library(dplyr)
# library(ggplot2)
# gt %>% arrange(gene_set, -pvalue) %>% group_by(gene_set) %>% top_n(10, -pvalue) %>% ungroup() %>%
#   mutate(Term=factor(V1, unique(as.character(V1)))) %>%  
#   ggplot(aes(x = Term, y = -log10(pvalue), fill=gene_set)) + 
#   geom_bar(stat = "identity") + 
#   coord_flip() + 
#   theme_bw() + theme_classic() +
#   xlab("GO terms") + ylab("-log10(Fisher's test P value)") + 
#   ggtitle(paste("Top 10 enriched terms (p <0.01) in"), subtitle="Merge_circexplorer_BC197.filtered.enriched.annotation.2genes.stat.ORA.xls")
# ggsave("Merge_circexp_norm_filtered_and_enriched.2genes.stat.ORA.pdf", width=10, height=8)

### tryout 4:  --- USED
## GOseq of circRNA-hosting genes vs. all expressed genes (only in brain)
## ------------------------
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
Merge_circexplorer_BC_annotation_per_cell=read.table("../results/Merge_circexplorer_BC197.annotation_per_cell.xls", sep="\t", header=T, stringsAsFactors = F); head(Merge_circexplorer_BC_annotation_per_cell)
# only circRNAs in SNDA and PY
Merge_circexplorer_BC_annotation_per_cell = filter(Merge_circexplorer_BC_annotation_per_cell, celltype3 != "NN", circType == "circRNA"); dim(Merge_circexplorer_BC_annotation_per_cell) # 6798   16

# clip <- pipe("pbcopy", "w")                       
# dump(unique(Merge_circexplorer_BC_annotation_per_cell$geneName), file=clip)                               
# close(clip)

# limit to lincRNA and protein-coding
genes_expressed_symbol=subset(genes_annotation, EnsID %in% rownames(genes_expressed) & type %in% c('protein_coding', 'lincRNA'),select=symbol, drop =T)
length(genes_expressed_symbol)

source("~/neurogen/pipeline/RNAseq/bin/lib.R")
#ORA(inputGenes=unique(Merge_circexplorer_BC_annotation_per_cell$geneName), allGenes=genes_expressed_symbol, output="Merge_circexp_norm_filtered_and_enriched.circRNA.expressedinNeuron")
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