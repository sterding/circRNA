## R script to test if circRNA can protect mRNA from degradation
## TODO: Limit to only single-exon circRNAs

library(RCurl)
library(tidyverse)
setwd("~/projects/circRNA/data/") 

covarianceTableURL="https://docs.google.com/spreadsheets/d/1I8nRImE9eJCCuZwpjfrrj-Uwx9bLebnO6o-ph7u6n8s/pub?gid=195725118&single=true&output=tsv"  # for all 140 samples
covs=read.delim(textConnection(getURL(covarianceTableURL)), stringsAsFactors =F)
covs=subset(covs, BRAINCODE.final.selection==1)

## RIN distribution
covs %>% filter(BRAINCODE.final.selection==1) %>% ggplot(aes(x=RIN)) + geom_histogram(na.rm = T, breaks=seq(5.5, 10, by=.5), fill='lightblue',col='white') + geom_vline(xintercept=c(7,8.5), col='red', lwd=1.5,lty=2) + theme_classic()
ggsave("../results/BC.RIN.hist.pdf", width = 3, height = 3)
## RIN vs. PMI
covs %>% filter(BRAINCODE.final.selection<2) %>% ggplot(aes(x=RIN,y=PMI)) + geom_point(na.rm = T) + theme_classic()
ggsave("../results/BC.RIN.vs.PMI.pdf", width = 3, height = 3)

# read coverage data for selected samples
geneBodyCoverage=data.frame();
for(i in (covs %>% filter(BRAINCODE.final.selection==1, RIN!='NA'))$sampleName){
  message(i)
  x=read.table(paste("~/neurogen/rnaseq_PD/run_output",i,"uniq/accepted_hits.bam.non-rRNA-mt.bam.geneBodyCoverage.txt",sep="/"), header = T, row.names = 1)
  rownames(x) = i
  geneBodyCoverage=rbind(geneBodyCoverage, x)
}
df=geneBodyCoverage 
df = (df-apply(df, 1, min))/(apply(df, 1, max)-apply(df,1,min))
df = df %>% rownames_to_column(var = 'sampleID') %>% mutate(binX50=X50, RIN=covs$RIN[match(sampleID, covs$sampleName)], cellType=gsub(".*_.*_(.*)_.*_.*","\\1",sampleID))
head(df); str(df)

## 3' bias is not exactly correlated with RIN. This is unexpected. Possibly because our RIN is already too good to separate the degradation?
# use 5'--3' bias to indicate degradation level 
# Ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3954844/figure/pone-0091851-g005/
# Ref: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0922-z#Fig7
# Ref: https://openi.nlm.nih.gov/detailedresult.php?img=PMC3185437_gkr547f4&req=4

dd=filter(df, cellType %in% c("SNDA","MCPY","TCPY")) %>% 
  gather(starts_with("X"), key='bins',value = "percentage") %>% 
  mutate(x=as.numeric(gsub("X","",bins))) 
pdf("../results/BC.geneBodyCoverage.pdf", width = 5, height = 5)
p=ggplot(dd,aes(x = x, y = percentage, group = 1)) + 
  stat_summary(fun.y = 'mean', colour = 'blue', geom = 'line') + 
  #stat_summary(fun.data = 'mean_sdl', fill = 'blue',geom = 'ribbon', alpha = 0.2) +
  stat_summary(data=subset(dd, binX50<.4), fun.y = 'mean', colour = 'red', geom = 'line') +
  #stat_summary(data=subset(df, RIN<7), fun.data = 'mean_sdl', fill = 'red',geom = 'ribbon', alpha = 0.2) +
  stat_summary(data=subset(dd, binX50>.75), fun.y = 'mean', colour = 'green', geom = 'line') + 
  #stat_summary(data=subset(df, RIN>9), fun.data = 'mean_sdl', fill = 'green',geom = 'ribbon', alpha = 0.2) +
  theme_bw() + theme(legend.position="top") 
print(p)
  
#ggplot(dd,aes(x = x, y = percentage, group = sampleID, col=cut(binX50, breaks=c(0,.4,.75,1)), alpha=RIN)) +
p=ggplot(dd,aes(x = x, y = percentage, group = sampleID, alpha=cut(binX50, breaks=c(0,.4,.75,1)))) +
  stat_summary(fun.y = 'mean', geom = 'line') + 
  labs(col="Grouping on P50", x="Gene body percentage (5' --> 3')", alpha="binP50") +
  theme_bw() + theme(legend.position="top") 
print(p)

ggplot(dd,aes(x = x, y = percentage, group = cut(binX50, breaks=c(0,.4,.75,1)))) +
  stat_summary(aes(alpha = cut(binX50, breaks=c(0,.4,.75,1))), fun.y = 'mean', geom = 'line') + 
  #stat_summary(aes(alpha = cut(binX50, breaks=c(0,.4,.75,1))), fun.data = 'mean_sdl', geom = 'ribbon', alpha = 0.2, show.legend = F) +
  labs(col="Grouping on P50", x="Gene body percentage (5' --> 3')", alpha="binP50") +
  theme_bw()   + theme(legend.position='top')
print(p)
dev.off()

table(cut(df$binX50, breaks=c(0,.4,.75,1)))

## S3/S5 in samples with different 3' bias level
# circRNAs: nS3 nS5
load(file = "Merge_circexplorer_BC.annotation.bed14.nS3nS5nU3nU5.RData")  # filtered_enriched_annotation in 106 samples
# controls: nS3 nS5
cnS=read.table("Merge_circexplorer_BC.annotation.bed14.matched.s3s5.gz", col.names=c("ID","type","sample","reads","count"), check.names = F, stringsAsFactors=F)
cnS3=filter(cnS,type=='s3') %>% select(ID, sample, reads) %>% spread(key=sample, value=reads, fill = 0) %>% column_to_rownames(var = "ID")
cnS5=filter(cnS,type=='s5') %>% select(ID, sample, reads) %>% spread(key=sample, value=reads, fill = 0) %>% column_to_rownames(var = "ID")
# matched ID
matched=read.table("Merge_circexplorer_BC.annotation.bed14.matched", stringsAsFactors = F); head(matched)
rows=matched$V4[match(rownames(nS3), matched$V15)]; 
cnS5=cnS5[rows, colnames(nS5)]; cnS3=cnS3[rows, colnames(nS3)]; dim(cnS3)

nS3S5=rbind(
  inner_join(x=rownames_to_column(nS3, var = 'ID') %>% gather(key="sample",value=nS3, contains("_")),
             y=rownames_to_column(nS5, var = 'ID') %>% gather(key="sample",value=nS5, contains("_")),
             by=c("ID","sample")) %>% 
    filter(nS5>0 & nS3>0) %>% mutate(celltype=gsub(".*_.*_(.*)_.*_.*","\\1",sample), rS3S5=(1+(nS3-nS5)/(nS3+nS5))/2, type='circRNA'),
  inner_join(x=rownames_to_column(cnS3, var = 'ID') %>% gather(key="sample",value=nS3, contains("_")),
             y=rownames_to_column(cnS5, var = 'ID') %>% gather(key="sample",value=nS5, contains("_")),
             by=c("ID","sample")) %>% 
    filter(nS5>0 & nS3>0) %>% mutate(celltype=gsub(".*_.*_(.*)_.*_.*","\\1",sample), rS3S5=(1+(nS3-nS5)/(nS3+nS5))/2, type='control')
)
head(nS3S5)
# limit to neuron, add X50, RIN
nS3S5 = filter(nS3S5, sample %in% df$sampleID, celltype %in% c("SNDA","MCPY","TCPY")) %>% 
  mutate(binX50=df$binX50[match(sample, df$sampleID)], RIN=df$RIN[match(sample, df$sampleID)]) %>% 
  separate(ID,c(NA,"start","end"),sep="_", convert=T)

library(plyr)

pdf("../results/BC.3p.bias.pdf", width = 5, height = 5)
p=ggplot(nS3S5, aes(x=nS3,y=nS5,col=celltype)) + 
  geom_point() + scale_y_log10() + scale_x_log10() + 
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(cols = vars(type))
print(p)

p=ggplot(nS3S5, aes(cut(RIN, breaks=c(0,7,8.5,10)),rS3S5)) +
  scale_y_continuous(breaks = seq(0, 1, by = .20)) + 
  geom_bar(aes(fill=type), position = "dodge", stat = "summary", fun.y = "mean") + 
  stat_summary(aes(col=type),position = position_dodge(width = 0.9),fun.data=mean_sdl, geom="errorbar", width=0.25) +
  labs(x="Grouped samples based on RIN", y="Ratio between S3 and S5 reads", title="No difference between different RIN") +
  theme_bw() +
  theme(legend.position='top')
print(p)

p=ggplot(nS3S5, aes(x=cut(binX50, breaks=c(0,.4,.75,1)),100*(rS3S5*2-1), alpha=cut(binX50, breaks=c(0,.4,.75,1)))) +
  scale_y_continuous(breaks = seq(-100, 100, by = 2), limits = c(-6,6), oob = rescale_none) + 
  geom_bar(aes(fill=type), position = "dodge", stat = "summary", fun.y = "mean") + 
  #stat_summary(aes(col=type),position = position_dodge(width = 0.9),fun.data=mean_sdl, geom="errorbar", width=0.25) +
  labs(x="Grouped samples based on 3' bias rate", y="Relative percentage betweeb S3 and S5 reads", alpha="P50",title="Difference between different 3' ratio") +
  theme_bw() +
  theme(legend.position='top')
print(p)

p=ggplot(nS3S5,aes(x=cut(binX50, breaks=c(0,.4,.75,1)),y=nS3-nS5, alpha=cut(binX50, breaks=c(0,.4,.75,1)))) +
  #scale_y_continuous(breaks = seq(0, 1, by = .20)) + 
  geom_bar(aes(fill=type),position = "dodge", stat = "summary", fun.y = "mean") + 
  #stat_summary(aes(col=type),position = position_dodge(width = 0.9),fun.data=mean_sdl, geom="errorbar", width=0.25) +
  labs(x="Grouped samples based on 3' bias rate", y="S3 - S5", alpha="P50",title="Difference between different 3' ratio") +
  theme_bw() + 
  theme(legend.position='top')
print(p)

# different length
p %+% subset(nS3S5, (end-start)>1500)
print(p)

dev.off()


#### update on single-exon circRNAs
#### =============================
setwd("~/projects/circRNA/data/") 

# circRNAs: nS3 nS5
load(file = "Merge_circexplorer_BC.annotation.bed14.nS3nS5nU3nU5.RData")  # filtered_enriched_annotation in 106 samples
single_exon_circRNAs=readRDS(file="Merge_circexplorer_BC106.filtered.enriched.annotation.bed14.rds") %>% 
  filter(exonCount==1, circType=='circRNA')
nS3=nS3[as.character(single_exon_circRNAs$ID),]; nS5=nS5[as.character(single_exon_circRNAs$ID),]

# controls: nS3 nS5
cnS=read.table("exons.internal.meta.bed.s3s5.gz", col.names=c("ID","type","sample","reads","count"), check.names = F, stringsAsFactors=F) %>% filter(grepl("protein_coding",ID)); dim(cnS)
cnS3=filter(cnS,type=='s3') %>% select(ID, sample, reads) %>% spread(key=sample, value=reads, fill = 0) %>% column_to_rownames(var = "ID")
cnS5=filter(cnS,type=='s5') %>% select(ID, sample, reads) %>% spread(key=sample, value=reads, fill = 0) %>% column_to_rownames(var = "ID")
# matched ID
matched=read.table("Merge_circexplorer_BC.annotation.bed14.matched", stringsAsFactors = F); head(matched)
rows=matched$V4[match(rownames(nS3), matched$V15)]; 
cnS5=cnS5[rows, colnames(nS5)]; cnS3=cnS3[rows, colnames(nS3)]; dim(cnS3)

nS3S5=rbind(
  inner_join(x=rownames_to_column(nS3, var = 'ID') %>% gather(key="sample",value=nS3, contains("_")),
             y=rownames_to_column(nS5, var = 'ID') %>% gather(key="sample",value=nS5, contains("_")),
             by=c("ID","sample")) %>% 
    filter(nS5>0 & nS3>0) %>% mutate(celltype=gsub(".*_.*_(.*)_.*_.*","\\1",sample), rS3S5=(1+(nS3-nS5)/(nS3+nS5))/2, type='circRNA'),
  inner_join(x=rownames_to_column(cnS3, var = 'ID') %>% gather(key="sample",value=nS3, contains("_")),
             y=rownames_to_column(cnS5, var = 'ID') %>% gather(key="sample",value=nS5, contains("_")),
             by=c("ID","sample")) %>% 
    filter(nS5>0 & nS3>0) %>% mutate(celltype=gsub(".*_.*_(.*)_.*_.*","\\1",sample), rS3S5=(1+(nS3-nS5)/(nS3+nS5))/2, type='control')
)
head(nS3S5)
with(nS3S5, table(type))
# limit to neuron, add X50, RIN
nS3S5 = filter(nS3S5, sample %in% df$sampleID, celltype %in% c("SNDA","MCPY","TCPY")) %>% 
  mutate(binX50=df$binX50[match(sample, df$sampleID)], RIN=df$RIN[match(sample, df$sampleID)]) %>% 
  separate(ID,c(NA,"start","end"),sep="_", convert=T)

library(plyr)

pdf("../results/BC.3p.bias.pdf", width = 5, height = 5)
p=ggplot(nS3S5, aes(x=nS3,y=nS5,col=celltype)) + 
  geom_point() + scale_y_log10() + scale_x_log10() + 
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(cols = vars(type))
print(p)

ggplot(nS3S5, aes(x=type, y=nS3-nS5)) + geom_boxplot(outlier.shape = NA) + ylim(c(-10,10))

p=ggplot(nS3S5, aes(cut(RIN, breaks=c(0,7,8.5,10)),rS3S5)) +
  scale_y_continuous(breaks = seq(0, 1, by = .20)) + 
  geom_bar(aes(fill=type), position = "dodge", stat = "summary", fun.y = "mean") + 
  stat_summary(aes(col=type),position = position_dodge(width = 0.9),fun.data=mean_sdl, geom="errorbar", width=0.25) +
  labs(x="Grouped samples based on RIN", y="Ratio between S3 and S5 reads", title="No difference between different RIN") +
  theme_bw() +
  theme(legend.position='top')
print(p)

p=ggplot(nS3S5, aes(x=cut(binX50, breaks=c(0,.4,.75,1)),100*(rS3S5*2-1), alpha=cut(binX50, breaks=c(0,.4,.75,1)))) +
  scale_y_continuous(breaks = seq(-100, 100, by = 2), limits = c(-6,6), oob = rescale_none) + 
  geom_bar(aes(fill=type), position = "dodge", stat = "summary", fun.y = "mean") + 
  #stat_summary(aes(col=type),position = position_dodge(width = 0.9),fun.data=mean_sdl, geom="errorbar", width=0.25) +
  labs(x="Grouped samples based on 3' bias rate", y="Relative percentage betweeb S3 and S5 reads", alpha="P50",title="Difference between different 3' ratio") +
  theme_bw() +
  theme(legend.position='top')
print(p)

p=ggplot(nS3S5,aes(x=cut(binX50, breaks=c(0,.4,.75,1)),y=nS3-nS5, alpha=cut(binX50, breaks=c(0,.4,.75,1)))) +
  #scale_y_continuous(breaks = seq(0, 1, by = .20)) + 
  geom_bar(aes(fill=type),position = "dodge", stat = "summary", fun.y = "mean") + 
  #stat_summary(aes(col=type),position = position_dodge(width = 0.9),fun.data=mean_sdl, geom="errorbar", width=0.25) +
  labs(x="Grouped samples based on 3' bias rate", y="S3 - S5", alpha="P50",title="Difference between different 3' ratio") +
  theme_bw() + 
  theme(legend.position='top')
print(p)

# different length
p %+% subset(nS3S5, (end-start)>1500)
print(p)

dev.off()