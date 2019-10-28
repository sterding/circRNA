# prepare for eQTL
require(RCurl)
library(tidyr)
library(dplyr)

# n=84 HCILB_SNDA samples
sampleID=scan("~/neurogen/rnaseq_PD/results/merged/samplelist.HCILB_SNDA84",character())

# cov
covarianceTableURL="https://docs.google.com/spreadsheets/d/1I8nRImE9eJCCuZwpjfrrj-Uwx9bLebnO6o-ph7u6n8s/pub?gid=195725118&single=true&output=tsv"  # for all 140 samples
covs=read.delim(textConnection(getURL(covarianceTableURL)), stringsAsFactors =F)
covs=subset(covs, sampleName %in% sampleID, select = c(subjectID, batch, RIN, sex, age, PMI))
# convert batch to batchxxx
covs$batch=paste0("batch", covs$batch)
write.table(t(covs), "~/projects/circRNA/data/QTL/covs.fastqtl.txt", col.names = F, quote = F, sep="\t", row.names = T)

## expression
# genes
df1=read.table("~/neurogen/rnaseq_PD/results/merged/genes.loci.txt", header=T)
df2=read.table("~/neurogen/rnaseq_PD/results/merged/genes.fpkm.cuffnorm.allSamples.uniq.xls", header=T, check.names = F)
df2=df2[,c('tracking_id', sampleID)]
# filter: >0.1 RPKM in at least 10 individuals (32813 out of 57816 remained)
df2=df2[rowSums(df2[,sampleID]>0.1)>=10,]
df = merge(df1,df2,by='tracking_id')[,c(2:4,1,5:88)]
df$chr=gsub("chr","", df$chr)
df=df[with(df, order(chr, s1, s2)),]
colnames(df)=gsub(".*_(.*)_SNDA_.*","\\1",colnames(df))
colnames(df)[1]=paste0("#",colnames(df)[1]) # vi to add # at the beginning
write.table(df,"~/projects/circRNA/data/QTL/phenotype/genes.expression.bed", row.names=F, quote=F, sep="\t")

annotation=read.table(file="~/projects/circRNA/results/Merge_circexp_norm_filtered.celltype.annotation.xls", sep="\t", header=T, stringsAsFactors=F)
circularRNA_SNDA = subset(annotation, subset = grepl("SNDA", celltype3), select = gene, drop =T)
circRNA_SNDA = subset(annotation, subset = grepl("SNDA", celltype3) & circType=='circRNA', select = gene, drop =T)

# circRNA
load("~/projects/circRNA/data/Merge_circexp.BC.Rdata")
df = Merge_circexp_norm_filtered
df = df[,sampleID]
df = df[rowSums(df>0) >= 5,]  # >0 in at least 5 samples
df = df %>% mutate(ID=sub("chr","",rownames(.))) %>% 
  separate(ID, c("Chr","Start","End"), remove=F) %>%
  select(Chr, Start, End, ID, dplyr::contains("_")) %>% 
  arrange(Chr, as.numeric(Start), as.numeric(End)) %>%
  setNames(gsub(".*_(.*)_SNDA_.*","\\1",names(.)))
# only circRNA or all circularRNA
colnames(df)[1]=paste0("#",colnames(df)[1]) # vi to add # at the beginning
write.table(df[paste0("chr", df$ID) %in% circularRNA_SNDA,],"~/projects/circRNA/data/QTL/phenotype/circularRNA.expression.bed", row.names=F, quote=F, sep="\t")
write.table(df[paste0("chr", df$ID) %in% circRNA_SNDA,],"~/projects/circRNA/data/QTL/phenotype/circRNA.expression.bed", row.names=F, quote=F, sep="\t")

## circularization
df = readRDS("~/projects/circRNA/data/Merge_circexplorer_BC.annotation.bed14.cRatio.circRNA_norm_filtered.rds")
df = df[,sampleID]; dim(df)
#df = df[rowSums(df>0) >= 5,]; dim(df)  # >0 in at least 5 samples
df = as.data.frame(df) %>% mutate(ID=sub("chr","",rownames(df))) %>% 
  separate(ID, c("Chr","Start","End"), remove=F) %>%
  select(Chr, Start, End, ID, dplyr::contains("_")) %>% 
  arrange(Chr, as.numeric(Start), as.numeric(End)) %>%
  setNames(gsub(".*_(.*)_SNDA_.*","\\1",names(.)))
# only circRNA or all circularRNA
colnames(df)[1]=paste0("#",colnames(df)[1]) # vi to add # at the beginning
write.table(df[paste0("chr", df$ID) %in% circularRNA_SNDA,],"~/projects/circRNA/data/QTL/phenotype/circularRNA.PCI.bed", row.names=F, quote=F, sep="\t")
write.table(df[paste0("chr", df$ID) %in% circRNA_SNDA,],"~/projects/circRNA/data/QTL/phenotype/circRNA.PCI.bed", row.names=F, quote=F, sep="\t")

## splicing
# bash
# for i in `seq 1 22`; do echo $i;  cat leafcutter/sqtl_tmpdir/leafcutter_perind.counts.gz.qqnorm_chr$i | sed 's/.accepted_hits.bam//g;s/HC_//g;s/ILB_//g;s/_SNDA_[1-9]_rep[12]//g' | awk '{OFS="\t"; s=0; for(i=5;i<=NF;i++) if($i>0) s++; if(NR==1 || s>=1) print}' > circRNA.PSI.chr$i.bed & done
# for i in *bed; do echo $i; bgzip $i && tabix -p bed $i.gz; done