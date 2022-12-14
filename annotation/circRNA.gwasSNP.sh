## GWAS snps distribution around circularized exons, comparing to the control

cd ~/projects/circRNA/data/
  
circRNA_annotation=Merge_circexplorer_BC.annotation.bed14  # BRAINCODE only
circRNA_annotation_control=Merge_circexplorer_BC.annotation.bed14.matched2

#GWAS
GENOME=~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19
snps_in_LD=$GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.SNAP.LD_w250.r2_0.8.bed  # in the final figure we used SNAP
## extract all autosomal.associations
[ -e $snps_in_LD.autosomal.associations.bed ] || awk 'BEGIN{FS="\t"; OFS="\t";}{split($8,a,"|");  n=split(a[2],b,";"); print $1,$2,$3,$7; for(i=1;i<n;i++) print $1,b[i]-1,b[i],$7;}' $snps_in_LD | grep -v chrX | grep -v chrY | sort -u > $snps_in_LD.autosomal.associations.bed

echo "## overlapped SNPs with each dataset"
### ##################
## [-1000,+1000] around the circRNA, then overlap with SNPs, the relateive location of SNP is then converted to a relative coordinate: [1-1000],[1001-2000],[2001-3000]
awk 'NR>1{OFS="\t"; print $1,$2-1000,$3+1000,$4,($13=="ciRNA")?0:1,$6;}' $circRNA_annotation | intersectBed -a - -b $snps_in_LD.autosomal.associations.bed -wo | awk '{FS="\t"; dis=$9-$2;len=$3-$2-2000;if($6=="+") coor=(dis<1000)?dis:((dis<(1000+len))?(1000+1000*(dis-1000)/len):(dis-1000-len+2000)); else coor=(dis<1000)?(2000+1000-dis):((dis<(1000+len))?(1000+1000*(len+1000-dis)/len):(1000-(dis-1000-len))); OFS="\t"; print $4,$5,$10,coor;}' > $circRNA_annotation.gwasSNP.relativeCoor.txt
paste <(sort -k15,15 $circRNA_annotation_control) <(awk 'NR>1' $circRNA_annotation| sort -k4,4) | cut -f1-15,28 | awk 'NR>1{OFS="\t"; print $1,$2-1000,$3+1000,$15,($16=="ciRNA")?0:1,$6;}' | intersectBed -a - -b $snps_in_LD.autosomal.associations.bed -wo | awk '{FS="\t"; dis=$9-$2;len=$3-$2-2000;if($6=="+") coor=(dis<1000)?dis:((dis<(1000+len))?(1000+1000*(dis-1000)/len):(dis-1000-len+2000)); else coor=(dis<1000)?(2000+1000-dis):((dis<(1000+len))?(1000+1000*(len+1000-dis)/len):(1000-(dis-1000-len))); OFS="\t"; print $4,$5,$10,coor;}' > $circRNA_annotation_control.gwasSNP.relativeCoor.txt

## [-1000,+1000] around the circRNA, and +/100bp into the circRNAs, then overlap with SNPs, the relateive location of SNP is then converted to a relative coordinate: [-1000,100] around both end of circRNAs
awk 'NR>1{OFS="\t"; print $1,$2-1000,$3+1000,$4,($13=="ciRNA")?0:1,$6;}' $circRNA_annotation | \
intersectBed -a - -b $snps_in_LD.autosomal.associations.bed -wo | \
awk '{FS="\t"; dis=$9-$2;len=$3-$2-2000;if($6=="+") coor=(dis<1100)?dis:((dis<(1000+len-100))?-1:(dis-1000-len+100+1100)); else coor=(dis<1100)?(1100+1100-dis):((dis<(1000+len-100))?-1:(2000+len-dis)); OFS="\t"; print $4,$5,$10,coor;}' > $circRNA_annotation.gwasSNP.relativeCoor1100.txt
paste <(sort -k15,15 $circRNA_annotation_control) <(awk 'NR>1' $circRNA_annotation| sort -k4,4) | cut -f1-15,28 | \
awk 'NR>1{OFS="\t"; print $1,$2-1000,$3+1000,$15,($16=="ciRNA")?0:1,$6;}' | \
intersectBed -a - -b $snps_in_LD.autosomal.associations.bed -wo | \
awk '{FS="\t"; dis=$9-$2;len=$3-$2-2000;if($6=="+") coor=(dis<1100)?dis:((dis<(1000+len-100))?-1:(dis-1000-len+100+1100)); else coor=(dis<1100)?(1100+1100-dis):((dis<(1000+len-100))?-1:(2000+len-dis)); OFS="\t"; print $4,$5,$10,coor;}' > $circRNA_annotation_control.gwasSNP.relativeCoor1100.txt

## plot in R
setwd("~/projects/circRNA/data/")
df1=read.delim("Merge_circexplorer_BC.annotation.bed14.gwasSNP.relativeCoor.txt", header = F, sep = "\t", quote = "", stringsAsFactors = F, col.names = c("circRNA","cirType","traits","coordinates"))
df2=read.delim("Merge_circexplorer_BC.annotation.bed14.matched2.gwasSNP.relativeCoor.txt", header = F, sep = "\t", quote = "", stringsAsFactors = F, col.names = c("circRNA","cirType","traits","coordinates"))
df=rbind(cbind(df1,type="circularized exon"),cbind(df2,type="uncircularized exon"))

library(tidyverse)
head(df)
ggplot(df, aes(coordinates, col = type)) + geom_histogram(aes(y=..density..), fill=NA,position = "identity", alpha=0.6,binwidth = 20) + geom_density(alpha=0.6)+ scale_color_brewer(palette="Dark2") + theme_minimal()+theme_classic()+theme(legend.position="top")

ggplot(subset(df, cirType==1), aes(coordinates, color = type)) + geom_histogram(aes(y=..density..), fill=NA,position = "identity", alpha=0.6,binwidth = 20) + geom_density(alpha=0.6)+ scale_color_brewer(palette="Dark2") + theme_minimal()+theme_classic()+theme(legend.position="top")

# total circular RNAs from BC106
BC106=readRDS("Merge_circexplorer_BC106.annotation.bed14.rds") # n=189128
filter(df, circRNA %in% BC106$ID) %>% ggplot(aes(coordinates, color = type)) + geom_histogram(aes(y=..density..), fill=NA,position = "identity", alpha=0.6,binwidth = 20) + geom_density(alpha=0.6)+ scale_color_brewer(palette="Dark2") + theme_minimal()+theme_classic()+theme(legend.position="top")

# For those colocalized with GWAS SNP, lenth distribution comparsion between called circRNAs and controls
filter(df, circRNA %in% BC106$ID) %>% separate(circRNA, c(NA,"start","end"), sep = "_", convert=T) %>% ggplot(aes(x=end-start, fill=type)) + geom_histogram(position = 'identity',bins = 50, alpha=.6)  + scale_colour_manual(values=c("red",'gray'))  + theme_minimal()+theme_classic()+theme(legend.position="top")+ scale_x_log10()  + ggsave("../results/Merge_circexplorer_BC106.annotation.lengthDistribution.+controls.pdf", width = 4,height = 3)
# lenth distribution for called circRNAs
filter(df, circRNA %in% BC106$ID, type=="circularized exon") %>% separate(circRNA, c(NA,"start","end"), sep = "_", convert=T) %>% ggplot(aes(x=end-start, fill=as.factor(cirType))) + geom_histogram(position = 'stack',bins = 60)  + scale_color_brewer(palette="Dark2") + theme_minimal()+theme_classic()+theme(legend.position="top")+ scale_x_log10() + scale_fill_manual(values=c("orange", "red")) + ggsave("../results/Merge_circexplorer_BC106.annotation.lengthDistribution.pdf", width = 4,height = 3)

filter(df, circRNA %in% BC106$ID) %>% ggplot(aes(coordinates, fill = type)) + 
  #geom_histogram(aes(y=..density..),fill=NA,position = "identity", alpha=0.6,binwidth = 20) + geom_density(alpha=0.6)+ 
  geom_histogram(position = "identity", alpha=0.8,binwidth = 20) + 
  facet_grid(rows = vars(cirType), scales = "free") +
  theme_minimal()+theme_classic()+theme(legend.position="top") + scale_fill_manual(values=c("red",'gray')) +
  ggsave("../results/Merge_circexplorer_BC106.annotation.GWASsnps.pdf", width = 5,height = 5)

# high-confident ones
high_confident_circRNAs=readRDS("Merge_circexplorer_BC106.filtered.enriched.annotation.bed14.rds") # n=10017
filter(df, circRNA %in% high_confident_circRNAs$ID, type=="circularized exon") %>% separate(circRNA, c(NA,"start","end"), sep = "_", convert=T) %>% ggplot(aes(x=end-start, fill=as.factor(cirType))) + geom_histogram(position = 'stack',bins = 60)  + scale_color_brewer(palette="Dark2") + theme_minimal()+theme_classic()+theme(legend.position="top")+ scale_x_log10() + scale_fill_manual(values=c("orange", "red"))+ ggsave("../results/Merge_circexplorer_BC106.filtered.enriched.annotation.lengthDistribution.pdf", width = 4,height = 3)
filter(df, circRNA %in% high_confident_circRNAs$ID) %>% ggplot(aes(coordinates, color = type)) + geom_histogram(aes(y=..density..), fill=NA,position = "identity", alpha=0.6,binwidth = 20) + geom_density(alpha=0.6)+ scale_color_brewer(palette="Dark2") + theme_minimal()+theme_classic()+theme(legend.position="top")

# exon length ==> mean > 100bp
mutate(high_confident_circRNAs,beginExon=gsub(",.*","",exonSizes),lastExon=gsub(".*,","",exonSizes)) %>% ggplot(aes(x=as.numeric(beginExon)+1, color=strand)) + geom_histogram(fill=NA,position = "identity", alpha=0.6) + scale_x_log10()


#### [-1000,+1000] around the circRNA, and +/100bp into the circRNAs, then overlap with SNPs, the relateive location of SNP is then converted to a relative coordinate: [-1000,100] around both end of circRNAs
df1=read.delim("Merge_circexplorer_BC.annotation.bed14.gwasSNP.relativeCoor1100.txt", header = F, sep = "\t", quote = "", stringsAsFactors = F, col.names = c("circRNA","cirType","traits","coordinates"))
df2=read.delim("Merge_circexplorer_BC.annotation.bed14.matched.gwasSNP.relativeCoor1100.txt", header = F, sep = "\t", quote = "", stringsAsFactors = F, col.names = c("circRNA","cirType","traits","coordinates"))
df=rbind(cbind(df1,type="circularized exon"),cbind(df2,type="uncircularized exon"))
filter(df, circRNA %in% BC106$ID, cirType==1, coordinates>0) %>% ggplot(aes(coordinates, color = type)) + geom_histogram(fill=NA,position = "identity", alpha=0.6,binwidth = 20) +  scale_color_brewer(palette="Dark2") + theme_minimal()+theme_classic()+theme(legend.position="top")
