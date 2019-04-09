#!/usr/bin/env Rscript
###########################################
# Rscript to make manhattan plot for QTL result
# author: Xianjun Dong
# email: xdong@rics.bwh.harvard.edu
# date: 5/21/2018
# version: 1.0
# Usage: $0 -i eQTL.nominal.txt.gz -p 0.01 -f fastQTL
###########################################

# install packages
require('qqman', quietly =T, warn.conflicts=F) || install.packages('qqman', repo='http://cran.revolutionanalytics.com') && require('qqman'); 
require('tidyverse',quietly=T, warn.conflicts=F) || install.packages('tidyverse', repo='http://cran.revolutionanalytics.com') && require('tidyverse');
require("optparse",quietly=T, warn.conflicts=F) || install.packages('optparse', repo='http://cran.revolutionanalytics.com') && require('optparse');

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Path for eQTL result file (zipped or unzipped)", metavar="character"),
  make_option(c("-p", "--pcutoff"), type="double", default=0.01, 
              help="P-value cutoff for display [default=%default]", metavar="double"),
  make_option(c("-P", "--Psignifiance"), type="double", default=0.001, 
              help="Where to draw a genome-wide sigificant line [default=%default]", metavar="double"),
  make_option(c("-f", "--format"), type="character", default='fastQTL', 
              help="eQTL result format [default=%default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

nominalPzip=opt$input 
p_cutoff=opt$pcutoff
Psignifiance=opt$Psignifiance
input_format=opt$format

#debug
# setwd("~/projects/circRNA/data/QTL_BC"); nominalPzip="eQTL.nominal.txt.gz"; p_cutoff=.01; Psignifiance=0.0009176055; input_format='fastQTL'
# setwd("~/neurogen/targetSeq_BRAINcode"); nominalPzip="eQTLgene.nominal.txt.gz"; p_cutoff=1; Psignifiance=0.005290778; input_format='fastQTL'

message(paste(nominalPzip, p_cutoff))
if(input_format=="fastQTL"){
  df = read.table(nominalPzip, header = F, col.names=c("gene","SNP","distance","P","beta"))
  df2 = df %>% filter(P <= p_cutoff) %>%
    separate(SNP, c("CHROM","BP"), sep=":|_", remove=F, convert =T) %>% 
    mutate(CHROM=gsub("chr","",CHROM)) %>%
    mutate(CHR=as.numeric(ifelse(CHROM=='X',23,ifelse(CHROM=='Y',24,ifelse(CHROM=='MT',25,CHROM)))), P) %>% select(gene, CHR, BP, P) %>% distinct()
}
if(input_format=='matrixEQTL' || input_format=='matrixQTL' || input_format=='matrixeQTL'){
  df = read.table(nominalPzip, header = T)  # header: SNP	gene	beta	t.stat	p.value	FDR	SNP.pos	SNP.chr
  df2 = df %>% rename(P = p.value, BP=SNP.pos, CHROM=SNP.chr) %>%
    filter(P <= p_cutoff) %>%
    mutate(CHROM=gsub("chr","",CHROM)) %>%
    mutate(CHR=as.numeric(ifelse(CHROM=='X',23,ifelse(CHROM=='Y',24,ifelse(CHROM=='MT',25,CHROM)))), P) %>% select(gene, CHR, BP, P) %>% distinct()
}
range(df2$P)
# down the resolution of the data to reduce the size of pdf
dim(df2)
df2 = mutate(df2, BP=as.integer(BP/1000000)*1000000, P=signif(P,2)) %>% distinct()
dim(df2)

pdf(file=paste0(nominalPzip, ".manhattan.pdf"), width = 7, height = 5, useDingbats=T)
manhattan(df2, suggestiveline = F, genomewideline =-log10(Psignifiance), chrlabs=c(1:22, "X"), main = paste("Manhattan plot for", nominalPzip), cex = 0.5, cex.axis = .5, ylim=c(0, 1.05*max(-log10(range(df2$P)))))

# for subset of circRNA eQTL 
if(nominalPzip=="eQTL.nominal.txt.gz"){
  annotation= readRDS("~/projects/circRNA/data/Merge_circexplorer_BC106.annotation.bed14.rds")
  df3 = filter(df2, gene %in% annotation$ID[annotation$circType=='circRNA']); dim(df3)
  manhattan(df3, suggestiveline = F, genomewideline =-log10(Psignifiance), chrlabs=c(1:22, "X"), main = paste("Manhattan plot for circRNA subset (n=752) in", nominalPzip), cex = 0.5, cex.axis = .5, ylim=c(0, 1.05*max(-log10(range(df2$P)))))
  
  df3 = filter(df2, gene %in% annotation$ID[annotation$circType=='ciRNA']); dim(df3)
  manhattan(df3, suggestiveline = F, genomewideline =-log10(Psignifiance), chrlabs=c(1:22, "X"), main = paste("Manhattan plot for ciRNA subset (n=302) in", nominalPzip), cex = 0.5, cex.axis = .5, ylim=c(0, 1.05*max(-log10(range(df2$P)))))
  
  annotation_filtered= readRDS("~/projects/circRNA/data/Merge_circexplorer_BC106.filtered.enriched.annotation.bed14.rds")
  df3 = filter(df2, gene %in% annotation_filtered$ID)  ; dim(df3)
  manhattan(df3, suggestiveline = F, genomewideline =-log10(Psignifiance), chrlabs=c(1:22, "X"), main = paste("Manhattan plot for filtered&enriched subset(n=601) in", nominalPzip), cex = 0.5, cex.axis = .5, ylim=c(0, 1.05*max(-log10(range(df2$P)))))
  
  df3 = filter(df2, gene %in% annotation_filtered$ID[annotation_filtered$circType=='circRNA'])  ; dim(df3)
  manhattan(df3, suggestiveline = F, genomewideline =-log10(Psignifiance), chrlabs=c(1:22, "X"), main = paste("Manhattan plot for filtered&enriched circRNA subset (n=566) in", nominalPzip), cex = 0.5, cex.axis = .5, ylim=c(0, 1.05*max(-log10(range(df2$P)))))

}
dev.off()
