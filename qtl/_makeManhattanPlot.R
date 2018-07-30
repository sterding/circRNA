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
              help="Path for eQTL result file", metavar="character"),
  make_option(c("-p", "--pcutoff"), type="double", default=0.01, 
              help="P-value cutoff for display [default=%default]", metavar="double"),
  make_option(c("-f", "--format"), type="character", default='fastQTL', 
              help="eQTL result format [default=%default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

nominalPzip=opt$input #nominalPzip="QTL_BC/eQTL.nominal.txt.gz"
p_cutoff=opt$pcutoff
input_format=opt$format

message(paste(nominalPzip, p_cutoff))
if(input_format=="fastQTL"){
  df = read.table(gzfile(nominalPzip), header = F, col.names=c("gene","SNP","distance","P","beta"))
  df2 = df %>% filter(P <= p_cutoff) %>%
    separate(SNP, c("CHROM","BP"), sep=":", remove=F, convert =T) %>% 
    mutate(CHR=as.numeric(ifelse(CHROM=='X',23,ifelse(CHROM=='Y',24,ifelse(CHROM=='MT',25,CHROM)))), P) %>% select(CHR, BP, P) %>% distinct()
}
if(input_format=='matrixEQTL' || input_format=='matrixQTL'){
  df = read.table(gzfile(nominalPzip), header = T)  # header: SNP	gene	beta	t.stat	p.value	FDR	SNP.pos	SNP.chr
  df2 = df %>% rename(P = p.value, BP=SNP.pos, CHROM=SNP.chr) %>%
    mutate(CHROM=gsub("chr","",CHROM)) %>%
    filter(P <= p_cutoff) %>%
    mutate(CHR=as.numeric(ifelse(CHROM=='X',23,ifelse(CHROM=='Y',24,ifelse(CHROM=='MT',25,CHROM)))), P) %>% select(CHR, BP, P) %>% distinct()
}
range(df2$P)
# down the resolution of the data to reduce the size of pdf
dim(df2)
df2 = mutate(df2, BP=as.integer(BP/10000)*10000, P=signif(P,2)) %>% distinct()
dim(df2)

pdf(file=paste0(nominalPzip, ".manhattan.pdf"), width = 7, height = 5, useDingbats=T)
manhattan(df2, suggestiveline = F, genomewideline =F, main = paste("Manhattan plot for", nominalPzip), cex = 0.5, cex.axis = 0.8)
dev.off()