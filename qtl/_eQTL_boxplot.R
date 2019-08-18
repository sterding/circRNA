###############################################
## Rscript for generating boxplot and sashimi plot of expression vs. genotype for gene-SNP pairs (e.g. from the eQTL output)
## Author: Xianjun Dong
## Date: 2018-7-2
## Version: 0.0
## Usage: 
## Rscript ~/projects/circRNA/src/qtl/_eQTL_boxplot.R -i eQTLcircularRNAmerged.nominal.txt.top.list -p ~/projects/circRNA/data/QTL_BC -g phenotype/circularRNAmerged.expression.bed.gz -e eQTLcircularRNAmerged.nominal.txt.gz -f fastQTL

###############################################
require(MatrixEQTL)
require(tidyverse)
library(vcfR) # install.packages('vcfR')

require("optparse",quietly=T, warn.conflicts=F) || install.packages('optparse', repo='http://cran.revolutionanalytics.com') && require('optparse');

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Path for gene SNP list file", metavar="character"),
  make_option(c("-p", "--path"), type="character", default=getwd(), 
              help="working path", metavar="character"),
  make_option(c("-g", "--geneexp"), type="character", default="expression.postSVA.xls", 
              help="Gene expression matrix [default=%default]", metavar="character"),
  make_option(c("-e", "--eqtlfile"), type="character", default='final.cis.eQTL.d1e6.p1e-2.xls', 
              help="eQTL result file [default=%default]", metavar="character"),
  make_option(c("-f", "--format"), type="character", default='matrixQTL', 
              help="eQTL result format [default=%default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

GSfile=opt$input 
path=opt$path
eqtlfile=opt$eqtlfile
input_format=opt$format
geneexpfile=opt$geneexp

setwd(path); 

# setwd("~/projects/circRNA/results/eQTL/HCILB_SNDA"); eqtlfile = "final.cis.eQTL.d1e6.p1e-2.xls"; geneexpfile="expression.postSVA.xls"; input_format="matrixQTL"; G="chr10_32761408_32807433"; S="exm818678_G:T"; 
# setwd("~/projects/circRNA/data/QTL_BC"); eqtlfile="eQTLcircularRNAmerged.nominal.txt.gz"; geneexpfile="phenotype/circularRNAmerged.expression.bed.gz"; input_format="fastQTL"; GSfile = "eQTLcircularRNAmerged.nominal.txt.top.list"; 

## example
# GS=subset(d, ID=='chr2_40655612_40657444', select = c("ID","SNP")); # where d is from end of eQTL.main.R

message("# load data...")
######################

## GS list
GS=read.table(GSfile, header = F, stringsAsFactors = F)

## SNP loci info
snpspos_rds="~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID.rds"
if(file.exists(snpspos_rds)) snpspos=readRDS(snpspos_rds) else {
  snpspos = read.table("~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID", header = TRUE, stringsAsFactors = FALSE); # header: snp	chr	pos	impute2
  snpspos = snpspos %>% unite(chrpos,chr,pos,sep = ":", remove = F) %>% separate(snp, c("RS","GT"), sep = "_", remove = F, convert = T) %>% separate(GT, c("REF","ALT"),sep = ":", remove = T, convert = T)
  saveRDS(snpspos, snpspos_rds)
}

## eQTL, gene expression
if(input_format=="fastQTL"){
  genesnp = read.table(eqtlfile, header = F, col.names=c("gene","SNP","distance","p.value","beta"), stringsAsFactors = F)
  GS = genesnp %>% filter(gene %in% GS[,1], SNP %in% GS[,2]) %>% left_join(y=snpspos, by=c("SNP" = "chrpos"))
  genes = read.table(geneexpfile, check.names = F, header = T, comment.char = "")[,-c(1:3)] %>% column_to_rownames(var = "ID")
}
if(input_format=='matrixEQTL' || input_format=='matrixQTL' || input_format=='matrixeQTL'){
  genesnp = read.table(eqtlfile, header = T, stringsAsFactors = F)  # header: SNP	gene	beta	t.stat	p.value	FDR	SNP.pos	SNP.chr
  GS = genesnp %>% filter(gene %in% GS[,1], SNP %in% GS[,2]) %>% left_join(y=snpspos, by=c("SNP" = "snp"))
  genes = read.table(geneexpfile, check.names=F, header = T)
}

#sort by chr in order to save time in loading genotype for different chr
GS = arrange(GS, chr);
CHR="";gt=NULL;
######################
apply(GS, 1, function(gs) {
  G=gs['gene']; S=gs['SNP'];
  message(paste0("# making eQTL plot for ",S," and ",G," ..."))

  RS=sub("([^:]*):.*","\\1", gs['RS']) ## the part before the first ":" in the SNP ID
  REF = gs['REF']  ## get the REF allele
  ALT = gs['ALT']  ## get the ALT allele
  p = signif(as.numeric(gs['p.value']), 3)
  beta = gs['beta']
  snp=gs['snp']
  chr=gs['chr']
  
  ## read genotype from vcf
  if(chr != CHR){
    vcf_file <- file.path("genotype", paste0(chr,".dos.postQC.vcf.gz"))
    vcf <- read.vcfR( vcf_file, verbose = FALSE )
    # remove duplicated SNPs ## some multialleic SNPs are save into multiple lines of bialleic SNPs in our case, which should have been remove ahead using "bcftools norm -m +snps"
    vcf = vcf[!duplicated(vcf@fix[,'ID']),]
    gt <- extract.gt(vcf, element = 'GT')    
  }

  # extract the SNP and recode using the number of REF (same way as the BRAINcode All.Matrix.txt by Ganqiang)
  gt0 = sapply(gt[S,], str_count, pattern = "0") # takes long time...
  
  ## create data.frame including expression and SNP
  samplenames = intersect(names(gt0), colnames(genes)); length(samplenames)
  df=data.frame(expression=as.numeric(genes[G,samplenames]), SNP=as.numeric(gt0[samplenames]), row.names = samplenames)
  ## write data to txt file
  write.table(df, file=paste("eQTLboxplot",G,gsub(":","..",S),"xls",sep="."), col.names = NA,quote=F, sep="\t")
  

  df$SNP1=factor(df$SNP, levels=2:0)  ## in the current All.Matrix.txt, the number is number of REF allele (since we use REF in the --recode-allele) -- Ganqiang
  if(is.na(p) || length(p)==0) {
    test=aov(expression~SNP1, df)
    p=signif(summary(test)[[1]][["Pr(>F)"]][[1]],3)  # in case SNP:eQTL pair is not in the eqtl result file (e.g. those not or less significant pairs)
  }
  df$SNP2=ifelse(as.numeric(as.character(df$SNP))==0,0,1)  ## 0: without REF; 1:with REF
  df$SNP3=ifelse(as.numeric(as.character(df$SNP))==2,0,1)  ## 0: without ALT; 1:with ALT
  
  ## boxplot
  
  pdf(file=paste("eQTLboxplot",G,gsub(":","..",S),"pdf",sep="."), width=6, height=4)
  par(mfrow=c(1,3), mar=c(4,4,2,1), oma = c(0, 0, 2, 0))
  bp=boxplot(expression~SNP1, data=df, ylab="Normalized expression", xaxt='n', main="",  col='lightgreen', outpch=NA)
  stripchart(expression~SNP1, data=df, vertical = TRUE, method = "jitter", pch = 1, col = "darkred", cex=1, add = TRUE) 
  title(main=paste0("additive effect (p = ", p,")"), cex.main=1, line=0.5)
  mtext(c("Homo Ref","Het","Homo Alt"), cex=0.7,side=1,line=.5,at=1:3)
  mtext(paste0(c(REF,REF,ALT),c(REF,ALT,ALT)),  cex=0.7,side=1,line=1.5,at=1:3)
  mtext(paste0("N = ", bp$n),  cex=0.7,side=1,line=2.5,at=1:3)
  
  # with/without REF
  df$SNP2=factor(df$SNP2, levels=1:0)
  if(length(unique(df$SNP2))>1 & min(table(df$SNP2))>1) p0=signif(t.test(expression~SNP2, df)$p.value,3) else {p0=NA; warning("Grouping factor must have exactly 2 levels  AND at least 2 values in each group")}
  bp=boxplot(expression~SNP2, data=df, ylab="Normalized expression", xaxt='n', main="", col='lightblue', outpch=NA)
  stripchart(expression~SNP2, data=df, vertical = TRUE, method = "jitter", pch = 1, col = "darkred", cex=1, add = TRUE)
  title(main=paste0("dominant effect (p = ",p0,")"), cex.main=1, line=0.5)
  mtext(c("with REF","without REF"), cex=0.7,side=1,line=0.5,at=1:2)
  mtext(c(paste0(c(REF,REF),c(REF,ALT),collapse = "/"), paste0(ALT,ALT)),  cex=0.7,side=1,line=1.5,at=1:2)
  mtext(paste0("N = ", bp$n), cex=0.7,side=1,line=2.5,at=1:2)
  
  # with/without ALT
  df$SNP3=ifelse(as.numeric(as.character(df$SNP))==2,0,1)  ## 0: without ALT; 1:with ALT
  if(length(unique(df$SNP3))>1 & min(table(df$SNP3))>1) p0=signif(t.test(expression~SNP3, df)$p.value,3) else {p0=NA; warning("Grouping factor must have exactly 2 levels  AND at least 2 values in each group")}
  df$SNP3=factor(df$SNP3, levels=0:1)
  bp=boxplot(expression~SNP3, data=df, ylab="Normalized expression", xaxt='n', main="", col='lightblue', outpch=NA)
  stripchart(expression~SNP3, data=df, vertical = TRUE, method = "jitter", pch = 1, col = "darkred", cex=1, add = TRUE)
  title(main=paste0("dominant effect (p = ",p0,")"), cex.main=1, line=0.5)
  mtext(c("without ALT","with ALT"), cex=0.7,side=1,line=0.5,at=1:2)
  mtext(c(paste0(REF,REF), paste0(c(REF,ALT),c(ALT,ALT),collapse = "/")),  cex=0.7,side=1,line=1.5,at=1:2)
  mtext(paste0("N = ", bp$n), cex=0.7,side=1,line=2.5,at=1:2)
  
  mtext(paste("eQTL for",G,"and",S), outer = TRUE, cex = 1)     
  
  ## without jitter and outliers (optional)
  ##############################
  
  par(mfrow=c(1,3), mar=c(4,4,2,1), oma = c(0, 0, 2, 0))
  bp=boxplot(expression~SNP1, data=df, ylab="Normalized expression", xaxt='n', main="",  col='lightgreen', outline = F, outpch=NA)
  title(main=paste0("additive effect (p = ", p,")"), cex.main=1, line=0.5)
  mtext(c("Homo Ref","Het","Homo Alt"), cex=0.7,side=1,line=.5,at=1:3)
  mtext(paste0(c(REF,REF,ALT),c(REF,ALT,ALT)),  cex=0.7,side=1,line=1.5,at=1:3)
  mtext(paste0("N = ", bp$n),  cex=0.7,side=1,line=2.5,at=1:3)
  
  # with/without REF
  if(length(unique(df$SNP2))>1 & min(table(df$SNP2))>1) p0=signif(t.test(expression~SNP2, df)$p.value,3) else {p0=NA; warning("Grouping factor must have exactly 2 levels  AND at least 2 values in each group")}
  bp=boxplot(expression~SNP2, data=df, ylab="Normalized expression", xaxt='n', main="", col='lightblue', outline = F, outpch=NA)
  title(main=paste0("dominant effect (p = ",p0,")"), cex.main=1, line=0.5)
  mtext(c("with REF","without REF"), cex=0.7,side=1,line=0.5,at=1:2)
  mtext(c(paste0(c(REF,REF),c(REF,ALT),collapse = "/"), paste0(ALT,ALT)),  cex=0.7,side=1,line=1.5,at=1:2)
  mtext(paste0("N = ", bp$n), cex=0.7,side=1,line=2.5,at=1:2)
  
  # with/without ALT
  if(length(unique(df$SNP3))>1 & min(table(df$SNP3))>1) p0=signif(t.test(expression~SNP3, df)$p.value,3) else {p0=NA; warning("Grouping factor must have exactly 2 levels  AND at least 2 values in each group")}
  bp=boxplot(expression~SNP3, data=df, ylab="Normalized expression", xaxt='n', main="", col='lightblue', outline = F, outpch=NA)
  title(main=paste0("dominant effect (p = ",p0,")"), cex.main=1, line=0.5)
  mtext(c("without ALT","with ALT"), cex=0.7,side=1,line=0.5,at=1:2)
  mtext(c(paste0(REF,REF), paste0(c(REF,ALT),c(ALT,ALT),collapse = "/")),  cex=0.7,side=1,line=1.5,at=1:2)
  mtext(paste0("N = ", bp$n), cex=0.7,side=1,line=2.5,at=1:2)
  
  mtext(paste("eQTL for",G,"and",S), outer = TRUE, cex = 1)     
  
  dev.off() 
  
})