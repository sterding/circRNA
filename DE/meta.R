## Rscript for meta analysis

# modified based on invnorm function in the library('metaRNASeq')
invnorm <- function(indpval, nrep, BHth = 0.05) {
  listres = vector("list", 4)
  qnormpval= do.call(cbind, lapply(indpval, FUN=function(x) qnorm(1-x))) 
  nrepcorr=t(apply(qnormpval,1,FUN=function(x){
    nrep2=nrep
    nrep2[which(is.na(x))]=0
    nrep2
  }))
  nreptot=apply(nrepcorr,1,sum)
  weight=sqrt(nrepcorr/nreptot)
  wqnormp=weight*qnormpval
  statc=apply(wqnormp,1, FUN=function(x) sum(x,na.rm=TRUE))
  nan.index <- which(apply(wqnormp, 1, function(x) sum(is.nan(x))) == ncol(wqnormp))
  statc[nan.index] <- NA
  
  ## Update code for very small pvalues (from Umber Dube)
  ## Based on: https://stackoverflow.com/questions/5932276/adding-floating-point-precision-to-qnorm-pnorm-in-r
  rpvalc_gen <- function(x) { 
    y <- statc[[x]]
    if (y >= 8) {
      rpvalc = pnorm(y, lower.tail = FALSE)
    } else {
      rpvalc = 1 - pnorm(y)
    } 
  }
  rpvalc <- unlist(lapply(1:length(statc), rpvalc_gen))
  
  res = which(p.adjust(rpvalc, method = "BH") <= BHth) 
  listres[[1]] = res
  listres[[2]] = statc
  listres[[3]] = rpvalc
  listres[[4]] = p.adjust(rpvalc, method = "BH")
  names(listres) = c("DEindices", "TestStatistic", "rawpval", "adjpval") 
  return(listres)
}
library(tidyverse)
setwd("~/projects/circRNA/results/")


res_TCPY_AD = read.table("DE2gene_TCPY/DEresult.DE2gene_TCPY.CONDITION_AD_vs_HC.xls.gz", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)
head(res_TCPY_AD); dim(res_TCPY_AD)

Ns = list(BM10=154,BM22=135,BM36=132,BM44=129)
i=2
for(i in 1:4)
{
  res_NN_AD = read.table(paste0("../data/Dube2019NN/meta_both_condition_",names(Ns)[i],"_table_Box.txt"), skip = 3, stringsAsFactors = F, sep = "\t", row.names = 1, 
                         col.names = c("circID", "chr",NA,"log2FC_parietal","pvalue_parietal",NA,"log2FC_MSBB","pvalue_MSBB",NA,"pvalue_meta")) %>% dplyr::select(-contains("NA"))
  #head(res_NN_AD); dim(res_NN_AD)
  # trim
  res_TCPY_AD_trim = res_TCPY_AD[intersect(rownames(res_TCPY_AD),rownames(res_NN_AD)),]
  res_NN_AD_trim = res_NN_AD[intersect(rownames(res_TCPY_AD),rownames(res_NN_AD)),]
  dim(res_NN_AD_trim); dim(res_TCPY_AD_trim)
  rawpval_AD = list("pval1"=res_TCPY_AD_trim[["pvalue"]],"pval2"=res_NN_AD_trim[["pvalue_MSBB"]])
  # how many significant ones in discovery also signicant in replication?
  replicated = cbind(res_TCPY_AD_trim[rawpval_AD$pval1<0.05 & rawpval_AD$pval2<0.05, 1:12], res_NN_AD_trim[rawpval_AD$pval1<0.05 & rawpval_AD$pval2<0.05,])
  cor(replicated$log2FoldChange, replicated$log2FC_MSBB)
  
  invnormcomb_AD = invnorm(rawpval_AD,nrep=c(83,as.numeric(Ns)[i]), BHth = 0.05)
  invnormcomb_AD$DEindices <- NULL
  invnormcomb_AD = as.data.frame(invnormcomb_AD)
  rownames(invnormcomb_AD) <- rownames(res_NN_AD_trim)
  #head(invnormcomb_AD)
  message(paste(names(Ns)[i], sum(invnormcomb_AD$adjpval<0.05)))
  options(scipen=2)
  options(digits=3)
  select(res_TCPY_AD_trim, geneDescription, pvalue_TCPY=pvalue, log2FC_TCPY = log2FoldChange) %>% 
    rownames_to_column() %>%
    bind_cols(select(res_NN_AD_trim,Chr=chr,log2FC_MSBB,pvalue_MSBB), select(invnormcomb_AD, pvalue_meta=rawpval, padj_meta = adjpval)) %>% 
    filter(padj_meta <=0.05) %>% arrange(padj_meta) %>%
    mutate(pvalue_TCPY=format(pvalue_TCPY, scientific=T), pvalue_MSBB=format(pvalue_MSBB,scientific=T), 
           pvalue_meta=format(pvalue_meta, scientific=T), padj_meta=format(padj_meta, scientific=T),
           log2FC_TCPY = round(log2FC_TCPY,2), log2FC_MSBB=round(log2FC_MSBB, 2)) %>%
    column_to_rownames() %>% select(Chr, pvalue_TCPY, log2FC_TCPY, pvalue_MSBB, log2FC_MSBB, pvalue_meta, padj_meta, geneDescription) %>%
    write.table(paste0("DE2gene_TCPY/summaryTable.",names(Ns)[i],".meta.padj0.05.xls"), sep="\t", quote =F, na="", row.names=T, col.names = NA)
  
  # barplot
  df=read.table(paste0("DE2gene_TCPY/summaryTable.",names(Ns)[i],".meta.padj0.05.xls"), sep="\t", header = T, stringsAsFactors = F, row.names = 1)
  df = subset(df, pvalue_TCPY<=0.05 & pvalue_MSBB<=0.05, select=c("pvalue_meta","pvalue_MSBB","pvalue_TCPY")); dim(df)
  pdf(paste0("DE2gene_TCPY/summaryTable.",names(Ns)[i],".meta.padj0.05.pdf"), width=5, height=3)
  par(mar=c(4,7,2,3));
  barplot(as.matrix(t(-log10(df))), beside=T, horiz=TRUE, las=1, xlab=expression('-log'[10]*"(pvalue)"), main=paste0("DE2gene_TCPY/summaryTable.",names(Ns)[i],".meta.padj0.05")); 
  abline(v=c(-log10(0.05)))
  dev.off()
}
# BM10 19
# BM22 21
# BM36 30
# BM44 46


## meta between BC dopamine samples vs. Bennett midbrain samples
setwd("~/projects/circRNA/results/")

res1 = read.table("DE2gene_VMB/DEresult.DE2gene_VMB.CONDITION_PD_vs_HC.xls.gz", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)
head(res1); dim(res1)

res2 = read.table("DE2gene_SNDA/DEresult.DE2gene_SNDA.CONDITION2_PD_vs_HC.xls", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)
head(res2); dim(res2)

# trim
res1_trim = res1[intersect(rownames(res1),rownames(res2)),]
res2_trim = res2[intersect(rownames(res1),rownames(res2)),]
dim(res1); dim(res2); dim(res1_trim); dim(res2_trim)
rawpval = list("pval1"=res1_trim[["pvalue"]],"pval2"=res2_trim[["pvalue"]])

invnormcomb = invnorm(rawpval,nrep=c(23,104), BHth = 0.05)
invnormcomb$DEindices <- NULL
invnormcomb = as.data.frame(invnormcomb)
rownames(invnormcomb) <- rownames(res2_trim)
#head(invnormcomb_AD)
message(paste(sum(invnormcomb$adjpval<0.05)))
options(scipen=2)
options(digits=3)
select(res1_trim, geneDescription, pvalue_1=pvalue, log2FC_1 = log2FoldChange) %>% 
  rownames_to_column() %>%
  bind_cols(select(res2_trim,pvalue_2=pvalue,log2FC_2=log2FoldChange), select(invnormcomb, pvalue_meta=rawpval, padj_meta = adjpval)) %>% 
  filter(pvalue_meta <=0.05) %>% arrange(pvalue_meta) %>%
  column_to_rownames() %>% select(pvalue_VMB = pvalue_1, log2FC_VMB = log2FC_1, pvalue_SNDA=pvalue_2, log2FC_SNDA=log2FC_2, pvalue_meta, padj_meta, geneDescription) %>%
  write.table(paste0("DE2gene_SNDA/summaryTable.VMB.meta.p0.05.xls"), sep="\t", quote =F, na="", row.names=T, col.names = NA)

# barplot
df=read.table(paste0("DE2gene_SNDA/summaryTable.VMB.meta.p0.05.xls"), sep="\t", header = T, stringsAsFactors = F, row.names = 1)
head(df)
df = subset(df, pvalue_SNDA<=0.05 & pvalue_VMB<=0.1, select=c("pvalue_meta","pvalue_VMB","pvalue_SNDA")) 
pdf("DE2gene_SNDA/summaryTable.VMB.meta.p0.05.pdf", width=5, height=3)
par(mar=c(4,7,2,3));barplot(as.matrix(t(-log10(df))), beside=T, horiz=TRUE, las=1, xlab=expression('-log'[10]*"(pvalue)")); abline(v=c(-log10(0.05)))
dev.off()

## AD
res1 = read.table("DE2gene_FCX/DEresult.DE2gene_FCX.CONDITION_AD_vs_HC.xls.gz", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)
head(res1); dim(res1)

res2 = read.table("DE2gene_TCPY/DEresult.DE2gene_TCPY.CONDITION_AD_vs_HC.xls", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)
head(res2); dim(res2)

# trim
res1_trim = res1[intersect(rownames(res1),rownames(res2)),]
res2_trim = res2[intersect(rownames(res1),rownames(res2)),]
dim(res1); dim(res2); dim(res1_trim); dim(res2_trim)
rawpval = list("pval1"=res1_trim[["pvalue"]],"pval2"=res2_trim[["pvalue"]])

invnormcomb = invnorm(rawpval,nrep=c(19,83), BHth = 0.05)
invnormcomb$DEindices <- NULL
invnormcomb = as.data.frame(invnormcomb)
rownames(invnormcomb) <- rownames(res2_trim)
#head(invnormcomb_AD)
message(paste(sum(invnormcomb$adjpval<0.05)))
options(scipen=2)
options(digits=3)
select(res1_trim, geneDescription, pvalue_1=pvalue, log2FC_1 = log2FoldChange) %>% 
  rownames_to_column() %>%
  bind_cols(select(res2_trim,pvalue_2=pvalue,log2FC_2=log2FoldChange), select(invnormcomb, pvalue_meta=rawpval, padj_meta = adjpval)) %>% 
  filter(pvalue_meta <=0.05) %>% arrange(pvalue_meta) %>%
  column_to_rownames() %>% select(pvalue_FCX = pvalue_1, log2FC_FCX = log2FC_1, pvalue_TCPY=pvalue_2, log2FC_TCPY=log2FC_2, pvalue_meta, padj_meta, geneDescription) %>%
  write.table(paste0("DE2gene_TCPY/summaryTable.FCX.meta.p0.05.xls"), sep="\t", quote =F, na="", row.names=T, col.names = NA)



