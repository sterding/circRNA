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
library(dplyr)
setwd("~/projects/circRNA/results/")
res_TCPY_AD = read.table("DE2gene_TCPY/DEresult.DE2gene_TCPY.CONDITION_AD_vs_HC.xls", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)
head(res_TCPY_AD); dim(res_TCPY_AD)

Ns = list(BM10=154,BM22=135,BM36=132,BM44=129)
i=4
for(i in 1:4)
{
  res_NN_AD = read.table(paste0("../data/Dube2019NN/meta_both_condition_",names(Ns)[i],"_table_Box.txt"), skip = 3, stringsAsFactors = F, sep = "\t", row.names = 1, 
                         col.names = c("circID", "chr",NA,"log2FC_parietal","pvalue_parietal",NA,"log2FC_MSBB","pvalue_MSBB",NA,"pvalue_meta")) %>% dplyr::select(-contains("NA"))
  #head(res_NN_AD); dim(res_NN_AD)
  # trim
  res_TCPY_AD_trim = res_TCPY_AD[intersect(rownames(res_TCPY_AD),rownames(res_NN_AD)),]
  res_NN_AD_trim = res_NN_AD[intersect(rownames(res_TCPY_AD),rownames(res_NN_AD)),]
  #dim(res_NN_AD_trim); dim(res_TCPY_AD_trim)
  rawpval_AD = list("pval1"=res_TCPY_AD_trim[["pvalue"]],"pval2"=res_NN_AD_trim[["pvalue_MSBB"]])
  
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
}
# BM10 19
# BM22 21
# BM36 30
# BM44 46


