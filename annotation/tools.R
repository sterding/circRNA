###############################################################################
# modified based on https://github.com/Bioconductor-mirror/cummeRbund/blob/master/R/tools.R
# Author: lgoff
###############################################################################

JSdistFromP<-function(mat,q){
  #row_js<-apply(mat,MARGIN=1,shannon.entropy)
  res<-apply(mat,MARGIN=1,function(p) {
    JSdistVec(p,q)
  }
  )
  res
}
JSdistVec<-function(p,q){
  JSdiv<-shannon.entropy((p+q)/2)-(shannon.entropy(p)+shannon.entropy(q))*0.5
  JSdist<-sqrt(JSdiv)
  JSdist
}

shannon.entropy <- function(p) {
  if (min(p) < 0 || sum(p) <=0)
    return(Inf)
  p.norm<-p[p>0]/sum(p)
  -sum( log2(p.norm)*p.norm)
}

# normalize each column of a matrix to its sum
makeprobs<-function(a){
  colSums<-apply(a,2,sum)
  b<-t(t(a)/colSums)
  b[is.na(b)] = 0
  b
}

# mat is a genes x samples matrix
specificity<-function(mat,logMode=T,pseudocount=1,relative=FALSE,...){
  if(logMode){
    mat<-log10(mat+pseudocount) # make sure min(mat)>=0
  }
  mat<-t(makeprobs(t(mat)))  # normalize each row to its sum e.g. sum=1 for each gene
  
  d<-diag(ncol(mat)) # for case with one sample per cell type
  
  res<-apply(d,MARGIN=1,function(q){
    JSdistFromP(mat,q)
  })
  colnames(res)<-paste(colnames(mat),"_spec",sep="")
  
  if(relative){
    res<-res/max(res)
  }
  1-res
}

## example code
# mat=matrix(rnorm(64, 0,0.5),8,8)+diag(1,8,8)
# mat[mat<0]=0
# colnames(mat)=paste0("sample",1:8)
# rownames(mat)=paste0("gene",1:8)
# specificity(mat, logMode = F)

groupmean_to_specificity <- function(groupmean){
  # The input groupmean is a matrix of FPKM, where the first column as gene names and the rest columns as sampleID
  rownames(groupmean)=groupmean[,1]; groupmean=groupmean[,-1]
  groupmean[groupmean<1e-10]=0  # limit to 1e-10
  groupmean = groupmean[rowMeans(groupmean)>0,]
  overmean = groupmean %>% mutate(gene=rownames(groupmean)) %>% melt(id.vars="gene", variable.name = "sampleID", value.name='fpkm') %>% 
    group_by(gene) %>%
    summarise(overall.mean=mean(fpkm), overall.sd=sd(fpkm)) %>% data.frame()
  rownames(overmean)=overmean[,1]; overmean=overmean[,-1]
  # specificty score
  groupmean_s = specificity(groupmean*1000, logMode = T)
  # Being specific: specificity score S>=0.5 AND mean expression > mean+s.d. of overall expression (Zheng et al., doi:10.1038/ncomms11215)
  celltypes=colnames(groupmean) # c("SNDA","MCPY","TCPY","FB","PBMC")
  df = data.frame(gene=rownames(groupmean_s), S=apply(groupmean_s,1,max), celltype=apply(groupmean_s,1,which.max))
  df = df %>% mutate(mean=groupmean[cbind(1:nrow(df), celltype)], m2sd=with(overmean,overall.mean+1*overall.sd)) %>% mutate(celltype=celltypes[celltype], Private_or_not=ifelse(S>=0.5 & mean>m2sd, 1, 0))
  return(cbind(groupmean_s, df))
}