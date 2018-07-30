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

mat=matrix(rnorm(64, 0,0.5),8,8)+diag(1,8,8)
mat[mat<0]=0
colnames(mat)=paste0("sample",1:8)
rownames(mat)=paste0("gene",1:8)
specificity(mat, logMode = F)
