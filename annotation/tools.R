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

## rewrite the disease_enrichment function in the disgenet2r package

disease_enrichment_v2 <- function (genes, vocabulary = "HGNC", verbose = TRUE, gdas = disgenet_CURATED, warnings = TRUE) 
{
  if (length(genes) != length(unique(genes))) {
    genes <- unique(genes)
    warning("Removing duplicates from input genes list.")
  }
  type <- "symbol"
  if (vocabulary == "HGNC") {
  }
  else if (vocabulary == "ENTREZ") {
    type <- "entrez"
  }
  else if (class(genes[1]) == "factor") {
    message("Your input genes are in factor format.\n    
            Please, revise your input genes and save them as numeric or as character.")
    stop()
  }
  else {
    message(paste0("Your input genes are a wrong vocabulary ", 
                   vocabulary, "Please, revise your input vocabulary. Remember that genes should be identified\n                   
                   using HGNC gene symbols or NCBI Entrez identifiers"))
    stop()
  }
  diseases <- unique(gdas[c("diseaseId", "diseaseName")])
  if (type == "symbol") {
    genes <- unique(subset(gdas, geneSymbol %in% genes)$geneId) # Note: here only the genes in gdas$geneSymbol are taken into test
    universe <- unique(as.character(gdas$geneSymbol))
  }
  else {
    genes <- intersect(genes, gdas$geneId)
    universe <- unique(as.character(gdas$geneId))
  }
  if (length(genes) > 0) {
    message(paste0("A total of ", length(genes), " from the initial list are annotated in DisGENET"))
  }
  else {
    stop("The list does not contain genes annotated in DisGENET \n\n
         remember that the genes should be identified as entrez gene identifiers")
  }
  if (length(universe) == 1) {
    message(paste0("A total of ", length(universe), " extracted from DisGeNET are being used as universe \n"))
  }
  else {
    message(paste0("A total of ", length(universe), " are being used as universe \n"))
  }
  data <- data.frame(ID = character(), Description = character(), 
                     GeneRatio = character(), BgRatio = character(), OR=numeric(), geneID = character(), 
                     pvalue = numeric(), Count = integer(), stringsAsFactors = FALSE)
  i <- 1
  diseaselist <- as.character(unique(gdas[gdas$geneId %in% genes, ]$diseaseId))  # test only those diseases with at least 1 associated genes in the input gene set
  for (dd in diseaselist) {
    aa <- subset(gdas, diseaseId == dd)$geneId
    inter <- length(intersect(aa, genes))
    t <- matrix(c(inter, length(aa) - inter, length(genes) - 
                    inter, length(universe) + inter - length(aa) - length(genes)), 
                nrow = 2, dimnames = list(module = c("in", "out"), 
                                          Pheno = c("phen", "nophen")))
    test <- fisher.test(t, alternative = "greater")
    pv <- test$p.value
    OR <- as.numeric(test$estimate) ## Odds ratio
    GeneRatio <- paste(inter, "/", length(genes), sep = "")
    BgRatio <- paste(length(aa), "/", length(universe), sep = "")
    geneID <- paste(intersect(aa, genes), collapse = "/") ## can be changed to geneSymbol 
    data[i, ] <- c(as.character(dd), as.character(subset(diseases, diseaseId == dd)$diseaseName), 
                   GeneRatio, BgRatio, OR, geneID, pv, inter)
    i <- i + 1
  }
  data$pvalue <- as.numeric(data$pvalue)
  data$FDR <- p.adjust(data$pvalue, method = "BH")
  data$Count <- as.numeric(as.character(data$Count))
  data$gg <- data$Count/length(genes)
  data$OR <- as.numeric(data$OR)
  data <- data[order(as.numeric(as.character(data$FDR))), ]
  return(data)
}
# revised plot_enrichment @ disgenet2r package
plot_enrichment <- function( input, cutoff , count, title="DisGeNET enrichment", limit ) {
  input<- subset(input, FDR < cutoff & Count > count  )
  if ( dim( input )[ 1 ] > limit ){
    input <- input[ 1:limit ,]
    show( paste0("warning: dataframe of ", nrow(input), " rows has been reduced to ", limit, " rows." ))
  }
  
  idx <- order(input$OR, decreasing = T)
  input$Description<-factor(input$Description, levels=rev(unique(input$Description[idx])))
  p <- ggplot2::ggplot(input, ggplot2::aes_string(x= input$OR, y="Description", size=input$Count, color=input$FDR)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_continuous(low="red", high="blue", name = "FDR", guide=ggplot2::guide_colorbar(reverse=TRUE)) +
    ggplot2::xlab("Odds Ratio") + ggplot2::ylab(NULL) +
    ggplot2::ggtitle(title) + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(colour = "black",  size = 14, vjust = 1),
                   axis.text.y = ggplot2::element_text(colour = "black", size = 14, hjust = 1),
                   axis.title = ggplot2::element_text(margin = ggplot2::margin(10, 5, 0, 0), color = "black", size = 12),
                   axis.title.y = ggplot2::element_text(angle = 90)) + ggplot2::scale_size(range=c(3, 8))
  print(p)
}
