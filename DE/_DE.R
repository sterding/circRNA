###########################################
# A general framework in R for running differntial expression analysis using DESeq2
# author: Xianjun Dong
# email: xdong@rics.bwh.harvard.edu
# date: 9/3/2019
# version: 2.0
# Usage: Rscript --vanilla ../src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -C CONDITION:PD:HC
# Requirement:
# 1. The headers in covariate file should be in CAPITAL format for SAMPLE_ID (required), SUBJECT_ID (if any), CELLTYPE (optional), CONDITION (required)
###########################################
library("optparse")
options(stringsAsFactors=F)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input expression dataset file name", metavar="character"),
  make_option(c("-c", "--covariate"), type="character", default=NULL, 
              help="Covariate file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="./", 
              help="output directory path [default=%default]", metavar="character"),
  make_option(c("-O", "--output_addition_columns"), type="character", default=NULL, 
              help="output additioal columns than DEseq2 default", metavar="character"),
  make_option(c("-C", "--comparison"), type="character", default=NULL, 
              help="Comparison in format of Variable:X:Y", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default='hg19', 
              help="Genome assembly [default=%default]", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", default='~/projects/circRNA/data/Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds', 
              help="circRNA annotation [default=%default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
input_expression_filename=opt$input
input_covariance_filename=opt$covariate
output_dir=opt$out
output_additonal_columns=opt$output_addition_columns
index=opt$genome
annotation_path=opt$annotation
comparison=opt$comparison

# debug
# setwd("~/neurogen/rnaseq_PD/results/merged"); input_expression_filename="genes.htseqcount.cufflinks.TCPY.uniq.xls"; input_covariance_filename="TCPY.covariates.txt"; output_dir="TCPY"; index="hg19"
# setwd("~/neurogen/rnaseq_Rot/results/merged"); input_expression_filename="genes.htseqcount.cufflinks.mm9.uniq.xls"; input_covariance_filename="covariate.mm9.txt"; output_dir="DE.mm9"; index="mm9"; output_additonal_columns="mi"
# setwd("~/projects/circRNA/results/DE_SNDA"); input_expression_filename="../data/Merge_circexplorer_BC106.filtered.enriched.rawcount.rds"; input_covariance_filename="Table.PD.SNDA.pathology.covariates.xls"; output_dir="DE_SNDA"; index="hg19"; comparison="CONDITION:PD:HC"
# setwd("~/projects/circRNA/results/DE_SNDA"); input_expression_filename="../data/Merge_circexplorer_BC106.filtered.enriched.rawcount.rds"; input_covariance_filename="Table.PD.SNDA.pathology.covariates.xls"; output_dir="DE_SNDA"; index="hg19"; comparison="PD.pathology.group:early:no"
# setwd("~/projects/circRNA/results"); input_expression_filename="../data/Merge_circexplorer_CSF87.filtered.enriched.rawcount.rds"; input_covariance_filename="Table.PD.CSF.pathology.covariates.xls"; output_dir="DE_CSF"; index="hg19"; comparison="CONDITION:PD:HC"; annotation_path="~/projects/circRNA/data/Merge_circexplorer_CSF87.filtered.annotation.bed14.rds"
# setwd("~/projects/circRNA/results/DE_TCPY"); input_expression_filename="../data/Merge_circexplorer_BC106.filtered.enriched.rawcount.rds"; input_covariance_filename="Table.AD.TCPY.pathology.covariates.xls"; output_dir="DE_TCPY"; index="hg19"; comparison="CONDITION:AD:HC"; annotation_path='~/projects/circRNA/data/Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds'

if(is.null(comparison)){
  print_help(opt_parser)
  stop("Comparison for DE analysis must be supplied, e.g. CONDITION:PD:HC", call.=FALSE)
}

if(is.null(input_expression_filename) || is.null(input_covariance_filename)){
  print_help(opt_parser)
  stop("Both input expression matrix and covariate file must be supplied", call.=FALSE)
}

str_comparison = unlist(strsplit(comparison, "[:]"))
variable=str_comparison[1];
variable_ALT=str_comparison[2];
variable_REF=str_comparison[3];

# check input
if(!file.exists(input_expression_filename)) {stop(paste(input_expression_filename, "doesn't exist. Exit!"), call.=FALSE);}
if(!file.exists(input_covariance_filename)) {stop(paste(input_covariance_filename, "doesn't exist. Exit!"), call.=FALSE);}
genome_name=switch(index, "hg19" = "Homo_sapiens", "mm9" = "Mus_musculus", "rn6" = "Rattus_norvegicus")

# Create folder if the directory doesn't exist
dir.create(file.path(output_dir,'report/figures'), recursive =T, showWarnings = FALSE)

pwd=getwd()
setwd(output_dir)

# install packages
library('tidyverse',quietly=T, warn.conflicts=F, logical.return=T) || install.packages('tidyverse', repo='http://cran.revolutionanalytics.com');
library('RCurl',quietly=T, warn.conflicts=F, logical.return=T) || install.packages('RCurl', repo='http://cran.revolutionanalytics.com');
library('hexbin',quietly=T, warn.conflicts=F, logical.return=T) || install.packages('hexbin', repo='http://cran.revolutionanalytics.com');
library('pheatmap',quietly=T, warn.conflicts=F, logical.return=T) || install.packages('pheatmap', repo='http://cran.revolutionanalytics.com');
library('RColorBrewer',quietly=T, warn.conflicts=F, logical.return=T) || install.packages('RColorBrewer', repo='http://cran.revolutionanalytics.com');
library('hwriter',quietly=T, warn.conflicts=F, logical.return=T) || install.packages('hwriter', repo='http://cran.revolutionanalytics.com');
source("https://bioconductor.org/biocLite.R"); 
library('vsn',quietly=T, warn.conflicts=F, logical.return=T) || biocLite('vsn');
library('DESeq2',quietly=T, warn.conflicts=F, logical.return=T) || biocLite('DESeq2');
#library('ReportingTools',quietly=T, warn.conflicts=F, logical.return=T) || biocLite('ReportingTools');
library('BiocParallel',quietly=T, warn.conflicts=F, logical.return=T) || biocLite('BiocParallel');
library('limma',quietly=T, warn.conflicts=F, logical.return=T) || biocLite('limma');


###########################################
message("#step1: load data...")
###########################################
if(file.exists(file.path("DESeq2.RData"))) load(file.path("DESeq2.RData")) else {
  
  # debug
  # genome_name="Rattus_norvegicus"; index='rn6'; input_expression_filename='genes.htseqcount.cufflinks.allSamples.uniq.xls'; input_covariance_filename='covariate.txt'
  
  #annotation (Note: this is downloaded from biomart, so EnsID format is ENSGxxxxx, no version number tailing)
  genes_annotation = read.table(file.path("~/neurogen/referenceGenome",genome_name,"UCSC",index,"Annotation/Genes/annotation.genes.bed6+3"), sep="\t", quote="", header = F, stringsAsFactors = F, col.names = c("chr","start","end","geneID","score","strand","geneName","geneType","geneDescription"));

  ## for circRNAs
  if(tolower(tools::file_ext(annotation_path)) == "rds") {
    circRNA_annotation=readRDS(file.path(annotation_path))
  }else circRNA_annotation=read.delim(file.path(input_expression_filename), header = T, row.names = 1,check.names =F)
  circRNA_annotation$geneDescription = genes_annotation$geneDescription[match(sub("\\..*","", circRNA_annotation$geneID), genes_annotation$geneID)]
  circRNA_annotation$geneDescription = sub(" \\[Source:.*","", circRNA_annotation$geneDescription)
  genes_annotation = circRNA_annotation
  
  # raw reads count
  if(tolower(tools::file_ext(input_expression_filename)) == "rds") {
    cts=readRDS(file.path(pwd,input_expression_filename))
  }else cts=read.delim(file.path(pwd,input_expression_filename), row.names = 1,check.names =F)
  # remove those non-geneID rows, e.g. __no_feature (pre-mRNA reads) and __ambiguous (see http://htseq.readthedocs.io/en/master/count.html )
  dim(cts); cts=cts[grep("^__", rownames(cts), invert = T),]; dim(cts);
  
  # covariance table
  covarianceTable = read.table(file.path(pwd,input_covariance_filename), sep="\t", header = T,check.names =F, stringsAsFactors = F)
  head(covarianceTable)
  rownames(covarianceTable) = covarianceTable$SAMPLE_ID;
  
  # varialbe must be one of the columns in covariance table
  if(!(variable %in% names(covarianceTable))) stop("Variable of interest is not in your covariate file! Please double check it", call. =F);
  
  # subset and re-order
  all(rownames(covarianceTable) %in% colnames(cts))
  dim(cts); cts = cts[, intersect(colnames(cts), rownames(covarianceTable))]; dim(cts)
  dim(covarianceTable); covarianceTable = covarianceTable[colnames(cts), ]; dim(covarianceTable);
  covarianceTable[] <- lapply(covarianceTable, function(x) if(is.factor(x)) factor(x) else x) # drop levels  after subsetting
  all(rownames(covarianceTable) == colnames(cts))
  
  # check the structure of corariance table, to see if it's necessary to factorize some of the columns
  str(covarianceTable)
  
  # factorize batch, age, sex, RIN, PMI, condition, if exist
  if("CONDITION" %in% colnames(covarianceTable)) covarianceTable$CONDITION=factor(covarianceTable$CONDITION)
  if(variable %in% colnames(covarianceTable) && length(str_comparison)==3) covarianceTable[[variable]]=factor(covarianceTable[[variable]])
  if("CELLTYPE" %in% colnames(covarianceTable)) covarianceTable$CELLTYPE=factor(covarianceTable$CELLTYPE)
  if("BATCH" %in% colnames(covarianceTable)) covarianceTable$BATCH=factor(covarianceTable$BATCH)
  if("SEX" %in% colnames(covarianceTable)) covarianceTable$SEX=factor(covarianceTable$SEX, levels = c("F","M"))
  if("AGE" %in% colnames(covarianceTable)) {x=cut(covarianceTable$AGE, breaks = c(20,80,110)); levels(x)=gsub("\\((.*),(.*)\\]","\\1_\\2",levels(x)); covarianceTable$AGE=x;}
  if("RIN" %in% colnames(covarianceTable)) {x=cut(covarianceTable$RIN, breaks = c(5,8,10)); levels(x)=gsub("\\((.*),(.*)\\]","\\1_\\2",levels(x)); covarianceTable$RIN=x;}
  if("PMI" %in% colnames(covarianceTable)) {x=cut(covarianceTable$PMI, breaks = c(0,4,100)); levels(x)=gsub("\\((.*),(.*)\\]","\\1_\\2",levels(x)); covarianceTable$PMI=x}
  # TODO: add CDR, Braak etc.
  
  str(covarianceTable)
  
  ###########################################
  message("#step2: load data to DEseq")
  ###########################################
  
  # Note: With no arguments to results, the results will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the first level.
  # Ref: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
  
  # making formula
  #fmla <- names(covarianceTable)[!names(covarianceTable) %in% c("SAMPLE_ID", "SUBJECT_ID")]
  #fmla <- as.formula(paste(" ~ ", paste(fmla, collapse= "+")))
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = covarianceTable,
                                design= ~ 1)
  
  ## pre-filtering
  #keep <- rowSums(counts(dds)) >= 10
  keep <- rowSums(counts(dds) >= 5) >= 4  # at least 4 samples with more than 5 reads
  dim(dds); dds <- dds[keep,]; dim(dds);
  #head(sort(rowMeans(counts(dds)), decreasing = T))
  # head(counts(dds)[order(rowMeans(counts(dds)), decreasing =T),])
  
  ###########################################
  message("#step3: QA of the data [optional]")
  ###########################################
  # Note: This part is not necessary for DEseq, but important for data QA
  
  ##--------------------------------------
  ## 3.1: compare different vairance stablization methods
  ##--------------------------------------
  
  ntd <- normTransform(dds) # log2(x+1)
  vsd <- varianceStabilizingTransformation(dds, blind=T) # Note: blind to the design, equal to design = ~ 1
  
  # # using limma to remove covariates, it returns adjusted values in log2 scale
  vsd_adjusted_log2 <- removeBatchEffect(assay(vsd), batch=vsd$BATCH, batch2=vsd$SEX, covariates = colData(vsd)[,c('AGE','PMI','RIN')])
  
  pdf("diagnosis.pdf")
  par(mar=c(10,5,3,3));boxplot(assay(ntd), las=2, cex.axis=.7, ylab="log2(x+1) transform")
  msd <- meanSdPlot(counts(dds), ranks = FALSE); msd$gg + ggtitle("no transformation")
  msd <- meanSdPlot(assay(ntd), ranks = FALSE); msd$gg + ggtitle("log2(x+1) transform")
  msd <- meanSdPlot(assay(vsd), ranks = FALSE); msd$gg + ggtitle("VST")
  msd <- meanSdPlot(vsd_adjusted_log2, ranks = FALSE); msd$gg + ggtitle("vsd_adjusted_log2")
  dev.off()
  
  ## PCA
  if(length(str_comparison)==3) {
    pcaData <- plotPCA(vsd, intgroup = c(variable), returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData, aes_string(x = "PC1", y = "PC2", color = variable)) +
      geom_point(size =3) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      ggtitle(paste0(output_dir,": PCA")) + 
      coord_fixed()  
    ggsave(file.path(paste("PCA", variable ,"pdf", sep = ".")))
  }
  
  ##--------------------------------------
  ## 3.2: save normalized reads count
  ##--------------------------------------
  
  ## save the raw reads count 
  write.table(counts(dds), 
              file=paste0(basename(input_expression_filename), ".filtered.raw.xls"), 
              sep="\t", quote = F, 
              col.names = NA, row.names = TRUE)
  ## save the variance-stabilized data
  write.table(assay(vsd), 
              file=paste0(basename(input_expression_filename), ".filtered.vsd.xls"), 
              sep="\t", quote = F, 
              col.names = NA, row.names = TRUE)
  ## save the vsd adjusted
  write.table(vsd_adjusted_log2,
              file=paste0(basename(input_expression_filename), ".filtered.vsd_adjusted_log2.xls"),
              sep="\t", quote = F,
              col.names = NA, row.names = TRUE)
  

  save(genes_annotation, covarianceTable, dds, vsd, vsd_adjusted_log2, file="DESeq2.RData")

}
# script to generate html report page
makeNewImages <- function(df,variable,...){
  imagename <- c()
  tablename <- c()
  for (i in 1:nrow(df)){
    ensId <- rownames(df)[i]
    symbol <- df$symbol[i]    
    imagename[i] <- paste('plot', ensId, symbol, variable, 'pdf', sep = ".")
    tablename[i] <- paste('plot', ensId, symbol, variable, 'txt', sep = ".")
    
    d <- data.frame(samples=colnames(assay(ntd)), 
                    expression_ntd=assay(ntd)[ensId,], 
                    expression_raw=assay(dds)[ensId,], 
                    condition=colData(dds)[[variable]])
  
    if(!file.exists(file.path('report/figures',tablename[i]))) {
      write.table(d,file.path('report/figures',tablename[i]),sep="\t", quote =F, row.names=F, col.names = T)
    }
    if(!file.exists(file.path('report/figures',imagename[i]))) 
    {
      N=length(levels(colData(dds)[[variable]]))
      p=ggplot(d, aes_string(x=variable, y="expression_ntd")) + 
        geom_boxplot(position=position_dodge(.8), width=.5, outlier.shape = NA) +
        geom_jitter(size=1.5, position = position_jitter(width=.15)) +
        theme_bw() +
        xlab(variable) + ylab("log2(counts+1)") + ggtitle(symbol,subtitle =ensId)
      #png(file.path('report/figures', imagename[i]), height = 250, width = 600)
      pdf(file.path('report/figures', imagename[i]), height = 4, width = 3*N)
      print(p)
      dev.off()
    }
  }
  ## Using the following code to show thumb figures. It's slow to display if many
  # df$Boxplot <- hwriteImage(paste0('figures/', imagename), 
  #                           link=paste0('figures/', imagename), 
  #                           table=FALSE, width=100)
  df$Boxplot <- hwrite('boxplot', link = paste0('figures/', imagename), table=F)
  df$Rawdata <- hwrite("data", link = paste0('figures/', tablename), table=F)
  df$symbol <- hwrite(as.character(df$symbol), 
                      link = paste0("http://useast.ensembl.org/",genome_name,"/Gene/Summary?db=core;g=",as.character(rownames(df))), 
                      table=F)
  return(df)
}

###########################################
message("#step4: Run DE")
###########################################
if(length(str_comparison)==1){  # e.g. --comparison="PMI"
  message(paste("# comparison will be performed on a quantitative / continuous variable:", variable));  

  # subsetting
  dds=dds[, !is.na(dds[[variable]])]
  covariances=c();
  for(i in c('AGE', 'RIN', 'SEX', 'PMI', 'BATCH', 'CELLTYPE')) if(i %in% colnames(covarianceTable) & !(i %in% covariances)) {dds[[i]] <- droplevels(dds[[i]]); if(length(levels(dds[[i]]))>1) covariances=c(covariances,i);}
  keep <- rowSums(counts(dds) >= 5) >= 4  # at least 4 samples with more than 5 reads
  dim(dds); dds <- dds[keep,]; dim(dds);
  
  #fmla <- names(covarianceTable)[!names(covarianceTable) %in% c("SAMPLE_ID", "SUBJECT_ID")]
  design(dds) <- as.formula(paste(" ~ ", paste(c(covariances, variable), collapse= " + ")))
  
  # In versions >=1.16, the default is set to FALSE, and shrunken LFCs are obtained afterwards using lfcShrink.
  dds <- DESeq(dds, betaPrior=T, parallel=TRUE, BPPARAM=MulticoreParam(4))  
  resultsNames(dds)
  
  # If the variable is continuous then the results can be extracted using the name argument to results, where the name is one of elements returned by resultsNames(dds).
  res <- results(dds, name = variable,   
                 alpha = 0.1, 
                 cooksCutoff=FALSE, # regardless of outlier. See https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
                 parallel=TRUE, BPPARAM=MulticoreParam(4))
  
  com_name = comparison
  
  summary(res)
  head(res); dim(res)
  # decimal value of Fold-change
  res$FoldChange <- 2**res$log2FoldChange
  
  # add annotation
  res <- cbind(res, genes_annotation[match(sub("\\..*","",row.names(res)), genes_annotation$ID), c('circType','geneName','geneType','geneDescription')])
  
  # add additional columns in the output
  if(!is.null(output_additonal_columns) && grepl("i", output_additonal_columns)){ # individual raw expression values of each sample 
    individual <- counts(dds,normalized=FALSE)
    colnames(individual) = paste0("ind_raw.", colnames(individual))
    res = cbind(res, individual)
  }
  if(!is.null(output_additonal_columns) && grepl("I", output_additonal_columns)){ # individual normalized expression values of each sample
    individual <- counts(dds,normalized=TRUE)
    colnames(individual) = paste0("ind_norm.", colnames(individual))
    res = cbind(res, individual)
  }
  res <- res[order(res$pvalue),]
  head(res); dim(res)
  
  # remove NA
  res <- na.omit(res)
  dim(res)
  
  # write to xls
  write.table(as.data.frame(res), 
              file=file.path(paste("DEresult",output_dir,com_name ,"xls",sep=".")), 
              sep="\t", quote =F, na="", row.names=T, col.names = NA)
  
  DE = subset(res, !is.na(pvalue) & pvalue<=0.05 & abs(log2FoldChange)>=1)
  NDE = subset(res, pvalue>0.05 | abs(log2FoldChange)<1)
  # ## Note: 20% of non-significant(NS) genes are randomly selected in order to increase the performance of the generating graph
  if(nrow(NDE)>5000){
    n_NS=nrow(NDE)
    NDE=NDE[sample(n_NS,round(n_NS * .20)),]
    dim(res); res=DESeqResults(rbind(NDE, DE)); dim(res);
  }
  
  ## MAKING PLOTS
  pdf(file.path(paste("DEresult",output_dir, com_name ,"pdf", sep = ".")), paper = 'USr')
  ## ==============================
  # MA plot
  ## ==============================
  DESeq2::plotMA(DESeqResults(res), alpha = 0.05, colNonSig = "gray", main=paste0(output_dir,": ",com_name))
  
  RES=as.data.frame(res)
  DE =as.data.frame(DE)

  ## ==============================
  # vocano plot
  ## ==============================
  with(RES, plot(log2FoldChange, -log10(pvalue), 
                 pch=20, cex=0.5, main=paste0(output_dir,": ",com_name), col='gray',
                 xlab=bquote(~Log[2]~fold~change), 
                 ylab=bquote(~-log[10]~pvalue)))
  if(nrow(DE)>0) {
    with(DE, points(log2FoldChange, -log10(pvalue), pch=20, col="red", cex=1))
    with(DE, text(log2FoldChange, -log10(pvalue), labels=geneName, cex=0.5, pos=1, offset=0.2))
  }
  abline(v=c(0,-1, 1), lty=c(3,4,4), lwd=c(1,2,2))
  abline(h=-log10(0.05), col="black", lty=4, lwd=2.0)
  
  
  ## ==============================
  message("# heatmap for top 10 DE genes")
  ## ==============================
  RES=RES[order(RES$pvalue),] # sort by padj in increasing order
  DE10=rbind(head(subset(RES, log2FoldChange<0),10),head(subset(RES, log2FoldChange>0),10)) # top 10 down-regulated and top 10 up-regulated
  topDE=assay(vsd)[rownames(DE10),rownames(colData(dds))]
  rownames(topDE) = paste(DE10$geneName, rownames(topDE), sep=".")
  annotation_row = dplyr::select(DE10, circType, geneType, log2FoldChange, geneName) %>% 
    mutate(updown=ifelse(log2FoldChange>0,"up","down")) %>% select(updown, circType, geneType)
  rownames(annotation_row) = rownames(topDE) # DE10$geneName
  annotation_col = dplyr::select(as.data.frame(colData(dds)), variable)
  ann_colors = list(
    updown = c(up = "red", down = "blue"),
    circType = c(circRNA = "red", ciRNA = "orange"),
    geneType = c(protein_coding = "darkblue", lincRNA = "orange", Mt_rRNA='pink', 
                 antisense='yellow', pseudogene='gray', processed_transcript='lightgray')
  )
  #tmp=c("green","red"); names(tmp)=c(variable_REF, variable_ALT);
  #ann_colors=c(list(tmp), ann_colors); names(ann_colors)[1]=variable
  
  # sort column of topDE by variable
  topDE=topDE[,rownames(annotation_col[order(annotation_col[[variable]]),, drop=F])]
  
  ## Scale/center each genes (by rows)
  topDE=t(scale(t(as.matrix(topDE))))
  ## trim max and min
  MAX=2; MIN=-2; topDE[topDE>MAX]=MAX;topDE[topDE<MIN]=MIN
  
  ## add noise to avoid SD=0 cases
  topDE=topDE + matrix(rnorm(prod(dim(topDE)),0,0.01),nrow = nrow(topDE),ncol = ncol(topDE))
  
  par(cex=0.5, mar=c(5, 8, 4, 1))
  pheatmap(topDE,
           fontsize = 5,
           main =paste0(output_dir,": heatmap for top 10 DE genes"),
           #width=8, height=5,filename=paste0("DEresult.padj_05.", com_name ,".heatmap.pdf"),
           border_color = NA,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           annotation_row = annotation_row,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           drop_levels = TRUE,
           scale = "none", 
           clustering_method = 'ward.D', 
           cluster_rows = TRUE,
           clustering_distance_rows = "correlation",
           cutree_rows = 2,
           #cutree_cols=2,
           cluster_cols = FALSE,
           clustering_distance_cols = "correlation")
  
  ## ==============================
  message("# heatmap for top 20 DE genes")
  ## ==============================
  RES=RES[order(RES$pvalue),] # sort by padj in increasing order
  DE20=rbind(head(subset(RES, log2FoldChange<0),20),head(subset(RES, log2FoldChange>0),20)) # top 10 down-regulated and top 10 up-regulated
  topDE=assay(vsd)[rownames(DE20),rownames(colData(dds))]; 
  rownames(topDE) = paste(DE10$geneName, rownames(topDE), sep=".")
  annotation_row = dplyr::select(DE20, circType,geneType, log2FoldChange, geneName) %>% 
    mutate(updown=ifelse(log2FoldChange>0,"up","down")) %>% select(updown, circType,geneType)
  rownames(annotation_row) = rownames(topDE) # DE20$geneName
  annotation_col = dplyr::select(as.data.frame(colData(dds)), variable)
  ann_colors = list(
    updown = c(up = "red", down = "blue"),
    circType = c(circRNA = "red", ciRNA = "orange"),
    geneType = c(protein_coding = "darkblue", lincRNA = "orange", Mt_rRNA='pink', 
                 antisense='yellow', pseudogene='gray', processed_transcript='lightgray')
  )
  #tmp=c("green","red"); names(tmp)=c(variable_REF, variable_ALT);
  #ann_colors=c(list(tmp), ann_colors); names(ann_colors)[1]=variable
  
  topDE=topDE[,rownames(annotation_col[order(annotation_col[[variable]]),, drop=F])]
  
  ## Scale/center each genes (by rows)
  topDE=t(scale(t(as.matrix(topDE))))
  ## trim max and min
  MAX=2; MIN=-2; topDE[topDE>MAX]=MAX;topDE[topDE<MIN]=MIN
  
  ## add noise to avoid SD=0 cases
  topDE=topDE + matrix(rnorm(prod(dim(topDE)),0,0.01),nrow = nrow(topDE),ncol = ncol(topDE))
  
  par(cex=0.5, mar=c(5, 8, 4, 1))
  pheatmap(topDE,
           main =paste0(output_dir,": heatmap for top 20 DE genes"),
           fontsize = 5,
           border_color = NA,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           annotation_row = annotation_row,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           drop_levels = TRUE,
           scale = "none", 
           clustering_method = 'ward.D', 
           cluster_rows = TRUE,
           clustering_distance_rows = "correlation",
           cutree_rows = 2,
           #cutree_cols=2,
           cluster_cols = FALSE,
           clustering_distance_cols = "correlation")
  
  dev.off() 
  
  
}else if(length(str_comparison)==3){  # e.g. --comparison="DX:PD:HC"  TODO: can be multi comparisons e.g. DX:PD:HC;DX:ILB:HC

  message(paste("# comparison will be performed on a factor variable:", variable, "between level", variable_ALT, "and", variable_REF));  
  
  ## subsetting
  dim(dds)
  dds=dds[, dds[[variable]] %in% c(variable_ALT, variable_REF)]
  covariances=c();
  for(i in c('AGE', 'RIN', 'SEX', 'PMI', 'BATCH', 'CELLTYPE', variable)) if(i %in% colnames(covarianceTable) & !(i %in% covariances)) {dds[[i]] <- droplevels(dds[[i]]); if(length(levels(dds[[i]]))>1) covariances=c(covariances,i);}
  keep <- rowSums(counts(dds) >= 5) >= 4  # at least 4 samples with more than 5 reads
  dim(dds); dds <- dds[keep,]; dim(dds);

  #fmla <- names(covarianceTable)[!names(covarianceTable) %in% c("SAMPLE_ID", "SUBJECT_ID")]
  design(dds) <- as.formula(paste(" ~ ", paste(covariances, collapse= " + ")))
  
  # In versions >=1.16, the default is set to FALSE, and shrunken LFCs are obtained afterwards using lfcShrink.
  dds <- DESeq(dds, betaPrior=T, parallel=TRUE, BPPARAM=MulticoreParam(4))  
  resultsNames(dds)
  
  com_name= paste(variable, variable_ALT, "vs", variable_REF, sep="_")
  message(paste("processing comparison:",com_name))
  res <- results(dds, contrast = c(variable, variable_ALT, variable_REF),   # contract format: factor name, numerator in the fold change, denominator in the fold change
                 alpha = 0.1, 
                 cooksCutoff=FALSE, # regardless of outlier. See https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
                 parallel=TRUE, BPPARAM=MulticoreParam(4))
  ## Shrink log2 fold changes [required for version >=1.16]
  #res2 <- lfcShrink(dds, contrast = c(x[1], x[2], x[3]), res=res, type='normal', parallel=TRUE, BPPARAM=MulticoreParam(4))
  ## You can get the shrunken LFC either with lfcShrink like above or with betaPrior=TRUE. It will be the same shrunken LFC and the same as previously calculated in DESeq2. The difference is that betaPrior=TRUE will give you a p-value for the shrunken LFC, while lfcShrink (at the moment) is only giving you the LFC, and is keeping the p-value for the test of the MLE LFC. 
  ## see https://support.bioconductor.org/p/95695/ and https://support.bioconductor.org/p/98833/#98843
  
  summary(res)
  head(res); dim(res)
  # decimal value of Fold-change
  res$FoldChange <- 2**res$log2FoldChange
  
  # add annotation
  res <- cbind(res, genes_annotation[match(sub("\\..*","",row.names(res)), genes_annotation$ID), c('circType','geneName','geneType','geneDescription')])
  
  # add additional columns in the output
  if(!is.null(output_additonal_columns) && grepl("m", output_additonal_columns)){ # mean of raw expression values for each group, 
    baseMeanPerLvl <- sapply( levels(dds[[variable]]), function(lvl) rowMeans( counts(dds,normalized=FALSE)[,dds[[variable]] == lvl] ) )
    colnames(baseMeanPerLvl) = paste0("baseMean_raw.", colnames(baseMeanPerLvl))
    res = cbind(res, baseMeanPerLvl)
  }
  if(!is.null(output_additonal_columns) && grepl("M", output_additonal_columns)){ # mean of normalized expression values for each group, 
    baseMeanPerLvl <- sapply( levels(dds[[variable]]), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds[[variable]] == lvl] ) )
    colnames(baseMeanPerLvl) = paste0("baseMean_norm.", colnames(baseMeanPerLvl))
    res = cbind(res, baseMeanPerLvl)
  }
  if(!is.null(output_additonal_columns) && grepl("i", output_additonal_columns)){ # individual raw expression values of each sample 
    individual <- counts(dds,normalized=FALSE)
    colnames(individual) = paste0("ind_raw.", colnames(individual))
    res = cbind(res, individual)
  }
  if(!is.null(output_additonal_columns) && grepl("I", output_additonal_columns)){ # individual normalized expression values of each sample
    individual <- counts(dds,normalized=TRUE)
    colnames(individual) = paste0("ind_norm.", colnames(individual))
    res = cbind(res, individual)
  }
  res <- res[order(res$pvalue),]
  head(res); dim(res)
  
  # remove NA
  res <- na.omit(res)
  dim(res)
  
  # write to xls
  write.table(as.data.frame(res), 
              file=file.path(paste("DEresult",output_dir,com_name ,"xls",sep=".")), 
              sep="\t", quote =F, na="", row.names=T, col.names = NA)
  
  DE = subset(res, !is.na(pvalue) & pvalue<=0.05 & abs(log2FoldChange)>=1)
  NDE = subset(res, pvalue>0.05 | abs(log2FoldChange)<1)
  # ## Note: 20% of non-significant(NS) genes are randomly selected in order to increase the performance of the generating graph
  if(nrow(NDE)>5000){
    n_NS=nrow(NDE)
    NDE=NDE[sample(n_NS,round(n_NS * .20)),]
    dim(res); res=DESeqResults(rbind(NDE, DE)); dim(res);
  }
  
  ## MAKING PLOTS
  pdf(file.path(paste("DEresult",output_dir, com_name ,"pdf", sep = ".")), paper = 'USr')
  ## ==============================
  # MA plot
  ## ==============================
  DESeq2::plotMA(DESeqResults(res), alpha = 0.05, colNonSig = "gray", main=paste0(output_dir,": ",com_name))
  
  RES=as.data.frame(res)
  DE =as.data.frame(DE)
  # scatter plot
  if(!is.null(output_additonal_columns) && grepl("m", output_additonal_columns)){ # mean of raw expression values for each group, 
    par(pty="s"); 
    REF=paste0("baseMean_raw.", variable_REF); ALT=paste0("baseMean_raw.", variable_ALT)
    plot(jitter(RES[[REF]])+1, jitter(RES[[ALT]])+1, 
         log='xy', pch=19, cex=1, cex.axis=2, cex.lab=2,
         xlab=REF,ylab=ALT,
         xlim=range(RES[[REF]], RES[[ALT]])+1, 
         ylim=range(RES[[REF]], RES[[ALT]])+1, 
         col=ifelse(RES$circType=="ciRNA",'orange','red'))
    RES=RES[RES[[REF]]>=5 | RES[[ALT]]>=5,]
    text(x=1+RES[[REF]], y=1+RES[[ALT]], labels = RES$geneName, cex=.5, adj=c(0.5,0), pos=3);
    abline(a=0,b=1, lty=2, col='black')
  }
  
  ## ==============================
  # vocano plot
  ## ==============================
  with(RES, plot(log2FoldChange, -log10(pvalue), 
                  pch=20, cex=0.5, main=paste0(output_dir,": ",com_name), col='gray',
                  xlab=bquote(~Log[2]~fold~change), 
                  ylab=bquote(~-log[10]~pvalue)))
  if(nrow(DE)>0) {
    with(DE, points(log2FoldChange, -log10(pvalue), pch=20, col="red", cex=1))
    with(DE, text(log2FoldChange, -log10(pvalue), labels=geneName, cex=0.5, pos=1, offset=0.2))
  }
  abline(v=c(0,-1, 1), lty=c(3,4,4), lwd=c(1,2,2))
  abline(h=-log10(0.05), col="black", lty=4, lwd=2.0)
  
  
  ## ==============================
  message("# heatmap for top 10 DE genes")
  ## ==============================
  #DE=DE[order(-DE$baseMean),] # sort by baseMean in decreasing order
  #DE=DE[order(-abs(DE$log2FoldChange)),] # sort by abs(log2FoldChange) in decreasing order
  RES=RES[order(RES$pvalue),] # sort by padj in increasing order
  DE10=rbind(head(subset(RES, log2FoldChange<0),10),head(subset(RES, log2FoldChange>0),10)) # top 10 down-regulated and top 10 up-regulated
  topDE=assay(vsd)[rownames(DE10),rownames(colData(dds))]
  rownames(topDE) = paste(DE10$geneName, rownames(topDE), sep=".")
  annotation_row = dplyr::select(DE10, circType, geneType, log2FoldChange, geneName) %>% 
    mutate(updown=ifelse(log2FoldChange>0,"up","down")) %>% select(updown, circType, geneType)
  rownames(annotation_row) = rownames(topDE) # DE10$geneName
  annotation_col = dplyr::select(as.data.frame(colData(dds)), variable)
  ann_colors = list(
    updown = c(up = "red", down = "blue"),
    circType = c(circRNA = "red", ciRNA = "orange"),
    geneType = c(protein_coding = "darkblue", lincRNA = "orange", Mt_rRNA='pink', 
                 antisense='yellow', pseudogene='gray', processed_transcript='lightgray')
  )
  tmp=c("green","red"); names(tmp)=c(variable_REF, variable_ALT);
  ann_colors=c(list(tmp), ann_colors); names(ann_colors)[1]=variable
  
  ## Scale/center each genes (by rows)
  topDE=t(scale(t(as.matrix(topDE))))
  ## trim max and min
  MAX=2; MIN=-2; topDE[topDE>MAX]=MAX;topDE[topDE<MIN]=MIN
  
  ## add noise to avoid SD=0 cases
  topDE=topDE + matrix(rnorm(prod(dim(topDE)),0,0.01),nrow = nrow(topDE),ncol = ncol(topDE))
  
  par(cex=0.5, mar=c(5, 8, 4, 1))
  pheatmap(topDE,
           fontsize = 5,
           main =paste0(output_dir,": heatmap for top 10 DE genes"),
           #width=8, height=5,filename=paste0("DEresult.padj_05.", com_name ,".heatmap.pdf"),
           border_color = NA,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           annotation_row = annotation_row,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           drop_levels = TRUE,
           scale = "none", 
           clustering_method = 'ward.D', 
           cluster_rows = TRUE,
           clustering_distance_rows = "correlation",
           cutree_rows = 2,cutree_cols=2,
           cluster_cols = TRUE,
           clustering_distance_cols = "correlation")
  
  ## ==============================
  message("# heatmap for top 20 DE genes")
  ## ==============================
  RES=RES[order(RES$pvalue),] # sort by padj in increasing order
  DE20=rbind(head(subset(RES, log2FoldChange<0),20),head(subset(RES, log2FoldChange>0),20)) # top 10 down-regulated and top 10 up-regulated
  topDE=assay(vsd)[rownames(DE20),rownames(colData(dds))]; 
  rownames(topDE) = paste(DE10$geneName, rownames(topDE), sep=".")
  annotation_row = dplyr::select(DE20, circType,geneType, log2FoldChange, geneName) %>% 
    mutate(updown=ifelse(log2FoldChange>0,"up","down")) %>% select(updown, circType,geneType)
  rownames(annotation_row) = rownames(topDE) # DE20$geneName
  annotation_col = dplyr::select(as.data.frame(colData(dds)), variable)
  ann_colors = list(
    updown = c(up = "red", down = "blue"),
    circType = c(circRNA = "red", ciRNA = "orange"),
    geneType = c(protein_coding = "darkblue", lincRNA = "orange", Mt_rRNA='pink', 
                 antisense='yellow', pseudogene='gray', processed_transcript='lightgray')
  )
  tmp=c("green","red"); names(tmp)=c(variable_REF, variable_ALT);
  ann_colors=c(list(tmp), ann_colors); names(ann_colors)[1]=variable
  
  ## Scale/center each genes (by rows)
  topDE=t(scale(t(as.matrix(topDE))))
  ## trim max and min
  MAX=2; MIN=-2; topDE[topDE>MAX]=MAX;topDE[topDE<MIN]=MIN
  
  ## add noise to avoid SD=0 cases
  topDE=topDE + matrix(rnorm(prod(dim(topDE)),0,0.01),nrow = nrow(topDE),ncol = ncol(topDE))
  
  par(cex=0.5, mar=c(5, 8, 4, 1))
  pheatmap(topDE,
           main =paste0(output_dir,": heatmap for top 20 DE genes"),
           fontsize = 5,
           border_color = NA,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           annotation_row = annotation_row,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           drop_levels = TRUE,
           scale = "none", 
           clustering_method = 'ward.D', 
           cluster_rows = TRUE,
           clustering_distance_rows = "correlation",
           cutree_rows = 2,cutree_cols=2,
           cluster_cols = TRUE,
           clustering_distance_cols = "correlation")
  
  dev.off() 
  
  # message("## Exporting results to HTML")
  # htmlRep <- HTMLReport(shortName=com_name, title=com_name,
  #                       reportDirectory="report")
  # publish(makeNewImages(head(as.data.frame(subset(res, abs(log2FoldChange)>=1 & padj<0.05, select=!grepl("^ind_", colnames(res)))),1000), variable), htmlRep)
  # finish(htmlRep)
  
} else stop("Unrecognized comparison format. Please use style like --comparison=DX:PD:HC", call. = F);

# ###########################################
# message("#step5: scp html result to web server")
# ###########################################
# system(paste0("rsync -a --chmod=u+rwx,g+rwx,o+rwx report/* xd010@panda.dipr.partners.org:~/public_html/DE_reports/",output_dir))
# message(paste0("Visit the HTML result at: http://panda.partners.org/~xd010/DE_reports/",output_dir))

message("Done!")