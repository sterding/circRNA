###########################################
# A general framework in R for running differntial expression analysis using DESeq2
# author: Xianjun Dong
# email: xdong@rics.bwh.harvard.edu
# date: 9/3/2019
# version: 2.0
# Usage: Rscript --vanilla ../src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -C CONDITION:PD:HC
# Requirement:
# 1. The headers in covariate file should be in CAPITAL format for SAMPLE_ID (required), SUBJECT_ID (if any), CELLTYPE (optional), CONDITION (required)
# Log:
# 9/10/2019: Update to R 3.6.1 and bioC 3.9 
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
  make_option(c("-k", "--collapse"), action="store_true", default=FALSE,  
              help="Collapse circRNAs to gene level"),
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
COLLAPSE=opt$collapse

# debug
# setwd("~/projects/circRNA/results/"); input_expression_filename="../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds"; input_covariance_filename="Table.PD.SNDA.pathology.covariates.xls "; output_additonal_columns='mi'; output_dir="DE_SNDA"; index="hg19"; comparison="CONDITION:PD:HC"
# setwd("~/projects/circRNA/results/"); input_expression_filename="../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds"; input_covariance_filename="Table.PD.SNDA.pathology.covariates.xls"; output_additonal_columns='mi'; output_dir="DE_SNDA"; index="hg19"; comparison="PD.pathology.group:early:no"
# setwd("~/projects/circRNA/results/"); input_expression_filename="../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds"; input_covariance_filename="Table.PD.SNDA.pathology.covariates.xls"; output_additonal_columns='i'; output_dir="DE_SNDA"; index="hg19"; comparison="MUSS"; annotation_path='~/projects/circRNA/data/Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds'
# setwd("~/projects/circRNA/results/"); input_expression_filename="../data/Merge_circexplorer_CSF87.filtered.enriched.rawcount.rds"; input_covariance_filename="Table.PD.CSF.pathology.covariates.xls"; output_additonal_columns='mi'; output_dir="DE_CSF"; index="hg19"; comparison="CONDITION:PD:HC"; annotation_path='~/projects/circRNA/data/Merge_circexplorer_CSF87.filtered.enriched.annotation.bed14.rds'

# setwd("~/projects/circRNA/results/"); input_expression_filename="../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds"; input_covariance_filename="Table.PD.SNDA.pathology.covariates.xls"; output_additonal_columns='Mi'; output_dir="DE_SNDA"; index="hg19"; comparison="CONDITION2:PD:HC"; annotation_path='~/projects/circRNA/data/Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds'
# setwd("~/projects/circRNA/results/"); input_expression_filename="../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds"; input_covariance_filename="Table.AD.TCPY.pathology.covariates.xls"; output_additonal_columns='Mi'; output_dir="DE_TCPY"; index="hg19"; comparison="CONDITION:AD:HC"; annotation_path='~/projects/circRNA/data/Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds'
# setwd("~/projects/circRNA/results/"); input_expression_filename="../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds"; input_covariance_filename="Table.PD.SNDA.pathology.covariates.xls"; output_additonal_columns='Mi'; output_dir="DE_SNDA2gene"; index="hg19"; comparison="CONDITION2:PD:HC"; annotation_path='~/projects/circRNA/data/Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds'; COLLAPSE=TRUE

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
genome_name=switch(index, "hg19" = "Homo_sapiens", "mm9" = "Mus_musculus", "rn6" = "Rattus_norvegicus")

# check input
if(!file.exists(input_expression_filename)) {stop(paste(input_expression_filename, "doesn't exist. Exit!"), call.=FALSE);}
if(!file.exists(input_covariance_filename)) {stop(paste(input_covariance_filename, "doesn't exist. Exit!"), call.=FALSE);}

# Create folder if the directory doesn't exist
dir.create(file.path(output_dir,'report/figures'), recursive =T, showWarnings = FALSE)

pwd=getwd()
setwd(output_dir)

# install packages
suppressPackageStartupMessages(library('tidyverse',logical.return=T) || install.packages('tidyverse', repo='http://cran.revolutionanalytics.com'))
suppressPackageStartupMessages(library('RCurl',logical.return=T) || install.packages('RCurl', repo='http://cran.revolutionanalytics.com'))
suppressPackageStartupMessages(library('hexbin',logical.return=T) || install.packages('hexbin', repo='http://cran.revolutionanalytics.com'))
suppressPackageStartupMessages(library('pheatmap',logical.return=T) || install.packages('pheatmap', repo='http://cran.revolutionanalytics.com'))
suppressPackageStartupMessages(library('RColorBrewer',logical.return=T) || install.packages('RColorBrewer', repo='http://cran.revolutionanalytics.com'))
suppressPackageStartupMessages(library('hwriter',logical.return=T) || install.packages('hwriter', repo='http://cran.revolutionanalytics.com'))
suppressPackageStartupMessages(library('ggforce',logical.return=T) || install.packages('ggforce', repo='http://cran.revolutionanalytics.com'))
# source("https://bioconductor.org/biocLite.R"); 
if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager") # change to use bioC 3.9
suppressPackageStartupMessages(library('vsn',logical.return=T) || BiocManager::install('vsn'));
suppressPackageStartupMessages(library('DESeq2',logical.return=T) || BiocManager::install('DESeq2'))
suppressPackageStartupMessages(library('ReportingTools',logical.return=T) || BiocManager::install('ReportingTools'))
suppressPackageStartupMessages(library('BiocParallel',logical.return=T) || BiocManager::install('BiocParallel'))
suppressPackageStartupMessages(library('limma',logical.return=T) || BiocManager::install('limma'))
suppressPackageStartupMessages(library('IHW',logical.return=T) || BiocManager::install('IHW'))

###########################################
message("#step1: load data...")
###########################################
if(file.exists(file.path("DESeq2.RData"))) load(file.path("DESeq2.RData")) else {
  
  #annotation (Note: this is downloaded from biomart, so EnsID format is ENSGxxxxx, no version number tailing)
  genes_annotation = read.table(file.path("~/neurogen/referenceGenome",genome_name,"UCSC",index,"Annotation/Genes/annotation.genes.bed6+3"), sep="\t", quote="", header = F, stringsAsFactors = F, col.names = c("chr","start","end","geneID","score","strand","geneName","geneType","geneDescription"));
  
  ## for circRNAs
  if(tolower(tools::file_ext(annotation_path)) == "rds") {
    circRNA_annotation=readRDS(file.path(annotation_path))
  }else circRNA_annotation=read.delim(file.path(annotation_path), header = T, row.names = 1,check.names =F)
  circRNA_annotation$geneDescription = genes_annotation$geneDescription[match(sub("\\..*","", circRNA_annotation$geneID), genes_annotation$geneID)]
  circRNA_annotation$geneDescription = sub(" \\[Source:.*","", circRNA_annotation$geneDescription)
  genes_annotation = select(circRNA_annotation, ID, circType, geneID, geneName, geneType, geneDescription) %>% mutate(circID=paste0(sub("RNA","",as.character(circType)),as.character(geneName)))
  
  # raw reads count
  if(tolower(tools::file_ext(input_expression_filename)) == "rds") {
    cts=readRDS(file.path(pwd,input_expression_filename))
  }else cts=read.delim(file.path(pwd,input_expression_filename), row.names = 1,check.names =F)
  # remove those non-geneID rows, e.g. __no_feature (pre-mRNA reads) and __ambiguous (see http://htseq.readthedocs.io/en/master/count.html )
  dim(cts); cts=cts[grep("^__", rownames(cts), invert = T),]; dim(cts);
  
  ## collapse circRNA to gene, e.g. circDNAJC6 (https://stat.ethz.ch/pipermail/bioconductor/2012-February/043410.html)
  if(COLLAPSE){
    dim(cts);  cts = rownames_to_column(cts) %>% mutate(rowname=genes_annotation$circID[match(rowname, genes_annotation$ID)]) %>% group_by(rowname) %>% summarise_all(sum) %>% column_to_rownames(); dim(cts)
    genes_annotation = mutate(genes_annotation,ID=circID) %>% select(-circID) %>% distinct()
  }
 
  # covariance table
  covarianceTable = read.table(file.path(pwd,input_covariance_filename), sep="\t", header = T,check.names =F, stringsAsFactors = F)
  covarianceTable[covarianceTable==""]=NA
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
  if("CELLTYPE" %in% colnames(covarianceTable)) covarianceTable$CELLTYPE=factor(covarianceTable$CELLTYPE)
  if("BATCH" %in% colnames(covarianceTable)) covarianceTable$BATCH=factor(covarianceTable$BATCH)
  if("SEX" %in% colnames(covarianceTable)) covarianceTable$SEX=factor(covarianceTable$SEX, levels = c("F","M"))
  #if("AGE" %in% colnames(covarianceTable)) {x=cut(covarianceTable$AGE, breaks = c(20,80,110)); levels(x)=gsub("\\((.*),(.*)\\]","\\1_\\2",levels(x)); covarianceTable$AGE=x;}
  #if("RIN" %in% colnames(covarianceTable)) {x=cut(covarianceTable$RIN, breaks = c(5,8,10)); levels(x)=gsub("\\((.*),(.*)\\]","\\1_\\2",levels(x)); covarianceTable$RIN=x;}
  #if("PMI" %in% colnames(covarianceTable)) {x=cut(covarianceTable$PMI, breaks = c(0,4,100)); levels(x)=gsub("\\((.*),(.*)\\]","\\1_\\2",levels(x)); covarianceTable$PMI=x}
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
  
  if(sum(apply(cts, 1, min))==0) cts=cts+1; # to avoid 0 geometric mean (https://help.galaxyproject.org/t/error-with-deseq2-every-gene-contains-at-least-one-zero/564/2)
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = covarianceTable,
                                design= ~ 1)
  
  ## pre-filtering
  #keep <- rowSums(counts(dds)) >= 10
  keep <- rowSums(counts(dds) >= 5) >= 4  # at least 4 samples with more than 5 reads
  dim(dds); dds <- dds[keep,]; dim(dds);
  #head(sort(rowMeans(counts(dds)), decreasing = T))
  # head(counts(dds)[order(rowMeans(counts(dds)), decreasing =T),])
  
  if(nrow(dds)<5) stop("Too few genes left after filtering. Quit!", call. =F);
  
  ###########################################
  message("#step3: QA of the data [optional]")
  ###########################################
  # Note: This part is not necessary for DEseq, but important for data QA
  
  ##--------------------------------------
  message("## 3.1: compare different vairance stablization methods")
  ##--------------------------------------
  
  ntd <- normTransform(dds) # log2(x+1)
  vsd <- varianceStabilizingTransformation(dds, blind=T) # Note: blind to the design, equal to design = ~ 1
  
    # # using limma to remove covariates, it returns adjusted values in log2 scale
  vsd_adjusted_log2 <- removeBatchEffect(assay(vsd), batch=vsd$BATCH, batch2=vsd$SEX, covariates = colData(vsd)[,colnames(colData(vsd)) %in% c('AGE','PMI','RIN')])

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
  message("## 3.2: save normalized reads count")
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
  
  DDS=dds
  save(genes_annotation, covarianceTable, DDS, vsd, vsd_adjusted_log2, file="DESeq2.RData")
  
}
# script to generate html report page
makeNewImages <- function(df,...){
  imagename <- c()
  tablename <- c()
  for (i in 1:nrow(df)){
    ensId <- df$circID[i]
    geneDescription <- df$geneDescription[i]    
    pvalue <- df$pvalue[i]
    imagename[i] <- paste('plot', ensId, df$geneName[i], variable, 'pdf', sep = ".")
    tablename[i] <- paste('plot', ensId, df$geneName[i], variable, 'txt', sep = ".")
    
    message(paste(" # processing", ensId))
    
    d <- data.frame(samples=colnames(assay(vsd)), 
                    expression_vsd=assay(vsd)[ensId,], 
                    expression_log2vsd=vsd_adjusted_log2[ensId,],
                    expression_raw=assay(DDS)[ensId,]) 
    d[[variable]]=colData(vsd)[[variable]]
    
    # filter
    if(variable=="MUSS") d=d[d$MUSS %in% c(0:4),]
    d=d[!(is.na(d[[variable]]) | d[[variable]]=="NA"),]
    
    if(!file.exists(file.path('report/figures',tablename[i]))) {
      write.table(d,file.path('report/figures',tablename[i]),sep="\t", quote =F, row.names=F, col.names = T)
    }
    if(!file.exists(file.path('report/figures',imagename[i]))) 
    {
      N=length(levels(colData(dds)[[variable]]))  # here has to be dds, as variable will be set to factor in dds (not in DDS or vsd) if it's a factor
      if(N>0){  # non factor
        p=ggplot(d, aes_string(x=variable, y="expression_log2vsd")) + 
          geom_boxplot(position=position_dodge(.8), width=.5, outlier.shape = NA) +
          geom_jitter(size=1.5, position = position_jitter(width=.15)) +
          theme_bw() +
          xlab(variable) + ylab("log2(normalized adjusted expression)") + ggtitle(ensId, subtitle = geneDescription)
      }else{
        N=3; # random number to make a nearly square figure
        p=ggplot(d, aes_string(x=variable, y="expression_log2vsd")) + 
          geom_jitter(size=1.5, position = position_jitter(width=.15)) +
          theme_bw() +
          xlab(variable) + ylab("log2(normalized adjusted expression)") + ggtitle(ensId, subtitle = geneDescription)
      }
      #png(file.path('report/figures', imagename[i]), height = 250, width = 600)
      pdf(file.path('report/figures', imagename[i]), height = 4, width = 1.5*N)
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
  df$geneName <- hwrite(as.character(df$geneName), 
                        #link = paste0("http://useast.ensembl.org/",genome_name,"/Gene/Summary?db=core;g=",as.character(df$geneName)), 
                        link = paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",as.character(df$geneName)), 
                        table=F)
  return(df)
}

###########################################
message("#step4: Run DE")
###########################################
if(length(str_comparison)==1){  # e.g. --comparison="PMI" 
  message(paste("# comparison will be performed on a quantitative / continuous variable:", variable));  
  
  # subsetting
  dds=DDS[, !(is.na(DDS[[variable]]) | DDS[[variable]]=="")]
  if(variable=="MUSS") dds=dds[,dds$MUSS %in% c(0:4)]
  keep <- rowSums(counts(dds) >= 5) >= 4  # at least 4 samples with more than 5 reads
  dim(dds); dds <- dds[keep,]; dim(dds);
  
  # factorize
  covariances=c();
  for(i in c('AGE', 'RIN', 'PMI')) if(i %in% colnames(covarianceTable) & !(i %in% covariances)) {covariances=c(covariances,i);}
  for(i in c('SEX', 'BATCH', 'CELLTYPE')) if(i %in% colnames(covarianceTable) & !(i %in% covariances)) {dds[[i]] <- droplevels(dds[[i]]); if(length(levels(dds[[i]]))>1) covariances=c(covariances,i);}
  
  #message(paste(covariances))
  
  design(dds) <- as.formula(paste(" ~ ", paste(c(covariances, variable), collapse= " + ")))
  
  # Use the similar setting for single-cell RNAseq, since circRNA data here is also in ZINB distribution. See http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#recommendations-for-single-cell-analysis
  dds <- DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
  
  # We tried LRT test (similar to ANOVA) below. See http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#likelihood-ratio-test
  # It turned out LRT produced less number of DE genes. So we sticked back to Wald test (to test the null hypothesis of log2FC==0)
  #dds <- DESeq(dds, test="LRT", reduced=reduced, sfType="poscounts", minmu=1e-6, minReplicatesForReplace=Inf)
  
  resultsNames(dds)
  #plotDispEsts(dds)
  
  # If the variable is continuous then the results can be extracted using the name argument to results, where the name is one of elements returned by resultsNames(dds).
  res <- lfcShrink(dds, coef=variable,
                   cooksCutoff=FALSE, # regardless of outlier. See https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
                   parallel=TRUE, BPPARAM=MulticoreParam(4),
                   type = "normal", 
                   quiet =T)
  
  com_name = comparison
  
  summary(res)
  head(res); dim(res)
  # decimal value of Fold-change
  res$FoldChange <- 2**res$log2FoldChange
  
  # add annotation
  res <- cbind(res, genes_annotation[match(row.names(res), genes_annotation$ID), -1])
  
  
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
  res <- res[!is.na(res$pvalue),]
  dim(res)
  
  # write to xls
  write.table(as.data.frame(res), 
              file=file.path(paste("DEresult",output_dir,com_name ,"xls",sep=".")), 
              sep="\t", quote =F, na="", row.names=T, col.names = NA)
  
  DE = subset(res, !is.na(pvalue) & pvalue<=0.05 & abs(log2FoldChange)>=1)
  DE2 = subset(res, !is.na(pvalue) & pvalue<=0.05 & abs(log2FoldChange)<1)
  NDE = subset(res, pvalue>0.05)
  # ## Note: 20% of non-significant(NS) genes are randomly selected in order to increase the performance of the generating graph
  if(nrow(NDE)>5000){
    n_NS=nrow(NDE)
    NDE=NDE[sample(n_NS,round(n_NS * .20)),]
    dim(res); res=DESeqResults(rbind(NDE, DE2, DE)); dim(res);
  }
  
  ## MAKING PLOTS
  pdf(file.path(paste("DEresult",output_dir, com_name ,"pdf", sep = ".")), paper = 'USr')
  ## ==============================
  # MA plot
  ## ==============================
  #DESeq2::plotMA(DESeqResults(res), alpha = 0.05, colNonSig = "gray", main=paste0(output_dir,": ",com_name))
  
  RES=as.data.frame(res)
  DE =as.data.frame(DE)
  DE2 =as.data.frame(DE2)
  
  ## ==============================
  # vocano plot
  ## ==============================
  with(RES, plot(log2FoldChange, -log10(pvalue), 
                 pch=20, cex=0.5, main=paste0(output_dir,": ",com_name), col='gray',
                 xlab=bquote(~Log[2]~fold~change), 
                 ylab=bquote(~-log[10]~pvalue)))
  if(nrow(DE)>0) {
    with(DE, points(log2FoldChange, -log10(pvalue), pch=20, col="red", cex=1))
    with(DE, text(log2FoldChange, -log10(pvalue), labels=geneName, cex=0.5, pos=3, offset = .3, adj=c(0.5,0)))
  }
  if(nrow(DE2)>0) {
    with(DE2, points(log2FoldChange, -log10(pvalue), pch=20, col="orange", cex=.8))
    with(DE2, text(log2FoldChange, -log10(pvalue), labels=geneName, cex=0.5, pos=3, offset = .3, adj=c(0.5,0)))
  }
  abline(v=c(0,-1, 1), lty=c(3,4,4), lwd=1)
  abline(h=-log10(0.05), col="black", lty=4, lwd=1)
  
  
  ## ==============================
  message("# heatmap for top 10 DE genes")
  ## ==============================
  RES=RES[order(RES$pvalue),] # sort by padj in increasing order
  DE10=rbind(head(subset(RES, log2FoldChange<0),20),head(subset(RES, log2FoldChange>0),20)) # top 10 down-regulated and top 10 up-regulated
  
  topDE=vsd_adjusted_log2[rownames(DE10),rownames(colData(dds))]; 
  if(!COLLAPSE) rownames(topDE) = paste(DE10$geneName, rownames(topDE), sep=".")
  
  annotation_row = dplyr::select(DE10, circType,geneType, log2FoldChange, baseMean, geneName) %>% rownames_to_column() %>% 
    rowwise() %>% mutate(updown=ifelse(log2FoldChange>0,"up","down"), rowname=ifelse(COLLAPSE, rowname, paste(geneName, rowname, sep="."))) %>% 
    arrange(updown, log2FoldChange) %>% select(rowname, updown) %>% column_to_rownames()

  annotation_col = dplyr::select(as.data.frame(colData(dds)), variable) %>% rownames_to_column() %>% 
    arrange(!!as.name(variable)) %>% column_to_rownames()
  
  # reorder topDE
  topDE=topDE[rownames(annotation_row), rownames(annotation_col)]
  
  ann_colors = list(
    updown = c(up = "red", down = "blue"),
    circType = c(circRNA = "red", ciRNA = "orange"),
    geneType = c(protein_coding = "darkblue", lincRNA = "orange", Mt_rRNA='pink', 
                 antisense='yellow', pseudogene='gray', processed_transcript='lightgray')
  )
  #tmp=c("green","red"); names(tmp)=c(variable_REF, variable_ALT);
  #ann_colors=c(list(tmp), ann_colors); names(ann_colors)[1]=variable
  
  ## Scale/center each genes (by rows)
  topDE=t(scale(t(as.matrix(topDE))))
  ## trim max and min
  MAX=2; MIN=-2; topDE[topDE>MAX]=MAX;topDE[topDE<MIN]=MIN
  
  ## add noise to avoid SD=0 cases
  #topDE=topDE + matrix(rnorm(prod(dim(topDE)),0,0.01),nrow = nrow(topDE),ncol = ncol(topDE))
  
  par(cex=0.5, mar=c(5, 8, 4, 1))
  pheatmap(topDE,
           fontsize = 5,
           main =paste0(output_dir,": heatmap for top 20 DE genes"),
           border_color = NA,
           color = colorRampPalette(c("lightblue", "white", "red"))(50),
           annotation_row = annotation_row,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           drop_levels = TRUE,
           scale = "none", 
           clustering_method = 'ward.D', 
           cluster_rows = F,
           clustering_distance_rows = "correlation",
           gaps_row = as.vector(table(annotation_row$updown))[1],
           #cutree_rows = 2,
           #cutree_cols=2,
           cluster_cols = F,
           clustering_distance_cols = "correlation")
  
  if(nrow(DE)>1){
    ## ==============================
    message("# heatmap for all DE genes")
    ## ==============================
    topDE=vsd_adjusted_log2[rownames(DE),rownames(colData(dds))]; 
    if(!COLLAPSE) rownames(topDE) = paste(DE10$geneName, rownames(topDE), sep=".")
    
    annotation_row = dplyr::select(DE10, circType,geneType, log2FoldChange, baseMean, geneName) %>% rownames_to_column() %>% 
      rowwise() %>% mutate(updown=ifelse(log2FoldChange>0,"up","down"), rowname=ifelse(COLLAPSE, rowname, paste(geneName, rowname, sep="."))) %>% 
      arrange(updown, log2FoldChange) %>% select(rowname, updown) %>% column_to_rownames()
    
    annotation_col = dplyr::select(as.data.frame(colData(dds)), variable) %>% rownames_to_column() %>% 
      arrange(!!as.name(variable)) %>% column_to_rownames()
    
    # reorder topDE
    topDE=topDE[rownames(annotation_row), rownames(annotation_col)]
    
    ann_colors = list(
      updown = c(up = "red", down = "blue"),
      circType = c(circRNA = "red", ciRNA = "orange"),
      geneType = c(protein_coding = "darkblue", lincRNA = "orange", Mt_rRNA='pink', 
                   antisense='yellow', pseudogene='gray', processed_transcript='lightgray')
    )
    
    ## Scale/center each genes (by rows)
    topDE=t(scale(t(as.matrix(topDE))))
    
    ## trim max and min
    MAX=2; MIN=-2; topDE[topDE>MAX]=MAX;topDE[topDE<MIN]=MIN
    
    ## add noise to avoid SD=0 cases
    #topDE=topDE + matrix(rnorm(prod(dim(topDE)),0,0.01),nrow = nrow(topDE),ncol = ncol(topDE))
    
    par(cex=0.5, mar=c(5, 8, 4, 1))
    pheatmap(topDE,
             main =paste0(output_dir,": heatmap for all DE genes"),
             fontsize = 5,
             border_color = NA,
             color = colorRampPalette(c("lightblue", "white", "red"))(50),
             annotation_row = annotation_row,
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             drop_levels = TRUE,
             scale = "none", 
             clustering_method = 'ward.D', 
             cluster_rows = F,
             clustering_distance_rows = "correlation",
             gaps_row = as.vector(table(annotation_row$updown))[1],
             #cutree_rows = 2,
             #cutree_cols=2,
             cluster_cols = F,
             clustering_distance_cols = "correlation")
  }
  
  message("# violin plots between the groups")
  
  if(nrow(rbind(DE,DE2))>0){
    individual <- counts(dds,normalized=TRUE)[rownames(rbind(DE,DE2)),]
    colnames(individual) = dds[[variable]]
    
    df = rownames_to_column(as.data.frame(individual)) %>% gather(key=variable, value = "norm_expressed", -1) %>% 
      mutate(!!variable := gsub("(.*)\\..*","\\1",variable)) # %>% # filter(rowname=='chr13_78293666_78320990')
    
    # remove NA
    df=df[!is.na(df[[variable]]),]
    df=df[!is.na(df$norm_expressed),]
    df=df[!is.na(df$rowname),]
    
    for(i in 1:ceiling(length(unique(df$rowname))/16)){
      p = ggplot(df, aes_string(x=variable, y="norm_expressed")) + 
        geom_violin(trim=FALSE, fill="gray")+
        geom_boxplot(width=0.1)+
        scale_y_log10() + 
        facet_wrap_paginate(~rowname, nrow = 4, ncol = 4, page = i, scales='free') + 
        labs(title="Plot of expression by case and control",x="Groups", y = "Normalized expression")+
        theme_classic() 
      class(p) <- c('gg_multiple', class(p))
      print(p)
    }
  }
  
  dev.off() 
  
  message("## Exporting results to HTML")
  htmlRep <- HTMLReport(shortName=com_name, title=com_name,
                        reportDirectory="report")
  publish(makeNewImages(rownames_to_column(RES, var="circID") %>% dplyr::select(starts_with("circ"), starts_with("baseMean"), log2FoldChange, pvalue, padj, geneName, geneType, geneDescription)), htmlRep)
  finish(htmlRep)
  
  
}else if(length(str_comparison)==3){  # e.g. --comparison="DX:PD:HC"  TODO: can be multi comparisons e.g. DX:PD:HC;DX:ILB:HC
  
  message(paste("# comparison will be performed on a factor variable:", variable, "between level", variable_ALT, "and", variable_REF));  
  
  ## subsetting
  dim(DDS)
  dds=DDS[, DDS[[variable]] %in% c(variable_ALT, variable_REF)]
  keep <- rowSums(counts(dds) >= 5) >= 4  # at least 4 samples with more than 5 reads
  dim(dds); dds <- dds[keep,]; dim(dds);
  
  dds[[variable]]=factor(dds[[variable]])
  levels(dds[[variable]]) = union(variable_REF, levels(dds[[variable]])) # put the REF in the first in the levels
  
  # factorize
  covariances=c();
  for(i in c('AGE', 'RIN', 'PMI')) if(i %in% colnames(covarianceTable) & !(i %in% covariances)) {covariances=c(covariances,i);}
  for(i in c('SEX', 'BATCH', 'CELLTYPE', variable)) if(i %in% colnames(covarianceTable) & !(i %in% covariances)) {dds[[i]] <- droplevels(dds[[i]]); if(length(levels(dds[[i]]))>1) covariances=c(covariances,i);}
  
  str(colData(dds))
  
  design(dds) <- as.formula(paste(" ~ ", paste(covariances, collapse= " + ")))
  
  # Zero-inflated NB distribution, Wald test
  dds <- DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6, parallel=TRUE, BPPARAM=MulticoreParam(4))
  
  com_name= paste(variable, variable_ALT, "vs", variable_REF, sep="_")
  message(paste("processing comparison:",com_name))
  
  res <- lfcShrink(dds, contrast = c(variable, variable_ALT, variable_REF),   # contract format: factor name, numerator in the fold change, denominator in the fold change
                   type='normal',
                   cooksCutoff=FALSE, # regardless of outlier. See https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
                   parallel=TRUE, BPPARAM=MulticoreParam(4))
  ## Shrink log2 fold changes [required for version >=1.16]
  #res2 <- lfcShrink(dds, contrast = c(x[1], x[2], x[3]), res=res, type='normal', parallel=TRUE, BPPARAM=MulticoreParam(4))
  ## You can get the shrunken LFC either with lfcShrink like above or with betaPrior=TRUE. It will be the same shrunken LFC and the same as previously calculated in DESeq2. 
  ## The difference is that betaPrior=TRUE will give you a p-value for the shrunken LFC, while lfcShrink (at the moment) is only giving you the LFC, and is keeping the p-value for the test of the MLE LFC. 
  ## see https://support.bioconductor.org/p/95695/ and https://support.bioconductor.org/p/98833/#98843
  
  summary(res)
  head(res); dim(res)
  # decimal value of Fold-change
  res$FoldChange <- 2**res$log2FoldChange
  
  # add annotation
  res <- cbind(res, genes_annotation[match(row.names(res), genes_annotation$ID), -1])
  
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
  
  # remove lines with NA in p.value
  res <- res[!is.na(res$pvalue),]
  dim(res)
  
  # write to xls
  write.table(as.data.frame(res), 
              file=file.path(paste("DEresult",output_dir,com_name ,"xls",sep=".")), 
              sep="\t", quote =F, na="", row.names=T, col.names = NA)
  
  DE = subset(res, !is.na(pvalue) & pvalue<=0.05 & abs(log2FoldChange)>=1)
  DE2 = subset(res, !is.na(pvalue) & pvalue<=0.05 & abs(log2FoldChange)<1)  
  NDE = subset(res, pvalue>0.05)
  # ## Note: 20% of non-significant(NS) genes are randomly selected in order to increase the performance of the generating graph
  if(nrow(NDE)>5000){
    n_NS=nrow(NDE)
    NDE=NDE[sample(n_NS,round(n_NS * .20)),]
    dim(res); res=DESeqResults(rbind(NDE, DE2, DE)); dim(res);
  }
  
  ## MAKING PLOTS
  pdf(file.path(paste("DEresult",output_dir, com_name ,"pdf", sep = ".")), paper = 'USr')
  ## ==============================
  # MA plot
  ## ==============================
  #DESeq2::plotMA(DESeqResults(res), alpha = 0.05, colNonSig = "gray", main=paste0(output_dir,": ",com_name))
  
  RES=as.data.frame(res)
  DE =as.data.frame(DE)
  DE2 =as.data.frame(DE2)
  # scatter plot
  if(!is.null(output_additonal_columns) && grepl("M", output_additonal_columns)){ # mean of raw expression values for each group, 
    par(pty="s"); 
    REF=paste0("baseMean_norm.", variable_REF); ALT=paste0("baseMean_norm.", variable_ALT)
    plot(jitter(RES[[REF]])+1, jitter(RES[[ALT]])+1, 
         log='xy', pch=19, cex=1, cex.axis=2, cex.lab=2,
         xlab=REF,ylab=ALT,main= paste0(output_dir,": ",com_name),
         xlim=range(RES[[REF]], RES[[ALT]])+1, 
         ylim=range(RES[[REF]], RES[[ALT]])+1, 
         col=ifelse(RES$circType=="circRNA",'red','orange'))
    RES0=RES[RES[[REF]]>=max(RES[[REF]], RES[[ALT]])/2 | RES[[ALT]]>=max(RES[[REF]], RES[[ALT]])/2,]
    text(x=1+RES0[[REF]], y=1+RES0[[ALT]], labels = RES0$geneName, cex=.5, adj=c(0.5,0), pos=3, offset = .3);
    abline(a=0,b=1, lty=2, col='black')
  }
  
  ## ==============================
  # vocano plot
  ## ==============================
  with(RES, plot(log2FoldChange, -log10(pvalue), 
                 pch=20, cex=1, main=paste0(output_dir,": ",com_name), col='darkgray',
                 xlab=bquote(~Log[2]~fold~change), 
                 ylab=bquote(~-log[10]~pvalue)))
  if(nrow(DE)>0) {
    with(DE, points(log2FoldChange, -log10(pvalue), pch=20, col="red", cex=1))
    with(DE, text(log2FoldChange, -log10(pvalue), labels=geneName, cex=0.6, pos=3, offset = .3, adj=c(0.5,0)))
  }
  if(nrow(DE2)>0) {
    with(DE2, points(log2FoldChange, -log10(pvalue), pch=20, col="orange", cex=.8))
    with(DE2, text(log2FoldChange, -log10(pvalue), labels=geneName, cex=0.5, pos=3, offset = .3, adj=c(0.5,0)))
  }
  abline(v=c(0,-1, 1), lty=c(3,4,4), lwd=c(1,1,1))
  abline(h=-log10(0.05), col="black", lty=4, lwd=1)
  
  ## ==============================
  message("# heatmap for top 20 DE genes")
  ## ==============================
  #DE=DE[order(-DE$baseMean),] # sort by baseMean in decreasing order
  #DE=DE[order(-abs(DE$log2FoldChange)),] # sort by abs(log2FoldChange) in decreasing order
  RES=RES[order(RES$pvalue),] # sort by padj in increasing order
  DE10=rbind(head(subset(RES, log2FoldChange<0),20),head(subset(RES, log2FoldChange>0),20)) # top 10 down-regulated and top 10 up-regulated
  
  topDE=vsd_adjusted_log2[rownames(DE10),rownames(colData(dds))]; 
  #topDE=assay(vsd)[rownames(DE10),rownames(colData(dds))]
  if(!COLLAPSE) rownames(topDE) = paste(DE10$geneName, rownames(topDE), sep=".")
  
  annotation_row = dplyr::select(DE10, circType,geneType, log2FoldChange, baseMean, geneName) %>% rownames_to_column() %>% 
    rowwise() %>% mutate(updown=ifelse(log2FoldChange>0,"up","down"), rowname=ifelse(COLLAPSE, rowname, paste(geneName, rowname, sep="."))) %>% 
    arrange(updown, log2FoldChange) %>% select(rowname, updown) %>% column_to_rownames()
  
  if(variable=="PDpathologygroup"){
    annotation_col = rownames_to_column(as.data.frame(colData(dds))) %>% arrange(!!as.name(variable), PD.pathology.group) %>% dplyr::select(rowname,variable, PD.pathology.group) %>% column_to_rownames()
  } else if(variable=="CONDITION2") {
    annotation_col = rownames_to_column(as.data.frame(colData(dds))) %>% arrange(!!as.name(variable), CONDITION) %>% dplyr::select(rowname,variable, CONDITION) %>% column_to_rownames()
  } else
    annotation_col = rownames_to_column(as.data.frame(colData(dds))) %>% arrange(!!as.name(variable)) %>% dplyr::select(rowname,variable) %>% column_to_rownames()
  
  # reorder topDE
  topDE=topDE[rownames(annotation_row), rownames(annotation_col)]
  
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
  #topDE=topDE + matrix(rnorm(prod(dim(topDE)),0,0.01),nrow = nrow(topDE),ncol = ncol(topDE))
  
  par(cex=0.5, mar=c(5, 8, 4, 1))
  pheatmap(topDE,
           fontsize = 5,
           main =paste0(output_dir,": heatmap for top 20 DE genes"),
           #width=8, height=5,filename=paste0("DEresult.padj_05.", com_name ,".heatmap.pdf"),
           border_color = NA,
           color = colorRampPalette(c("lightblue", "white", "red"))(50),
           annotation_row = annotation_row,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           drop_levels = TRUE,
           scale = "none", 
           clustering_method = 'ward.D', 
           cluster_rows = F,
           gaps_row = as.vector(table(annotation_row$updown))[1],
           gaps_col = as.vector(table(annotation_col[[variable]]))[1],
           clustering_distance_rows = "correlation",
           cutree_rows = 2,cutree_cols=2,
           cluster_cols = F,
           clustering_distance_cols = "correlation")
  
  if(nrow(DE)>1){
    ## ==============================
    message("# heatmap for all DE genes")
    ## ==============================
    #topDE=assay(vsd)[rownames(DE),rownames(colData(dds))]; 
    topDE=vsd_adjusted_log2[rownames(DE),rownames(colData(dds))]; 
    if(!COLLAPSE) rownames(topDE) = paste(DE$geneName, rownames(topDE), sep=".") 
    
    annotation_row = dplyr::select(DE, circType,geneType, log2FoldChange, baseMean, geneName) %>% rownames_to_column() %>% 
      rowwise() %>% mutate(updown=ifelse(log2FoldChange>0,"up","down"), rowname=ifelse(COLLAPSE, rowname, paste(geneName, rowname, sep="."))) %>% 
      arrange(updown, -baseMean) %>% select(rowname, updown) %>% column_to_rownames()
    
    annotation_col = dplyr::select(as.data.frame(colData(dds)), variable) %>% rownames_to_column() %>% 
      arrange(!!as.name(variable)) %>% column_to_rownames()
    
    # reorder topDE
    topDE=topDE[rownames(annotation_row), rownames(annotation_col)]
    
    ann_colors = list(
      updown = c(up = "red", down = "blue"),
      circType = c(circRNA = "red", ciRNA = "orange"),
      geneType = c(protein_coding = "darkblue", lincRNA = "orange", Mt_rRNA='pink', 
                   antisense='yellow', pseudogene='gray', processed_transcript='lightgray')
    )
    tmp=c("green","red"); names(tmp)=c(variable_REF, variable_ALT);
    ann_colors=c(list(tmp), ann_colors); names(ann_colors)[1]=variable
    
    ## Scale/center each genes (by rows), similar to z-score
    topDE=t(scale(t(as.matrix(topDE)), center = T, scale = T))
    ## trim max and min
    MAX=2; MIN=-2; topDE[topDE>MAX]=MAX;topDE[topDE<MIN]=MIN
    
    ## add noise to avoid SD=0 cases
    #topDE=topDE + matrix(rnorm(prod(dim(topDE)),0,0.01),nrow = nrow(topDE),ncol = ncol(topDE))
    
    par(cex=0.5, mar=c(5, 8, 4, 1))
    pheatmap(topDE,
             main =paste0(output_dir,": heatmap for all DE genes"),
             fontsize = 5,
             border_color = NA,
             color = colorRampPalette(c("lightblue","white","red"))(50),
             annotation_row = annotation_row,
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             drop_levels = TRUE,
             scale = "none", 
             clustering_method = 'ward.D', 
             cluster_rows = F,
             gaps_row = as.vector(table(annotation_row$updown))[1],
             gaps_col = as.vector(table(annotation_col[[variable]]))[1],
             clustering_distance_rows = "correlation",
             cutree_rows = 2,cutree_cols=2,
             cluster_cols = F,
             clustering_distance_cols = "correlation")
  }
  
  message("# violin plots between the groups")
  
  if(nrow(rbind(DE,DE2))>0){
    individual <- counts(dds,normalized=TRUE)[rownames(rbind(DE,DE2)),]
    colnames(individual) = dds[[variable]]
    
    df = rownames_to_column(as.data.frame(individual)) %>% gather(key=variable, value = "norm_expressed", -1) %>% 
      mutate(!!variable := gsub("(.*)\\..*","\\1",variable))
    
    # remove NA
    df=df[!is.na(df[[variable]]),]
    df=df[!is.na(df$norm_expressed),]
    df=df[!is.na(df$rowname),]
    
    for(i in 1:ceiling(length(unique(df$rowname))/16)){
      p = ggplot(df, aes_string(x=variable, y="norm_expressed")) + 
        geom_violin(trim=FALSE, fill="gray")+
        geom_boxplot(width=0.1)+
        scale_y_log10() + 
        facet_wrap_paginate(~rowname, nrow = 4, ncol = 4, page = i, scales='free') + 
        labs(title="Plot of expression by case and control",x="Groups", y = "Normalized expression")+
        theme_classic() 
      class(p) <- c('gg_multiple', class(p))
      print(p)
    }
  }
  
  dev.off() 
  
  message("## Exporting results to HTML")
  
  htmlRep <- HTMLReport(shortName=com_name, title=com_name,
                        reportDirectory="report")
  publish(makeNewImages(rownames_to_column(RES, var="circID") %>% dplyr::select(starts_with("circ"), starts_with("baseMean"), log2FoldChange, pvalue, padj, geneName, geneType, geneDescription)), htmlRep)
  finish(htmlRep)
  
} else stop("Unrecognized comparison format. Please use style like --comparison=DX:PD:HC", call. = F);

###########################################
message("#step5: scp html result to web server")
###########################################
## Note: prerequisites: 1. generate public key via "ssh-keygen" if not already (~/.ssh/id_rsa.pub); 2. add ~/.ssh/id_rsa.pub to remote server's ~/.ssh/authorized_keys 
system(paste0("rsync -a --chmod=u+rwx,g+rwx,o+rwx report/* xd010@panda.dipr.partners.org:~/public_html/DE_reports/",output_dir))
message(paste0("Visit the HTML result at: http://panda.partners.org/~xd010/DE_reports/",output_dir))

message("Done!")