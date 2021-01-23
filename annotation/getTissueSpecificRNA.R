###########################################
# R script to get cell type specific RNAs/genes
# Usage: Rscript $0 Merge_circexplorer_BC.filtered.normRPM.rds
# Input: Merge_circexplorer_BC.filtered.normRPM.rds
# Output: ../results/Merge_circexplorer_BC.annotation_per_cell.xls, ../results/Merge_circexplorer_BC.cellspecific_heatmap.pdf, and ../results/Merge_circexplorer_BC.cellspecific_heatmap+.pdf
# Author: Xianjun Dong
# Date: 2017-10-28
###########################################
library(dplyr)
library(tidyverse)
library(reshape2)

setwd("~/projects/circRNA/data")

# calculate specificity score (S) as cummeRbund (https://github.com/Bioconductor-mirror/cummeRbund/blob/master/R/tools.R)
source("~/projects/circRNA/src/annotation/tools.R")

Merge_circexp_norm_filtered_and_enriched=readRDS("Merge_circexplorer_BC197.filtered.enriched.normRPM.rds")
annotation = readRDS("Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds")

# Numbers for Fig. 1b
Merge_circexp_norm_filtered_and_enriched %>% rownames_to_column('gene') %>% gather(key = "sampleID", value='fpkm', -gene) %>%
  separate(sampleID, c("Dx", "subjectID","cellType5","batch","rep"), sep = "_", remove=FALSE) %>% 
  mutate(cellType3=ifelse(cellType5 %in% c("TCPY","MCPY"), "PY", ifelse(cellType5=="SNDA","SNDA","NN"))) %>% 
  filter(fpkm>0) %>% dplyr::select(gene, cellType3, cellType5) %>% arrange(gene, cellType3, cellType5) %>% distinct() %>% 
  select(gene, cellType5) %>% group_by(cellType5) %>% summarise(n=n_distinct(gene))
# cellType5     n
# 1 FB         2193
# 2 MCPY       1263
# 3 PBMC       5985
# 4 SNDA       6699
# 5 TCPY       8381

## only 109 HC samples (59 SNDA + 43 PY + 7 NN)
Merge_circexp_norm_filtered_and_enriched = select(Merge_circexp_norm_filtered_and_enriched, starts_with("HC_"));
Merge_circexp_norm_filtered_and_enriched = Merge_circexp_norm_filtered_and_enriched[rowMeans(Merge_circexp_norm_filtered_and_enriched)>0,]
dim(Merge_circexp_norm_filtered_and_enriched); dim(annotation)
# 11636   109

# two validated circRNAs from ERC1 (Note: mean of non-zero expression goes to Fig. 2e)
Merge_circexp_raw_filtered_and_enriched=select(readRDS("Merge_circexplorer_BC197.filtered.enriched.rawcount.rds"), starts_with("HC_"))
Merge_circexp_raw_filtered_and_enriched['chr12_1480998_1519619',]  # mean of non-zero expression
###########################################
## circRNA expressed in each cell #########
###########################################
#definition: “Expressed” circular RNAs: expressed overall (e.g. >=2 unique reads overall) and in at least 1 sample in a specific cell group 
df = Merge_circexp_norm_filtered_and_enriched
#df = df %>% mutate(gene=rownames(df)) %>% melt(id.vars="gene", variable.name = "sampleID", value.name='fpkm') %>%
df = df %>% mutate(gene=rownames(df)) %>% gather(key = "sampleID", value='fpkm', -gene) %>%
  separate(sampleID, c("Dx", "subjectID","cellType5","batch","rep"), sep = "_", remove=FALSE) %>% 
  mutate(cellType3=ifelse(cellType5 %in% c("TCPY","MCPY"), "PY", ifelse(cellType5=="SNDA","SNDA","NN"))) %>% 
  filter(fpkm>0) %>% 
  dplyr::select(gene, cellType3, cellType5) %>% arrange(gene, cellType3, cellType5) %>% distinct()

## In individual circType
df %>% merge(annotation, by.x='gene',by.y='ID', all.x=T, all.y=F) %>% select(gene, cellType3, cellType5, circType) %>% group_by(cellType5, circType) %>% summarise(n=n_distinct(gene)) %>% dcast(cellType5 ~ circType, value.var="n") 
#   cellType5 circRNA ciRNA
# 1        FB    1929   264
# 2      MCPY    1246    17
# 3      PBMC    5817   168
# 4      SNDA    4830   106
# 5      TCPY    5846   101
df %>% merge(annotation, by.x='gene',by.y='ID', all.x=T, all.y=F) %>% select(gene, cellType3, cellType5, circType) %>% group_by(cellType3, circType) %>% summarise(n=n_distinct(gene)) %>% dcast(cellType3 ~ circType, value.var="n") #%>% write.table(file=clip, sep="\t", row.names=F)
#   cellType3 circRNA ciRNA
# 1        NN    6483   309
# 2        PY    6246   108
# 3      SNDA    4830   106

## annotation per cell type
df %>% group_by(gene) %>% summarize(celltype3 = paste(unique(cellType3), collapse = ','), celltype5 = paste(unique(cellType5), collapse = ',')) %>% ungroup() %>% merge(annotation, by.x='gene',by.y='ID', all.x=T, all.y=F) %>% write.table(file="../results/Merge_circexplorer_BC109.annotation_per_cell.xls", sep="\t", quote=F, row.names=F)

celltype3_n = df %>% group_by(gene) %>% summarize(celltype3 = paste(unique(cellType3), collapse = ','), celltype5 = paste(unique(cellType5), collapse = ',')) %>% group_by(celltype3) %>% summarise(n=n())
#  celltype3      n
# 1 NN          3149
# 2 NN,PY       1223
# 3 NN,PY,SNDA  1666
# 4 NN,SNDA      754
# 5 PY          2328
# 6 PY,SNDA     1137
# 7 SNDA        1379

## call eulerAPE to generate Merge_circexplorer_BC109.venn_celltype3.svg

## Q: how many dopamine neuron circRNAs also detected in blood? (5/9/2020: for AMP PD grant)
celltype5_n = df %>% group_by(gene) %>% summarize(celltype3 = paste(unique(cellType3), collapse = ','), celltype5 = paste(unique(cellType5), collapse = ',')) %>% group_by(celltype5) %>% summarise(n=n())
print(celltype5_n, n=100)
filter(celltype5_n, grepl("PBMC", celltype5), grepl("SNDA", celltype5)) %>% select(n) %>% sum()
# 2747
filter(celltype5_n, grepl("SNDA", celltype5)) %>% select(n) %>% sum()
# 6699
# 41%

# 5-ways venn diagram (not area-proportional)
library(gplots); # install.packages('gplots')
pdf("../results/Merge_circexplorer_BC109.venn_celltype5.pdf")
venn(list(SNDA = unique(df$gene[df$cellType5=='SNDA']),
          TCPY = unique(df$gene[df$cellType5=='TCPY']),
          MCPY = unique(df$gene[df$cellType5=='MCPY']),
          PBMC = unique(df$gene[df$cellType5=='PBMC']),
          FB   = unique(df$gene[df$cellType5=='FB'])
))
dev.off()

###########################################
############ cell specificity   ###########
###########################################
df = Merge_circexp_norm_filtered_and_enriched

# 5 groups mean
group5mean = df %>% mutate(gene=rownames(df)) %>% melt(id.vars="gene", variable.name = "sampleID", value.name='fpkm') %>%
  separate(sampleID, c("Dx", "subjectID","cellType","batch","rep"), sep = "_", remove=FALSE) %>%
  group_by(gene, cellType) %>%
  summarise(meanFPKM=mean(fpkm)) %>% 
  mutate(cellType=factor(cellType, levels = c("SNDA","MCPY","TCPY","FB","PBMC"))) %>% # redefine factor of cellType
  dcast(gene ~ cellType, value.var='meanFPKM')

# 3 groups mean
group3mean = df %>% mutate(gene=rownames(df)) %>% melt(id.vars="gene", variable.name = "sampleID", value.name='fpkm') %>%
  separate(sampleID, c("Dx", "subjectID","cellType","batch","rep"), sep = "_", remove=FALSE) %>%
  mutate(cellType=ifelse(cellType %in% c("TCPY","MCPY"), "PY", ifelse(cellType=="SNDA","SNDA","NN"))) %>%
  group_by(gene, cellType) %>%
  summarise(meanFPKM=mean(fpkm)) %>% 
  mutate(cellType=factor(cellType, levels = c("SNDA","PY","NN"))) %>% # redefine factor of cellType
  dcast(gene ~ cellType, value.var='meanFPKM')

## cell-specific circRNAs based on 3 major groups
groupmean_s3 = groupmean_to_specificity(group3mean)
df3 = groupmean_s3 %>% filter(Private_or_not==1) ## n = 9694
group_by(df3, celltype) %>% summarise(count=n())
#   celltype count
# 1 NN        4860
# 2 PY        3308
# 3 SNDA      1526

groupmean_s5 = groupmean_to_specificity(group5mean)
df5 = groupmean_s5 %>% filter(Private_or_not==1) 
group_by(df5, celltype) %>% summarise(count=n())
#   celltype count
# 1 FB         793
# 2 MCPY       462
# 3 PBMC      3761
# 4 SNDA      1505
# 5 TCPY      2761

## df3 and df5 are very similar
rbind(data.frame(gene=df3$gene, DF="df3"), data.frame(gene=df5$gene, DF="df5")) %>% group_by(gene) %>% summarize(DFs = paste(unique(DF), collapse = ',')) %>% group_by(DFs) %>% summarise(n=n())
#   DFs         n
# 1 df3       859
# 2 df3,df5  8835
# 3 df5       447

# Consistence of the qPCR validated ones in cell-specificity
library(googlesheets4)
qPCR_validated = read_sheet("16xwR1MjCJKzkhfEk3yfjpf9JEOPUWnkLuL0HRdGT6XY", sheet="Table S2. qPCR.circRNAs") # https://docs.google.com/spreadsheets/d/16xwR1MjCJKzkhfEk3yfjpf9JEOPUWnkLuL0HRdGT6XY/edit
# join group3mean and groupmean_s3
write_sheet(left_join(qPCR_validated, rename_if(group3mean, is.numeric, list(~str_c("mean.lcRNAseq.", .))), by=c("circRNA_ID" = "gene")) %>% 
              left_join(y=groupmean_s3, by = c("circRNA_ID" = "gene")), ss="16xwR1MjCJKzkhfEk3yfjpf9JEOPUWnkLuL0HRdGT6XY", sheet="Table S2. qPCR.circRNAs")


###########################################
# Clustering 109 samples based on cell specific circRNAs, using tSNE
###########################################
# input: df3 (the cell-specific circRNAs baesd on major groups)
res=t(Merge_circexp_norm_filtered_and_enriched[df3$gene,])  # sample on the row, gene on the columns
#dim(res)
#res=log10(res*1000+1)
##Cell Type
cellType = do.call("rbind",strsplit(rownames(res),"_"))[,3]
cellType = factor(cellType, levels = c("SNDA", "TCPY", "MCPY", "PBMC", "FB"))

## PCA
# # note: According to the R help, SVD has slightly better numerical accuracy. Therefore, the function prcomp() is preferred compared to princomp().
# # see ref: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/
# res.pca = prcomp(res, scale = F)
# library(factoextra) # install.packages("factoextra")
# fviz_eig(res.pca) 
# fviz_pca_ind(res.pca,
#              col.ind = "cos2", # Color by the quality of representation
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     # Avoid text overlapping
# )

## t-SNE
set.seed(1) # for reproducibility
library(Rtsne) # install.packages("Rtsne")

pdf("../results/Merge_circexplorer_BC109.tSNE.PCA.plot.pdf", width=10, height=10)
xy <- Rtsne(res, dims = 2, perplexity=30, max_iter = 500)$Y
# visualizing
colors = c("#F22A7B", "#0029DB", "#CC00FF", "#00FF66FF", "#CCFF00FF")
names(colors) = unique(cellType)

plot(xy, t='p', main="t-SNE", col="black", pch = 21, bg = colors[cellType], 
     cex.axis=1.5, cex.lab=1.5, cex=1.5, 
     xlim=range(xy[,1])*1.2, ylim=range(xy[,2])*1.2)
legend("topleft", cex = 1, bty = "n", 
       legend = c("SNDA", "TCPY", "MCPY", "PBMC", "FB"), 
       text.col = "black", col = "black" ,pch = 21, 
       pt.bg = c("#F22A7B", "#0029DB", "#CC00FF", "#00FF66FF", "#CCFF00FF")
)
dev.off()

## UMAP
# library(umap) #install.packages('umap')
# xy = umap(res)$layout
# # plot umap
# plot(xy, t='p', main="UMAP", col="black", pch = 21, bg = colors[cellType], 
#      cex.axis=1.5, cex.lab=1.5, cex=1.5, 
#      xlim=range(xy[,1])*1.2, ylim=range(xy[,2])*1.2)
# legend("topleft", cex = 1, bty = "n", 
#        legend = c("SNDA", "TCPY", "MCPY", "PBMC", "FB"), 
#        text.col = "black",col = "black", pch = 21,
#        pt.bg = c("#F22A7B", "#0029DB", "#CC00FF", "#00FF66FF", "#CCFF00FF")
# )


###########################################
# heatmpa of cell specific circRNAs
###########################################

#=====
## cell-specific circRNAs in the 5 minor groups
# df=log10(groupmean[df5$gene,]*1000+1)
# dissimilarity <- 1 - cor(df)
# distance <- as.dist(dissimilarity)
# hc=hclust(distance)
# hccol=as.dendrogram(hc)
# weights.dd <- order(match(c("PBMC","FB","SNDA","TCPY","MCPY"), celltypes)) # the trick is to call order() on the specific index of target dendragram
# plot(reorder(hccol, wts = weights.dd, agglo.FUN = mean))
# 
# hm.parameters <- list(df, 
#                       color = colorRampPalette(c("#999999","#ffffff", "#ca0020"))(100),
#                       scale = "row",
#                       treeheight_row = 50,
#                       treeheight_col = 30,
#                       kmeans_k = NA,
#                       show_rownames = F, show_colnames = T,
#                       clustering_method = "average",
#                       cluster_rows = TRUE, cluster_cols = as.hclust(reorder(hccol, wts = weights.dd, agglo.FUN = mean)), 
#                       clustering_distance_rows = 'correlation', 
#                       clustering_distance_cols = 'correlation')
# do.call("pheatmap", hm.parameters)
# do.call("pheatmap", c(hm.parameters, width=3, height=5, filename="../results/Merge_circexplorer_BC.cellspecific_heatmap.pdf"))

# ## cell-specific circRNAs in both 3 major and 5 minor groups
# df=log10(groupmean[unique(df3$gene,df5$gene),]*1000+1)
# dissimilarity <- 1 - cor(df)
# distance <- as.dist(dissimilarity)
# hc=hclust(distance)
# hccol=as.dendrogram(hc)
# weights.dd <- order(match(c("PBMC","FB","SNDA","MCPY","TCPY"), celltypes)) # the trick is to call order() on the specific index of target dendragram
# plot(reorder(hccol, wts = weights.dd))
# 
# dissimilarity <- 1 - cor(t(df))
# distance <- as.dist(dissimilarity)
# hcrow=hclust(distance, method = 'average')
# roworders = hcrow$order # order of the new rows in the old matrix
# library(RColorBrewer) 
# hm.parameters <- list(df, 
#                       color = colorRampPalette(c("#999999","#ffffff", "#ca0020"))(100),
#                       scale = "row",
#                       treeheight_row = 50,
#                       treeheight_col = 30,
#                       kmeans_k = NA,
#                       show_rownames = F, show_colnames = T,
#                       clustering_method = "average",
#                       cluster_rows = hcrow, cluster_cols = as.hclust(reorder(hccol, wts = weights.dd, agglo.FUN = mean)), 
#                       clustering_distance_rows = 'correlation', 
#                       clustering_distance_cols = 'correlation')
# hm=do.call("pheatmap", hm.parameters)
# do.call("pheatmap", c(hm.parameters, width=3, height=5, filename="../results/Merge_circexplorer_BC.cellspecific_heatmap+.pdf"))
#=====

## cell-specific circRNAs in 3 major groups, but show expression of group5mean
groupmean = group5mean; rownames(groupmean)=groupmean[,1]; groupmean=groupmean[,-1]
df=log10(groupmean[df3$gene,]*1000+1)
celltypes=colnames(df)
# cluster columns
dissimilarity <- 1 - cor(df)
distance <- as.dist(dissimilarity)
hc=hclust(distance)
hccol=as.dendrogram(hc)
weights.dd <- order(match(c("SNDA","TCPY","MCPY","PBMC","FB"), celltypes)) # the trick is to call order() on the specific index of target dendragram (see https://www.biostars.org/p/237067/#260551)
hccol=as.hclust(reorder(hccol, wts = weights.dd, agglo.FUN = max))

# cluster rows
dissimilarity <- 1 - cor(t(df))
distance <- as.dist(dissimilarity)
hcrow=hclust(distance, method = 'average')
## reorder to make sure in the same order as columns, e.g. SNDA-TCPY-MCPY-PBMC-FB
hcrow=as.dendrogram(hcrow)
neworder = c(labels(hcrow)[4861:nrow(df)],labels(hcrow)[1:4860])  # 4860 is the total number of NN-specific genes
weights.rows <- order(match(neworder,rownames(df)))
hcrow=as.hclust(reorder(hcrow, wts = weights.rows, agglo.FUN=max))
#plot(dendsort(hcrow))
roworders = hcrow$order # order of the new rows in the old matrix
library(RColorBrewer) 
library('pheatmap')
hm.parameters <- list(df, 
                      color = colorRampPalette(c("#999999","#ffffff", "#ca0020"))(100),
                      scale = "row",
                      treeheight_row = 50,
                      treeheight_col = 30,
                      kmeans_k = NA,
                      show_rownames = F, show_colnames = T,
                      clustering_method = "average",
                      cluster_rows = hcrow, cluster_cols = hccol,
                      #cluster_rows = T, cluster_cols = T,
                      cutree_cols = 5,#cutree_rows =5, # to allow small gap space between columns
                      clustering_distance_rows = 'correlation', 
                      clustering_distance_cols = 'correlation')
do.call("pheatmap", hm.parameters)
do.call("pheatmap", c(hm.parameters, width=3, height=5, filename="../results/Merge_circexplorer_BC109.cellspecific_heatmap++.pdf"))

###########################################
# host genes of cell-specific circRNAs
###########################################
# circRNAs in the order of presence in above figure
df_roworders = df3$gene[roworders]; head(df_roworders); length(df_roworders)  # n = 9694
DF3=groupmean_s3[df_roworders,]; dim(DF3); head(DF3); table(DF3$Private_or_not)

## add host gene of the cell-specific circRNAs
filtered_enriched_annotation=readRDS("Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds"); head(filtered_enriched_annotation); dim(filtered_enriched_annotation)
DF3$hostgene=filtered_enriched_annotation$geneName[match(DF3$gene, filtered_enriched_annotation$ID)]
DF3$hostgeneID=filtered_enriched_annotation$geneID[match(DF3$gene, filtered_enriched_annotation$ID)] 
dim(DF3); DF3=filter(DF3, !is.na(hostgeneID)); dim(DF3) # 9694

## How many host genes per celltype
DF3 %>% dplyr::select(hostgene, celltype) %>% distinct() %>% group_by(celltype) %>% summarise(n=n())
# 1 NN        2188
# 2 PY        1814
# 3 SNDA      1102

## save the host gene list into table for downstream functional enrichment analysis
DF3 %>% dplyr::select(celltype, hostgene) %>% distinct() %>% group_by(celltype) %>% mutate(n=n(), hostgene = paste0(unique(hostgene), collapse = ";")) %>% distinct() %>% write.table(file="../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.txt", sep = "\t", col.names = T, quote=F, row.names = F)
DF3 %>% dplyr::select(celltype, hostgeneID) %>% mutate(hostgeneID=sub("\\..*","",hostgeneID)) %>% distinct() %>% group_by(celltype) %>% mutate(n=n(), hostgeneID = paste0(unique(hostgeneID), collapse = ";")) %>% distinct() %>% write.table(file="../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.EnsID.txt", sep = "\t", col.names = T, quote=F, row.names = F)
## 1. GO to DAVID website to run KEGG pathway analysis, then go to getTissueSpecificRNA.makingGObarplot.R to make the figures
## 2. run DisGeNET enrichment analysis for the host genes of private circRNAs
source("../src/annotation/tools.R")  # rewrite some of DisGeNET functions
gsd=readRDS("~/neurogen/external_download/externalData/others/DisGeNET.curated_gene_disease_associations.RDS") %>% select(geneId, geneSymbol, diseaseId, diseaseName) %>% distinct()
for(i in c("NN","PY","SNDA")) {
  disease_enrichment_v2(genes = unique(subset(DF3, celltype==i, select = 'hostgene',drop = T)), gdas = gsd) %>% 
    filter(Count>3, FDR<0.05) %>% 
    select("ID","Description", "pvalue", "FDR", "GeneRatio",  "BgRatio","OR") %>% 
    write.table(paste0("../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.",i,".DisGeNet.xls"),sep="\t", na="", row.names=F) 
}


## divide into venn diagram --> ## call eulerAPE to generate Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.venn_final.pdf
DF3 %>% arrange(celltype, hostgene)%>% group_by(hostgene) %>% summarize(celltype3 = paste(unique(celltype), collapse = ',')) %>% group_by(celltype3) %>% summarise(n=n(), gene3 = paste(unique(hostgene), collapse = ';'))
DF3 %>% arrange(celltype, hostgene)%>% group_by(hostgene) %>% summarize(celltype3 = paste(unique(celltype), collapse = ',')) %>% group_by(celltype3) %>% summarise(n=n(), gene3 = paste(unique(hostgene), collapse = ';')) %>% write.table(file="../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.overlapping.txt", sep = "\t", col.names = T, quote=F, row.names = F)
# # celltype3      n gene3                                                                                                                                                                            
# 1 NN          1191 AAGAB;AASS;AATF;ABCB10;ABCC1;ABCC3;ABCC4;ABHD2;AC004381.6;AC068533.7;AC093838.4;ACBD6;ACPL2;ACTR10;ACTR2;ACTR3;ACYP2;ADAM23;ADAMTS12;ADAMTS6;ADAMTSL1;AFG3L2;AGFG1;AGGF1;AGO4;AGPAT5;AGPAT6;AHSA2;AKAP7;AKIRIN1;AKT2;ALAS1;ALG6;ALKBH8;ALS2CR12;AM…
# 2 NN,PY        509 ABI1;ACADM;ADAL;ADD1;ADK;AEBP2;AFF1;AFF3;AFF4;AFTPH;AGO3;AHCTF1;AKAP10;ALG9;ALMS1;ALS2;AMN1;ANAPC1;ANKHD1;ANKRD12;ANKRD28;ANKRD42;ANO10;APC;APLP2;ARFGEF2;ARHGAP10;ARHGAP19-SLIT1;ARHGAP21;ARID2;ARIH1;ARSG;ASAP1;ASCC1;ASH1L;ASH2L;ATAD2;ATAD2B;A…
# 3 NN,PY,SNDA   306 AAK1;ABI2;ACSL3;ACSL4;ADAM32;ADD3;ADRBK2;AGTPBP1;AHI1;AK5;AKAP13;AKAP6;AKT3;ALG8;AMPH;ANKAR;ANKRD13C;ANKRD17;ANKRD27;AQR;ARHGAP26;ARHGEF12;ARNT2;ASPH;ATF2;ATF6;ATP2B4;ATP8A1;ATP8B4;ATP9B;ATXN1;ATXN10;ATXN7;BARD1;BBS9;BCAS3;BIRC6;BIVM;BPTF;BRA…
# 4 NN,SNDA      182 ACACA;ADAM10;ADAM17;AK9;AKAP9;ALS2CR11;AMOT;ANKRD32;ANO6;ARFGEF1;ARHGEF9;ARID4B;ASXL2;ATP6V0A2;B3GALNT2;BAZ2B;BRWD1;C10orf32-ASMT;C11orf80;CACNA1D;CANX;CARF;CASC3;CBFA2T2;CCDC91;CCP110;CDH8;CEP170;CEP192;CHMP5;CIRH1A;CNOT4;CNTRL;COG3;CSDE1;CS…
# 5 PY           730 AACS;ABCA5;ABCC9;ABCD3;ABL2;AC004076.9;AC034228.7;AC066593.1;ACADSB;ACAT1;ACBD5;ACTN4;ACTR3B;ADAM22;ADAM9;ADAMTS17;ADAMTS19;ADCK4;ADPGK;AGBL4;AGK;AGO2;AHR;AKTIP;ANKH;ANKMY2;ANKRD19P;ANKRD24;ANKRD26;ANKRD30B;ANKS1A;ANO3;ANP32B;ANTXR1;AP000304.…
# 6 PY,SNDA      269 ABCA3;ABLIM1;ABLIM2;ACAD11;ADARB1;AGAP1;AGL;AMBRA1;ANK2;ANK3;ANKIB1;ANKS1B;ANO4;ANO5;APBB2;APP;ARHGAP32;ARHGEF28;ARHGEF7;ARMC9;ARNTL2;ATG7;ATP10B;ATP1A1OS;ATP2C1;ATP6V0A1;ATP6V1H;ATP8A2;ATP8B1;ATP9A;BAI3;BTBD9;BTRC;C18orf8;C6orf89;C7orf55-LUC…
# 7 SNDA         345 AC004893.11;AC012307.3;AC067956.1;ACACB;ACP6;ADAMTSL3;ADCY2;AF127936.7;AGAP3;AGPAT3;AGPS;ALDH3A2;ANKRD20A7P;AP2B1;AP2M1;ARAP1;ARFGAP3;ARHGEF10L;ARHGEF11;ARHGEF4;ASTN2;ASXL3;ATF6B;ATG5;BACE2;BEND5;BEND7;BRINP1;BSCL2;BUD13;BZW2;C10orf90;C14orf1…

## ===========
## cell-specificity of the host genes (in celltype3)
## ===========
gene_fpkm=read.table("~/neurogen/rnaseq_PD/results/merged/genes.fpkm.cufflinks.allSamples.BCv2.uniq.xls", header = T, row.names = 1, stringsAsFactors = F, check.names = F); 
gene_fpkm=gene_fpkm[,colnames(Merge_circexp_norm_filtered_and_enriched)]; head(gene_fpkm); dim(gene_fpkm);
matrix_fpkm_to_tpm <- function(fpkm) t(t(fpkm) / colSums(fpkm)) * 1e6  # Note: matrix / vector is by rows, transpose it if by columns
gene_tpm = as.data.frame(matrix_fpkm_to_tpm(gene_fpkm))  ## convert FPKM to TPM, in order to have equal sum of expression for each sample
gene_group5mean = gene_tpm %>% rownames_to_column(var = 'gene') %>% melt(id.vars="gene", variable.name = "sampleID", value.name='fpkm') %>%
  separate(sampleID, c("Dx", "subjectID","cellType","batch","rep"), sep = "_", remove=FALSE) %>%
  group_by(gene, cellType) %>%
  summarise(meanFPKM=mean(fpkm)) %>%
  mutate(cellType=factor(cellType, levels = c("SNDA","MCPY","TCPY","FB","PBMC"))) %>% # redefine factor of cellType
  dcast(gene ~ cellType, value.var='meanFPKM')
head(gene_group5mean)

gene_group3mean = gene_tpm %>% rownames_to_column(var = 'gene') %>% melt(id.vars="gene", variable.name = "sampleID", value.name='fpkm') %>%
  separate(sampleID, c("Dx", "subjectID","cellType","batch","rep"), sep = "_", remove=FALSE) %>%
  mutate(cellType=ifelse(cellType %in% c("TCPY","MCPY"), "PY", ifelse(cellType=="SNDA","SNDA","NN"))) %>%
  group_by(gene, cellType) %>%
  summarise(meanFPKM=mean(fpkm)) %>% 
  mutate(cellType=factor(cellType, levels = c("SNDA","PY","NN"))) %>% # redefine factor of cellType
  dcast(gene ~ cellType, value.var='meanFPKM') 
head(gene_group3mean)

## background: cell-specific genes
groupmean_to_specificity(gene_group3mean) %>% filter(Private_or_not==1) %>% group_by(celltype) %>% summarise(n=n())
# 1 NN         523
# 2 PY        2066
# 3 SNDA      5431

gene_group3mean_DF3 = gene_group3mean[match(unique(DF3$hostgeneID), gene_group3mean$gene),]
gene_groupmean_s3_DF3 = groupmean_to_specificity(gene_group3mean_DF3)

## cell-specificity of circRNA VS. cell-specificity of hostgene
DF3 = inner_join(DF3, gene_groupmean_s3_DF3,by = c("hostgeneID" = "gene"), suffix=c(".circRNA",".gene")); head(DF3)
DF3 = mutate(DF3, celltype.specific.gene = ifelse(Private_or_not.gene==1, celltype.gene, "NS"))
dim(DF3)

## save the cell specificty table
write.table(DF3, file = "../results/Merge_circexplorer_BC109.cellspecificity.circRNA3.genes3.xls",sep="\t", na="", row.names=F)

## Any host genes for cell-specific circRNAs are private genes? None!
for(i in c("NN","PY","SNDA")) {print(i); filter(DF3, celltype.circRNA==i, Private_or_not.gene==1) %>% dplyr::select(hostgene, celltype.gene) %>% distinct()%>% group_by(celltype.gene) %>% summarise(n=n()) %>% print()}
# or table(dplyr::select(DF3, Private_or_not.circRNA, Private_or_not.gene))
# zero

for(i in c("NN","PY","SNDA")) {print(i); filter(DF3, celltype.circRNA==i) %>% dplyr::select(hostgene, celltype.gene) %>% distinct() %>% group_by(celltype.gene) %>% summarise(n=n()) %>% print()}
ggplot(DF3, aes(x=factor(celltype.gene), fill=factor(celltype.circRNA))) +
  geom_bar(stat="count", width=0.7, fill="steelblue")+ 
  coord_flip() +
  facet_grid(cols=vars(celltype.circRNA), scales ='free') +
  theme_minimal()
ggsave("../results/Merge_circexplorer_BC109.cellspecific_barplot++.hostgene.pdf", width = 10, height = 2)
# [1] "NN"
# celltype.gene     n
# 1 NN             1094
# 2 PY              693
# 3 SNDA            346

# [1] "PY"
# celltype.gene     n
# 1 NN              631
# 2 PY             1104
# 3 SNDA            498

# [1] "SNDA"
# celltype.gene     n
# 1 NN              336
# 2 PY              570
# 3 SNDA            383

## expression level of parental genes of specific circRNAs in each cell type


## spec.circRNA vs. spec.gene for the same cell type
xx=DF3 %>% mutate(x=match(paste0(celltype.circRNA,"_spec.gene"), colnames(DF3))) 
xx = xx %>% mutate(S.circRNA.gene = DF3[cbind(seq_along(x), x)]) %>% dplyr::select(celltype.circRNA, S.circRNA, S.circRNA.gene) %>% gather("circRNA_or_gene", "specificity_score", -1, convert = T) %>% mutate(circRNA_or_gene=as.factor(circRNA_or_gene), celltype.circRNA=as.factor(celltype.circRNA), specificity_score=as.numeric(specificity_score))
ggplot(xx, aes(x=celltype.circRNA, y=specificity_score, fill=circRNA_or_gene)) + 
  geom_violin(trim = T, scale = "width", position = position_dodge(width = .51)) + 
  coord_flip() +
  theme_bw()
  #geom_boxplot(col='white',width=0.05, outlier.shape = NA, position = position_dodge(width = .51)) 
  #geom_jitter(height = 0, width = 0.1)
  #geom_boxplot() +
  #geom_point(aes(fill = circRNA_or_gene), size = 5, shape = 21, position = position_jitterdodge())
ggsave("../results/Merge_circexplorer_BC109.cellspecific_boxplot++.specificityscore.pdf", width=6, height = 2)

## heatmap of host genes of cell-specific circRNAs
X=gene_group5mean[match(DF3$hostgeneID, gene_group5mean$gene),]; head(X); dim(X); 
X=log10(X[,-1]+0.01)
dissimilarity <- 1 - cor(X)
distance <- as.dist(dissimilarity)
hc=hclust(distance)
hccol=as.dendrogram(hc)
weights.dd <- order(match(c("SNDA","TCPY","MCPY","PBMC","FB"), colnames(X))) # the trick is to call order() on the specific index of target dendragram
#plot(reorder(hccol, wts = weights.dd, agglo.FUN = mean))

my_gene_col <- data.frame(celltype=DF3$celltype.circRNA, row.names = rownames(X))
# Specify colors
ann_colors = list(
  celltype = c(SNDA="#F22A7B",PY="#3182bd", NN="#513931")
)

hm.parameters <- list(X, 
                      color = colorRampPalette(c("#999999","#ffffff", "#ca0020"))(100),
                      scale = "row",
                      treeheight_row = 50,
                      treeheight_col = 30,
                      kmeans_k = NA,
                      show_rownames = F, show_colnames = T,
                      annotation_row = my_gene_col, annotation_colors = ann_colors,
                      cutree_cols = 5, # to allow small gap space between columns
                      cluster_rows = F, cluster_cols = as.hclust(reorder(hccol, wts = weights.dd, agglo.FUN = mean)))
library(pheatmap);
do.call("pheatmap", hm.parameters)
do.call("pheatmap", c(hm.parameters, width=3, height=5, filename="../results/Merge_circexplorer_BC109.cellspecific_heatmap++.hostgene.pdf"))

#############################################
# cell-specific isoform of circRNAs, esp. those from the same host gene (esp. the 575 shared genes btw SNDA and PY)
#############################################
Merge_circexp_raw_filtered_and_enriched = readRDS(file="Merge_circexplorer_BC197.filtered.enriched.rawcount.rds")
Merge_circexp_raw_filtered_and_enriched = Merge_circexp_raw_filtered_and_enriched[rownames(Merge_circexp_norm_filtered_and_enriched),colnames(Merge_circexp_norm_filtered_and_enriched)]; dim(Merge_circexp_raw_filtered_and_enriched)
# 3 groups max
group3max = Merge_circexp_raw_filtered_and_enriched %>% mutate(gene=rownames(Merge_circexp_raw_filtered_and_enriched)) %>% 
  melt(id.vars="gene", variable.name = "sampleID", value.name='fpkm') %>%
  separate(sampleID, c("Dx", "subjectID","cellType","batch","rep"), sep = "_", remove=FALSE) %>%
  mutate(cellType=ifelse(cellType %in% c("TCPY","MCPY"), "PY", ifelse(cellType=="SNDA","SNDA","NN"))) %>%
  group_by(gene, cellType) %>%
  summarise(meanFPKM=mean(fpkm)) %>% 
  mutate(cellType=factor(cellType, levels = c("SNDA","PY","NN"))) %>% # redefine factor of cellType
  dcast(gene ~ cellType, value.var='meanFPKM')
group3nonzero = Merge_circexp_raw_filtered_and_enriched %>% mutate(gene=rownames(Merge_circexp_raw_filtered_and_enriched)) %>%
  melt(id.vars="gene", variable.name = "sampleID", value.name='fpkm') %>%
  separate(sampleID, c("Dx", "subjectID","cellType","batch","rep"), sep = "_", remove=FALSE) %>%
  mutate(cellType=ifelse(cellType %in% c("TCPY","MCPY"), "PY", ifelse(cellType=="SNDA","SNDA","NN"))) %>%
  group_by(gene, cellType) %>%
  summarise(nNONZERO=sum(fpkm>0)) %>% 
  mutate(cellType=factor(cellType, levels = c("SNDA","PY","NN"))) %>% # redefine factor of cellType
  dcast(gene ~ cellType, value.var='nNONZERO')

x=DF3 %>% filter(Private_or_not.circRNA==1) %>% arrange(celltype.circRNA) %>% 
  inner_join(y=group3max, by='gene') %>% mutate(SNDA=signif(SNDA,2), PY=signif(PY,2), NN=signif(NN,2), max3=pmax(SNDA,PY,NN)) %>%
  #filter(max3>.001) %>%
  unite('RPM_SNDA_PY_NN', c("SNDA","PY","NN")) %>%
  inner_join(y=group3nonzero, by='gene') %>% mutate(max3nonzero=pmax(SNDA,PY,NN)) %>%
  #filter(max3nonzero>3) %>%
  group_by(hostgene) %>% 
  summarise(celltype3s.circRNA = paste((celltype.circRNA), collapse = ','), 
            celltype3n.circRNA = length(unique(celltype.circRNA)), 
            genes = paste((gene), collapse = ','),
            geneN = length((gene)),
            RPMs_SNDA_PY_NN=paste(RPM_SNDA_PY_NN, collapse = ','),
            max3nonzeros = paste((max3nonzero), collapse = ','),
            max3mean=max(max3))  %>%
  arrange(-max3mean) %>% 
  filter(celltype3n.circRNA>1, hostgene %in% c('PHF21A','SETD2','ERC1','MAN2A1','PSEN1','ATXN10'))
