###########################################
# R script to get cell type specific RNAs/genes
# Usage: Rscript $0 Merge_circexplorer_BC197.filtered.enriched.normRPM.rds
# Input: Merge_circexplorer_BC197.filtered.enriched.normRPM.rds
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

# 3 groups: SNDA, PY, and NN
Merge_circexp_norm_filtered_and_enriched %>% rownames_to_column('gene') %>% gather(key = "sampleID", value='fpkm', -gene) %>%
  separate(sampleID, c("Dx", "subjectID","cellType5","batch","rep"), sep = "_", remove=FALSE) %>% 
  mutate(cellType3=ifelse(cellType5 %in% c("TCPY","MCPY"), "PY", ifelse(cellType5=="SNDA","SNDA","NN"))) %>% 
  filter(fpkm>0) %>% dplyr::select(gene, cellType3, cellType5) %>% arrange(gene, cellType3, cellType5) %>% distinct() %>% 
  select(gene, cellType3) %>% group_by(cellType3) %>% summarise(n=n_distinct(gene))
# 1 NN         6792
# 2 PY         8658
# 3 SNDA       6699

## only 109 HC samples (59 SNDA + 43 PY + 7 NN)
Merge_circexp_norm_filtered_and_enriched = select(Merge_circexp_norm_filtered_and_enriched, starts_with("HC_"));
Merge_circexp_norm_filtered_and_enriched = Merge_circexp_norm_filtered_and_enriched[rowMeans(Merge_circexp_norm_filtered_and_enriched)>0,]
dim(Merge_circexp_norm_filtered_and_enriched); dim(annotation)
# 11636   109

# 3 groups: SNDA, PY, and NN
Merge_circexp_norm_filtered_and_enriched %>% rownames_to_column('gene') %>% gather(key = "sampleID", value='fpkm', -gene) %>%
  separate(sampleID, c("Dx", "subjectID","cellType5","batch","rep"), sep = "_", remove=FALSE) %>% 
  mutate(cellType3=ifelse(cellType5 %in% c("TCPY","MCPY"), "PY", ifelse(cellType5=="SNDA","SNDA","NN"))) %>% 
  filter(fpkm>0) %>% dplyr::select(gene, cellType3, cellType5) %>% arrange(gene, cellType3, cellType5) %>% distinct() %>% 
  select(gene, cellType3) %>% group_by(cellType3) %>% summarise(n=n_distinct(gene))
# 1 NN         6792
# 2 PY         6354
# 3 SNDA       4936

# two validated circRNAs from ERC1 (Note: mean of non-zero expression goes to Fig. 2e)
Merge_circexp_raw_filtered_and_enriched=select(readRDS("Merge_circexplorer_BC197.filtered.enriched.rawcount.rds"), starts_with("HC_")) # only 109 HC remained
write.table(Merge_circexp_raw_filtered_and_enriched, file="Merge_circexplorer_BC109.filtered.enriched.raw.xls", sep = "\t", quote = F, row.names = T, col.names = NA) 

group3meancount = Merge_circexp_raw_filtered_and_enriched %>% mutate(gene=rownames(Merge_circexp_raw_filtered_and_enriched)) %>% melt(id.vars="gene", variable.name = "sampleID", value.name='readcount') %>%
  separate(sampleID, c("Dx", "subjectID","cellType","batch","rep"), sep = "_", remove=FALSE) %>%
  mutate(cellType=ifelse(cellType %in% c("TCPY","MCPY"), "PY", ifelse(cellType=="SNDA","SNDA","NN"))) %>%
  #filter(readcount>0) %>% # mean of non-zero expression only
  group_by(gene, cellType) %>%
  summarise(meancount=mean(readcount)) %>% 
  mutate(cellType=factor(cellType, levels = c("SNDA","PY","NN"))) %>% # redefine factor of cellType
  dcast(gene ~ cellType, value.var='meancount')

# For fig 2e ERC1 and fig S6 APP genes
filter(group3meancount, gene %in% c('chr12_1480998_1519619', 'chr21_27326903_27354790','chr21_27347382_27354790'))
filter(group3mean, gene %in% c('chr12_1480998_1519619', 'chr21_27326903_27354790','chr21_27347382_27354790'))

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
df %>% group_by(gene) %>% summarize(celltype3 = paste(unique(cellType3), collapse = ','), celltype5 = paste(unique(cellType5), collapse = ',')) %>% ungroup() %>% merge(annotation, by.x='gene',by.y='ID', all.x=T, all.y=F) %>% 
  write.table(file="../results/Merge_circexplorer_BC109.annotation_per_cell.xls", sep="\t", quote=F, row.names=F)

# gene    celltype3       celltype5       chrom   start   end     score   strand  thickStart      thickEnd        itemRgb exonCount       exonSizes       exonOffsets     circType        geneID  geneName        geneType
# chr1_10032075_10041228  NN,PY,SNDA      PBMC,MCPY,SNDA  chr1    10032075        10041228        0       +       10032075        10032075        0,0,0   3       171,184,140     0,3574,9013     circRNA ENSG00000173614.9       NMNAT1  protein_coding
# chr1_100335955_100343384        PY,SNDA TCPY,SNDA       chr1    100335955       100343384       0       +       100335955       100335955       0,0,0   7       182,112,124,103,98,140,188      0,358,4287,4754,4958,6058,7241  circRNA ENSG00000162688.11      AGL      protein_coding
# chr1_100340242_100343384        SNDA    SNDA    chr1    100340242       100343384       0       +       100340242       100340242       0,0,0   5       124,103,98,140,188      0,467,671,1771,2954     circRNA ENSG00000162688.11      AGL     protein_coding
# chr1_100515464_100535241        NN      PBMC    chr1    100515464       100535241       0       +       100515464       100515464       0,0,0   7       96,63,128,115,219,114,72        0,8757,9969,11926,18068,18564,19705     circRNA ENSG00000156875.9       HIAT1    protein_coding

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
## add host gene of the cell-specific circRNAs
filtered_enriched_annotation=readRDS("Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds"); head(filtered_enriched_annotation); dim(filtered_enriched_annotation)
groupmean_s3$hostgene=filtered_enriched_annotation$geneName[match(groupmean_s3$gene, filtered_enriched_annotation$ID)]
groupmean_s3$hostgeneID=filtered_enriched_annotation$geneID[match(groupmean_s3$gene, filtered_enriched_annotation$ID)] 
dim(groupmean_s3); groupmean_s3=filter(groupmean_s3, !is.na(hostgeneID)); dim(groupmean_s3) # 11636

# circRNAs in the order of presence in above figure
df_roworders = df3$gene[roworders]; head(df_roworders); length(df_roworders)  # n = 9694
DF3=groupmean_s3[df_roworders,]; dim(DF3); head(DF3); table(DF3$Private_or_not)

# head(DF3)  ==> data frame for cell-specific circRNAs and their host genes
#                           SNDA_spec   PY_spec   NN_spec                     gene         S  celltype         mean         m2sd Private_or_not hostgene         hostgeneID
# chr20_35288739_35299800  5.584548e-02 0.3399192 0.5199240  chr20_35288739_35299800 0.5199240       NN 0.0048869395 0.0046959633              1    NDRG3 ENSG00000101079.16
# chr15_76146726_76165909  1.632801e-01 0.2455835 0.5367755  chr15_76146726_76165909 0.5367755       NN 0.0086112068 0.0079411967              1   UBE2Q2  ENSG00000140367.7
# chr3_171965322_172016577 7.056216e-02 0.2479733 0.6029572 chr3_171965322_172016577 0.6029572       NN 0.0065904238 0.0060734101              1   FNDC3B  ENSG00000075420.8
# chr1_213251037_213290752 0.000000e+00 0.2992331 0.5869380 chr1_213251037_213290752 0.5869380       NN 0.0011578848 0.0011063056              1  RPS6KC1  ENSG00000136643.7
# chr6_125330319_125379252 1.110223e-16 0.3062132 0.5797380 chr6_125330319_125379252 0.5797380       NN 0.0006942521 0.0006712448              1   RNF217 ENSG00000146373.12
# chr1_201823722_201828122 0.000000e+00 0.1836853 0.7098156 chr1_201823722_201828122 0.7098156       NN 0.0020827564 0.0019134237              1     IPO9  ENSG00000198700.5

# save DF3
saveRDS(DF3, "Merge_circexplorer_BC109.cellspecific_heatmap.circRNA3.rds")
saveRDS(groupmean_s3, "Merge_circexplorer_BC109.filtered.enriched.groupmean_s3.rds")
#DF3=readRDS("Merge_circexplorer_BC109.cellspecific_heatmap.circRNA3.rds")
#groupmean_s3=readRDS("Merge_circexplorer_BC109.filtered.enriched.groupmean_s3.rds")

## How many host genes per celltype
DF3 %>% dplyr::select(hostgene, celltype) %>% distinct() %>% group_by(celltype) %>% summarise(n=n())
# 1 NN        2188
# 2 PY        1814
# 3 SNDA      1102

## 1. GO/pathway enrichment analysis for the host genes of cell-specific circRNAs
## -----------
## Method1: save the host gene list into table and use DAVID for KEGG pathway analysis  -- (NOT USED FINALLY)
# DF3 %>% dplyr::select(celltype, hostgene) %>% distinct() %>% group_by(celltype) %>% mutate(n=n(), hostgene = paste0(unique(hostgene), collapse = ";")) %>% distinct() %>% write.table(file="../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.txt", sep = "\t", col.names = T, quote=F, row.names = F)
# DF3 %>% dplyr::select(celltype, hostgeneID) %>% mutate(hostgeneID=sub("\\..*","",hostgeneID)) %>% distinct() %>% group_by(celltype) %>% mutate(n=n(), hostgeneID = paste0(unique(hostgeneID), collapse = ";")) %>% distinct() %>% write.table(file="../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.EnsID.txt", sep = "\t", col.names = T, quote=F, row.names = F)
## GO to DAVID website to run KEGG pathway analysis, then go to getTissueSpecificRNA.makingGObarplot.R to make the figures

## Method2: Use ORA() for Fisher test with all genes as background -- (used in the final version)
source("~/neurogen/pipeline/RNAseq/bin/lib.R")
genes_annotation=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed", header = F, col.names = c('chr','start','end','EnsID','score','strand','symbol','type'), check.names = F, stringsAsFactors=F)
for(i in c("NN","PY","SNDA")) {
  message(paste("processing cell type",i,"..."));
  ORA(inputGenes = unique(subset(DF3, celltype==i, select = 'hostgene',drop = T)), allGenes=genes_annotation$symbol, topN=10, output=paste0("../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.",i))
}

# Q: We found that xxx (xx%), xxx (xx%), and xxx (xx%) cell-specific circRNAs annotated to endocytosis in each of the three cell types, respectively
for(i in c("NN","PY","SNDA")) {
  ora=read.table(paste0("../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.",i,".ORA.significant.xls"), sep = "\t", header = T, row.names = 1)
  x=subset(DF3, celltype==i, select = 'hostgene',drop = T) %in% strsplit(ora$genelist[ora$V1=="KEGG_ENDOCYTOSIS"], ",")[[1]]
  message(paste(i,sum(x), 100*mean(x)))
}

# Q: We found that xxx (xx%), xxx (xx%), and xxx (xx%) cell-specific circRNAs annotated to endocytosis/synapse/intersection in each of the three cell types, respectively
## synaptic endocytosis
## non-synaptic endocytosis
## non-endocytosis synaptic
GO_df = read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/msigdb_v7.2/MSigDB.C5.GO_terms_annotation.txt", sep="\t",header = TRUE)
endocytosis_genes=strsplit(GO_df %>% filter(Sub_Catgory=="GO:BP", Name=="GO_ENDOCYTOSIS") %>% pull(geneSymbols), ",")[[1]]
synapse_genes=strsplit(GO_df %>% filter(Sub_Catgory=="GO:CC", Name=="GO_SYNAPSE") %>% pull(geneSymbols), ",")[[1]]
endocytosis_AND_synapse=intersect(endocytosis_genes, synapse_genes) 
endocytosis_NOT_synapse=setdiff(endocytosis_genes, synapse_genes) 
synapse_NOT_endocytosis=setdiff(synapse_genes, endocytosis_genes) 
private_circRNA_percent = data.frame()
for(i in c("NN","PY","SNDA")) {
  ii=subset(DF3, celltype==i, select = 'hostgene',drop = T) # number of circRNA. For host genes, use unique(ii)
  message(paste(i,length(ii), sum(ii %in% endocytosis_NOT_synapse), sum(ii %in% endocytosis_AND_synapse), sum(ii %in% synapse_NOT_endocytosis)))
  private_circRNA_percent=rbind(private_circRNA_percent, c(sum(ii %in% endocytosis_NOT_synapse), sum(ii %in% endocytosis_AND_synapse), sum(ii %in% synapse_NOT_endocytosis))/length(ii))
}
# NN 4860 292 78 295
# PY 3308 104 69 520
# SNDA 1526 48 36 209
colnames(private_circRNA_percent)=c('endocytosis_NOT_synapse', 'endocytosis_AND_synapse', 'synapse_NOT_endocytosis')
rownames(private_circRNA_percent)=c("NN","PY","SNDA")

private_circRNA_percent %>% rownames_to_column('celltype') %>% pivot_longer(contains("_")) %>% 
  ggplot(aes(x=factor(celltype, levels = c("SNDA", "PY", "NN")), y=value, fill=celltype))+
  geom_bar(stat='identity')+
  facet_wrap(~factor(name, levels = c('endocytosis_NOT_synapse', 'endocytosis_AND_synapse', 'synapse_NOT_endocytosis')), scales="free_y") +
  scale_fill_manual(values=c("SNDA" = "#F22A7B", "PY" = "#3182bd", "NN" = "#513931")) +
  theme_classic()
ggsave("../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.endocytosis_OR_synapse.pdf", width=8, height = 6)
  
# combined heatmap ==> Now go to Fig 3
# (like Fig5 of https://www.cell.com/action/showPdf?pii=S2211-1247%2819%2931155-6)
ORA=cbind(set="all",read.delim("../results/Merge_circexplorer_BC190.filtered.enriched.circRNA.all.ORA.all.xls", row.names = 1))
#ORA=data.frame()
for(i in c("NN","PY","SNDA")) {
  message(paste("processing cell type",i,"..."));
  ora=read.delim(paste0("../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.",i,".ORA.all.xls"), row.names = 1)
  ora_slimmed=ORA_slim(ora)
  ORA=rbind(ORA, cbind(set=paste0(i,".specific"),ora_slimmed))
}
head(ORA); dim(ORA)

# TODO: save slimmed significant GO terms into Supplementary Table
saveRDS(ORA, file = "Merge_circexplorer_BC109.all_and_cellspecific_hostgenes.ORA.rds")
ORA %>% filter(FDR<0.05) %>% # dim()
  arrange(set, gene_set, FDR) %>%
  write.table(file = "../results/Merge_circexplorer_BC109.all_and_cellspecific_hostgenes.ORA.GO_KEGG.xls", sep = "\t", row.names = F)

library("RColorBrewer")
# include only the top ones and don't show p-values for the terms in the non-significant groups
ORA.pvalues = ORA %>% group_by(set, gene_set) %>% slice_min(order_by = pvalue, n = 10) %>%
  mutate(log10p=-log10(pvalue)) %>% # convert p to -log10(p)
  pivot_wider(id_cols = V1, names_from = set, values_from = log10p, values_fill = 0) %>% # convert to wide table
  column_to_rownames(var="V1") # GO term to rownames

# # include only the top ones but also show p-values for the terms in the non-significant groups
# significantTerms = ORA %>% group_by(set, gene_set) %>% slice_min(order_by = pvalue, n = 10) %>% pull(V1) %>% unique() # take top 10 in each set 
# ORA.pvalues = ORA %>% mutate(log10p=-log10(pvalue)) %>% # convert p to -log10(p)
#   pivot_wider(id_cols = V1, names_from = set, values_from = log10p, values_fill = 0) %>% # convert to wide table
#   filter(V1 %in% significantTerms) %>% 
#   column_to_rownames(var="V1") # GO term to rownames

ORA.pvalues[ORA.pvalues>100]=100
pheatmap::pheatmap(ORA.pvalues, color = colorRampPalette(c("#FFFFFF",rev(heat.colors(7))))(100))#, clustering_method = "ward.D")
pheatmap::pheatmap(ORA.pvalues, color = colorRampPalette(c("#FFFFFF",rev(heat.colors(7))))(100), filename="../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.ORA.heatmap.pdf", width=8, height = 10)
#ORA.pvalues[ORA.pvalues>20]=20
# separate GO and KEGG
# GO only
pheatmap::pheatmap(ORA.pvalues[grep("GO_", rownames(ORA.pvalues)),], color = colorRampPalette(c("#FFFFFF",rev(heat.colors(7))))(100))
pheatmap::pheatmap(ORA.pvalues[grep("GO_", rownames(ORA.pvalues)),], color = colorRampPalette(c("#FFFFFF",rev(heat.colors(7))))(100), clustering_method = "ward.D", filename="../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.ORA.GO.heatmap.pdf", width=8, height = 8)

# KEGG only
pheatmap::pheatmap(ORA.pvalues[grep("KEGG_", rownames(ORA.pvalues)),], color = colorRampPalette(c("#FFFFFF",rev(heat.colors(7))))(100))
pheatmap::pheatmap(ORA.pvalues[grep("KEGG_", rownames(ORA.pvalues)),], color = colorRampPalette(c("#FFFFFF",rev(heat.colors(7))))(100), clustering_method = "ward.D", filename="../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.ORA.KEGG.heatmap.pdf", width=8, height = 4)

# # GO only
# # include only the significant GO terms but also show p-values for the terms in the non-significant groups
# significantTerms = ORA %>% filter(grepl("c5.go", gene_set), bonferroni<1e-4, OR>1, n>15) %>% pull(V1) %>% unique() # terms with bonferroni<0.05 in any of GO categories
# length(significantTerms)
# ORA.pvalues = ORA %>% mutate(log10p=-log10(pvalue)) %>% # convert p to -log10(p)
#   pivot_wider(id_cols = V1, names_from = set, values_from = log10p, values_fill = 0) %>% # convert to wide table
#   filter(V1 %in% significantTerms) %>%
#   column_to_rownames(var="V1") # GO term to rownames
# ORA.pvalues = ORA %>% filter(grepl("c5.go", gene_set), bonferroni<1e-4, OR>1, n>15) %>%
#   mutate(log10p=-log10(pvalue)) %>% # convert p to -log10(p)
#   pivot_wider(id_cols = V1, names_from = set, values_from = log10p, values_fill = 0) %>% # convert to wide table
#   column_to_rownames(var="V1")
# pheatmap::pheatmap(ORA.pvalues[grep("GO_", rownames(ORA.pvalues)),], color = colorRampPalette(c("#FFFFFF",rev(heat.colors(7))))(100))
# 
# # Create a scatter plot
# library("plot3D")
# with(ORA.pvalues, scatter3D(NN.specific, PY.specific, SNDA.specific, phi = 0, bty = "g", pch = 20, cex = 0.5))
# # Add text
# text3D(ORA.pvalues$NN.specific, ORA.pvalues$PY.specific, ORA.pvalues$SNDA.specific,  labels = rownames(ORA.pvalues),
#        add = TRUE, colkey = FALSE, cex = 0.5)

## 2. run DisGeNET disease enrichment analysis for the host genes of private circRNAs
source("../src/annotation/tools.R")  # rewrite some of DisGeNET functions
gsd=readRDS("~/neurogen/external_download/externalData/others/DisGeNET.curated_gene_disease_associations.RDS") %>% select(geneId, geneSymbol, diseaseId, diseaseName) %>% distinct()
results=data.frame()
for(i in c("NN","PY","SNDA")) {
  message(paste("processing cell type",i,"..."));
  res_enrich = disease_enrichment_v2(genes = unique(subset(DF3, celltype==i, select = 'hostgene',drop = T)), gdas = gsd)
  results = rbind(results, cbind(celltype=i, res_enrich))
  filter(res_enrich, Count>3, FDR<0.05) %>% 
    select("ID","Description", "pvalue", "FDR", "GeneRatio",  "BgRatio","OR") %>% 
    write.table(paste0("../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.",i,".DisGeNet.xls"),sep="\t", na="", row.names=F) 
}
write.table(results, paste0("../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.DisGeNet.ALL.xls"),sep="\t", na="", row.names=F) 
results = read.table(paste0("../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.DisGeNet.ALL.xls"), sep = "\t", header = T)
## long to wide
# select diseases that one or more p values lower than the Bonferroni threshold
N=results %>% select("ID","Description", "pvalue", "celltype") %>% pivot_wider(names_from = celltype, values_from = pvalue) %>% nrow()
CUTOFF=0.05/N # n = 4169 diseases tested in total
results %>% select("ID","Description", "pvalue", "celltype") %>% 
  pivot_wider(names_from = celltype, values_from = pvalue) %>% 
  filter(NN<=CUTOFF | PY<=CUTOFF | SNDA<=CUTOFF) %>%
  write.table("../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.DisGeNet.SIG.Bonferroni_0.05_4169.xls",sep="\t", na="", row.names=F)

# ## draw barplot as Clemens's proposed (on 3/12/2021) -- sort by NN
# results %>% select("ID","Description", "pvalue", "celltype") %>% 
#   pivot_wider(names_from = celltype, values_from = pvalue) %>% 
#   filter(NN<=CUTOFF | PY<=CUTOFF | SNDA<=CUTOFF) %>%
#   arrange(desc(NN)) %>%
#   mutate(Description = factor(Description, unique(as.character(Description)))) %>% 
#   pivot_longer(NN:SNDA, names_to="celltype", values_to="pvalue") %>%
#   ggplot(aes(x=Description, y=-log10(pvalue), fill=celltype, color=celltype, ymax=max(-log10(pvalue))*1.1)) + 
#   geom_bar(width=.2, position = position_dodge(width=1), stat="identity") + 
#   geom_point(position = position_dodge(width=1), stat="identity") +
#   geom_hline(yintercept=-log10(CUTOFF), size=.5,linetype = 2) +  ## Bonferroni correction cutoff
#   scale_fill_manual(values=c("SNDA" = "#F22A7B", "PY" = "#3182bd", "NN" = "#513931")) +
#   scale_color_manual(values=c("SNDA" = "#F22A7B", "PY" = "#3182bd", "NN" = "#513931")) +
#   theme_bw() +
#   theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8), legend.justification=c(0,1), legend.position=c(0,1)) +
#   ggtitle("Disease/traits enriched in celltype-specific circRNAs' hostgene (sorted by pvalue)") +
#   ylim(-2, 15)
# 
# x=c("Bipolar Disorder", "Autistic Disorder", "Atrial Fibrillation", "Paroxysmal atrial fibrillation", "Persistent atrial fibrillation", "familial atrial fibrillation", 
#     "Substance abuse problem", "Drug abuse", "Drug habituation", "Drug Use Disorders", "Organic Mental Disorders, Substance-Induced", "Substance Dependence", "Substance Use Disorders", "Substance-Related Disorders", "Drug Dependence", "Prescription Drug Abuse", 
#     "Profound Mental Retardation", "Mental Retardation, Psychosocial", "Mental deficiency","Global developmental delay",  "Intellectual Disability",  "Neurodevelopmental Disorders", 
#     "Leukemia, Myelocytic, Acute", "Adenocarcinoma of large intestine")
# 
# results %>% select("ID","Description", "pvalue", "OR", "celltype") %>% 
#   pivot_wider(names_from = celltype, values_from = c(pvalue, OR)) %>% 
#   filter(pvalue_NN<=CUTOFF | pvalue_PY<=CUTOFF | pvalue_SNDA<=CUTOFF) %>%
#   #arrange(desc(pvalue_NN)) %>% mutate(Description = factor(Description, unique(as.character(Description)))) %>% 
#   mutate(Description = factor(Description, levels = x)) %>%
#   pivot_longer(contains("_"), names_to=c(".value", "celltype"), names_pattern = "(.*)_(.*)") %>%
#   ggplot(aes(x=Description, y=-log10(pvalue), fill=celltype, color=celltype, alpha=pvalue<=CUTOFF, ymax=max(-log10(pvalue))*1.1)) + 
#   #geom_bar(colour =NA, width=.1, position = position_dodge(width=1), stat="identity") + 
#   #geom_point(stroke = 0, shape=16, aes(size=OR), position = position_dodge(width=1)) +
#   geom_point(stroke = 0, shape=16, aes(size=OR)) +
#   geom_hline(yintercept=-log10(CUTOFF), size=.5,linetype = 2) +  ## Bonferroni correction cutoff
#   scale_fill_manual(values=c("SNDA" = "#F22A7B", "PY" = "#3182bd", "NN" = "#513931")) +
#   scale_color_manual(values=c("SNDA" = "#F22A7B", "PY" = "#3182bd", "NN" = "#513931")) +
#   theme_bw() +
#   theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8)) +
#   ggtitle("Disease/traits enriched in hostgenes of celltype-specific circRNA") 
# ggsave("../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.DisGeNet.SIG.Bonferroni_0.05_4169.pdf", width = 8, height = 6)

CUTOFF=0.05/length(unique(gsd$diseaseId)) # n = all 11181 diseaese in DisGeNET
results %>% select("ID","Description", "pvalue", "celltype") %>% 
  pivot_wider(names_from = celltype, values_from = pvalue) %>% 
  filter((as.numeric(NN)<=CUTOFF | as.numeric(PY)<=CUTOFF | as.numeric(SNDA)<=CUTOFF)) %>%
  write.table("../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.DisGeNet.SIG.Bonferroni_0.05_11181.xls",sep="\t", na="", row.names=F)

## Add PD/AD-gwas gene set
ADgenes=unique(scan("~/neurogen/external_download/externalData/GWAS/AD/AD_gwas_associatedGene.txt", character())); length(ADgenes)
#PDgenes=unique(system("grep -v Candidate ~/neurogen/external_download/externalData/GWAS/PD/PD.GWAS.Chang2017.table1n2.xls | cut -f3 | sed 's/,/\\\n/g' | sort -u", intern = TRUE)); length(PDgenes)  # Chang et al. 2017
PDgenes=unique(system("grep -v Nearest ~/neurogen/external_download/externalData/GWAS/PD/PD.GWAS.Nalls2019.TableS2.txt | cut -f4,5 | sed 's/\\\t$//;s/\\\t/\\\n/g' | sort -u", intern = TRUE)); length(PDgenes) # Nalls et al. 2019
annotation=readRDS(file="Merge_circexplorer_BC190.annotation.bed14.rds")
GENCODEv19=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v19.annotation.genes.bed", 
                      col.names = c("chr","start","end","geneName","geneID","geneType","strand"), stringsAsFactors = F, header = F)
gwas_list = list(AD_risk_genes=ADgenes, PD_risk_genes=PDgenes)
universe=unique(GENCODEv19$geneName)
DF3=readRDS("Merge_circexplorer_BC109.cellspecific_heatmap.circRNA3.rds")

results_gwas=data.frame()
for(i in c("NN","PY","SNDA")) {
  message(paste("processing cell type",i,"..."));
  genes = unique(subset(DF3, celltype==i, select = 'hostgene',drop = T))
  
  data <- data.frame(ID = character(), Description = character(), 
                     GeneRatio = character(), BgRatio = character(), OR=numeric(), geneID = character(), 
                     pvalue = numeric(), Count = integer(), gg=numeric(), stringsAsFactors = FALSE)
  for(j in 1:length(gwas_list)){
    aa <- gwas_list[[j]]
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
    data[j, ] <- c(as.character(names(gwas_list)[j]), as.character(names(gwas_list)[j]), GeneRatio, BgRatio, OR, geneID, pv, inter, inter/length(genes))
  }
  results_gwas = rbind(results_gwas, cbind(celltype=i, data))
}
results_gwas$pvalue <- as.numeric(results_gwas$pvalue)
results_gwas$Count <- as.numeric(as.character(results_gwas$Count))
results_gwas$OR <- as.numeric(results_gwas$OR)
results_gwas$gg <- as.numeric(results_gwas$gg)

results_disgenet = read.table(paste0("../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.DisGeNet.ALL.xls"), sep = "\t", header = T)

results = rbind(results_gwas %>% select("celltype", "ID", "Description", "GeneRatio", "BgRatio", "OR", "geneID", "pvalue", "Count", "gg"), 
       results_disgenet %>% select("celltype", "ID", "Description", "GeneRatio", "BgRatio", "OR", "geneID", "pvalue", "Count", "gg"))

write.table(results, "../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.DisGeNet+GWAS.ALL.xls",sep="\t", na="", row.names=F)

# N=results %>% select("ID","Description", "pvalue", "celltype") %>% pivot_wider(names_from = celltype, values_from = pvalue) %>% nrow()
# CUTOFF=0.05/N # n = 4169+2 diseases tested in total
# 
# x=c("Bipolar Disorder", "Autistic Disorder", "Atrial Fibrillation", "Paroxysmal atrial fibrillation", "Persistent atrial fibrillation", "familial atrial fibrillation", 
#     "Substance abuse problem", "Drug abuse", "Drug habituation", "Drug Use Disorders", "Organic Mental Disorders, Substance-Induced", "Substance Dependence", "Substance Use Disorders", "Substance-Related Disorders", "Drug Dependence", "Prescription Drug Abuse", 
#     "Profound Mental Retardation", "Mental Retardation, Psychosocial", "Mental deficiency","Global developmental delay",  "Intellectual Disability",  "Neurodevelopmental Disorders", 
#     "Leukemia, Myelocytic, Acute", "Adenocarcinoma of large intestine",
#     "AD_risk_genes","PD_risk_genes")
# 
# results %>% pivot_wider(names_from = celltype, values_from = c(pvalue, OR)) %>% 
#   filter(Description %in% x) %>%
#   #arrange(desc(pvalue_NN)) %>% mutate(Description = factor(Description, unique(as.character(Description)))) %>% 
#   mutate(Description = factor(Description, levels = x)) %>%
#   pivot_longer(contains("_"), names_to=c(".value", "celltype"), names_pattern = "(.*)_(.*)") %>%
#   ggplot(aes(x=Description, y=-log10(pvalue), fill=celltype, color=celltype, alpha=pvalue<=CUTOFF, ymax=max(-log10(pvalue))*1.1)) + 
#   #geom_bar(colour =NA, width=.1, position = position_dodge(width=1), stat="identity") + 
#   #geom_point(stroke = 0, shape=16, aes(size=OR), position = position_dodge(width=1)) +
#   geom_point(stroke = 0, shape=16, aes(size=OR)) +
#   geom_hline(yintercept=-log10(CUTOFF), size=.5,linetype = 2) +  ## Bonferroni correction cutoff
#   scale_fill_manual(values=c("SNDA" = "#F22A7B", "PY" = "#3182bd", "NN" = "#513931")) +
#   scale_color_manual(values=c("SNDA" = "#F22A7B", "PY" = "#3182bd", "NN" = "#513931")) +
#   theme_bw() +
#   theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8)) +
#   ggtitle("Disease/traits enriched in hostgenes of celltype-specific circRNA") 
# ggsave("../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.DisGeNet+GWAS.SIG.Bonferroni_0.05_4171.pdf", width = 8, height = 6)

##all+cell-specific combined
library("scales")
all = read.table("../results/Merge_circexplorer_BC190.filtered.enriched.circRNA.expressedinNeuron.DisGeNet+gwas.ALL.xls", header = T) 
cell_specific = read.table("../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.DisGeNet+GWAS.ALL.xls", header = T)
combined = rbind(select(all, celltype, ID, Description, OR, pvalue,Count), select(cell_specific, celltype, ID, Description, OR, pvalue, Count))
CUTOFF=0.05/length(unique(combined$ID)) # N = 4638; CUTOFF = 1.078051e-05

# combined %>% pivot_wider(names_from = celltype, values_from = c(pvalue, OR, Count)) %>% 
#   filter(pvalue_all <= CUTOFF | pvalue_NN <= CUTOFF | pvalue_PY <= CUTOFF | pvalue_SNDA <= CUTOFF) %>%
#   arrange(desc(-OR_all)) %>% mutate(Description = factor(Description, unique(as.character(Description)))) %>% 
#   pivot_longer(contains("_"), names_to=c(".value", "celltype"), names_pattern = "(.*)_(.*)") %>%
#   ggplot(aes(x=OR, y=Description, size=Count, color=-log10(pvalue))) + 
#   geom_point() +
#   facet_wrap(vars(celltype), nrow = 1) + 
#   scale_colour_gradientn(colours = c("white","black","red"), name = "-log10(pvalue)", values = rescale(c(0, -log10(CUTOFF), 5.5)), breaks =c(0, -log10(CUTOFF), 5.5), limits=c(0,5.5), oob=scales::squish) +
#   theme_bw() +
#   theme(legend.position="bottom", axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8)) +
#   ggtitle("Disease/traits enriched in hostgenes of all & celltype-specific circRNAs") 

# different color scales for significant and non-significant ones (see Ref: https://eliocamp.github.io/ggnewscale/)
library(ggnewscale)
XX=combined %>% pivot_wider(names_from = celltype, values_from = c(pvalue, OR, Count)) %>% 
  filter(pvalue_all <= CUTOFF | pvalue_NN <= CUTOFF | pvalue_PY <= CUTOFF | pvalue_SNDA <= CUTOFF) %>%
  arrange(desc(-OR_all)) %>% mutate(Description = factor(Description, unique(as.character(Description)))) %>% 
  pivot_longer(contains("_"), names_to=c(".value", "celltype"), names_pattern = "(.*)_(.*)") 
ggplot() + 
  geom_point(data=subset(XX,pvalue > CUTOFF), shape=16, aes(x=OR, y=Description, size=Count,  color=-log10(pvalue))) +
  facet_wrap(vars(factor(celltype, levels = c("all","SNDA","PY","NN"))), nrow = 1) + 
  scale_y_discrete(drop=FALSE) + 
  scale_colour_gradientn(colours = c("#eeeeee","black"), name = "-log10P") +
  new_scale_color() +
  geom_point(data=subset(XX,pvalue <= CUTOFF), shape=16, aes(x=OR, y=Description, size=Count, color=-log10(pvalue))) +
  scale_colour_gradientn(colours = c("#ff000033","#ff0000ff"), name = "", values = rescale(c(-log10(CUTOFF), 6)), limits=c(-log10(CUTOFF),6), oob=scales::squish) +
  theme_bw() + 
  theme(legend.position="bottom", axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8)) +
  ggtitle("Disease/traits enriched in hostgenes of all & celltype-specific circRNAs") 

ggsave("../results/Merge_circexplorer_BC109.all+cellspecific.combined.DisGeNet+GWAS.SIG.Bonferroni_0.05.pdf", width = 7.5, height = 6)


rbind(select(all, celltype, ID, Description, OR, pvalue,GeneRatio, BgRatio), select(cell_specific, celltype, ID, Description, OR, pvalue, GeneRatio, BgRatio)) %>% 
  pivot_wider(names_from = celltype, values_from = c(pvalue, OR, GeneRatio, BgRatio)) %>% 
  filter(pvalue_all <= CUTOFF | pvalue_NN <= CUTOFF | pvalue_PY <= CUTOFF | pvalue_SNDA <= CUTOFF) %>%
  arrange(desc(OR_all)) %>%
  write.table("../results/Merge_circexplorer_BC109.all+cellspecific.combined.DisGeNet+GWAS.SIG.Bonferroni_0.05.xls", sep = "\t", row.names = F)


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

## Clemens: for the 306 super-genes (which host cell-specific circRNAs in three cell types), how many private circRNAs they host in each cell type?
DF3 %>% arrange(celltype, hostgene)%>% group_by(hostgene) %>% summarize(celltype3_n = n_distinct(celltype)) %>% filter(celltype3_n==3) %>% #head() # dim()
  inner_join(DF3, by='hostgene') %>% # pull(hostgene) %>% unique() %>% length()
  group_by(celltype) %>% summarise(n=n()) %>% ggplot(aes(y=1,x=1:3,label=paste(celltype,":",n), size=n)) + 
  geom_point() + scale_size(range=c(5,10)) + geom_text(size=2, vjust = 0, nudge_y = 0.2) + ylim(c(0.5,1.5)) + theme_classic() + ggsave("../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.306supergene_privateCircRNAs.pdf")

# find a key illustrating example of one locus that produces different circRNA isoforms in the 3 cell types.We want to pick a locus that is ideally related to endocytosis AND to PD or AD.e.g. RIMS1, RIMS2, DNAJC6, PICALM, VPS13C could be options.
DF3 %>% arrange(celltype, hostgene)%>% group_by(hostgene) %>% summarize(celltype3_n = n_distinct(celltype)) %>% filter(celltype3_n==3) %>% #head() # dim()
  inner_join(DF3, by='hostgene') %>% filter(hostgene %in% c('AAK1','RIMS1', 'RIMS2', 'DNAJC6', 'PICALM', 'VPS13C','ERC1')) %>% print(n=Inf)
# hostgene celltype3_n SNDA_spec PY_spec NN_spec gene                       S celltype     mean     m2sd Private_or_not hostgeneID        
# <chr>          <int>     <dbl>   <dbl>   <dbl> <chr>                  <dbl> <chr>       <dbl>    <dbl>          <dbl> <fct>             
#   1 AAK1               3     1       0       0     chr2_69746085_69752244 1     SNDA     0.000167 0.000152              1 ENSG00000115977.14
# 2 AAK1               3     0       1       0     chr2_69741602_69754451 1     PY       0.000686 0.000625              1 ENSG00000115977.14
# 3 AAK1               3     0       1       0     chr2_69723116_69736592 1     PY       0.00107  0.000970              1 ENSG00000115977.14
# 4 AAK1               3     0       1       0     chr2_69723116_69754451 1     PY       0.00165  0.00150               1 ENSG00000115977.14
# 5 AAK1               3     0.322   0.564   0     chr2_69757139_69759294 0.564 PY       0.00189  0.00182               1 ENSG00000115977.14
# 6 AAK1               3     0       0       1     chr2_69732700_69736592 1     NN       0.00334  0.00304               1 ENSG00000115977.14
# 7 AAK1               3     0       0       1     chr2_69736362_69754451 1     NN       0.0145   0.0132                1 ENSG00000115977.14
# 8 ERC1               3     1       0       0     chr12_1289705_1346070  1     SNDA     0.000566 0.000515              1 ENSG00000082805.15
# 9 ERC1               3     1       0       0     chr12_1372199_1399178  1     SNDA     0.000482 0.000439              1 ENSG00000082805.15
# 10 ERC1               3     1       0       0     chr12_1225031_1346070  1     SNDA     0.000311 0.000283              1 ENSG00000082805.15
# 11 ERC1               3     1       0       0     chr12_1136913_1481143  1     SNDA     0.000553 0.000503              1 ENSG00000082805.15
# 12 ERC1               3     1       0       0     chr12_1299024_1399178  1     SNDA     0.000694 0.000632              1 ENSG00000082805.15
# 13 ERC1               3     1       0       0     chr12_1250785_1299218  1     SNDA     0.00317  0.00289               1 ENSG00000082805.15
# 14 ERC1               3     1       0       0     chr12_1372199_1553916  1     SNDA     0.00172  0.00157               1 ENSG00000082805.15
# 15 ERC1               3     1       0       0     chr12_1372199_1481143  1     SNDA     0.000555 0.000505              1 ENSG00000082805.15
# 16 ERC1               3     0.272   0.615   0     chr12_1399017_1481143  0.615 PY       0.00203  0.00190               1 ENSG00000082805.15
# 17 ERC1               3     0.350   0.535   0     chr12_1289705_1399178  0.535 PY       0.00113  0.00113               1 ENSG00000082805.15
# 18 ERC1               3     0       1       0     chr12_1219357_1225199  1     PY       0.00639  0.00582               1 ENSG00000082805.15
# 19 ERC1               3     0       1       0     chr12_1399017_1399178  1     PY       0.000389 0.000354              1 ENSG00000082805.15
# 20 ERC1               3     0.165   0.730   0     chr12_1399017_1519619  0.730 PY       0.00790  0.00722               1 ENSG00000082805.15
# 21 ERC1               3     0       1       0     chr12_1399017_1517413  1     PY       0.00121  0.00110               1 ENSG00000082805.15
# 22 ERC1               3     0       0       1     chr12_1213915_1225199  1     NN       0.00167  0.00152               1 ENSG00000082805.15
# 23 ERC1               3     0       0.273   0.614 chr12_1225031_1299218  0.614 NN       0.00694  0.00642               1 ENSG00000082805.15
# 24 ERC1               3     0       0       1     chr12_1192329_1192746  1     NN       0.00463  0.00422               1 ENSG00000082805.15

## any in PD familiar genes (see Table 1 of Cornelis et al. Lancet Neurology 2020)
DF3 %>% filter(celltype == "SNDA", hostgene %in% c('SNCA', 'PRKN', 'UCHL1', 'PARK7', 'LRRK2', 'PINK1', 'POLG', 'HTRA2', 'ATP13A2', 'FBXO7', 'GIGYF2', 'GBA', 'PLA2G6', 'EIF4G1', 'VPS35', 'DNAJC6', 'SYNJ1', 'DNAJC13', 'TMEM230', 'VPS13C', 'LRP10'))
#                          SNDA_spec PY_spec NN_spec                     gene S celltype         mean         m2sd Private_or_not hostgene         hostgeneID
# chr22_32874967_32881196          1       0       0  chr22_32874967_32881196 1     SNDA 0.0010858469 0.0009888630              1    FBXO7 ENSG00000100225.13
# chr2_233612324_233626146         1       0       0 chr2_233612324_233626146 1     SNDA 0.0005296589 0.0004823517              1   GIGYF2 ENSG00000204120.10
# chr1_65830317_65871816           1       0       0   chr1_65830317_65871816 1     SNDA 0.0005846243 0.0005324077              1   DNAJC6 ENSG00000116675.11

# circRNAs in PD genes in all three cell types
groupmean_s3 %>% filter(hostgene %in% PDgenes)

DF3 %>% arrange(celltype, hostgene)%>% group_by(hostgene) %>% summarize(celltype3_n = n_distinct(celltype)) %>% filter(celltype3_n==3) %>% #head() # dim()
  inner_join(DF3, by='hostgene') %>% filter(hostgene %in% PDgenes) %>% print(n=Inf)
  
# circRNAs in  leukemia or adenoCA
groupmean_s3 %>% filter(celltype == "NN",  Private_or_not==1, hostgene %in% gsd$geneSymbol[gsd$diseaseName %in% c("Leukemia, Myelocytic, Acute", "Adenocarcinoma of large intestine")])

# circRNAs in  addiction
groupmean_s3 %>% filter(celltype == "SNDA",  Private_or_not==1, hostgene %in% gsd$geneSymbol[gsd$diseaseName %in% c("Substance-Related Disorders", "Drug Dependence")])
            
# AD familiar genes
DF3 %>% filter(hostgene %in% c('APP', 'PSEN1', 'PSEN2', 'APOE')) 
#                            SNDA_spec   PY_spec   NN_spec                    gene         S celltype         mean         m2sd Private_or_not hostgene         hostgeneID
# chr14_73614502_73626846 1.000000e+00 0.0000000 0.0000000 chr14_73614502_73626846 1.0000000     SNDA 0.0001829260 0.0001665877              1    PSEN1 ENSG00000080815.14
# chr21_27423315_27425664 1.000000e+00 0.0000000 0.0000000 chr21_27423315_27425664 1.0000000     SNDA 0.0002560425 0.0002331737              1      APP ENSG00000142192.16
# chr21_27394155_27484463 1.000000e+00 0.0000000 0.0000000 chr21_27394155_27484463 1.0000000     SNDA 0.0001387411 0.0001263493              1      APP ENSG00000142192.16
# chr21_27394155_27394358 1.000000e+00 0.0000000 0.0000000 chr21_27394155_27394358 1.0000000     SNDA 0.0008005572 0.0007290543              1      APP ENSG00000142192.16
# chr21_27347382_27354790 2.854724e-01 0.6011852 0.0000000 chr21_27347382_27354790 0.6011852       PY 0.0003955135 0.0003793747              1      APP ENSG00000142192.16
# chr21_27326903_27354790 3.382194e-01 0.5469224 0.0000000 chr21_27326903_27354790 0.5469224       PY 0.0013489531 0.0013253998              1      APP ENSG00000142192.16
# chr21_27269884_27328069 0.000000e+00 1.0000000 0.0000000 chr21_27269884_27328069 1.0000000       PY 0.0025242630 0.0022988050              1      APP ENSG00000142192.16
# chr21_27347382_27372497 0.000000e+00 1.0000000 0.0000000 chr21_27347382_27372497 1.0000000       PY 0.0003610881 0.0003288370              1      APP ENSG00000142192.16
# chr14_73614502_73673180 0.000000e+00 1.0000000 0.0000000 chr14_73614502_73673180 1.0000000       PY 0.0006917505 0.0006299659              1    PSEN1 ENSG00000080815.14
# chr21_27394155_27425664 0.000000e+00 1.0000000 0.0000000 chr21_27394155_27425664 1.0000000       PY 0.0029208270 0.0026599492              1      APP ENSG00000142192.16
# chr21_27326903_27372497 1.736703e-01 0.7209670 0.0000000 chr21_27326903_27372497 0.7209670       PY 0.0026699587 0.0024481137              1      APP ENSG00000142192.16
# chr14_73614502_73664837 1.301553e-01 0.1961118 0.6093113 chr14_73614502_73664837 0.6093113       NN 0.0075090841 0.0068866944              1    PSEN1 ENSG00000080815.14
# chr14_73614502_73653628 0.000000e+00 0.0000000 1.0000000 chr14_73614502_73653628 1.0000000       NN 0.0016686853 0.0015196444              1    PSEN1 ENSG00000080815.14
# chr14_73614502_73640415 1.110223e-16 0.3721196 0.5124613 chr14_73614502_73640415 0.5124613       NN 0.0520719727 0.0488625855              1    PSEN1 ENSG00000080815.14
# chr14_73614502_73659572 0.000000e+00 0.0000000 1.0000000 chr14_73614502_73659572 1.0000000       NN 0.0060097613 0.0054729911              1    PSEN1 ENSG00000080815.14

# Schizophrenia genes (https://medlineplus.gov/genetics/condition/schizophrenia/#causes)
DF3 %>% filter(hostgene %in% c('AKT1', 'COMT', 'YWHAE', 'ABCA13', 'C4A', 'DGCR2', 'DGCR8', 'DRD2', 'MIR137', 'NOS1AP', 'NRXN1', 'OLIG2', 'RTN4R', 'SYN2', 'TOP3B', 'ZDHHC8')) 
#                        SNDA_spec   PY_spec NN_spec                   gene         S celltype         mean         m2sd Private_or_not hostgene         hostgeneID
# chr2_50692579_50733755 1.0000000 0.0000000 0.0000000 chr2_50692579_50733755 1.0000000     SNDA 1.857447e-04 1.691546e-04              1    NRXN1 ENSG00000179915.16
# chr2_50723042_50765774 1.0000000 0.0000000 0.0000000 chr2_50723042_50765774 1.0000000     SNDA 1.877574e-04 1.709876e-04              1    NRXN1 ENSG00000179915.16
# chr2_50692579_50780163 1.0000000 0.0000000 0.0000000 chr2_50692579_50780163 1.0000000     SNDA 3.833756e-04 3.491339e-04              1    NRXN1 ENSG00000179915.16
# chr2_50758364_50847321 1.0000000 0.0000000 0.0000000 chr2_50758364_50847321 1.0000000     SNDA 9.287234e-05 8.457732e-05              1    NRXN1 ENSG00000179915.16
# chr2_50692579_50847321 1.0000000 0.0000000 0.0000000 chr2_50692579_50847321 1.0000000     SNDA 2.668736e-03 2.430374e-03              1    NRXN1 ENSG00000179915.16
# chr2_50280408_50733755 1.0000000 0.0000000 0.0000000 chr2_50280408_50733755 1.0000000     SNDA 1.060800e-03 9.660530e-04              1    NRXN1 ENSG00000179915.16
# chr2_50692579_50765774 0.2786282 0.6082995 0.0000000 chr2_50692579_50765774 0.6082995       PY 1.277369e-02 1.176156e-02              1    NRXN1 ENSG00000179915.16
# chr2_51149019_51149795 0.3533957 0.5314616 0.0000000 chr2_51149019_51149795 0.5314616       PY 9.021483e-03 8.602454e-03              1    NRXN1 ENSG00000179915.16
# chr3_12182150_12211406 0.0000000 1.0000000 0.0000000 chr3_12182150_12211406 1.0000000       PY 7.653028e-04 6.969487e-04              1     SYN2 ENSG00000157152.12
# chr2_51149021_51149795 0.0000000 1.0000000 0.0000000 chr2_51149021_51149795 1.0000000       PY 2.099061e-03 1.911580e-03              1    NRXN1 ENSG00000179915.16
# chr17_1257504_1265302  0.0000000 1.0000000 0.0000000  chr17_1257504_1265302 1.0000000       PY 2.566912e-03 2.337644e-03              1    YWHAE ENSG00000108953.12
# chr2_50779724_50850753 0.0000000 1.0000000 0.0000000 chr2_50779724_50850753 1.0000000       PY 2.145934e-03 1.954267e-03              1    NRXN1 ENSG00000179915.16
# chr2_50692579_50850753 0.0000000 1.0000000 0.0000000 chr2_50692579_50850753 1.0000000       PY 1.031438e-03 9.393134e-04              1    NRXN1 ENSG00000179915.16
# chr2_50723042_50758568 0.0000000 1.0000000 0.0000000 chr2_50723042_50758568 1.0000000       PY 8.097933e-04 7.374654e-04              1    NRXN1 ENSG00000179915.16
# chr2_50723042_50850753 0.0000000 1.0000000 0.0000000 chr2_50723042_50850753 1.0000000       PY 3.817297e-03 3.476350e-03              1    NRXN1 ENSG00000179915.16
# chr3_12182150_12224872 0.0000000 1.0000000 0.0000000 chr3_12182150_12224872 1.0000000       PY 7.668848e-04 6.983894e-04              1     SYN2 ENSG00000157152.12
# chr2_50699435_50847321 0.1009630 0.8067966 0.0000000 chr2_50699435_50847321 0.8067966       PY 2.105432e-03 1.922437e-03              1    NRXN1 ENSG00000179915.16
# chr2_50692579_50758568 0.2196339 0.6706193 0.0000000 chr2_50692579_50758568 0.6706193       PY 1.447698e-02 1.324333e-02              1    NRXN1 ENSG00000179915.16
# chr2_50318460_50733755 0.2917177 0.5947101 0.0000000 chr2_50318460_50733755 0.5947101       PY 2.755785e-03 2.589169e-03              1    NRXN1 ENSG00000179915.16
# chr2_50463926_50733755 0.1447586 0.7539036 0.0000000 chr2_50463926_50733755 0.7539036       PY 8.779922e-03 8.014362e-03              1    NRXN1 ENSG00000179915.16
# chr2_50847159_50850753 0.2763917 0.6106287 0.0000000 chr2_50847159_50850753 0.6106287       PY 9.161212e-03 8.453001e-03              1    NRXN1 ENSG00000179915.16
# chr17_1264385_1268352  0.0000000 0.2736297 0.6135081  chr17_1264385_1268352 0.6135081       NN 1.531227e-03 1.439078e-03              1    YWHAE ENSG00000108953.12
# chr17_1268152_1273035  0.0000000 0.0000000 1.0000000  chr17_1268152_1273035 1.0000000       NN 3.337371e-03 3.039289e-03              1    YWHAE ENSG00000108953.12
# chr17_1264385_1265302  0.0000000 0.0000000 1.0000000  chr17_1264385_1265302 1.0000000       NN 1.668685e-03 1.519644e-03              1    YWHAE ENSG00000108953.12

# ADHD (https://www.omim.org/entry/143465)
DF3 %>% filter(hostgene %in% c('DRD5', 'ADHD4', 'ADHD3', 'DRD4', 'ADHD1', 'ADHD2', 'DAT1', 'SLC6A3', 'HTR1B', 'ADRA2A', 'SCN8A', 'SNAP25', 'COMT'))
# SNDA_spec PY_spec NN_spec                    gene S celltype         mean         m2sd Private_or_not hostgene        hostgeneID
# chr12_52139686_52145377         0       1       0 chr12_52139686_52145377 1       PY 0.0007748661 0.0007056578              1    SCN8A ENSG00000196876.9
# chr12_52077957_52082880         0       1       0 chr12_52077957_52082880 1       PY 0.0006587036 0.0005998706              1    SCN8A ENSG00000196876.9
# chr12_52180325_52188425         0       0       1 chr12_52180325_52188425 1       NN 0.0018380480 0.0016738801              1    SCN8A ENSG00000196876.9

# addiction (Table S3 of https://www.medrxiv.org/content/10.1101/2022.01.06.22268753v1.full-text)
DF3 %>% filter(hostgene %in% c("PDE4B", "GTF3C2", "IFT172", "GCKR", "C2orf16", "PLCL2", "PRKAR2A", "SLC25A20", "ARIH2", "P4HTM", "WDR6", "QRICH1", 
                               "QARS", "USP19", "C3orf84", "CCDC36", "RP11-3B7.1", "USP4", "RHOA", "NICN1", "DAG1", "ADD1", "NOP14", "BANK1", "PPP6C", 
                               "SOX6", "BDNF", "MTCH2", "ZDHHC5", "TMX2", "TMX2-CTNND1", "C11orf31", "RP11-691N7.6", "CTNND1", "ANKK1", "DRD2", "HS6ST3", 
                               "ARID4A", "PPP1R13B", "SEMA6D", "FTO", "C20orf112",
                               "ALDH2","ADH1B", "ADH1C","ADH4", "GCKR", "SIX3","SLC39A8" ,"DRD4", "CHRNA5", "COMT","OPRM1"
))

# Q: we need to show that circRNAs from some neuropsych disease loci show cell-type specific expression OR cell type specific isoforms across the 3 cell types
sum(DF3$hostgene %in% ADgenes)

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
# Weighted sankey / alluvial diagram from cell-specific circRNA to their host genes
#############################################
DF3=readRDS("Merge_circexplorer_BC109.cellspecific_heatmap.circRNA3.rds")

# generate data.frame in ggalluvial format, see ref below
# https://stackoverflow.com/a/48133004

d=DF3 %>% select(celltype, hostgeneID) %>% mutate(celltype=factor(celltype, levels = c("SNDA","PY","NN")))
d1=group_by(d, celltype, hostgeneID) %>% summarise(n_circ=n())
d2=group_by(d1, hostgeneID) %>% summarise(celltype_gene=paste(unique(celltype), collapse = ',')) 
dx=inner_join(d1,d2,by="hostgeneID") %>% rename(celltype_circ=celltype) %>% 
  group_by(celltype_circ, celltype_gene) %>%
  summarise(N_circ=sum(n_circ), N_gene=n_distinct(hostgeneID)) %>%
  tibble::rownames_to_column() %>%
  pivot_longer(cols = contains('_'), names_to =c(".value", "source"), names_sep = "_") 
  
library(ggalluvial)
ggplot(
  data = dx,
  aes(
    x = factor(source, levels = c("gene", "circ")),
    stratum = factor(celltype, levels=c("SNDA","PY","NN","SNDA,PY,NN", "SNDA,PY","SNDA,NN","PY,NN")),
    alluvium = rowname,
    y = N
  )
) +
  geom_stratum(aes(fill = factor(celltype, levels=c("SNDA","PY","NN","SNDA,PY,NN", "SNDA,PY","SNDA,NN","PY,NN")))) +
  geom_flow() +
  guides(fill=guide_legend(title="celltype")) + theme_classic()
ggsave("../results/Merge_circexplorer_BC109.cellspecific_circRNA2hostgene.alluvial.pdf", height = 5, width = 8)



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
