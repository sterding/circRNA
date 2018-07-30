###########################################
# R script to get cell type specific RNAs/genes
# Usage: Rscript $0 Merge_circexplorer_BC.filtered.normRPM.rds
# Input: Merge_circexplorer_BC.filtered.normRPM.rds
# Output: ../results/Merge_circexplorer_BC.annotation_per_cell.xls, ../results/Merge_circexplorer_BC.cellspecific_heatmap.pdf, and ../results/Merge_circexplorer_BC.cellspecific_heatmap+.pdf
# Author: Xianjun Dong
# Date: 2017-10-28
###########################################
library(dplyr)
library(tidyr)

args<-commandArgs(TRUE)
normRPM_file=ifelse(is.na(args[1]),"Merge_circexplorer_BC.filtered.enriched.normRPM.rds",args[1])  
# normRPM_file = "Merge_circexplorer_BC.filtered.enriched.normRPM.rds"

setwd("~/projects/circRNA/data")
# calculate specificity score (S) as cummeRbund (https://github.com/Bioconductor-mirror/cummeRbund/blob/master/R/tools.R)
source("~/projects/circRNA/src/annotation/tools.R")

Merge_circexp_norm_filtered_and_enriched=readRDS(normRPM_file); dim(Merge_circexp_norm_filtered_and_enriched)

###########################################
## circRNA expressed in each cell #########
###########################################
#definition: “Expressed” circular RNAs: expressed overall (e.g. >=2 unique reads overall) and in at least 1 sample in a specific cell group 
df = Merge_circexp_norm_filtered_and_enriched
df = df %>% mutate(gene=rownames(df)) %>% melt(id.vars="gene", variable.name = "sampleID", value.name='fpkm') %>%
  separate(sampleID, c("Dx", "subjectID","cellType5","batch","rep"), sep = "_", remove=FALSE) %>% 
  mutate(cellType3=ifelse(cellType5 %in% c("TCPY","MCPY"), "PY", ifelse(cellType5=="SNDA","SNDA","NN"))) %>% 
  filter(fpkm>0) %>% 
  select(gene, cellType3, cellType5) %>% arrange(gene, cellType3, cellType5) %>% distinct()

df %>% group_by(gene) %>% summarize(celltype3 = paste(unique(cellType3), collapse = ','), celltype5 = paste(unique(cellType5), collapse = ',')) %>% ungroup() %>% merge(annotation, by.x='gene',by.y='ID', all.x=T, all.y=F) %>% write.table(file="../results/Merge_circexplorer_BC.annotation_per_cell.xls", sep="\t", quote=F, row.names=F)

df %>% group_by(gene) %>% summarize(celltype3 = paste(unique(cellType3), collapse = ','), celltype5 = paste(unique(cellType5), collapse = ',')) %>% group_by(celltype3) %>% summarise(n=n())
# # A tibble: 7 × 2
# celltype3     n
# <chr> <int>
# 1         NN  3119
# 2      NN,PY   516
# 3 NN,PY,SNDA  1236
# 4    NN,SNDA  1450
# 5         PY   963
# 6    PY,SNDA   794
# 7       SNDA  1936

## call eulerAPE to generate Merge_circexplorer_BC.venn_celltype3.pdf


clip <- pipe("pbcopy", "w");
df %>% group_by(gene) %>% summarize(celltype3 = paste(unique(cellType3), collapse = ','), celltype5 = paste(unique(cellType5), collapse = ',')) %>% group_by(celltype5) %>% summarise(n=n()) %>% as.data.frame() %>% write.table(file=clip, sep="\t", row.names=F)
close(clip)

df %>% merge(annotation, by.x='gene',by.y='ID', all.x=T, all.y=F) %>% select(gene, cellType3, cellType5, circType) %>% group_by(cellType5, circType) %>% summarise(n=n_distinct(gene)) %>% dcast(cellType5 ~ circType, value.var="n") #%>% write.table(file=clip, sep="\t", row.names=F)
# cellType5 circRNA
# 1        FB    2090
# 2      MCPY    1197
# 3      PBMC    5569
# 4      SNDA    5416
# 5      TCPY    2884
df %>% merge(annotation, by.x='gene',by.y='ID', all.x=T, all.y=F) %>% select(gene, cellType3, cellType5, circType) %>% group_by(cellType3, circType) %>% summarise(n=n_distinct(gene)) %>% dcast(cellType3 ~ circType, value.var="n") #%>% write.table(file=clip, sep="\t", row.names=F)
# cellType3 circRNA
# 1        NN    6321
# 2        PY    3509
# 3      SNDA    5416


# 5-ways venn diagram (not area-proportional)
library(gplots)
pdf("../results/Merge_circexplorer_BC.venn_celltype5.pdf")
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

groupmean = group3mean
rownames(groupmean)=groupmean[,1]; groupmean=groupmean[,-1]
# overall mean and sd
overmean = groupmean %>% mutate(gene=rownames(groupmean)) %>% melt(id.vars="gene", variable.name = "sampleID", value.name='fpkm') %>% 
  group_by(gene) %>%
  summarise(overall.mean=mean(fpkm), overall.sd=sd(fpkm)) %>% data.frame()
# overmean = df %>% mutate(gene=rownames(df)) %>% melt(id.vars="gene", variable.name = "sampleID", value.name='fpkm') %>% 
#   group_by(gene) %>%
#   summarise(overall.mean=mean(fpkm), overall.sd=sd(fpkm)) %>% data.frame()
rownames(overmean)=overmean[,1]; overmean=overmean[,-1]
groupmean_s = specificity(groupmean*1000, logMode = T)
celltypes=c("SNDA","PY","NN"); #
df = data.frame(gene=rownames(groupmean_s), S=apply(groupmean_s,1,max), celltype=apply(groupmean_s,1,which.max))
df = df %>% mutate(mean=groupmean[cbind(1:nrow(df), celltype)], m2sd=with(overmean,overall.mean+1*overall.sd)) %>% mutate(celltype=celltypes[celltype])
df = df %>% filter(S>=0.5, mean>m2sd) 
df %>% group_by(celltype) %>% summarise(count=n())
# # A tibble: 3 × 2
# celltype count
# <chr> <int>
# 1       NN  4860
# 2       PY  1850
# 3     SNDA  2022
df3=df

groupmean = group5mean
rownames(groupmean)=groupmean[,1]; groupmean=groupmean[,-1]

overmean = groupmean %>% mutate(gene=rownames(groupmean)) %>% melt(id.vars="gene", variable.name = "sampleID", value.name='fpkm') %>% 
  group_by(gene) %>%
  summarise(overall.mean=mean(fpkm), overall.sd=sd(fpkm)) %>% data.frame()
rownames(overmean)=overmean[,1]; overmean=overmean[,-1]

# specificty score
groupmean_s = specificity(groupmean*1000, logMode = T)
# Being specific: specificity score S>=0.5 AND mean expression > mean+2s.d. of overall expression (Zheng et al., doi:10.1038/ncomms11215)
celltypes=c("SNDA","MCPY","TCPY","FB","PBMC")
df = data.frame(gene=rownames(groupmean_s), S=apply(groupmean_s,1,max), celltype=apply(groupmean_s,1,which.max))
df = df %>% mutate(mean=groupmean[cbind(1:nrow(df), celltype)], m2sd=with(overmean,overall.mean+1.5*overall.sd)) %>% mutate(celltype=celltypes[celltype])
df = df %>% filter(S>=0.5, mean>m2sd) 
df %>% group_by(celltype) %>% summarise(count=n())
# # A tibble: 5 × 2
# celltype count
# <chr> <int>
# 1       FB   795
# 2     MCPY   459
# 3     PBMC  3668
# 4     SNDA  1989
# 5     TCPY  1290
head(df)
df5=df

# heatmpa of cell specific genes
library('pheatmap')
df=log10(groupmean[df5$gene,]*1000+1)

dissimilarity <- 1 - cor(df)
distance <- as.dist(dissimilarity)
hc=hclust(distance)
hccol=as.dendrogram(hc)
weights.dd <- order(match(c("PBMC","FB","SNDA","TCPY","MCPY"), celltypes)) # the trick is to call order() on the specific index of target dendragram
plot(reorder(hccol, wts = weights.dd, agglo.FUN = mean))

hm.parameters <- list(df, 
                      color = colorRampPalette(c("green", "black", "red"))(50),
                      scale = "row",
                      treeheight_row = 50,
                      treeheight_col = 30,
                      kmeans_k = NA,
                      show_rownames = F, show_colnames = T,
                      clustering_method = "average",
                      cluster_rows = TRUE, cluster_cols = as.hclust(reorder(hccol, wts = weights.dd, agglo.FUN = mean)), 
                      clustering_distance_rows = 'correlation', 
                      clustering_distance_cols = 'correlation')
do.call("pheatmap", hm.parameters)
do.call("pheatmap", c(hm.parameters, width=3, height=5, filename="../results/Merge_circexplorer_BC.cellspecific_heatmap.pdf"))

## cell-specific circRNAs in both major and minor groups
df=log10(groupmean[unique(df3$gene,df5$gene),]*1000+1)
dissimilarity <- 1 - cor(df)
distance <- as.dist(dissimilarity)
hc=hclust(distance)
hccol=as.dendrogram(hc)
weights.dd <- order(match(c("PBMC","FB","SNDA","MCPY","TCPY"), celltypes)) # the trick is to call order() on the specific index of target dendragram
plot(reorder(hccol, wts = weights.dd))

hm.parameters <- list(df, 
                      color = colorRampPalette(c("green", "black", "red"))(50),
                      scale = "row",
                      treeheight_row = 50,
                      treeheight_col = 30,
                      kmeans_k = NA,
                      show_rownames = F, show_colnames = T,
                      clustering_method = "average",
                      cluster_rows = TRUE, cluster_cols = as.hclust(reorder(hccol, wts = weights.dd)), 
                      clustering_distance_rows = 'correlation', 
                      clustering_distance_cols = 'correlation')
do.call("pheatmap", hm.parameters)
do.call("pheatmap", c(hm.parameters, width=3, height=5, filename="../results/Merge_circexplorer_BC.cellspecific_heatmap+.pdf"))
