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
#df = df %>% mutate(gene=rownames(df)) %>% melt(id.vars="gene", variable.name = "sampleID", value.name='fpkm') %>%
df = df %>% mutate(gene=rownames(df)) %>% gather(key = "sampleID", value='fpkm', -gene) %>%
  separate(sampleID, c("Dx", "subjectID","cellType5","batch","rep"), sep = "_", remove=FALSE) %>% 
  mutate(cellType3=ifelse(cellType5 %in% c("TCPY","MCPY"), "PY", ifelse(cellType5=="SNDA","SNDA","NN"))) %>% 
  filter(fpkm>0) %>% 
  select(gene, cellType3, cellType5) %>% arrange(gene, cellType3, cellType5) %>% distinct()

## In individual circType
df %>% merge(annotation, by.x='gene',by.y='ID', all.x=T, all.y=F) %>% select(gene, cellType3, cellType5, circType) %>% group_by(cellType5, circType) %>% summarise(n=n_distinct(gene)) %>% dcast(cellType5 ~ circType, value.var="n") 
#    cellType5 circRNA ciRNA
# 1        FB    1839   252
# 2      MCPY    1185    14
# 3      PBMC    5420   150
# 4      SNDA    5330    89
# 5      TCPY    2854    31
df %>% merge(annotation, by.x='gene',by.y='ID', all.x=T, all.y=F) %>% select(gene, cellType3, cellType5, circType) %>% group_by(cellType3, circType) %>% summarise(n=n_distinct(gene)) %>% dcast(cellType3 ~ circType, value.var="n") #%>% write.table(file=clip, sep="\t", row.names=F)
#     cellType3 circRNA ciRNA
# 1        NN    6037   285
# 2        PY    3470    41
# 3      SNDA    5330    89

## annotation per cell type
df %>% group_by(gene) %>% summarize(celltype3 = paste(unique(cellType3), collapse = ','), celltype5 = paste(unique(cellType5), collapse = ',')) %>% ungroup() %>% merge(annotation, by.x='gene',by.y='ID', all.x=T, all.y=F) %>% write.table(file="../results/Merge_circexplorer_BC.annotation_per_cell.xls", sep="\t", quote=F, row.names=F)
#                       gene  celltype3      celltype5 chrom     start       end score strand thickStart  thickEnd itemRgb exonCount                   exonSizes                         exonOffsets circType geneName
# 1   chr1_10032075_10041228 NN,PY,SNDA PBMC,MCPY,SNDA  chr1  10032075  10041228     1      +   10032075  10032075   0,0,0         3                 171,184,140                         0,3574,9013  circRNA   NMNAT1
# 2 chr1_100340242_100343384       SNDA           SNDA  chr1 100340242 100343384     1      +  100340242 100340242   0,0,0         5          124,103,98,140,188                 0,467,671,1771,2954  circRNA      AGL
# 3 chr1_100342013_100347247         PY           TCPY  chr1 100342013 100347247     1      +  100342013 100342013   0,0,0         7 140,188,124,164,102,156,151     0,1183,3465,4174,4618,4834,5083  circRNA      AGL
# 4 chr1_100515464_100535241    NN,SNDA      PBMC,SNDA  chr1 100515464 100535241     1      +  100515464 100515464   0,0,0         7    96,63,128,115,219,114,72 0,8757,9969,11926,18068,18564,19705  circRNA    HIAT1
# 5 chr1_100889777_100908552         NN           PBMC  chr1 100889777 100908552     1      +  100889777 100889777   0,0,0         3                    80,67,63                       0,15710,18712  circRNA   CDC14A
# 6 chr1_100949847_100964818         NN           PBMC  chr1 100949847 100964818     1      +  100949847 100949847   0,0,0         5          160,113,48,123,334           0,10526,11710,13793,14637  circRNA   CDC14A

celltype3_n = df %>% group_by(gene) %>% summarize(celltype3 = paste(unique(cellType3), collapse = ','), celltype5 = paste(unique(cellType5), collapse = ',')) %>% group_by(celltype3) %>% summarise(n=n())
# # A tibble: 7 × 2
# celltype3     n
# <chr> <int>
# 1         NN  3119
# 2      NN,PY   516
# 3 NN,PY,SNDA  1237
# 4    NN,SNDA  1450
# 5         PY   963
# 6    PY,SNDA   795
# 7       SNDA  1937

## call eulerAPE to generate Merge_circexplorer_BC.venn_celltype3.pdf

# 5-ways venn diagram (not area-proportional)
library(gplots); # install.packages('gplots')
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

groupmean_to_specificity <- function(groupmean){
  # The input groupmean is a matrix of FPKM, where the first column as gene names and the rest columns as sampleID
  rownames(groupmean)=groupmean[,1]; groupmean=groupmean[,-1]
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

groupmean_s3 = groupmean_to_specificity(group3mean)
df3 = groupmean_s3 %>% filter(Private_or_not==1) 
group_by(df3, celltype) %>% summarise(count=n())
# # A tibble: 3 × 2
# celltype count
# <chr> <int>
# 1       NN  4860
# 2       PY  1850
# 3     SNDA  2023

groupmean_s5 = groupmean_to_specificity(group5mean)
df5 = groupmean_s5 %>% filter(Private_or_not==1) 
group_by(df5, celltype) %>% summarise(count=n())
# # A tibble: 5 × 2
# celltype count
# <chr> <int>
# 1 FB         795
# 2 MCPY       461
# 3 PBMC      3672
# 4 SNDA      1991
# 5 TCPY      1295

## df3 and df5 are very similar
rbind(data.frame(gene=df3$gene, DF="df3"), data.frame(gene=df5$gene, DF="df5")) %>% group_by(gene) %>% summarize(DFs = paste(unique(DF), collapse = ',')) %>% group_by(DFs) %>% summarise(n=n())
#   DFs         n
# 1 df3       760
# 2 df3,df5  7973
# 3 df5       230

###########################################
# heatmpa of cell specific circRNAs
###########################################
library('pheatmap')

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
dissimilarity <- 1 - cor(df)
distance <- as.dist(dissimilarity)
hc=hclust(distance)
hccol=as.dendrogram(hc)
weights.dd <- order(match(c("PBMC","FB","SNDA","MCPY","TCPY"), celltypes)) # the trick is to call order() on the specific index of target dendragram
#plot(reorder(hccol, wts = weights.dd))

dissimilarity <- 1 - cor(t(df))
distance <- as.dist(dissimilarity)
hcrow=hclust(distance, method = 'average')
roworders = hcrow$order # order of the new rows in the old matrix
library(RColorBrewer) 
hm.parameters <- list(df, 
                      color = colorRampPalette(c("#999999","#ffffff", "#ca0020"))(100),
                      scale = "row",
                      treeheight_row = 50,
                      treeheight_col = 30,
                      kmeans_k = NA,
                      show_rownames = F, show_colnames = T,
                      clustering_method = "average",
                      cluster_rows = hcrow, cluster_cols = as.hclust(reorder(hccol, wts = weights.dd, agglo.FUN = mean)), 
                      clustering_distance_rows = 'correlation', 
                      clustering_distance_cols = 'correlation')
hm=do.call("pheatmap", hm.parameters)
do.call("pheatmap", c(hm.parameters, width=3, height=5, filename="../results/Merge_circexplorer_BC.cellspecific_heatmap++.pdf"))

###########################################
# host genes of cell-specific circRNAs
###########################################
# circRNAs in the order of presence in above figure
df_roworders = df3$gene[roworders]; head(df_roworders); length(df_roworders)  # n = 8733
DF3=groupmean_s3[df_roworders,]; dim(DF3); head(DF3); table(DF3$Private_or_not)

gene_fpkm=read.table("~/neurogen/rnaseq_PD/results/merged/genes.fpkm.cuffnorm.allSamples.uniq.xls", header = T, row.names = 1, stringsAsFactors = F, check.names = F); head(gene_fpkm)
gene_fpkm=gene_fpkm[,colnames(Merge_circexp_norm_filtered_and_enriched)]; dim(gene_fpkm)
gene_group5mean = gene_fpkm %>% rownames_to_column(var = 'gene') %>% melt(id.vars="gene", variable.name = "sampleID", value.name='fpkm') %>%
  separate(sampleID, c("Dx", "subjectID","cellType","batch","rep"), sep = "_", remove=FALSE) %>%
  group_by(gene, cellType) %>%
  summarise(meanFPKM=mean(fpkm)) %>% 
  mutate(cellType=factor(cellType, levels = c("SNDA","MCPY","TCPY","FB","PBMC"))) %>% # redefine factor of cellType
  dcast(gene ~ cellType, value.var='meanFPKM') 
head(gene_group5mean)

gene_group3mean = gene_fpkm %>% rownames_to_column(var = 'gene') %>% melt(id.vars="gene", variable.name = "sampleID", value.name='fpkm') %>%
  separate(sampleID, c("Dx", "subjectID","cellType","batch","rep"), sep = "_", remove=FALSE) %>%
  mutate(cellType=ifelse(cellType %in% c("TCPY","MCPY"), "PY", ifelse(cellType=="SNDA","SNDA","NN"))) %>%
  group_by(gene, cellType) %>%
  summarise(meanFPKM=mean(fpkm)) %>% 
  mutate(cellType=factor(cellType, levels = c("SNDA","PY","NN"))) %>% # redefine factor of cellType
  dcast(gene ~ cellType, value.var='meanFPKM') 
head(gene_group3mean)

filtered_enriched_annotation=readRDS("Merge_circexplorer_BC106.filtered.enriched.annotation.bed14.rds"); head(filtered_enriched_annotation)
genes=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed", header = F, col.names = c("chr","start","end","EnsID","score","strand","geneName","geneType"), stringsAsFactors = F)

DF3$hostgene=filtered_enriched_annotation$geneName[match(DF3$gene, filtered_enriched_annotation$ID)]
DF3$hostgeneID=genes$EnsID[match(DF3$hostgene, genes$geneName)]  
# around 40 genes are discarded, as not found in the gencode table  (Note: in the future, circExplorer should use the same gene annotation as  GENCODE)
dim(DF3); DF3=filter(DF3, !is.na(hostgeneID)); dim(DF3) # 8733--> 8622

gene_group3mean_DF3 = gene_group3mean[match(unique(DF3$hostgeneID), gene_group3mean$gene),]
gene_groupmean_s3_DF3 = groupmean_to_specificity(gene_group3mean_DF3)

DF3 = inner_join(DF3, gene_groupmean_s3_DF3,by = c("hostgeneID" = "gene"), suffix=c(".circRNA",".gene")); head(DF3)
DF3 = mutate(DF3, celltype.specific.gene = ifelse(Private_or_not.gene==1, celltype.gene, "NS"))

table(DF3$celltype.circRNA, DF3$celltype.gene)
#       NN   PY SNDA
# NN   4357  148  282
# PY    981  492  356
# SNDA 1095  345  566
table(DF3$celltype.circRNA, DF3$celltype.specific.gene)
#        NN   NS
# NN     10 4777
# PY      0 1829
# SNDA    0 2006

DF3 %>% select(celltype, hostgene) %>% distinct() %>% group_by(celltype) %>% mutate(n=n(), hostgene = paste0(hostgene, collapse = ";")) %>% distinct() %>% write.table(file="../results/Merge_circexplorer_BC.cellspecific_heatmap5.genes.txt", sep = "\t", col.names = T, quote=F, row.names = F)

df5 %>% select(celltype, hostgene) %>% distinct() %>% mutate(celltype3=ifelse(celltype %in% c("TCPY","MCPY"), "PY", ifelse(celltype=="SNDA","SNDA","NN"))) %>% group_by(hostgene) %>% summarize(celltype3 = paste(unique(celltype3), collapse = ',')) %>% group_by(celltype3) %>% summarise(n=n(), gene3 = paste(unique(hostgene), collapse = ';')) %>% write.table(file="../results/Merge_circexplorer_BC.cellspecific_heatmap5.genes3.txt", sep = "\t", col.names = T, quote=F, row.names = F)

# How many host genes per celltype3
df5 %>% mutate(celltype3=ifelse(celltype %in% c("TCPY","MCPY"), "PY", ifelse(celltype=="SNDA","SNDA","NN"))) %>% select(hostgene, celltype3) %>% distinct() %>% group_by(celltype3) %>% summarise(n=n())

## heatmap of host gene for each group of cell-specific circRNAs
df5 = df5 %>% mutate(celltype3=ifelse(celltype %in% c("TCPY","MCPY"), "PY", ifelse(celltype=="SNDA","SNDA","NN")))
for(i in c("PBMC","FB","SNDA","TCPY","MCPY")){
  X=log10(gene_group5mean[unique(as.character(df5$hostgeneID[df5$celltype==i])),]*1000+1)
  # remove NA
  dim(X); X=X[!(rownames(X) =='NA'),]; dim(X)
  dissimilarity <- 1 - cor(X)
  distance <- as.dist(dissimilarity)
  hc=hclust(distance)
  hccol=as.dendrogram(hc)
  weights.dd <- order(match(c("PBMC","FB","SNDA","TCPY","MCPY"), celltypes)) # the trick is to call order() on the specific index of target dendragram
  #plot(reorder(hccol, wts = weights.dd, agglo.FUN = mean))
  
  hm.parameters <- list(X, 
                        color = colorRampPalette(c("#999999","#ffffff", "#ca0020"))(100),
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
  
}




# divide into venn diagram --> ## call eulerAPE to generate Merge_circexplorer_BC.cellspecific_heatmap5.genes3.venn_final.pdf
df5 %>% mutate(celltype3=ifelse(celltype %in% c("TCPY","MCPY"), "PY", ifelse(celltype=="SNDA","SNDA","NN"))) %>% arrange(celltype3, hostgene)%>% group_by(hostgene) %>% summarize(celltype3 = paste(unique(celltype3), collapse = ',')) %>% group_by(celltype3) %>% summarise(n=n(), gene3 = paste(unique(hostgene), collapse = ';'))

## run DAVID to get the enriched KEGG pathway for the SNDA-specific host genes
txt="Category	Term	Count	%	PValue	Genes	List Total	Pop Hits	Pop Total	Fold Enrichment	Bonferroni	Benjamini	FDR
KEGG_PATHWAY	hsa04721:Synaptic vesicle cycle	8	1.5873015873015872	0.0012399844447873723	SYT1, AP2B1, ATP6V1H, UNC13C, CLTCL1, UNC13B, AP2M1, CACNA1B	182	63	6879	4.799581371009943	0.24264930626033876	0.24264930626033876	1.5715045203099187
KEGG_PATHWAY	hsa05032:Morphine addiction	9	1.7857142857142856	0.0025988720185478507	ADCY2, KCNJ6, GABRB3, ADCY9, PDE10A, PRKCG, GABBR2, PDE8A, CACNA1B	182	91	6879	3.7381354908827436	0.44172583858068404	0.25282253686335265	3.267540392309576
KEGG_PATHWAY	hsa04961:Endocrine and other factor-regulated calcium reabsorption	6	1.1904761904761905	0.006215939205881658	AP2B1, ADCY9, ATP1A3, PRKCG, CLTCL1, AP2M1	182	45	6879	5.039560439560439	0.7525918417745279	0.3722240526486088	7.65163659074013
KEGG_PATHWAY	hsa04727:GABAergic synapse	8	1.5873015873015872	0.00682260021789801	PLCL1, ADCY2, KCNJ6, GABRB3, ADCY9, PRKCG, GABBR2, CACNA1B	182	85	6879	3.557336780866193	0.7842208062099019	0.3184426721502881	8.368752220142351
KEGG_PATHWAY	hsa04713:Circadian entrainment	8	1.5873015873015872	0.012278102794673108	RPS6KA5, ADCY2, KCNJ6, ADCY9, CAMK2G, RYR3, PER2, PRKCG	182	95	6879	3.18288027761712	0.9371697552915499	0.42504513881937067	14.590782790719492
KEGG_PATHWAY	hsa04915:Estrogen signaling pathway	8	1.5873015873015872	0.015161915605377912	ADCY2, KCNJ6, ADCY9, ATF6B, GABBR2, SHC3, GRM1, HSPA8	182	99	6879	3.0542790542790548	0.9673620420293375	0.4346893554149077	17.72014993982144
KEGG_PATHWAY	hsa05100:Bacterial invasion of epithelial cells	7	1.3888888888888888	0.01640351199278537	ARHGEF26, GAB1, SHC3, CLTCL1, CTNNA3, CTNNA2, ELMO1	182	78	6879	3.392011834319527	0.9753962591360783	0.410962396421404	19.034634905614855
KEGG_PATHWAY	hsa04723:Retrograde endocannabinoid signaling	8	1.5873015873015872	0.016771772517084036	ADCY2, KCNJ6, GABRB3, ADCY9, PRKCG, MAPK10, GRM1, CACNA1B	182	101	6879	2.993798280926994	0.9773758766924169	0.3772389103012044	19.420776103199533
KEGG_PATHWAY	hsa04728:Dopaminergic synapse	9	1.7857142857142856	0.019293437469790037	KCNJ6, CAMK2G, KIF5C, ATF6B, PRKCG, MAPK10, PPP2R2C, PPP2R2A, CACNA1B	182	128	6879	2.6575807005494507	0.987272173552954	0.3842319682313815	22.01961444882292
KEGG_PATHWAY	hsa04144:Endocytosis	13	2.579365079365079	0.024813933275234866	RAB7A, ARFGAP3, KIF5C, PSD3, RAB11FIP5, AP2B1, WIPF2, PDCD6IP, NEDD4L, CLTCL1, HSPA8, ARAP1, AP2M1	182	241	6879	2.0388263189093063	0.996405841641247	0.4304134561152373	27.44161703573994
KEGG_PATHWAY	hsa04261:Adrenergic signaling in cardiomyocytes	9	1.7857142857142856	0.028751765867603175	RPS6KA5, ADCY2, ADCY9, CAMK2G, ATF6B, ATP1A3, CACNB2, PPP2R2C, PPP2R2A	182	138	6879	2.4650023889154324	0.9985479744681158	0.44792553615139197	31.094470420416307
KEGG_PATHWAY	hsa04310:Wnt signaling pathway	9	1.7857142857142856	0.028751765867603175	CHD8, CSNK2A1, BTRC, CAMK2G, LRP6, LEF1, PRKCG, MAPK10, RBX1	182	138	6879	2.4650023889154324	0.9985479744681158	0.44792553615139197	31.094470420416307
KEGG_PATHWAY	hsa04724:Glutamatergic synapse	8	1.5873015873015872	0.030272428049298993	ADCY2, ADCY9, GRIK1, GRIK2, GRIK4, PRKCG, GRIN3A, GRM1	182	114	6879	2.6524002313476	0.9989777840746686	0.43662804957012	32.45913316872849
KEGG_PATHWAY	hsa04340:Hedgehog signaling pathway	4	0.7936507936507936	0.03297287102926383	BTRC, PTCH1, HHIP, GLI3	182	27	6879	5.5995115995116	0.999452657236921	0.4388268604804727	34.82132262200695
KEGG_PATHWAY	hsa04918:Thyroid hormone synthesis	6	1.1904761904761905	0.03651269160273091	GSR, ADCY2, ADCY9, ATF6B, ATP1A3, PRKCG	182	70	6879	3.2397174254317114	0.999759284747653	0.44851212043900646	37.8024650999936
KEGG_PATHWAY	hsa04930:Type II diabetes mellitus	5	0.992063492063492	0.036727907019762446	HK1, CACNA1E, MAPK10, INSR, CACNA1B	182	48	6879	3.937156593406593	0.9997710339021808	0.42810265699984007	37.979596254425005
KEGG_PATHWAY	hsa04913:Ovarian steroidogenesis	5	0.992063492063492	0.039199900719651105	CYP2J2, ADCY2, ADCY9, SCARB1, INSR	182	49	6879	3.856806458847275	0.9998712341210456	0.4287022024398337	39.981080244704756
KEGG_PATHWAY	hsa04066:HIF-1 signaling pathway	7	1.3888888888888888	0.04036335542540474	FLT1, CAMK2G, TEK, HK1, PRKCG, INSR, RBX1	182	96	6879	2.7560096153846154	0.9999018413738955	0.41892647758849155	40.902331618241625"
kegg=read.delim(text = txt, stringsAsFactors = F) %>% arrange(-PValue)
kegg$Term <- sub("hsa\\d+:(.*)","\\1",kegg$Term)
kegg$Term <- factor(kegg$Term, as.character(kegg$Term))

pdf("../results/Merge_circexplorer_BC.cellspecific_heatmap5.genes3.SNDA.KEGG.pdf", width=9, height=nrow(kegg)/5)
p = ggplot(kegg, aes(x = Term, y = -log10(PValue))) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  theme_bw() + theme_classic() +
  xlab("KEGG pathway") + ylab("-log10(EASE p value)") 
print(p)
dev.off() 

txt="Category	Term	Count	%	PValue	Genes	List Total	Pop Hits	Pop Total	Fold Enrichment	Bonferroni	Benjamini	FDR
KEGG_PATHWAY	hsa04724:Glutamatergic synapse	11	2.9729729729729732	5.64570739802734E-5	GLS2, GRM5, ADCY1, DLGAP1, GRM3, GNAO1, GRM8, GRIN1, GRIN2A, PRKACB, CACNA1A	132	114	6879	5.028508771929824	0.011395681234243615	0.011395681234243615	0.07091243871435449
KEGG_PATHWAY	hsa04727:GABAergic synapse	9	2.4324324324324325	1.962765869391299E-4	GLS2, ADCY1, GABRG3, GNAO1, TRAK2, GABRA4, GABRB2, PRKACB, CACNA1A	132	85	6879	5.517914438502673	0.03906456538858294	0.019726857140614507	0.24633237456836987
KEGG_PATHWAY	hsa04723:Retrograde endocannabinoid signaling	9	2.4324324324324325	6.390785063163251E-4	GRM5, ADCY1, GABRG3, GNAO1, GABRA4, GABRB2, MAPK9, PRKACB, CACNA1A	132	101	6879	4.6437893789378935	0.12170645649942813	0.04233584982438343	0.8000107593896133
KEGG_PATHWAY	hsa05033:Nicotine addiction	6	1.6216216216216217	8.941913852676818E-4	GABRG3, GABRA4, GABRB2, GRIN1, GRIN2A, CACNA1A	132	40	6879	7.817045454545456	0.16606686759699552	0.044385332070218	1.1177175697051678
KEGG_PATHWAY	hsa04720:Long-term potentiation	7	1.891891891891892	0.0015159882074382769	GRM5, ADCY1, CAMK4, GRIN1, GRIN2A, RAF1, PRKACB	132	66	6879	5.527203856749311	0.26506933054385495	0.05973715672358382	1.8881580262935516
KEGG_PATHWAY	hsa05032:Morphine addiction	8	2.1621621621621623	0.001659522515263193	ADCY1, GABRG3, GNAO1, GABRA4, GABRB2, ARRB1, PRKACB, CACNA1A	132	91	6879	4.581418581418581	0.2862074554867622	0.05464410960590338	2.0652182059219926
KEGG_PATHWAY	hsa04020:Calcium signaling pathway	11	2.9729729729729732	0.0021025770748828157	ATP2B1, GRM5, GNA14, ADCY1, CAMK4, GRIN1, PDGFRA, HTR4, GRIN2A, PRKACB, CACNA1A	132	179	6879	3.202513966480447	0.3477142136360851	0.05921338313187474	2.609905416925007
KEGG_PATHWAY	hsa04024:cAMP signaling pathway	11	2.9729729729729732	0.004360499218209016	ATP2B1, ADCY1, CAMK4, TIAM1, GRIN1, HTR4, CREB3L2, GRIN2A, RAF1, MAPK9, PRKACB	132	198	6879	2.8952020202020203	0.588158571601695	0.10496245402837834	5.342717404914088
KEGG_PATHWAY	hsa04540:Gap junction	7	1.891891891891892	0.006438557039101843	GRM5, ADCY1, PDGFRA, RAF1, PRKG2, PRKACB, LPAR1	132	88	6879	4.145402892561983	0.7305177926418227	0.13557848794759375	7.79528508698204
KEGG_PATHWAY	hsa05030:Cocaine addiction	5	1.3513513513513513	0.013736623466359347	GRM3, GRIN1, CREB3L2, GRIN2A, PRKACB	132	49	6879	5.317717996289425	0.9396662909722703	0.24481039660948256	15.95269984121389
KEGG_PATHWAY	hsa04080:Neuroactive ligand-receptor interaction	12	3.2432432432432434	0.016212477344186474	GRM5, GABRG3, GRM3, RXFP1, GABRA4, GRM8, GABRB2, GRIN1, HTR4, GRIN2A, CALCRL, LPAR1	132	277	6879	2.257630456186413	0.9637785398477117	0.26039999594918617	18.565509276080615
KEGG_PATHWAY	hsa04015:Rap1 signaling pathway	10	2.7027027027027026	0.01845904486929561	ADCY1, GNAO1, TIAM1, GRIN1, PDGFRA, GRIN2A, RAF1, LPAR1, DOCK4, FARP2	132	210	6879	2.4816017316017316	0.9772275185985022	0.2703450557309799	20.87143915471653
KEGG_PATHWAY	hsa04725:Cholinergic synapse	7	1.891891891891892	0.01890973086551148	KCNQ5, ADCY1, GNAO1, CAMK4, CREB3L2, PRKACB, CACNA1A	132	111	6879	3.286445536445537	0.9792546604038203	0.25778051760418785	21.32673324840424
KEGG_PATHWAY	hsa05230:Central carbon metabolism in cancer	5	1.3513513513513513	0.03316816758967316	GLS2, PDGFRA, RAF1, MTOR, PTEN	132	64	6879	4.071377840909091	0.9989377151234551	0.38681906580053016	34.54517950210243
KEGG_PATHWAY	hsa04713:Circadian entrainment	6	1.6216216216216217	0.034445796451056956	ADCY1, GNAO1, GRIN1, GRIN2A, PRKG2, PRKACB	132	95	6879	3.291387559808612	0.9991878022997882	0.3777320253544443	35.623692915696324
KEGG_PATHWAY	hsa04728:Dopaminergic synapse	7	1.891891891891892	0.034970716759053166	GNAO1, CREB3L2, GRIN2A, MAPK9, PRKACB, PPP2R2D, CACNA1A	132	128	6879	2.8499644886363638	0.9992726894851702	0.36341340137964095	36.06204579381457
KEGG_PATHWAY	hsa05031:Amphetamine addiction	5	1.3513513513513513	0.0365578442400563	CAMK4, GRIN1, CREB3L2, GRIN2A, PRKACB	132	66	6879	3.9480027548209367	0.9994792779640422	0.35899831091651113	37.37077604785078
KEGG_PATHWAY	hsa04916:Melanogenesis	6	1.6216216216216217	0.041561531015856724	DCT, ADCY1, GNAO1, CREB3L2, RAF1, PRKACB	132	100	6879	3.126818181818182	0.999819057720442	0.38043687911628066	41.33710052325436
KEGG_PATHWAY	hsa04141:Protein processing in endoplasmic reticulum	8	2.1621621621621623	0.04200364107431749	NGLY1, RAD23B, TUSC3, HSPA4L, MAPK9, SKP1, MARCH6, SEC24D	132	169	6879	2.4669176976869287	0.9998352356438934	0.3677518892411684	41.676192468752035
KEGG_PATHWAY	hsa04976:Bile secretion	5	1.3513513513513513	0.042004583973712116	ADCY1, NCEH1, LDLR, PRKACB, ABCG2	132	69	6879	3.7763504611330703	0.999835268560695	0.3530975072320862	41.67691372687984
"

kegg=read.delim(text = txt, stringsAsFactors = F) %>% arrange(-PValue)
kegg$Term <- sub("hsa\\d+:(.*)","\\1",kegg$Term)
kegg$Term <- factor(kegg$Term, as.character(kegg$Term))

pdf("../results/Merge_circexplorer_BC.cellspecific_heatmap5.genes3.PY.KEGG.pdf", width=9, height=nrow(kegg)/5)
p = ggplot(kegg, aes(x = Term, y = -log10(PValue))) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  theme_bw() + theme_classic() +
  xlab("KEGG pathway") + ylab("-log10(EASE p value)") 
print(p)
dev.off() 

txt="Category	Term	Count	%	PValue	Genes	List Total	Pop Hits	Pop Total	Fold Enrichment	Bonferroni	Benjamini	FDR
KEGG_PATHWAY	hsa04120:Ubiquitin mediated proteolysis	24	2.013422818791946	1.999388290751245E-5	ANAPC1, UBE2Z, UBE3B, DDB1, ANAPC4, UBE2F, SAE1, UBE2C, MID1, UBE2R2, FANCL, CBLB, FBXW7, CUL5, FBXW8, UBE2D2, CUL4A, UBE2K, MAP3K1, MDM2, FBXO4, PIAS1, ANAPC7, CUL1	447	137	6879	2.6959290648116396	0.004886571952388397	0.004886571952388397	0.025888307136945343
KEGG_PATHWAY	hsa04710:Circadian rhythm	9	0.7550335570469799	6.132477247086505E-4	CRY2, CREB1, PRKAB1, PRKAA1, BHLHE40, RORA, PER3, FBXL3, CUL1	447	31	6879	4.467850184022515	0.13954312458126727	0.0723918524405186	0.7912317454761841
KEGG_PATHWAY	hsa04110:Cell cycle	19	1.5939597315436242	9.816135383885326E-4	ANAPC1, CDC14A, DBF4, ANAPC4, RBL1, SMAD4, PRKDC, SMAD2, ATR, MCM3, YWHAE, ATM, CCNB1, MDM2, BUB1B, ANAPC7, ORC2, CUL1, STAG2	447	124	6879	2.358032041567439	0.21385453282567402	0.07707239035798463	1.263727569408779
KEGG_PATHWAY	hsa05130:Pathogenic Escherichia coli infection	11	0.9228187919463088	0.001363473251467702	ARPC1A, CDC42, NCK2, TUBA8, ROCK1, FYN, ROCK2, HCLS1, TLR5, ITGB1, CTNNB1	447	51	6879	3.3192525332280565	0.2841458891606605	0.08017304787577328	1.7513390143389262
KEGG_PATHWAY	hsa04144:Endocytosis	29	2.4328859060402683	0.0017595958540909907	USP8, SNX6, CYTH1, SNX5, CHMP4A, CAPZA1, ASAP2, SNX1, PIP5K1A, CDC42, IGF1R, DAB2, KIAA1033, SPG21, STAM, FAM21C, TGFBR2, SMAD2, HLA-E, ARPC1A, CBLB, RAB11FIP3, IGF2R, IST1, RAB5A, MDM2, RAB10, RAB11FIP1, DNM2	447	241	6879	1.8518198780250075	0.3504534088125151	0.08267746701592005	2.2548143489506844
KEGG_PATHWAY	hsa04520:Adherens junction	13	1.0906040268456376	0.0018280028806322531	PTPRJ, ERBB2, MET, TGFBR2, SMAD4, SMAD2, TCF7L2, CTNNB1, CDC42, IGF1R, MAPK1, FYN, YES1	447	71	6879	2.8177521504868133	0.36126814864941137	0.07198906089589019	2.3415188260386355
KEGG_PATHWAY	hsa05205:Proteoglycans in cancer	24	2.013422818791946	0.004998955668653182	ROCK1, ROCK2, HCLS1, ERBB2, MET, RPS6KB1, PDCD4, ITGB1, CTNNB1, PTPN11, CDC42, IGF1R, MAPK1, PDPK1, CBLB, HIF1A, CD44, PLCG2, SOS2, PPP1R12A, CAMK2D, MDM2, THBS1, FN1	447	200	6879	1.8467114093959731	0.707068240548993	0.16088056492982516	6.28361298996305
KEGG_PATHWAY	hsa04068:FoxO signaling pathway	18	1.5100671140939599	0.005683788394581946	USP7, SGK1, TGFBR2, SMAD4, PRKAB1, BNIP3, FOXO1, SMAD2, STK4, ATM, CCNB1, IGF1R, MAPK1, PDPK1, SOS2, MDM2, PRKAA1, EGF	447	134	6879	2.0672142642492233	0.7525385914031621	0.1601757029890878	7.115472053941996
KEGG_PATHWAY	hsa03022:Basal transcription factors	9	0.7550335570469799	0.0075427170393805735	TAF2, GTF2E1, TAF1, TAF12, GTF2A1, TAF5, GTF2B, GTF2H2, ERCC2	447	45	6879	3.077852348993288	0.8435424044038107	0.1862546556529191	9.339267314685829
KEGG_PATHWAY	hsa03013:RNA transport	21	1.761744966442953	0.007587335854862915	NCBP1, NUP133, NUP153, RGPD8, EIF5, NUP155, SMN2, SMN1, EIF4G2, NUP214, UPF3B, EIF3B, EIF3H, EIF3E, NUP50, CYFIP2, NUP107, PABPC1, TPR, THOC2, TGS1	447	172	6879	1.878921492117996	0.8452563184758807	0.17022319155969867	9.392034893973623
KEGG_PATHWAY	hsa04070:Phosphatidylinositol signaling system	14	1.174496644295302	0.010241552728108595	PIK3C2A, PIP5K1A, PI4KB, PI4K2B, TMEM55A, MTM1, MTMR1, DGKE, DGKD, DGKG, PLCG2, PIP4K2A, PIP4K2B, IP6K2	447	98	6879	2.198465963566635	0.9197106085653977	0.20489685801765212	12.480490036020742
KEGG_PATHWAY	hsa03018:RNA degradation	12	1.006711409395973	0.010287553401714264	EXOSC10, PATL1, CNOT8, CNOT6L, DCP2, DCP1A, CNOT3, DHX36, CNOT1, PABPC1, CNOT6, HSPA9	447	77	6879	2.398326505709056	0.9206196820643899	0.19032747420131169	12.533149646721842
KEGG_PATHWAY	hsa04114:Oocyte meiosis	15	1.2583892617449663	0.012103119039741617	PPP2R1B, ANAPC1, PPP2R5C, ANAPC4, AURKA, YWHAE, CCNB1, MAPK1, IGF1R, SLK, CAMK2D, PPP3CB, PPP3CA, ANAPC7, CUL1	447	111	6879	2.079629965536006	0.949377278089741	0.205060149724416	14.588332806137682
KEGG_PATHWAY	hsa05200:Pathways in cancer	38	3.1879194630872485	0.012956513253588481	ERBB2, FOXO1, ITGB1, TCF7L2, TPM3, ARNT, CTNNB1, IGF1R, CDC42, RASGRP1, SOS2, TPR, RUNX1, EGF, LAMB1, CSF2RA, AXIN1, FN1, ROCK1, MSH2, FLT3, RALBP1, ROCK2, MET, TGFBR2, SMAD4, SMAD2, STK4, MAPK1, VEGFC, CBLB, HIF1A, GNAQ, GNB1, ETS1, PLCG2, MDM2, GNB4	447	393	6879	1.488020219615076	0.9590371346053286	0.20405139951988338	15.538879386080673
KEGG_PATHWAY	hsa04360:Axon guidance	16	1.342281879194631	0.01693564254322581	ROCK1, ROCK2, MET, NTNG1, ITGB1, CDC42, MAPK1, NCK2, FYN, PPP3CB, SEMA4B, SEMA3C, PPP3CA, NFATC2, RASA1, SRGAP1	447	127	6879	1.9388046292871108	0.9847742800366913	0.2434484549722855	19.843494701275265
KEGG_PATHWAY	hsa04666:Fc gamma R-mediated phagocytosis	12	1.006711409395973	0.019073306293196472	ARPC1A, CDC42, MAPK1, PTPRC, MYO10, VAV3, GAB2, LYN, PLCG2, ASAP2, RPS6KB1, PIP5K1A	447	84	6879	2.1984659635666346	0.991067860552125	0.2553799631150948	22.07152098631042
KEGG_PATHWAY	hsa04611:Platelet activation	16	1.342281879194631	0.020564480506745177	ROCK1, LYN, ROCK2, PRKG1, APBB1IP, ITGB1, MAPK1, GNAQ, FYN, RASGRP1, PLCG2, PPP1R12A, GUCY1A3, GUCY1B3, COL1A1, SNAP23	447	130	6879	1.89406298399587	0.9938470527094566	0.25878180724157085	23.591738779425995
KEGG_PATHWAY	hsa04141:Protein processing in endoplasmic reticulum	19	1.5939597315436242	0.025083865916147495	MBTPS2, SEC24A, PDIA3, SEC63, STT3B, ATXN3, EIF2AK1, UBE2D2, STT3A, DNAJA1, ERN1, DNAJC5, NFE2L2, DNAJC3, UGGT1, TRAM1, CUL1, DNAJC1, SEC23B	447	169	6879	1.7301536873039198	0.9980185727527369	0.29232782746683994	28.03360789417555
KEGG_PATHWAY	hsa04922:Glucagon signaling pathway	13	1.0906040268456376	0.02570796113304237	CREB1, PRKAB1, PDE3B, FOXO1, CREB5, CPT1A, PDHB, PKM, GNAQ, CAMK2D, PPP3CB, PRKAA1, PPP3CA	447	99	6879	2.0208121483289267	0.9983062757904018	0.2852574857423088	28.627915555894024
KEGG_PATHWAY	hsa05215:Prostate cancer	12	1.006711409395973	0.02609929661045482	MAPK1, IGF1R, PDPK1, CREB1, ERBB2, SOS2, MDM2, FOXO1, CREB5, EGF, TCF7L2, CTNNB1	447	88	6879	2.098535692495424	0.9984650420691525	0.27672271858279895	28.998259177096376
KEGG_PATHWAY	hsa04022:cGMP-PKG signaling pathway	18	1.5100671140939599	0.02655432926537824	SLC8A3, MEF2A, ROCK1, ROCK2, CREB1, PDE3B, PDE3A, CREB5, PRKG1, MAPK1, GNAQ, PDE5A, PPP1R12A, PPP3CB, GUCY1A3, GUCY1B3, PPP3CA, NFATC2	447	158	6879	1.7532070342366832	0.9986311025686646	0.26947226322989726	29.426653067320263
KEGG_PATHWAY	hsa04152:AMPK signaling pathway	15	1.2583892617449663	0.027476835164073093	PPP2R1B, IGF1R, PDPK1, CD36, PFKFB4, CREB1, PPP2R5C, PRKAB1, FOXO1, CREB5, PRKAA1, RPS6KB1, RAB10, RPTOR, CPT1A	447	123	6879	1.8767392371910294	0.9989148521708198	0.26675434910741236	30.287843174220242
KEGG_PATHWAY	hsa04660:T cell receptor signaling pathway	13	1.0906040268456376	0.0275701147681274	PTPRC, VAV3, NFKBIE, CDC42, MAPK1, NCK2, PDPK1, FYN, RASGRP1, SOS2, PPP3CB, PPP3CA, NFATC2	447	100	6879	2.0006040268456373	0.998940056173864	0.2575546839328745	30.37438076078205
KEGG_PATHWAY	hsa03420:Nucleotide excision repair	8	0.6711409395973155	0.030332067016129872	XPA, ERCC6, POLE2, RFC1, CUL4A, DDB1, GTF2H2, ERCC2	447	47	6879	2.6194488076538627	0.9994719886563169	0.269797242150669	32.892227198776695
KEGG_PATHWAY	hsa04810:Regulation of actin cytoskeleton	22	1.8456375838926176	0.03062381722494732	ITGAL, VAV3, ROCK1, ROCK2, DIAPH2, ARHGEF6, DIAPH3, ITGA11, IQGAP2, PIP5K1A, ITGB1, NCKAP1, ARPC1A, CDC42, MAPK1, SOS2, CYFIP2, PPP1R12A, PIP4K2A, EGF, FN1, PIP4K2B	447	210	6879	1.6122083732821988	0.9995095163451669	0.26273268388278503	33.15322616451268
KEGG_PATHWAY	hsa05203:Viral carcinogenesis	21	1.761744966442953	0.04282466324261256	USP7, SP100, LYN, IL6ST, CREB1, DDB1, RBL1, CREB5, HLA-E, GTF2B, YWHAE, GTF2H2, PKM, CDC42, GTF2E1, MAPK1, REL, GTF2A1, SND1, MDM2, RASA2	447	205	6879	1.576460959240465	0.99997797492279	0.3379642531077439	43.26561230785367
KEGG_PATHWAY	hsa04510:Focal adhesion	21	1.761744966442953	0.044694798386530565	VAV3, ROCK1, ROCK2, ERBB2, MET, ITGA11, ITGB1, CTNNB1, CDC42, IGF1R, VEGFC, MAPK1, PDPK1, FYN, SOS2, PPP1R12A, COL1A1, THBS1, LAMB1, EGF, FN1	447	206	6879	1.5688082361373559	0.9999863596690871	0.33959807598523695	44.68441833195712
KEGG_PATHWAY	hsa03450:Non-homologous end-joining	4	0.33557046979865773	0.047567724069199485	XRCC5, XRCC4, XRCC6, PRKDC	447	13	6879	4.735157459989675	0.999993478379005	0.3471725418327013	46.80034288247007
KEGG_PATHWAY	hsa05213:Endometrial cancer	8	0.6711409395973155	0.04903164430683975	MAPK1, PDPK1, ERBB2, SOS2, EGF, TCF7L2, AXIN1, CTNNB1	447	52	6879	2.3675787299948374	0.9999955261019459	0.3460554880397838	47.849562135395075"

kegg=read.delim(text = txt, stringsAsFactors = F) %>% arrange(-PValue)
kegg$Term <- sub("hsa\\d+:(.*)","\\1",kegg$Term)
kegg$Term <- factor(kegg$Term, as.character(kegg$Term))

pdf("../results/Merge_circexplorer_BC.cellspecific_heatmap5.genes3.NN.KEGG.pdf", width=9, height=nrow(kegg)/5)
p = ggplot(kegg, aes(x = Term, y = -log10(PValue))) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  theme_bw() + theme_classic() +
  xlab("KEGG pathway") + ylab("-log10(EASE p value)") 
print(p)
dev.off() 

#############################################
# cell-specific isoform of circRNAs, esp. those from the same host gene (esp. the 488 shared genes btw SNDA and PY)
#############################################
df5_by_start = df5 %>% separate(gene, c("chr","start","end"), remove = F) %>% mutate(celltype3=ifelse(celltype %in% c("TCPY","MCPY"), "PY", ifelse(celltype=="SNDA","SNDA","NN"))) %>% arrange(celltype3) %>% group_by(chr, start) %>% summarize(celltype3s = paste(unique(celltype3), collapse = ','), celltype3n = length(unique(celltype3)), celltype5s = paste(unique(celltype), collapse = ','), celltype5n = length(unique(celltype)), hostgenes = paste(unique(hostgene), collapse = ','), genes = paste(gene, collapse = ','), Ss=paste(S, collapse = ','), means=paste(mean, collapse = ',')) 

df5_by_end = df5 %>% separate(gene, c("chr","start","end"), remove = F) %>% mutate(celltype3=ifelse(celltype %in% c("TCPY","MCPY"), "PY", ifelse(celltype=="SNDA","SNDA","NN"))) %>% arrange(celltype3) %>% group_by(chr, end) %>% summarize(celltype3s = paste(unique(celltype3), collapse = ','), celltype3n = length(unique(celltype3)), celltype5s = paste(unique(celltype), collapse = ','), celltype5n = length(unique(celltype)), hostgenes = paste(unique(hostgene), collapse = ','), genes = paste(gene, collapse = ','), Ss=paste(S, collapse = ','), means=paste(mean, collapse = ',')) 

## PY and SNDA
df5_by_end %>% ungroup() %>% filter(celltype3s=="PY,SNDA") %>% dplyr::rowwise() %>% mutate(minMean=min(as.numeric(strsplit(means,",")[[1]]))) %>% filter(minMean>0.001) %>% arrange(-minMean)
df5_by_start %>% ungroup() %>% filter(celltype3s=="PY,SNDA") %>% dplyr::rowwise() %>% mutate(minMean=min(as.numeric(strsplit(means,",")[[1]]))) %>% filter(minMean>0.001) %>% arrange(-minMean)

Merge_circexp_raw_filtered_and_enriched = readRDS(file="Merge_circexplorer_BC106.filtered.enriched.rawcount.rds")

df5_by_start_sub = df5 %>% filter(gene %in% rownames(Merge_circexp_raw_filtered_and_enriched)[rowSums(Merge_circexp_raw_filtered_and_enriched>0)>10]) %>% separate(gene, c("chr","start","end"), remove = F) %>% mutate(celltype3=ifelse(celltype %in% c("TCPY","MCPY"), "PY", ifelse(celltype=="SNDA","SNDA","NN"))) %>% arrange(celltype3) %>% group_by(chr, start) %>% summarize(celltype3s = paste(unique(celltype3), collapse = ','), celltype3n = length(unique(celltype3)), celltype5s = paste(unique(celltype), collapse = ','), celltype5n = length(unique(celltype)), hostgenes = paste(unique(hostgene), collapse = ','), genes = paste(gene, collapse = ','), Ss=paste(S, collapse = ','), means=paste(mean, collapse = ',')) %>% print(n = Inf)

df5_by_end_sub = df5 %>% separate(gene, c("chr","start","end"), remove = F) %>% mutate(celltype3=ifelse(celltype %in% c("TCPY","MCPY"), "PY", ifelse(celltype=="SNDA","SNDA","NN"))) %>% arrange(celltype3) %>% group_by(chr, end) %>% summarize(celltype3s = paste(unique(celltype3), collapse = ','), celltype3n = length(unique(celltype3)), celltype5s = paste(unique(celltype), collapse = ','), celltype5n = length(unique(celltype)), hostgenes = paste(unique(hostgene), collapse = ','), genes = paste(gene, collapse = ','), Ss=paste(S, collapse = ','), means=paste(mean, collapse = ',')) 

