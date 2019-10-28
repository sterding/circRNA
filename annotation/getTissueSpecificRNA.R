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

dim(Merge_circexp_norm_filtered_and_enriched); dim(annotation)

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
# 4      SNDA    6566   133
# 5      TCPY    8237   144
df %>% merge(annotation, by.x='gene',by.y='ID', all.x=T, all.y=F) %>% select(gene, cellType3, cellType5, circType) %>% group_by(cellType3, circType) %>% summarise(n=n_distinct(gene)) %>% dcast(cellType3 ~ circType, value.var="n") #%>% write.table(file=clip, sep="\t", row.names=F)
#   cellType3 circRNA ciRNA
# 1        NN    6483   309
# 2        PY    8511   147
# 3      SNDA    6566   133

## annotation per cell type
df %>% group_by(gene) %>% summarize(celltype3 = paste(unique(cellType3), collapse = ','), celltype5 = paste(unique(cellType5), collapse = ',')) %>% ungroup() %>% merge(annotation, by.x='gene',by.y='ID', all.x=T, all.y=F) %>% write.table(file="../results/Merge_circexplorer_BC197.annotation_per_cell.xls", sep="\t", quote=F, row.names=F)

celltype3_n = df %>% group_by(gene) %>% summarize(celltype3 = paste(unique(cellType3), collapse = ','), celltype5 = paste(unique(cellType5), collapse = ',')) %>% group_by(celltype3) %>% summarise(n=n())
#  celltype3      n
# 1 NN          2362
# 2 NN,PY       1362
# 3 NN,PY,SNDA  2305
# 4 NN,SNDA      763
# 5 PY          2978
# 6 PY,SNDA     2013
# 7 SNDA        1618

## call eulerAPE to generate Merge_circexplorer_BC197.venn_celltype3.svg

# 5-ways venn diagram (not area-proportional)
library(gplots); # install.packages('gplots')
pdf("../results/Merge_circexplorer_BC197.venn_celltype5.pdf")
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

groupmean_s3 = groupmean_to_specificity(group3mean)
df3 = groupmean_s3 %>% filter(Private_or_not==1) ## n = 10945
group_by(df3, celltype) %>% summarise(count=n())
#   celltype count
# 1 NN        4622
# 2 PY        4436
# 3 SNDA      1887

groupmean_s5 = groupmean_to_specificity(group5mean)
df5 = groupmean_s5 %>% filter(Private_or_not==1) 
group_by(df5, celltype) %>% summarise(count=n())
#   celltype count
# 1 FB         763
# 2 MCPY       426
# 3 PBMC      3627
# 4 SNDA      1904
# 5 TCPY      3968

## df3 and df5 are very similar
rbind(data.frame(gene=df3$gene, DF="df3"), data.frame(gene=df5$gene, DF="df5")) %>% group_by(gene) %>% summarize(DFs = paste(unique(DF), collapse = ',')) %>% group_by(DFs) %>% summarise(n=n())
#   DFs         n
# 1 df3       887
# 2 df3,df5 10058
# 3 df5       630

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
dissimilarity <- 1 - cor(df)
distance <- as.dist(dissimilarity)
hc=hclust(distance)
hccol=as.dendrogram(hc)
weights.dd <- order(match(c("PBMC","FB","SNDA","TCPY","MCPY"), celltypes)) # the trick is to call order() on the specific index of target dendragram
#plot(reorder(hccol, wts = weights.dd))

dissimilarity <- 1 - cor(t(df))
distance <- as.dist(dissimilarity)
hcrow=hclust(distance, method = 'average')
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
                      cluster_rows = hcrow, cluster_cols = as.hclust(reorder(hccol, wts = weights.dd, agglo.FUN = mean)), 
                      cutree_cols = 5, # to allow small gap space between columns
                      clustering_distance_rows = 'correlation', 
                      clustering_distance_cols = 'correlation')
hm=do.call("pheatmap", hm.parameters)
do.call("pheatmap", c(hm.parameters, width=3, height=5, filename="../results/Merge_circexplorer_BC197.cellspecific_heatmap++.pdf"))

###########################################
# host genes of cell-specific circRNAs
###########################################
# circRNAs in the order of presence in above figure
df_roworders = df3$gene[roworders]; head(df_roworders); length(df_roworders)  # n = 10945
DF3=groupmean_s3[df_roworders,]; dim(DF3); head(DF3); table(DF3$Private_or_not)

## add host gene of the cell-specific circRNAs
filtered_enriched_annotation=readRDS("Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds"); head(filtered_enriched_annotation); dim(filtered_enriched_annotation)
DF3$hostgene=filtered_enriched_annotation$geneName[match(DF3$gene, filtered_enriched_annotation$ID)]
DF3$hostgeneID=filtered_enriched_annotation$geneID[match(DF3$gene, filtered_enriched_annotation$ID)] 
dim(DF3); DF3=filter(DF3, !is.na(hostgeneID)); dim(DF3) # 10945

## How many host genes per celltype
DF3 %>% dplyr::select(hostgene, celltype) %>% distinct() %>% group_by(celltype) %>% summarise(n=n())
# 1 NN        2133
# 2 PY        2233
# 3 SNDA      1289

## save the list into table for downstream GO analysis
DF3 %>% dplyr::select(celltype, hostgene) %>% distinct() %>% group_by(celltype) %>% mutate(n=n(), hostgene = paste0(unique(hostgene), collapse = ";")) %>% distinct() %>% write.table(file="../results/Merge_circexplorer_BC197.cellspecific_heatmap5.genes3.txt", sep = "\t", col.names = T, quote=F, row.names = F)
DF3 %>% dplyr::select(celltype, hostgeneID) %>% mutate(hostgeneID=sub("\\..*","",hostgeneID)) %>% distinct() %>% group_by(celltype) %>% mutate(n=n(), hostgeneID = paste0(unique(hostgeneID), collapse = ";")) %>% distinct() %>% write.table(file="../results/Merge_circexplorer_BC197.cellspecific_heatmap5.genes3.EnsID.txt", sep = "\t", col.names = T, quote=F, row.names = F)
## GO to DAVID website to run KEGG pathway analysis, then go to getTissueSpecificRNA.makingGObarplot.R to make the figures

## divide into venn diagram --> ## call eulerAPE to generate Merge_circexplorer_BC197.cellspecific_heatmap5.genes3.venn_final.pdf
DF3 %>% arrange(celltype, hostgene)%>% group_by(hostgene) %>% summarize(celltype3 = paste(unique(celltype), collapse = ',')) %>% group_by(celltype3) %>% summarise(n=n(), gene3 = paste(unique(hostgene), collapse = ';'))
DF3 %>% arrange(celltype, hostgene)%>% group_by(hostgene) %>% summarize(celltype3 = paste(unique(celltype), collapse = ',')) %>% group_by(celltype3) %>% summarise(n=n(), gene3 = paste(unique(hostgene), collapse = ';')) %>% write.table(file="../results/Merge_circexplorer_BC.cellspecific_heatmap5.genes3.overlapping.txt", sep = "\t", col.names = T, quote=F, row.names = F)
# # celltype3      n gene3                                                                                                                                                                            
# 1 NN          1005 AAGAB;AASS;AATF;ABCB10;ABCC1;ABCC3;ABR;AC004381.6;AC068533.7;AC093838.4;ACTG1;ACTR10;ACTR2;ACTR3;ACYP2;ADAMTS12;ADAMTS6;ADAMTSL1;AFG3L2;AGFG1;AGGF1;AGPAT6;AHSA2;AK9;AKAP7;AKIRIN1;A…
# 2 NN,PY        581 ABI1;ACADM;ADAL;ADD1;ADK;AEBP2;AFF1;AFF3;AFF4;AGO3;ALG9;ALMS1;ALS2;ALS2CR12;AMN1;ANAPC1;ANKHD1;ANKRD12;ANKRD28;ANKRD32;APAF1;APBA1;APC;APLF;APLP2;ARAP2;ARHGAP12;ARHGAP18;ARHGAP19-S…
# 3 NN,PY,SNDA   373 AAK1;ABI2;ACSL4;ADAM10;ADAM32;ADD3;ADRBK2;AGTPBP1;AHCTF1;AHI1;AK5;AKAP10;AKAP13;AKT3;ALG8;ALS2CR11;AMPH;ANKAR;ANKRD13C;ANKRD17;ANKRD27;ANO10;AQR;ARFGEF1;ARFGEF2;ARHGAP10;ARHGAP26;A…
# 4 NN,SNDA      174 ABCC4;ADAM17;ADAM23;AGO4;AKAP6;AKAP9;AMOT;ANO6;ARHGAP6;ARSG;ATP6V0A2;ATRX;B4GALT6;BNC2;BTBD10;C10orf32-ASMT;C4orf29;CCP110;CDH8;CELF1;CEP350;CEP63;CERS6;CHMP5;CIRH1A;CNTRL;COG3;CRI…
# 5 PY           932 AACS;AASDH;ABAT;ABCA3;ABCA5;ABCC5;ABCC9;ABCD3;ABL2;ABLIM2;AC004076.9;AC034228.7;AC066593.1;AC067956.1;AC074391.1;ACADSB;ACAT1;ACBD5;ACBD6;ACER2;ACER3;ACOT9;ACP6;ACTN4;ACTR3B;ADAM19…
# 6 PY,SNDA      347 ABLIM1;ACAD11;ACAP2;ACSL3;ADAMTS19;ADARB1;AGAP1;AGL;AMBRA1;ANK3;ANKH;ANKIB1;ANKRD26;ANKS1A;ANKS1B;ANO4;ANO5;APBB2;ARHGAP32;ARHGEF26;ARHGEF28;ARHGEF7;ARID1B;ARMC2;ARMC9;ATG7;ATP10B;…
# 7 SNDA         395 ABCA1;AC004893.11;AC010642.1;AC012307.3;ACACA;ACACB;ACSS3;ADAMTSL3;ADCY2;AF127936.7;AGAP3;AGPS;AMMECR1L;ANK1;ANK2;ANKRD20A7P;AP2B1;AP2M1;APITD1-CORT;ARAP1;ARFGAP3;ARFIP1;ARHGEF11;A…

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
# 1 NN         319
# 2 PY        1513
# 3 SNDA      6377

gene_group3mean_DF3 = gene_group3mean[match(unique(DF3$hostgeneID), gene_group3mean$gene),]
gene_groupmean_s3_DF3 = groupmean_to_specificity(gene_group3mean_DF3)

## cell-specificity of circRNA VS. cell-specificity of hostgene
DF3 = inner_join(DF3, gene_groupmean_s3_DF3,by = c("hostgeneID" = "gene"), suffix=c(".circRNA",".gene")); head(DF3)
DF3 = mutate(DF3, celltype.specific.gene = ifelse(Private_or_not.gene==1, celltype.gene, "NS"))
dim(DF3)

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
ggsave("../results/Merge_circexplorer_BC197.cellspecific_barplot++.hostgene.pdf", width = 10, height = 2)
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
ggsave("../results/Merge_circexplorer_BC197.cellspecific_boxplot++.specificityscore.pdf", width=6, height = 2)

## heatmap of host genes of cell-specific circRNAs
X=gene_group5mean[match(DF3$hostgeneID, gene_group5mean$gene),]; head(X); dim(X); 
X=log10(X[,-1]+0.01)
dissimilarity <- 1 - cor(X)
distance <- as.dist(dissimilarity)
hc=hclust(distance)
hccol=as.dendrogram(hc)
weights.dd <- order(match(c("PBMC","FB","SNDA","TCPY","MCPY"), colnames(X))) # the trick is to call order() on the specific index of target dendragram
#plot(reorder(hccol, wts = weights.dd, agglo.FUN = mean))

my_gene_col <- data.frame(celltype=DF3$celltype.circRNA, row.names = rownames(X))

hm.parameters <- list(X, 
                      color = colorRampPalette(c("#999999","#ffffff", "#ca0020"))(100),
                      scale = "row",
                      treeheight_row = 50,
                      treeheight_col = 30,
                      kmeans_k = NA,
                      show_rownames = F, show_colnames = T,
                      annotation_row = my_gene_col,
                      cutree_cols = 5, # to allow small gap space between columns
                      cluster_rows = F, cluster_cols = as.hclust(reorder(hccol, wts = weights.dd, agglo.FUN = mean)))
library(pheatmap);
do.call("pheatmap", hm.parameters)
do.call("pheatmap", c(hm.parameters, width=3, height=5, filename="../results/Merge_circexplorer_BC197.cellspecific_heatmap++.hostgene.pdf"))

#############################################
# cell-specific isoform of circRNAs, esp. those from the same host gene (esp. the 720 shared genes btw SNDA and PY)
#############################################
Merge_circexp_raw_filtered_and_enriched = readRDS(file="Merge_circexplorer_BC197.filtered.enriched.rawcount.rds")
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