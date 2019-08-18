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

normRPM_file = "Merge_circexplorer_BC.filtered.enriched.normRPM.rds"

setwd("~/projects/circRNA/data")
# calculate specificity score (S) as cummeRbund (https://github.com/Bioconductor-mirror/cummeRbund/blob/master/R/tools.R)
source("~/projects/circRNA/src/annotation/tools.R")

Merge_circexp_norm_filtered_and_enriched=readRDS(normRPM_file); dim(Merge_circexp_norm_filtered_and_enriched)
annotation = readRDS("Merge_circexplorer_BC106.annotation.bed14.rds")
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
# celltype3     n
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

groupmean_s3 = groupmean_to_specificity(group3mean)
df3 = groupmean_s3 %>% filter(Private_or_not==1) 
group_by(df3, celltype) %>% summarise(count=n())
# celltype count
# 1       NN  4860
# 2       PY  1850
# 3     SNDA  2023

groupmean_s5 = groupmean_to_specificity(group5mean)
df5 = groupmean_s5 %>% filter(Private_or_not==1) 
group_by(df5, celltype) %>% summarise(count=n())
# celltype count
# 1 FB         795
# 2 MCPY       461
# 3 PBMC      3672
# 4 SNDA      1991
# 5 TCPY      1295

## df3 and df5 are very similar
rbind(data.frame(gene=df3$gene, DF="df3"), data.frame(gene=df5$gene, DF="df5")) %>% group_by(gene) %>% summarize(DFs = paste(unique(DF), collapse = ',')) %>% group_by(DFs) %>% summarise(n=n())
#   DFs         n
# 1 df3       758
# 2 df3,df5  7975
# 3 df5       239

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

## add host gene of the cell-specific circRNAs
filtered_enriched_annotation=readRDS("Merge_circexplorer_BC106.filtered.enriched.annotation.bed14.rds"); head(filtered_enriched_annotation)
genes=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed", header = F, col.names = c("chr","start","end","EnsID","score","strand","geneName","geneType"), stringsAsFactors = F)
DF3$hostgene=filtered_enriched_annotation$geneName[match(DF3$gene, filtered_enriched_annotation$ID)]
DF3$hostgeneID=genes$EnsID[match(DF3$hostgene, genes$geneName)]  
# In total 111 genes are discarded, as not found in the gencode table  (Note: in the future, circExplorer should use the same gene annotation as  GENCODE)
dim(DF3); DF3=filter(DF3, !is.na(hostgeneID)); dim(DF3) # 8733--> 8622

## How many host genes per celltype
DF3 %>% select(hostgene, celltype) %>% distinct() %>% group_by(celltype) %>% summarise(n=n())
# 1 NN        2123
# 2 PY        1155
# 3 SNDA      1316

## save the list into table for downstream GO analysis
DF3 %>% select(celltype, hostgene) %>% distinct() %>% group_by(celltype) %>% mutate(n=n(), hostgene = paste0(unique(hostgene), collapse = ";")) %>% distinct() %>% write.table(file="../results/Merge_circexplorer_BC.cellspecific_heatmap5.genes3.txt", sep = "\t", col.names = T, quote=F, row.names = F)
## GO to DAVID website to run KEGG pathway analysis, then go to getTissueSpecificRNA.makingGObarplot.R to make the figures

## divide into venn diagram --> ## call eulerAPE to generate Merge_circexplorer_BC.cellspecific_heatmap5.genes3.venn_final.pdf
DF3 %>% arrange(celltype, hostgene)%>% group_by(hostgene) %>% summarize(celltype3 = paste(unique(celltype), collapse = ',')) %>% group_by(celltype3) %>% summarise(n=n(), gene3 = paste(unique(hostgene), collapse = ';'))
DF3 %>% arrange(celltype, hostgene)%>% group_by(hostgene) %>% summarize(celltype3 = paste(unique(celltype), collapse = ',')) %>% group_by(celltype3) %>% summarise(n=n(), gene3 = paste(unique(hostgene), collapse = ';')) %>% write.table(file="../results/Merge_circexplorer_BC.cellspecific_heatmap5.genes3.overlapping.txt", sep = "\t", col.names = T, quote=F, row.names = F)
# celltype3      n gene3                                                                                                                                                                                                              
# 1 NN          1264 AASS;AATF;ABCB10;ABCC1;ABCC3;ABCC4;ABR;ACADM;ACTR10;ACTR2;ACTR3;ADAM23;ADAMTS12;ADAMTS6;ADAMTSL1;ADD1;ADK;AFF4;AFG3L2;AGFG1;AGO3;AGO4;AGPAT5;AGPAT6;AKAP10;AKAP7;AKIRIN1;AKT2;ALAS1;ALG6;ALG9;ALKBH8;AMN1;AMZ2;ANAPC…
# 2 NN,PY        268 ACPL2;ADAL;AEBP2;AFF1;AGGF1;AHCTF1;ALMS1;ALS2;ALS2CR12;AP3B1;ARAP2;ARHGAP19-SLIT1;ARHGEF9;ARID2;ARIH1;ASCC1;ASH2L;ATF7;ATP8B4;BDP1;BICD1;BMPR2;BRD4;BTBD7;C11orf30;C2CD3;C2CD5;C9orf41;CAPN2;CAPN7;CCAR1;CCDC149;CCD…
# 3 NN,PY,SNDA   276 AAK1;ABCA3;ACSL3;ADAM32;ADD3;ADRBK2;AGTPBP1;AK5;AKAP13;AKT3;AMPH;ANKRD13C;ANKRD17;ANKRD27;APC;AQR;ARFGEF2;ARHGAP10;ARHGAP21;ARHGAP26;ARHGEF12;ARNT2;ARSG;ASAP1;ASPH;ATE1;ATF2;ATF6;ATP6V0A1;ATP8A1;ATP9B;ATXN1;ATXN1…
# 4 NN,SNDA      315 ABHD2;ABI1;ABI2;ACSL4;ADAM10;ADAM17;AFF3;AFTPH;AHI1;AKAP6;AKAP9;ALG8;ALS2CR11;AMOT;ANKRD32;ANO10;ANO6;ARFGEF1;ARHGAP6;ARHGEF7;ARMC8;ASH1L;ASXL2;ATF1;ATP11B;ATP2B4;ATP6V0A2;ATRX;B3GALNT2;BARD1;BAZ2B;BBS7;BCL2L13;B…
# 5 PY           382 ABCA5;ABCG2;ABLIM2;ACER2;ACP6;ACTR3B;ADAM9;ADAMTS17;ADAMTS19;ADCK4;ADCY1;AGBL3;AGBL4;AGK;AGO2;AIMP2;ANKMY2;ANKRD24;ANTXR1;ARHGAP32;ARHGAP5;ARPP21;ARRB1;ASCC2;ASTN2;ATG14;ATG2B;ATP10D;ATP2B1;B3GALTL;BAZ1B;BEND3;BE…
# 6 PY,SNDA      229 ACAD11;ACAP2;AGL;AMBRA1;ANKIB1;ANKS1A;ANKS1B;ANO4;ANO5;APBB2;APP;ARHGEF28;ARID1B;ARMC2;ARMC9;ARNTL2;ATG7;ATP2C1;ATP8A2;ATP9A;ATRNL1;BAI3;BTBD9;C11orf74;C18orf8;C9orf3;CACNA2D1;CACNA2D3;CACNB2;CADPS;CADPS2;CAMSAP2…
# 7 SNDA         496 ABLIM1;ACACA;ACACB;ADAMTSL3;ADARB1;ADCY2;ADCY9;AGAP1;AGPS;ALDH3A2;ANK2;ANK3;ANKAR;AP2B1;AP2M1;ARAP1;ARFGAP3;ARGLU1;ARHGEF10L;ARHGEF11;ARHGEF26;ARHGEF4;ARID4B;ASXL3;ATF6B;ATG5;ATP10B;ATP1A1OS;ATP1A3;ATP6V1H;ATP8B1…

## ===========
## cell-specificity of the host genes (in celltype3)
## ===========
gene_fpkm=read.table("~/neurogen/rnaseq_PD/results/merged/genes.fpkm.cuffnorm.allSamples.uniq.xls", header = T, row.names = 1, stringsAsFactors = F, check.names = F); head(gene_fpkm)
gene_fpkm=gene_fpkm[,colnames(Merge_circexp_norm_filtered_and_enriched)]; dim(gene_fpkm)
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
# 1 NN         472
# 2 PY        1348
# 3 SNDA      7626

gene_group3mean_DF3 = gene_group3mean[match(unique(DF3$hostgeneID), gene_group3mean$gene),]
gene_groupmean_s3_DF3 = groupmean_to_specificity(gene_group3mean_DF3)

## cell-specificity of circRNA VS. cell-specificity of hostgene
DF3 = inner_join(DF3, gene_groupmean_s3_DF3,by = c("hostgeneID" = "gene"), suffix=c(".circRNA",".gene")); head(DF3)
DF3 = mutate(DF3, celltype.specific.gene = ifelse(Private_or_not.gene==1, celltype.gene, "NS"))
dim(DF3)

for(i in c("NN","PY","SNDA")) {print(i); filter(DF3, celltype.circRNA==i, Private_or_not.gene==1) %>% select(hostgene, celltype.gene) %>% distinct()%>% group_by(celltype.gene) %>% summarise(n=n()) %>% print()}
# zero

for(i in c("NN","PY","SNDA")) {print(i); filter(DF3, celltype.circRNA==i) %>% select(hostgene, celltype.gene) %>% distinct() %>% group_by(celltype.gene) %>% summarise(n=n()) %>% print()}
ggplot(DF3, aes(x=factor(celltype.gene), fill=factor(celltype.circRNA))) +
  geom_bar(stat="count", width=0.7, fill="steelblue")+ 
  coord_flip() +
  facet_grid(cols=vars(celltype.circRNA), scales ='free') +
  theme_minimal()
ggsave("../results/Merge_circexplorer_BC.cellspecific_barplot++.hostgene.pdf", width = 10, height = 2)
# [1] "NN"
# celltype.gene     n
# 1 NN             1100
# 2 PY              629
# 3 SNDA            394

# [1] "PY"
# celltype.gene     n
# 1 NN              301
# 2 PY              596
# 3 SNDA            258

# [1] "SNDA"
# celltype.gene     n
# 1 NN              346
# 2 PY              555
# 3 SNDA            415

## spec.circRNA vs. spec.gene for the same cell type
xx=DF3 %>% mutate(x=match(paste0(celltype.circRNA,"_spec.gene"), colnames(DF3))) 
xx = xx %>% mutate(S.circRNA.gene = DF3[cbind(seq_along(x), x)]) %>% select(celltype.circRNA, S.circRNA, S.circRNA.gene) %>% gather("circRNA_or_gene", "specificity_score", -1, convert = T) %>% mutate(circRNA_or_gene=as.factor(circRNA_or_gene), celltype.circRNA=as.factor(celltype.circRNA), specificity_score=as.numeric(specificity_score))
ggplot(xx, aes(x=celltype.circRNA, y=specificity_score, fill=circRNA_or_gene)) + 
  geom_violin(trim = T, scale = "width", position = position_dodge(width = .51)) + 
  coord_flip() +
  theme_bw()
  #geom_boxplot(col='white',width=0.05, outlier.shape = NA, position = position_dodge(width = .51)) 
  #geom_jitter(height = 0, width = 0.1)
  #geom_boxplot() +
  #geom_point(aes(fill = circRNA_or_gene), size = 5, shape = 21, position = position_jitterdodge())
ggsave("../results/Merge_circexplorer_BC.cellspecific_boxplot++.specificityscore.pdf")

## heatmap of host genes of cell-specific circRNAs
X=gene_group5mean[match(DF3$hostgeneID, gene_group5mean$gene),]; head(X); dim(X); 
X=log10(X[,-1]+0.01)
dissimilarity <- 1 - cor(X)
distance <- as.dist(dissimilarity)
hc=hclust(distance)
hccol=as.dendrogram(hc)
weights.dd <- order(match(c("PBMC","FB","SNDA","MCPY","TCPY"), colnames(X))) # the trick is to call order() on the specific index of target dendragram
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
                      cutree_cols = 5,
                      cluster_rows = F, cluster_cols = as.hclust(reorder(hccol, wts = weights.dd, agglo.FUN = mean)))
library(pheatmap);
do.call("pheatmap", hm.parameters)
do.call("pheatmap", c(hm.parameters, width=3, height=5, filename="../results/Merge_circexplorer_BC.cellspecific_heatmap++.hostgene.pdf"))

#############################################
# cell-specific isoform of circRNAs, esp. those from the same host gene (esp. the 505 shared genes btw SNDA and PY)
#############################################
Merge_circexp_raw_filtered_and_enriched = readRDS(file="Merge_circexplorer_BC106.filtered.enriched.rawcount.rds")
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
            genen = length((gene)),
            RPMs_SNDA_PY_NN=paste(RPM_SNDA_PY_NN, collapse = ','),
            max3nonzeros = paste((max3nonzero), collapse = ','),
            max3mean=max(max3))  %>%
  arrange(-max3mean) %>% 
  filter(celltype3n.circRNA>1, hostgene %in% c('PHF21A','SETD2','ERC1','MAN2A1','PSEN1','ATXN10')) 









