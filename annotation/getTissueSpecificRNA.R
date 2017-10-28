#R script to get cell type specific RNAs/genes

setwd("~/projects/circRNA/data")
# calculate specificity score (S) as cummeRbund (https://github.com/Bioconductor-mirror/cummeRbund/blob/master/R/tools.R)
source("~/projects/circRNA/src/tools.R")
load("Merge_circexp.BC.Rdata")

df=Merge_circexp_norm_filtered
head(df)

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
# celltype count
# <chr> <int>
# 1       NN  1623
# 2       PY   722
# 3     SNDA   458
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
# celltype count
# <chr> <int>
# 1       FB   238
# 2     MCPY   169
# 3     PBMC   974
# 4     SNDA   427
# 5     TCPY   441
head(df)
df5=df

# heatmpa of cell specific genes
library('pheatmap')
df=log10(groupmean[df5$gene,]*1000+1)

dissimilarity <- 1 - cor(df)
distance <- as.dist(dissimilarity)
hc=hclust(distance)
hccol=as.dendrogram(hc)
weights.dd <- order(match(c("FB","PBMC","SNDA","TCPY","MCPY"), celltypes)) # the trick is to call order() on the specific index of target dendragram
#plot(reorder(hccol, wts = weights.dd))

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
do.call("pheatmap", c(hm.parameters, width=3, height=5, filename="../results/cellspecific_heatmap.pdf"))

## cell-specific circRNAs in both major and minor groups
df=log10(groupmean[unique(df3$gene,df5$gene),]*1000+1)
dissimilarity <- 1 - cor(df)
distance <- as.dist(dissimilarity)
hc=hclust(distance)
hccol=as.dendrogram(hc)
weights.dd <- order(match(c("FB","PBMC","SNDA","TCPY","MCPY"), celltypes)) # the trick is to call order() on the specific index of target dendragram
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
do.call("pheatmap", c(hm.parameters, width=3, height=5, filename="../results/cellspecific_heatmap+.pdf"))
