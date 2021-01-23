## Main script to compare RNase R vs. mock treated samples
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)

## See ~/projects/circRNA/src/annotation/main.sh for how to generate the input files
setwd("~/projects/circRNA/data/")

## load
Merge_circexp_raw = read.table("Merge_circexplorer_RM12.rawcount.long.txt", sep="\t",check.names =F, col.names = c("ID","readsCount","sampleID")) %>%
  group_by(ID, sampleID) %>% summarise(readsCount=sum(readsCount)) %>%
  spread(key = sampleID, value = readsCount, fill = 0) %>%
  column_to_rownames("ID")
# to long format
Merge_circexp_raw_long = Merge_circexp_raw %>% rownames_to_column(var = "ID") %>% gather(key = 'sampleID',value = 'readsCount',starts_with("HC_")) %>% 
  separate(sampleID, c("Dx","subjectID","cellType","Treatment"), sep="_", extra='drop') %>%
  arrange(cellType) %>%
  unite(sampleID,Dx,subjectID,cellType) %>%
  dcast(ID + sampleID ~ Treatment, value.var='readsCount', fill = 0)
head(Merge_circexp_raw_long)

# normalized RPM for treated samples
library(googlesheets4) # install.packages("googlesheets4")
RM_RNAseq = read_sheet("16xwR1MjCJKzkhfEk3yfjpf9JEOPUWnkLuL0HRdGT6XY", sheet="RM.RNAseq") # google table: https://docs.google.com/spreadsheets/d/16xwR1MjCJKzkhfEk3yfjpf9JEOPUWnkLuL0HRdGT6XY/edit#gid=1071178719
readsNum_million <- as.numeric(gsub(",","",RM_RNAseq$total_reads))/10^6; names(readsNum_million) = RM_RNAseq$SampleID
readsNum_million = readsNum_million[colnames(Merge_circexp_raw)]
Merge_circexp_norm<-sweep(Merge_circexp_raw,2,readsNum_million,"/")

## to get group3mean, then groupmean_s3 specifcity score, from there to define cell-specific circRNAs 
# 3 groups mean # limit to R groups
group3mean = select(Merge_circexp_norm, contains("_R_")) %>% rownames_to_column(var='gene') %>% 
  gather(key = 'sampleID',value = 'RPM',starts_with("HC_")) %>%
  separate(sampleID, c("Dx","subjectID","cellType","Treatment"), sep="_", extra='drop') %>%
  filter(cellType!="MB") %>%
  mutate(cellType=ifelse(cellType %in% c("MB","SN"), "SN", ifelse(cellType=="TC","TC","NN"))) %>%
  group_by(gene, cellType) %>%
  summarise(meanRPM=mean(RPM)) %>% 
  mutate(cellType=factor(cellType, levels = c("SN","TC","NN"))) %>% # redefine factor of cellType
  dcast(gene ~ cellType, value.var='meanRPM')
head(group3mean)
source("~/projects/circRNA/src/annotation/tools.R")
groupmean_s3 = groupmean_to_specificity(group3mean)
head(groupmean_s3)
# Questions to ask:
# 1. Consistence of the qPCR validated ones in cell-specificity
library(googlesheets4)
qPCR_validated = read_sheet("16xwR1MjCJKzkhfEk3yfjpf9JEOPUWnkLuL0HRdGT6XY", sheet="Table S2. qPCR.circRNAs") # https://docs.google.com/spreadsheets/d/16xwR1MjCJKzkhfEk3yfjpf9JEOPUWnkLuL0HRdGT6XY/edit
write_sheet(left_join(qPCR_validated, rename_if(group3mean, is.numeric, list(~str_c("mean.RNaseR.", .))), by=c("circRNA_ID" = "gene")), ss="16xwR1MjCJKzkhfEk3yfjpf9JEOPUWnkLuL0HRdGT6XY", sheet="Table S2. qPCR.circRNAs")
# join group3mean and groupmean_s3
left_join(qPCR_validated, rename_if(group3mean, is.numeric, list(~str_c("mean.lcRNAseq.", .))), by=c("circRNA_ID" = "gene")) %>% 
  left_join(y=groupmean_s3, by = c("circRNA_ID" = "gene")) %>% print(width=Inf)


# 2. Consistence of cell-specificity from this RM data vs. cell-specificty from the lcRNAseq data?
DF3=read.table(file = "../results/Merge_circexplorer_BC109.cellspecificity.circRNA3.genes3.xls",sep="\t", header = T)
dim(DF3)
left_join(select(DF3, gene, celltype.circRNA), select(groupmean_s3, gene, celltype), by="gene") %>% select(celltype, celltype.circRNA) %>% table()
#         celltype.circRNA
# celltype   NN   PY SNDA
# NN 2014  292  219
# SN  141   69   45
# TC 2705 2947 1262

Merge_circexplorer_RM12.annotation = read.table("Merge_circexplorer_RM12.annotation.bed14", header = T, stringsAsFactors = F); 
head(Merge_circexplorer_RM12.annotation);
dim(Merge_circexplorer_RM12.annotation)
## OPTION3: using a hard filter with R/M >=2 & R >=20  (The cutoff is chosen based on visual inspection of the scatterplot)
# scatter plot: R vs. M
df = Merge_circexp_raw_long %>% mutate(circType=Merge_circexplorer_RM12.annotation$circType[match(ID, Merge_circexplorer_RM12.annotation$ID)]) %>% 
  mutate(col=ifelse(R>=2*M & R>=20,'enriched','non-enriched')) %>% mutate(M=M+1,R=R+1)

ggplot(select(df, sampleID, R, M, circType) %>% distinct(), aes(x=M, y=R, color=col)) + 
  geom_point(alpha=1/5, size=0.5) +
  geom_abline(slope=1, intercept=0, color = "black",size=.5, linetype = 2) + 
  geom_abline(slope=1, intercept=log10(2), color = "red",size=.5, linetype = 2) + 
  geom_hline(yintercept=20, color = "red",size=.5, linetype = 2) + 
  scale_y_log10() + scale_x_log10()+ 
  coord_fixed(xlim=range(Merge_circexp_raw_long$M, Merge_circexp_raw_long$R)+1, ylim=range(Merge_circexp_raw_long$M, Merge_circexp_raw_long$R)+1) +
  scale_color_manual(values=c("blue", "darkgray"))  + theme_bw() + theme(legend.position="top") + 
  ylab("raw reads with RNase R treatment") + xlab("raw reads without RNase R treatment") +
  #  annotate("text", x = 1, y = max(R+1), label = paste("median enrichment fold =", round(median((R+1)/(M+1))),2),parse = TRUE) +
  facet_wrap( ~ sampleID, nrow=2)
ggsave("../results/RM12.scatterplot.pdf", width = 8, height = 5)

## Are enriched circular RNAs more likely be to exon-derived circRNAs?
df = left_join(Merge_circexplorer_RM12.annotation, filter(df, col=='enriched') %>% select(ID, col) %>% distinct(), by = 'ID') %>% replace_na(list(col = "non-enriched"))
dim(df)
with(df, table(circType, col))
fisher.test(with(df, table(circType, col)))

Merge_circexp_raw_long_filterd = filter(Merge_circexp_raw_long, R >= 2*M & R>=20) %>% 
  filter(sampleID %in% c('HC_H02018_PBMC','HC_ND34770_FB','HC_TCKY1217_TC', 'HC_TCKY1247_TC','HC_MC3290_SN', 'HC_MD6326_MD')) %>% 
  mutate(sampleGroup = sub(".*_(.*)","\\1",sampleID)) %>%
  mutate(sampleGroup3 = ifelse(grepl('FB|PBMC', sampleGroup),'NN', ifelse(grepl('SN|MB', sampleGroup),'SN','PY'))) %>%
  #select(name, sampleGroup3) %>% 
  arrange(sampleGroup3) %>% distinct() 

Merge_circexp_raw_long_filterd %>% select(ID, sampleGroup3) %>% group_by(sampleGroup3) %>% summarise(n=n())
# # A tibble: 3 × 2
# sampleGroup3     n
# <chr> <int>
# 1           NN 24383
# 2           PY 51478
# 3           SN  3348

## how many enriched in brain?  # this number goes to Fig. 1a
filter(Merge_circexp_raw_long_filterd, sampleGroup3!="NN") %>% pull(ID) %>% as.character() %>% unique() %>% length()
# [1] 39875

saveRDS(Merge_circexp_raw_long_filterd, file="Merge_circexp_raw_long_filterd.RM.rds")
