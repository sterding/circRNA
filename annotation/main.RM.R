## Main script to compare RNase R vs. mock treated samples
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)

## See ~/projects/circRNA/README.Rmd for how to generate the input files
setwd("~/projects/circRNA/data/")

## validae in RM experiment

# load
annotation<- read.table("Merge_circexplorer_RM.annotation.bed14",header=T, check.names = F) # include both ciRNA and circRNA
head(annotation);dim(annotation)
Merge_circexp_raw <- read.table("Merge_circexplorer_RM.rawcount.txt", header = T, check.names = F, sep="\t", comment.char ='') 
colnames(Merge_circexp_raw)=gsub("#|_candidates.bed_circReads","",colnames(Merge_circexp_raw))

# wide --> long
Merge_circexp_raw_long = Merge_circexp_raw %>% 
  mutate(name=paste(chrom,start,end,sep="_")) %>% 
  select(name, dplyr::contains("_")) %>%
  melt(id.vars='name', value.name='raw_count') %>%
  separate(variable, c("Dx","subjectID","cellType","Treatment"), sep="_", extra='drop') %>%
  arrange(cellType) %>%
  unite(sampleID,Dx,subjectID,cellType) %>%
  dcast(name + sampleID ~ Treatment, value.var='raw_count')

# #OPTION1: We know HC_MGH1026_SN failed and informally we can use it as technical replicates to infer the background noise for other samples
# # [Later, Clemens disgards this idea as he thought it's a failed experiment]
# baseline = Merge_circexp_raw_long %>% filter(sampleID == 'HC_MGH1026_SN') %>% select(M, R) %>% rbind(c(max(Merge_circexp_raw_long$R),max(Merge_circexp_raw_long$R))) %>% distinct() %>% mutate(M=M+1,R=R+1) %>% group_by(M) %>% summarise(R=max(2*M,max(R))) %>% as.data.frame()
# base.line=with(baseline, loess(R~M))
# bl=data.frame(x=baseline$M, y=predict(base.line,baseline$M))
# 
# # scatter plot: R vs. M
# ggplot(filter(Merge_circexp_raw_long, R>0 | M>0) %>% select(sampleID, R, M) %>% distinct() %>% na.omit() %>% mutate(M=M+1,R=R+1) %>% mutate(col=ifelse(predict(base.line, M)<=R,'enriched','non-enriched')), aes(x=M, y=R, color=col)) + 
#   geom_point(alpha=1/5, size=0.5) +
#   #geom_smooth(aes(M,R), data=baseline, method='loess', se=F, color='blue') +
#   geom_line(mapping=aes(x,y), data=bl, color = "red",size=.5, linetype = 1) + 
#   geom_abline(slope=1, intercept=0, color = "black",size=.25, linetype = 2) + 
#   #geom_abline(slope=1, intercept=log10(2), color = "red",size=.5, linetype = 2) + 
#   scale_y_log10() + scale_x_log10()+ 
#   coord_fixed(xlim=range(Merge_circexp_raw_long$M, Merge_circexp_raw_long$R)+1, ylim=range(Merge_circexp_raw_long$M, Merge_circexp_raw_long$R)+1) +
#   scale_color_manual(values=c("blue", "darkgray")) + theme_bw() +
#   ylab("raw reads with RNase R treatment") + xlab("raw reads without RNase R treatment") +
# #  annotate("text", x = 1, y = max(R+1), label = paste("median enrichment fold =", round(median((R+1)/(M+1))),2),parse = TRUE) +
#   facet_wrap( ~ sampleID, nrow=2)
# 
# ggsave("../results/RM.scatterplot1.pdf", width = 8, height = 5)

# OPTION2: Using t-test from the 6 vs. 6 pairs, assuming the enrichment of R vs. M is indepedant from the cell types
Merge_circexp_raw_RM_ttest= ungroup(Merge_circexp_raw_long) %>%
  filter(sampleID != 'HC_MGH1026_SN') %>%
  arrange(name, sampleID) %>% group_by(name) %>%
  summarise(M = list(M), R=list(R)) %>%
  group_by(name) %>%
  mutate(pvalue=t.test(R[[1]], M[[1]], alternative ='greater', paired =T)$p.value,
         pvalue.u=wilcox.test(R[[1]], M[[1]], alternative ='greater', paired =T)$p.value,
         log2FoldChange=log2((sum(R[[1]])+1)/(sum(M[[1]])+1)),
         symbol=annotation$geneName[match(name, annotation$ID)],
         meanR=mean(R[[1]])) 
# save
saveRDS(Merge_circexp_raw_RM_ttest, file="Merge_circexp_raw_RM_ttest.rds")


# ggplot(filter(Merge_circexp_raw_long, R>0 | M>0) %>% na.omit() %>% mutate(col=ifelse(name %in% Merge_circexp_raw_long_filtered$name, 'enriched','non-enriched')) %>% select(sampleID, M,R,col) %>% distinct() %>% arrange(desc(col)) %>% mutate(M=M+1,R=R+1), aes(x=M, y=R, color=col)) + 
#   geom_point(alpha=1/5, size=0.5) +
#   geom_abline(slope=1, intercept=0, color = "black",size=.5, linetype = 2) + 
#   geom_abline(slope=1, intercept=log10(2), color = "red",size=.5, linetype = 2) + 
#   geom_hline(yintercept=50, color = "red",size=.5, linetype = 2) + 
#   scale_y_log10() + scale_x_log10()+ 
#   coord_fixed(xlim=range(Merge_circexp_raw_long$M, Merge_circexp_raw_long$R)+1, ylim=range(Merge_circexp_raw_long$M, Merge_circexp_raw_long$R)+1) +
#   scale_color_manual(values=c("blue", "darkgray")) + theme_bw() +
#   ylab("raw reads with RNase R treatment") + xlab("raw reads without RNase R treatment") +
#   facet_wrap( ~ sampleID, nrow=2)
# ggsave("../results/RM.scatterplot2.pdf", width = 8, height = 5)

## OPTION3: using a hard filter with R/M >=2 & R >=20  (The cutoff is chosen based on visual inspection of the scatterplot)
# scatter plot: R vs. M
ggplot(filter(Merge_circexp_raw_long, R>0 | M>0) %>% select(sampleID, R, M) %>% distinct() %>% na.omit() %>% mutate(col=ifelse(R>=2*M & R>=20,'enriched','non-enriched')) %>% mutate(M=M+1,R=R+1), aes(x=M, y=R, color=col)) + 
  geom_point(alpha=1/5, size=0.5) +
  geom_abline(slope=1, intercept=0, color = "black",size=.5, linetype = 2) + 
  geom_abline(slope=1, intercept=log10(2), color = "red",size=.5, linetype = 2) + 
  geom_hline(yintercept=20, color = "red",size=.5, linetype = 2) + 
  scale_y_log10() + scale_x_log10()+ 
  coord_fixed(xlim=range(Merge_circexp_raw_long$M, Merge_circexp_raw_long$R)+1, ylim=range(Merge_circexp_raw_long$M, Merge_circexp_raw_long$R)+1) +
  scale_color_manual(values=c("blue", "darkgray")) + theme_bw() +
  ylab("raw reads with RNase R treatment") + xlab("raw reads without RNase R treatment") +
  #  annotate("text", x = 1, y = max(R+1), label = paste("median enrichment fold =", round(median((R+1)/(M+1))),2),parse = TRUE) +
  facet_wrap( ~ sampleID, nrow=2)
ggsave("../results/RM.scatterplot3.pdf", width = 8, height = 5)

Merge_circexp_raw_long_filterd = filter(Merge_circexp_raw_long, R >= 2*M & R>=20) %>% 
  filter(sampleID %in% c('HC_H02018_PBMC','HC_ND34770_FB','HC_TCKY1217_TC', 'HC_TCKY1247_TC','HC_MC3290_SN', 'HC_MD6326_MD')) %>% 
  mutate(sampleGroup = sub(".*_(.*)","\\1",sampleID)) %>%
  mutate(sampleGroup3 = ifelse(grepl('FB|PBMC', sampleGroup),'NN', ifelse(grepl('SN|MB', sampleGroup),'SN','PY'))) %>%
  #select(name, sampleGroup3) %>% 
  arrange(sampleGroup3) %>% distinct() 

Merge_circexp_raw_long_filterd %>% select(name, sampleGroup3) %>% group_by(sampleGroup3) %>% summarise(n=n())
# # A tibble: 3 × 2
# sampleGroup3     n
# <chr> <int>
# 1           NN 19399
# 2           PY 37603
# 3           SN  3221

saveRDS(Merge_circexp_raw_long_filterd, file="Merge_circexp_raw_long_filterd.RM.rds")