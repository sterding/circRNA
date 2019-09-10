## Main script to compare RNase R vs. mock treated samples
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)

## See ~/projects/circRNA/README.Rmd for how to generate the input files
setwd("~/projects/circRNA/data/")

## validae in RM experiment

# load
Merge_circexp_raw_long <- read.table("Merge_circexplorer_RM12.rawcount.long.txt", sep="\t",check.names =F, col.names = c("ID","readsCount","sampleID")) %>% 
  group_by(ID, sampleID) %>% summarise(readsCount=sum(readsCount)) %>%
  separate(sampleID, c("Dx","subjectID","cellType","Treatment"), sep="_", extra='drop') %>%
  arrange(cellType) %>%
  unite(sampleID,Dx,subjectID,cellType) %>%
  dcast(ID + sampleID ~ Treatment, value.var='readsCount', fill = 0)
head(Merge_circexp_raw_long)

## OPTION3: using a hard filter with R/M >=2 & R >=20  (The cutoff is chosen based on visual inspection of the scatterplot)
# scatter plot: R vs. M
ggplot(filter(Merge_circexp_raw_long, R>0 | M>0) %>% select(sampleID, R, M) %>% distinct() %>% na.omit() %>% mutate(col=ifelse(R>=2*M & R>=20,'enriched','non-enriched')) %>% mutate(M=M+1,R=R+1), aes(x=M, y=R, color=col)) + 
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

saveRDS(Merge_circexp_raw_long_filterd, file="Merge_circexp_raw_long_filterd.RM.rds")
