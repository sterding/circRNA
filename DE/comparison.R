## AD vs. PD DE comparison
library(tidyverse)

setwd("~/projects/circRNA/results/")
AD=read.table("DE_TCPY/DEresult.DE_TCPY.CONDITION_AD_vs_HC.xls", stringsAsFactors = F, row.names = 1, header=T, sep = "\t")
head(AD)
PD=read.table("DE_SNDA/DEresult.DE_SNDA.CONDITION2_PD_vs_HC.xls", stringsAsFactors = F, row.names = 1, header=T, sep = "\t")
head(PD)
dim(AD); dim(PD)
library(gplots) # install.packages('gplots')
venn(list(AD = rownames(AD), PD = rownames(PD)))
# the direction of the 264 shared genes
ADPD = inner_join(rownames_to_column(AD) %>% filter(pvalue<=1) %>% select(1:8,13:14), 
           rownames_to_column(PD) %>% filter(pvalue<=1) %>% select(1:8,13:14,9:12), 
           by = 'rowname', suffix=c(".AD",".PD"))

# DE genes in either comparison
filter(ADPD, pvalue.AD<=0.05 | pvalue.PD<=0.05) %>% with(table(log2FoldChange.AD>0, log2FoldChange.PD>0))
df = filter(ADPD, pvalue.AD<=0.05 | pvalue.PD<=0.05) %>% mutate(log2FoldChangeAD = log2FoldChange.AD>0, log2FoldChangePD = log2FoldChange.PD>0) %>%
  arrange(log2FoldChangeAD, log2FoldChangePD, log2FoldChange.AD) %>% 
  #unite('genename', c(rowname, " (",geneName, ")"), sep = "") %>% 
  mutate(genename = paste0(rowname, " (",geneName, ")")) %>%
  select(genename, AD=log2FoldChange.AD, PD=log2FoldChange.PD) 
df$genename= factor(df$genename, levels = as.character(df$genename))
gather(df, key="ADPD", value="log2FoldChange", 2:3) %>%
  ggplot(aes(x=genename, y=log2FoldChange, fill=ADPD)) +
  geom_bar(stat="identity", position=position_dodge()) + coord_flip()  + 
  theme_minimal() + theme(legend.justification=c(1,0), legend.position=c(1,0))
ggsave(filename = "ADPD.DE.barplot.pdf", width = 6.5, height = 9)


# DE genes in both comparisons
venn(list(AD = rownames(subset(AD, pvalue<=0.05)), PD = rownames(subset(PD, pvalue<=0.05))))

# the direction of the 1 shared DE genes
filter(ADPD, pvalue.AD<=0.05, pvalue.PD<=0.05) 
# DNAJC6