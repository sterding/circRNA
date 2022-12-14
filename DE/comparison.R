## AD vs. PD DE comparison
library(tidyverse)
library(gplots) # install.packages('gplots')
suppressPackageStartupMessages(library('ggrepel',logical.return=T) || install.packages('ggrepel', repo='http://cran.revolutionanalytics.com'))

setwd("~/projects/circRNA/results/")
genes_annotation = read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v19.annotation.genes.bed", sep="\t", quote="", header = F, stringsAsFactors = F, col.names = c("chr","start","end","geneName","geneID","geneType","strand"));

### =======================
## DE_TCPY vs. DE_SNDA
### =======================
df1=read.table("DE_TCPY/DEresult.DE_TCPY.CONDITION_AD_vs_HC.xls.gz", stringsAsFactors = F, row.names = 1, header=T, sep = "\t")
df2=read.table("DE_SNDA/DEresult.DE_SNDA.CONDITION2_PD_vs_HC.xls.gz", stringsAsFactors = F, row.names = 1, header=T, sep = "\t")
df1df2 = inner_join(rownames_to_column(df1) %>% filter(pvalue<=1) %>% select(1:8), 
                    rownames_to_column(df2) %>% filter(pvalue<=1) %>% select(1:12), 
                    by = 'rowname') # %>% mutate(log2FoldChange.y=-1*log2FoldChange.y) # sign of fold change for continous variable?
head(df1df2)
# add circID
df1df2$circID=paste0(sub("RNA","",as.character(df1df2$circType)),as.character(df1df2$geneName))

pdf(file = "ADPD.DE.barplot.pdf", paper = 'US')

# Intersection of genes qualified for DE analysis
gplots::venn(list(AD = rownames(df1), PD = rownames(df2)))

# correlation between the them
cortest = cor.test(df1df2$log2FoldChange.x, df1df2$log2FoldChange.y,alternative = "two.sided", method = "pearson")

ggplot(df1df2, aes(x=log2FoldChange.x, y=log2FoldChange.y, label=paste(circID, rowname, sep = "\n"))) +
  geom_hline(yintercept = 0, color='black', size=.5, linetype = 2) + 
  geom_vline(xintercept = 0, color='black', size=.5, linetype = 2) +
  geom_text_repel(
    data          = subset(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05),
    nudge_y       = 0,
    vjust         = 1,
    segment.size  = 0.2,
    segment.color = "grey50") +
  geom_point(color = ifelse(df1df2$pvalue.x<=0.05 | df1df2$pvalue.y<=0.05, "red", "black"), size=ifelse(df1df2$pvalue.x<=0.05 | df1df2$pvalue.y<=0.05, 3, 1)) +
  labs(title="Correlation coeffecient analysis (AD case vs. PD case)", 
       x="log2FoldChange (AD vs. HC)", 
       y="log2FoldChange (PD vs. HC)",
       subtitle = paste("Pearson's r =", round(cortest$estimate, 2), "(p-value:", format.pval(cortest$p.value),")")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# significant in either comparison
filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% with(table(log2FoldChange.x>0, log2FoldChange.y>0))
# save to excel
filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% mutate(Chr=genes_annotation$chr[match(geneID, genes_annotation$geneID)]) %>% 
  select(CircRNA=rowname, Chr, log2FoldChange.AD = log2FoldChange.x, pvalue.AD = pvalue.x, log2FoldChange.PD = log2FoldChange.y, pvalue.PD = pvalue.y) %>%
  write.table(file="ADPD.DE.xls", sep="\t", na="", row.names=F)

# DE genes in both comparisons
gplots::venn(list(AD = rownames(subset(df1, pvalue<=0.05)), PD = rownames(subset(df2, pvalue<=0.05))))

# barplot
df = filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% mutate(log2FoldChangeX = log2FoldChange.x>0, log2FoldChangeY = log2FoldChange.y>0) %>%
  arrange(log2FoldChangeX, log2FoldChangeY, log2FoldChange.x) %>% 
  mutate(genename = paste0(rowname," (",geneName,")")) %>%
  select(genename, AD_vs_HC=log2FoldChange.x, PD_vs_HC=log2FoldChange.y) 
df$genename= factor(df$genename, levels = as.character(df$genename))
p= ggplot(gather(df, key="comparison", value="log2FoldChange", 2:3), aes(x=genename, y=log2FoldChange, fill=comparison)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.5) + coord_flip()  + 
  # geom_text(aes(label=paste0(observed," (",round(OR,1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5) +  # TODO: add significance *
  theme_minimal() + theme(legend.justification=c(1,0), legend.position=c(1,0))
print(p)

dev.off()

### =======================
## DE2gene_TCPY vs. DE2gene_SNDA
### =======================
df1=read.table("DE2gene_TCPY/DEresult.DE2gene_TCPY.CONDITION_AD_vs_HC.xls.gz", stringsAsFactors = F, row.names = 1, header=T, sep = "\t")
df2=read.table("DE2gene_SNDA/DEresult.DE2gene_SNDA.CONDITION2_PD_vs_HC.xls.gz", stringsAsFactors = F, row.names = 1, header=T, sep = "\t")
df1df2 = inner_join(rownames_to_column(df1) %>% filter(pvalue<=1) %>% select(1:8),
                    rownames_to_column(df2) %>% filter(pvalue<=1) %>% select(1:13),
                    by = 'rowname') # %>% mutate(log2FoldChange.y=-1*log2FoldChange.y) # sign of fold change for continous variable?
## Actually, if including circRNAs tested in either AD or PD, we should use full_join (below), which will be 820 circRNA genes in total. In this figure, we only include those tested in BOTH AD and PD. That would be 423. 
# df1df2 = full_join(rownames_to_column(df1) %>% filter(pvalue<=1) %>% select(1:8), 
#                     rownames_to_column(df2) %>% filter(pvalue<=1) %>% select(1:13), 
#                     by = 'rowname') %>% mutate(pvalue.x = ifelse(is.na(pvalue.x), 1, pvalue.x),
#                                                pvalue.y = ifelse(is.na(pvalue.y), 1, pvalue.y),
#                                                FoldChange.x = ifelse(is.na(FoldChange.x), 0, FoldChange.x),
#                                                FoldChange.y = ifelse(is.na(FoldChange.y), 0, FoldChange.y))
dim(df1df2)
# test if AD and PD associated circRNAs are significantly overlapped?
fisher.test(with(df1df2, table(pvalue.x<=0.05, pvalue.y<=0.05)))
# p-value =1 

pdf(file = "ADPD.DE2gene.barplot.pdf", paper = 'US')

# Intersection of genes qualified for DE analysis
gplots::venn(list(AD = rownames(df1), PD = rownames(df2)))

# correlation between the them
cortest = cor.test(df1df2$log2FoldChange.x, df1df2$log2FoldChange.y,alternative = "two.sided", method = "pearson")
plot(df1df2$log2FoldChange.x, df1df2$log2FoldChange.y, 
     xlab="log2FoldChange (AD vs. HC)", ylab="log2FoldChange (PD vs. HC)", 
     main="Correlation coeffecient analysis (AD case vs. PD case)")
mtext(paste("Pearson's r =", round(cortest$estimate, 2), "(p-value:", format.pval(cortest$p.value),")"))

# significant in either comparison
filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% with(table(log2FoldChange.x>0, log2FoldChange.y>0))
# save to excel
filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% mutate(Chr=genes_annotation$chr[match(geneID, genes_annotation$geneID)]) %>% 
  select(CircRNA=rowname, Chr, log2FoldChange.AD = log2FoldChange.x, pvalue.AD = pvalue.x, log2FoldChange.PD = log2FoldChange.y, pvalue.PD = pvalue.y) %>%
  write.table(file="ADPD.DE2gene.xls", sep="\t", na="", row.names=F)

# DE genes in both comparisons
gplots::venn(list(AD = rownames(subset(df1, pvalue<=0.05)), PD = rownames(subset(df2, pvalue<=0.05))))

# barplot
df = filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% mutate(log2FoldChangeX = log2FoldChange.x>=0, log2FoldChangeY = log2FoldChange.y>=0) %>%
  arrange(log2FoldChangeX, -log2FoldChangeY, log2FoldChange.x) %>% 
  mutate(genename = rowname) %>%
  select(genename, log2FoldChange.x, log2FoldChange.y, pvalue.x, pvalue.y) 
df$genename= factor(df$genename, levels = as.character(df$genename))

df = df %>% gather(key, value, -genename) %>% 
  extract(key, c("merit", "comparison"), "(.*)\\.(.*)") %>% 
  spread(merit, value) %>% 
  mutate(comparison = case_when(comparison=="x" ~ "AD_vs_HC", comparison=="y" ~ "PD_vs_HC"))

p= df %>% ggplot(aes(x=genename, y=log2FoldChange, fill=comparison, alpha=pvalue<0.05)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.5) + coord_flip()  + 
  # geom_text(aes(label=paste0(observed," (",round(OR,1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5) +  # TODO: add significance *
  theme_minimal() + theme(legend.justification=c(1,0), legend.position=c(1,0))
print(p)

dev.off()

### =======================
## DE2gene_TCPY: AD vs. Braak
### =======================
df1=read.table("DE2gene_TCPY/DEresult.DE2gene_TCPY.CONDITION_AD_vs_HC.xls.gz", stringsAsFactors = F, row.names = 1, header=T, sep = "\t")
df2=read.table("DE2gene_TCPY/DEresult.DE2gene_TCPY.Braak_Braak_stage.xls.gz", stringsAsFactors = F, row.names = 1, header=T, sep = "\t")
df1df2 = inner_join(rownames_to_column(df1) %>% filter(pvalue<=1) %>% select(1:8), 
                    rownames_to_column(df2) %>% filter(pvalue<=1) %>% select(1:13), 
                    by = 'rowname') %>%
  mutate(log2FoldChange.y=-1*log2FoldChange.y) # sign of fold change for continous variable?
head(df1df2)

pdf(file = "ADBraak.DE2gene.barplot.pdf", paper = 'US')

# Intersection of genes qualified for DE analysis
gplots::venn(list(AD = rownames(df1), Braak = rownames(df2)))

# correlation between the them
cortest = cor.test(df1df2$log2FoldChange.x, df1df2$log2FoldChange.y,alternative = "two.sided", method = "pearson")
ggplot(df1df2, aes(x=log2FoldChange.x, y=log2FoldChange.y, label=rowname)) +
  geom_point(shape=ifelse(df1df2$pvalue.x<=0.05 | df1df2$pvalue.y<=0.05, 16, 16),
             color = ifelse(df1df2$pvalue.x<=0.05 | df1df2$pvalue.y<=0.05, "red", "gray"), 
             size=ifelse(df1df2$pvalue.x<=0.05 | df1df2$pvalue.y<=0.05, 3, 3)) +
  labs(x="log2FoldChange (AD vs. HC)", y="log2FoldChange (Braak Braak stage)", 
       title = "Correlation coeffecient analysis (AD case vs. Braak score)",
       subtitle = paste("Pearson's r =", round(cortest$estimate, 2), "(p-value:", format.pval(cortest$p.value),")")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_repel(
    data = subset(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05),
    nudge_y       = 0,
    vjust         = 1,
    direction = 'both',
    segment.size  = 0.2,
    segment.color = "red")

# significant in either comparison
filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% with(table(log2FoldChange.x>0, log2FoldChange.y>0))
# save to excel
filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% mutate(Chr=genes_annotation$chr[match(geneID, genes_annotation$geneID)]) %>% 
  select(CircRNA=rowname, Chr, log2FoldChange.AD = log2FoldChange.x, pvalue.AD = pvalue.x, log2FoldChange.Braak = log2FoldChange.y, pvalue.Braak = pvalue.y) %>%
  write.table(file="ADBraak.DE2gene.xls", sep="\t", na="", row.names=F)

# DE genes in both comparisons
gplots::venn(list(AD = rownames(subset(df1, pvalue<=0.05)), Braak = rownames(subset(df2, pvalue<=0.05))))

# barplot
df = filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% mutate(log2FoldChangeX = log2FoldChange.x>0, log2FoldChangeY = log2FoldChange.y>0) %>%
  arrange(log2FoldChangeX, log2FoldChangeY, log2FoldChange.x) %>% 
  mutate(genename = rowname) %>%
  select(genename, AD_vs_HC=log2FoldChange.x, Braak_Braak_stage=log2FoldChange.y) 
df$genename= factor(df$genename, levels = as.character(df$genename))
p= ggplot(gather(df, key="comparison", value="log2FoldChange", 2:3), aes(x=genename, y=log2FoldChange, fill=comparison)) +
  geom_bar(stat="identity", position=position_dodge()) + coord_flip()  + 
  # geom_text(aes(label=paste0(observed," (",round(OR,1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5) +  # TODO: add significance *
  theme_minimal() + theme(legend.justification=c(1,0), legend.position=c(1,0))
print(p)

dev.off()

### =======================
## DE2gene_SNDA: PD vs. MUSS
### =======================
df1=read.table("DE2gene_SNDA/DEresult.DE2gene_SNDA.CONDITION_PD_vs_HC.xls.gz", stringsAsFactors = F, row.names = 1, header=T, sep = "\t")
df2=read.table("DE2gene_SNDA/DEresult.DE2gene_SNDA.MUSS.xls.gz", stringsAsFactors = F, row.names = 1, header=T, sep = "\t")
df1df2 = inner_join(rownames_to_column(df1) %>% filter(pvalue<=1) %>% select(1:8), 
                    rownames_to_column(df2) %>% filter(pvalue<=1) %>% select(1:13), 
                    by = 'rowname') # %>% mutate(log2FoldChange.y=-1*log2FoldChange.y) # sign of fold change for continous variable?
head(df1df2)

pdf(file = "PDMUSS.DE2gene.barplot.pdf", paper = 'US')

# Intersection of genes qualified for DE analysis
gplots::venn(list(PD = rownames(df1), MUSS = rownames(df2)))

# correlation between the them
cortest = cor.test(df1df2$log2FoldChange.x, df1df2$log2FoldChange.y,alternative = "two.sided", method = "pearson")
ggplot(df1df2, aes(x=log2FoldChange.x, y=log2FoldChange.y, label=rowname)) +
  geom_point(shape=ifelse(df1df2$pvalue.x<=0.05 | df1df2$pvalue.y<=0.05, 16, 16),
             color = ifelse(df1df2$pvalue.x<=0.05 | df1df2$pvalue.y<=0.05, "red", "gray"), 
             size=ifelse(df1df2$pvalue.x<=0.05 | df1df2$pvalue.y<=0.05, 3, 3)) +
  labs(title="Correlation coeffecient analysis (PD case vs. MUSS)", 
       x="log2FoldChange (PD vs. HC)", 
       y="log2FoldChange (MUSS)",
       subtitle = paste("Pearson's r =", round(cortest$estimate, 2), "(p-value:", format.pval(cortest$p.value),")")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_repel(
    data = subset(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05),
    nudge_y       = 0,
    vjust         = 1,
    direction = 'both',
    segment.size  = 0.2,
    segment.color = "red")

# significant in either comparison
filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% with(table(log2FoldChange.x>0, log2FoldChange.y>0))
# save to excel
filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% mutate(Chr=genes_annotation$chr[match(geneID, genes_annotation$geneID)]) %>% 
  select(CircRNA=rowname, Chr, log2FoldChange.PD = log2FoldChange.x, pvalue.PD = pvalue.x, log2FoldChange.MUSS = log2FoldChange.y, pvalue.MUSS = pvalue.y) %>%
  write.table(file="PDMUSS.DE2gene.xls", sep="\t", na="", row.names=F)

# DE genes in both comparisons
gplots::venn(list(PD = rownames(subset(df1, pvalue<=0.05)), MUSS = rownames(subset(df2, pvalue<=0.05))))

# barplot
df = filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% mutate(log2FoldChangeX = log2FoldChange.x>0, log2FoldChangeY = log2FoldChange.y>0) %>%
  arrange(log2FoldChangeX, log2FoldChangeY, log2FoldChange.x) %>% 
  mutate(genename = rowname) %>%
  select(genename, PD_vs_HC=log2FoldChange.x, MUSS=log2FoldChange.y) 
df$genename= factor(df$genename, levels = as.character(df$genename))
p= ggplot(gather(df, key="comparison", value="log2FoldChange", 2:3), aes(x=genename, y=log2FoldChange, fill=comparison)) +
  geom_bar(stat="identity", position=position_dodge()) + coord_flip()  + 
  # geom_text(aes(label=paste0(observed," (",round(OR,1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5) +  # TODO: add significance *
  theme_minimal() + theme(legend.justification=c(1,0), legend.position=c(1,0))
print(p)

dev.off()

### =======================
## DE2gene_SNDA: SNDA vs. CSF ( no overlap)
### =======================
df1=read.table("DE2gene_SNDA/DEresult.DE2gene_SNDA.CONDITION2_PD_vs_HC.xls", stringsAsFactors = F, row.names = 1, header=T, sep = "\t")
df2=read.table("DE2gene_CSF/DEresult.DE2gene_CSF.CONDITION_PD_vs_HC.xls", stringsAsFactors = F, row.names = 1, header=T, sep = "\t")
df1df2 = inner_join(rownames_to_column(df1) %>% filter(pvalue<=1) %>% select(1:8), 
                    rownames_to_column(df2) %>% filter(pvalue<=1) %>% select(1:13), 
                    by = 'rowname') # %>% mutate(log2FoldChange.y=-1*log2FoldChange.y) # sign of fold change for continous variable?
head(df1df2)

pdf(file = "PD.SNDACSF.DE2gene.barplot.pdf", paper = 'US')

# Intersection of genes qualified for DE analysis
gplots::venn(list(SNDA = rownames(df1), CSF = rownames(df2)))

# correlation between the them
cortest = cor.test(df1df2$log2FoldChange.x, df1df2$log2FoldChange.y,alternative = "two.sided", method = "pearson")
plot(df1df2$log2FoldChange.x, df1df2$log2FoldChange.y, 
     xlab="log2FoldChange (PD vs. HC in SNDA)", ylab="log2FoldChange (PD vs. HC in CSF)", 
     main="Correlation coeffecient analysis (PD case in SNDA vs. CSF)")
mtext(paste("Pearson's r =", round(cortest$estimate, 2), "(p-value:", format.pval(cortest$p.value),")"))

# significant in either comparison
filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% with(table(log2FoldChange.x>0, log2FoldChange.y>0))
# save to excel
filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% mutate(Chr=genes_annotation$chr[match(geneID, genes_annotation$geneID)]) %>% 
  select(CircRNA=rowname, Chr, log2FoldChange.SNDA = log2FoldChange.x, pvalue.SNDA = pvalue.x, log2FoldChange.CSF = log2FoldChange.y, pvalue.CSF = pvalue.y) %>%
  write.table(file="PD.SNDACSF.DE2gene.xls", sep="\t", na="", row.names=F)

# DE genes in both comparisons
gplots::venn(list(SNDA = rownames(subset(df1, pvalue<=0.05)), CSF = rownames(subset(df2, pvalue<=0.05))))

# barplot
df = filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% mutate(log2FoldChangeX = log2FoldChange.x>0, log2FoldChangeY = log2FoldChange.y>0) %>%
  arrange(log2FoldChangeX, log2FoldChangeY, log2FoldChange.x) %>% 
  mutate(genename = rowname) %>%
  select(genename, SNDA=log2FoldChange.x, CSF=log2FoldChange.y) 
df$genename= factor(df$genename, levels = as.character(df$genename))
p= ggplot(gather(df, key="comparison", value="log2FoldChange", 2:3), aes(x=genename, y=log2FoldChange, fill=comparison)) +
  geom_bar(stat="identity", position=position_dodge()) + coord_flip()  + 
  # geom_text(aes(label=paste0(observed," (",round(OR,1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5) +  # TODO: add significance *
  theme_minimal() + theme(legend.justification=c(1,0), legend.position=c(1,0))
print(p)

dev.off()

### =======================
## DE2gene_SNDA circRNA vs. mRNAs
### =======================
df1=read.table("~/projects/circRNA/results/DE2gene_SNDA/DEresult.DE2gene_SNDA.CONDITION2_PD_vs_HC.xls.gz", stringsAsFactors = F, row.names = 1, header=T, sep = "\t")
df2=read.table("~/neurogen/rnaseq_PD/results/DE_SNDA.gene/DEresult.DE_SNDA.gene.CONDITION2_PD_vs_HC.xls.gz", stringsAsFactors = F, row.names = 1, header=T, sep = "\t")
head(df1); head(df2)
df1df2 = inner_join(rownames_to_column(df1) %>% filter(pvalue<=1) %>% select(1:10) %>% mutate(geneID=gsub("\\..*","",geneID)), 
                    df2 %>% filter(pvalue<=1) %>% select(1:16), 
                    by = 'geneID') # %>% mutate(log2FoldChange.y=-1*log2FoldChange.y) # sign of fold change for continous variable?
head(df1df2)

pdf(file = "PD.SNDA.circRNAmRNA.DE2gene.barplot.pdf", paper = 'US')

# Intersection of genes qualified for DE analysis
gplots::venn(list(circRNA = df1$geneName, gene = df2$geneName))

# correlation between their foldchange
cortest = cor.test(df1df2$log2FoldChange.x, df1df2$log2FoldChange.y,alternative = "two.sided", method = "pearson")
plot(df1df2$log2FoldChange.x, df1df2$log2FoldChange.y, 
     xlab="log2FoldChange (PD vs. HC for circRNAs)", ylab="log2FoldChange (PD vs. HC for host gene)", 
     main="Correlation coeffecient analysis (PD-associated circRNAs vs. host genes)")
mtext(paste("Pearson's r =", round(cortest$estimate, 2), "(p-value:", format.pval(cortest$p.value),")"))

# correlation between their expression level
cortest = cor.test(df1df2$baseMean.x, df1df2$baseMean.y,alternative = "two.sided", method = "spearman")
plot(log10(0.1+df1df2$baseMean.x), log10(df1df2$baseMean.y+0.1), 
     xlab="baseMean (PD vs. HC for circRNAs)", ylab="baseMean (PD vs. HC for host gene)", 
     main="Correlation coeffecient analysis (PD-associated circRNAs vs. host genes)")
mtext(paste("Pearson's r =", round(cortest$estimate, 2), "(p-value:", format.pval(cortest$p.value),")"))

# significant in either comparison
filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% with(table(log2FoldChange.x>0, log2FoldChange.y>0))
# save to excel
filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% mutate(Chr=chr) %>% 
  select(CircRNA=rowname, Chr, log2FoldChange.circRNA = log2FoldChange.x, pvalue.circRNA = pvalue.x, log2FoldChange.hostgene = log2FoldChange.y, pvalue.hostgene = pvalue.y) %>%
  write.table(file="PD.SNDA.circRNAmRNA.DE2gene.xls", sep="\t", na="", row.names=F)

# DE genes in both comparisons
gplots::venn(list(circRNA = subset(df1, pvalue<=0.05, geneName, drop = T), gene = subset(df2, pvalue<=0.05, geneName, drop = T)))

# barplot
df = filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% mutate(log2FoldChangeX = log2FoldChange.x>0, log2FoldChangeY = log2FoldChange.y>0) %>%
  arrange(log2FoldChangeX, log2FoldChangeY, log2FoldChange.x) %>% 
  select(geneName, circRNA=log2FoldChange.x, genes=log2FoldChange.y) 
df$geneName= factor(df$geneName, levels = as.character(df$geneName))
p= ggplot(gather(df, key="comparison", value="log2FoldChange", 2:3), aes(x=geneName, y=log2FoldChange, fill=comparison)) +
  geom_bar(stat="identity", position=position_dodge()) + coord_flip()  + 
  # geom_text(aes(label=paste0(observed," (",round(OR,1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5) +  # TODO: add significance *
  theme_minimal() + theme(legend.justification=c(1,0), legend.position=c(1,0))
print(p)

dev.off()

### =======================
## DE2gene_TCPY circRNA vs. mRNAs
### =======================
df1=read.table("~/projects/circRNA/results/DE2gene_TCPY/DEresult.DE2gene_TCPY.CONDITION_AD_vs_HC.xls.gz", stringsAsFactors = F, row.names = 1, header=T, sep = "\t")
df2=read.table("~/neurogen/rnaseq_PD/results/DE_TCPY.gene/DEresult.DE_TCPY.gene.CONDITION_AD_vs_HC.xls.gz", stringsAsFactors = F, row.names = 1, header=T, sep = "\t")
head(df1); head(df2)
df1df2 = inner_join(rownames_to_column(df1) %>% filter(pvalue<=1) %>% select(1:10) %>% mutate(geneID=gsub("\\..*","",geneID)), 
                    df2 %>% filter(pvalue<=1) %>% select(1:16), 
                    by = 'geneID') # %>% mutate(log2FoldChange.y=-1*log2FoldChange.y) # sign of fold change for continous variable?
head(df1df2)

pdf(file = "AD.TCPY.circRNAmRNA.DE2gene.barplot.pdf", paper = 'US')

# Intersection of genes qualified for DE analysis
gplots::venn(list(circRNA = df1$geneName, gene = df2$geneName))

# correlation between their foldchange
cortest = cor.test(df1df2$log2FoldChange.x, df1df2$log2FoldChange.y,alternative = "two.sided", method = "pearson")
plot(df1df2$log2FoldChange.x, df1df2$log2FoldChange.y, 
     xlab="log2FoldChange (AD vs. HC for circRNAs)", ylab="log2FoldChange (AD vs. HC for host gene)", 
     main="Correlation coeffecient analysis (AD-associated circRNAs vs. host genes)")
mtext(paste("Pearson's r =", round(cortest$estimate, 2), "(p-value:", format.pval(cortest$p.value),")"))

# correlation between their expression level
cortest = cor.test(df1df2$baseMean.x, df1df2$baseMean.y,alternative = "two.sided", method = "spearman")
plot(log10(0.1+df1df2$baseMean.x), log10(df1df2$baseMean.y+0.1), 
     xlab="baseMean (AD vs. HC for circRNAs)", ylab="baseMean (AD vs. HC for host gene)", 
     main="Correlation coeffecient analysis (AD-associated circRNAs vs. host genes)")
mtext(paste("Pearson's r =", round(cortest$estimate, 2), "(p-value:", format.pval(cortest$p.value),")"))

# significant in either comparison
filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% with(table(log2FoldChange.x>0, log2FoldChange.y>0))
# save to excel
filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% mutate(Chr=chr) %>% 
  select(CircRNA=rowname, Chr, log2FoldChange.circRNA = log2FoldChange.x, pvalue.circRNA = pvalue.x, log2FoldChange.hostgene = log2FoldChange.y, pvalue.hostgene = pvalue.y) %>%
  write.table(file="AD.TCPY.circRNAmRNA.DE2gene.xls", sep="\t", na="", row.names=F)

# DE genes in both comparisons
gplots::venn(list(circRNA = subset(df1, pvalue<=0.05, geneName, drop = T), gene = subset(df2, pvalue<=0.05, geneName, drop = T)))

# barplot
df = filter(df1df2, pvalue.x<=0.05 | pvalue.y<=0.05) %>% mutate(log2FoldChangeX = log2FoldChange.x>0, log2FoldChangeY = log2FoldChange.y>0) %>%
  arrange(log2FoldChangeX, log2FoldChangeY, log2FoldChange.x) %>% 
  select(geneName, circRNA=log2FoldChange.x, genes=log2FoldChange.y) 
df$geneName= factor(df$geneName, levels = as.character(df$geneName))
p= ggplot(gather(df, key="comparison", value="log2FoldChange", 2:3), aes(x=geneName, y=log2FoldChange, fill=comparison)) +
  geom_bar(stat="identity", position=position_dodge()) + coord_flip()  + 
  # geom_text(aes(label=paste0(observed," (",round(OR,1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5) +  # TODO: add significance *
  theme_minimal() + theme(legend.justification=c(1,0), legend.position=c(1,0))
print(p)

dev.off()
