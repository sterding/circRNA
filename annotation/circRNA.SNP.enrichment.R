## Rscript to test the enrichment of GWAS SNPs in circRNAs
#
# The Fisher exact test was performed to test the odds ratio of trait association with enhancers vs non-enhancers taking into account the genomic distribution of SNPs.  
# OR = (A/B) / (C/D), where
# 
# A: number of trait-associated SNPs in enhancers
# B: number of non-trait-associated SNPs in enhancers (from dbSNP - GWAS SNPs for that trait)
# C: number of trait-associated SNPs distal to enhancers
# D: number of non-trait-associated SNPs distal to enhancers (from dbSNP - GWAS SNPs for that trait)
library(tidyverse)

type="SNAP"; setwd("~/projects/circRNA/data/"); output_prefix=paste0("circRNA.SNP.enrichments.",type)

s=read.table(paste0("SNP.",type,".counts.summary"), header=F,row.names=1); 
results=data.frame();
for(i in c("circRNA","circRNA-SNDA","circRNA-PY","circRNA-NN","circRNA-SNDA-private","circRNA-PY-private","circRNA-NN-private","circRNA-Neuron","circRNA-matched-exons","circRNA-matched-exons-same-host","circRNA-flanking-intron5","circRNA-flanking-intron3","circRNA-internal-introns","circRNA-matched-introns")){
  if(!file.exists(paste0("SNP.",type,".count.",i))) next;
  n1=s[i,1]; n2=s['all',1];
  message(i);
  all=read.table(paste0("SNP.",type,".count.all")); rownames(all)=all[,1]
  x=read.table(paste0("SNP.",type,".count.",i)); rownames(x)=x[,1]
  df=cbind(x, all[rownames(x),2]); # only the traits occurred
  #df=merge(all, x, by="V1", all=T); df[is.na(df)]=0; rownames(df)=df[,1]; # all traits (0 for nonoccurance)  --> not working for fisher.test
  df=df[,-1]; colnames(df)=c('observed','all')
  results=rbind(results, cbind(Disease_or_Trait=rownames(df), 
                               df, 
                               pvalue=apply(df, 1, function(x) fisher.test(matrix(c(x[1],x[2]-x[1], n1-x[1], n2-n1-x[2]+x[1]), nrow = 2), alternative='greater')$p.value), 
                               OR=apply(df, 1, function(x) fisher.test(matrix(c(x[1],x[2]-x[1], n1-x[1], n2-n1-x[2]+x[1]), nrow = 2), alternative='greater')$estimate), 
                               type=i))
  results$Disease_or_Trait = as.character(results$Disease_or_Trait)
}

results$Disease_or_Trait=gsub("_"," ", results$Disease_or_Trait)
results = results[with(results, order(type, -OR)), ]

# save all result
write.table(results, paste0(output_prefix,".full.xls"), sep="\t", col.names = T, row.names = F)

# convert to wide table
#results = read.table(paste0(output_prefix,".full.xls"), header = T)
results = filter(results, type %in% c("circRNA", "circRNA-SNDA", "circRNA-PY", "circRNA-NN", "circRNA-SNDA-private", "circRNA-PY-private", "circRNA-NN-private")) %>% 
  pivot_wider(names_from = type,names_glue = "{type}.{.value}",values_from = c(observed, pvalue, OR)) 
#results[, c("Disease_or_Trait", "all", sort(colnames(results)[grep("circRNA", colnames(results))]))]
write.table(results, paste0(output_prefix,".full.wide.xls"), sep="\t", col.names = T, row.names = F)


results=subset(results, OR>1 & pvalue<0.01/1037 & observed>3)
table(results$type)
# circRNA                    circRNA-SNDA                      circRNA-PY                      circRNA-NN            circRNA-SNDA-private              circRNA-PY-private 
# 93                              71                              61                              83                              49                              41 
# circRNA-NN-private                  circRNA-Neuron           circRNA-matched-exons circRNA-matched-exons-same-host        circRNA-flanking-intron5        circRNA-flanking-intron3 
# 64                              60                              37                              46                              93                              77 
# circRNA-internal-introns         circRNA-matched-introns 
# 92                             165 

#### main plot


pdf(paste0(output_prefix,".pdf"), width=20, height=12); 
# Note: Don't use ggsave() with Rscript, which will generate another Rplot.pdf unnecessarily. See http://stackoverflow.com/questions/19382384/ggplot2-overwrite-one-another-in-rplots-pdf
for(i in c("circRNA","circRNA-SNDA","circRNA-PY","circRNA-NN","circRNA-SNDA-private","circRNA-PY-private","circRNA-NN-private","circRNA-Neuron","circRNA-matched-exons","circRNA-matched-exons-same-host","circRNA-flanking-intron5","circRNA-flanking-intron3","circRNA-internal-introns","circRNA-matched-introns")){
  results = read.table(paste0(output_prefix,".full.xls"), sep="\t", header=T,stringsAsFactors = F)
  
  results=filter(results, OR>1, pvalue<0.01/1037, observed>3, type == i)
  results$pvalue[results$pvalue==0]=2.2e-16
  
  # re-order the levels in the order of appearance in the data.frame
  results = results[with(results, order(type, pvalue)), ]
  results$Disease_or_Trait2 <- factor(results$Disease_or_Trait, unique(as.character(results$Disease_or_Trait)))
  p = ggplot(results, aes(x=Disease_or_Trait2, y=-log10(pvalue), fill=type, ymax=max(-log10(pvalue))*1.1)) 
  p = p + geom_bar(width=.2, position = position_dodge(width=1), stat="identity") 
  p = p + geom_hline(yintercept=-log10(0.01/1037), size=.5,linetype = 2)  ## Bonferroni correction, where 1037 is the number of disease/traits in GWAS (wc -l SNP.$type.count.all)
  p = p + theme_bw() 
  p = p + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=12), legend.justification=c(1,1), legend.position=c(1,1)) 
  p = p + geom_text(aes(label=paste0(observed," (",round(OR,1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5) 
  p = p + ggtitle(paste0("GWAS SNP enrichments (LD from ",type,", sorted by pvalue)")) 
  
  print(p);
}

dev.off();


