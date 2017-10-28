# Rscript to make manhattan plot for QTL result

library(qqman, quietly =T) # install.packages("qqman")
library(tidyr, quietly =T)
library(dplyr, quietly =T)

args<-commandArgs(TRUE)

nominalPzip=args[1] #nominalPzip="QTL_BC/eQTL.nominal.txt.gz"
p_cutoff=ifelse(is.na(args[2]), 0.01, as.numeric(args[2])) # 0.01

message(paste(nominalPzip, p_cutoff))

df = read.table(gzfile(nominalPzip), header = F, col.names=c("Gene","SNP","DIS","P","BETA"))
range(df$P)

df2 = df %>% filter(P <= p_cutoff) %>%
  separate(SNP, c("CHROM","BP"), sep=":", remove=F, convert =T) %>% 
  mutate(CHR=as.numeric(ifelse(CHROM=='X',23,ifelse(CHROM=='Y',24,ifelse(CHROM=='MT',25,CHROM)))), P) %>% select(CHR, BP, P) %>% distinct()
range(df2$P)
pdf(file=paste0(nominalPzip, ".manhattan.pdf"), width = 7, height = 5, useDingbats=T)
manhattan(df2, suggestiveline = F, genomewideline =F, main = paste("Manhattan plot for", nominalPzip), cex = 0.5, cex.axis = 0.8)
dev.off()

range(df2$P)