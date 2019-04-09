####################################
### post-eQTL analysis
####################################
## Q1: Are the eQTL circRNAs more likely hosted in an eQTL gene?
## Q2: Are circRNA and their hostgene likely associated with the same set of SNPs?
## Q3: How many of eQTL circRNAs expressed in more than 20 samples?
library(tidyverse)

annotation_refilterd4eqtl=readRDS("~/projects/circRNA/data/Merge_circexplorer_BC106.refiltered4eqtl.annotation.bed14.rds")  # n=1054
table(annotation_refilterd4eqtl$circType)
# circRNA   ciRNA 
# 752     302 
annotation_gene=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed", header = F, col.names = c('chr','start','end','EnsID','score','strand','symbol','type'), check.names = F, stringsAsFactors=F)

eqtl = read.table("~/projects/circRNA/data/QTL_BC/eQTL.nominal.txt.gz", header = F, col.names=c("gene","SNP","distance","P","beta"))
Psignifiance=as.numeric(try(system("grep threshold ~/projects/circRNA/data/QTL_BC/eQTL.post_fastQTL.log | cut -f2 -d':'", intern = T)))
eqtl = filter(eqtl, P<=Psignifiance) %>% separate(SNP, c("CHROM","BP"), sep=":", remove=F, convert =T)

eqtl_gene = read.table("~/projects/circRNA/data/QTL_BC/eQTLgene.nominal.txt.gz", header = F, col.names=c("gene","SNP","distance","P","beta"))
Psignifiance=as.numeric(try(system("grep threshold ~/projects/circRNA/data/QTL_BC/eQTLgene.post_fastQTL.log | cut -f2 -d':'", intern = T)))
eqtl_gene = filter(eqtl_gene, P<=Psignifiance) %>% separate(SNP, c("CHROM","BP"), sep=":", remove=F, convert =T)
eqtl_gene$symbol = annotation_gene$symbol[match(eqtl_gene$gene, annotation_gene$EnsID)]

saveRDS(eqtl, "~/projects/circRNA/data/QTL_BC/eQTL.nominal.txt.gz.rds")
saveRDS(eqtl_gene, "~/projects/circRNA/data/QTL_BC/eQTLgene.nominal.txt.gz.rds")

eqtl=readRDS("~/projects/circRNA/data/QTL_BC/eQTL.nominal.txt.gz.rds")
eqtl_gene=readRDS("~/projects/circRNA/data/QTL_BC/eQTLgene.nominal.txt.gz.rds")
annotation_refilterd4eqtl = mutate(annotation_refilterd4eqtl, circRNA_is_eGene=(ID %in% eqtl$gene), hostgene_is_eGene = (geneName %in% eqtl_gene$symbol))

## ================
## annotate the top eQTL pair
## ================
eqtl %>% filter(P<=1e-6, CHROM=='chr18') %>% left_join(y=select(annotation_refilterd4eqtl, ID, exonCount, circType, geneName, circRNA_is_eGene, hostgene_is_eGene), by=c("gene"="ID")) %>% head()

## ================
## Q1: Are the eQTL circRNAs more likely hosted in an eQTL gene? 
## A1: No
## ================
table(select(annotation_refilterd4eqtl, circRNA_is_eGene, hostgene_is_eGene))
#                   hostgene_is_eGene
# circRNA_is_eGene FALSE TRUE
#           FALSE   454    5
#           TRUE    590    5
fisher.test(table(select(annotation_refilterd4eqtl, circRNA_is_eGene, hostgene_is_eGene)))
# 
# Fisher's Exact Test for Count Data
# 
# data:  table(select(annotation_refilterd4eqtl, circRNA_is_eGene, hostgene_is_eGene))
# p-value = 0.7545
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 0.1760195 3.3660344
# sample estimates:
# odds ratio 
# 0.7696913 

## ================
## Q2: Are circRNA and their hostgene likely associated with the same set of SNPs?
## A2: No
## ================
totalSNP=scan("~/projects/circRNA/data/QTL_BC/genotype.snp.pos.txt",character())
table(totalSNP %in% eqtl$SNP, totalSNP %in% eqtl_gene$SNP)
#         FALSE    TRUE
# FALSE 4296694    6724
# TRUE     5193       5
fisher.test(table(totalSNP %in% eqtl$SNP, totalSNP %in% eqtl_gene$SNP))
# Fisher's Exact Test for Count Data
# 
# data:  table(totalSNP %in% eqtl$SNP, totalSNP %in% eqtl_gene$SNP)
# p-value = 0.376
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.199624 1.437445
# sample estimates:
# odds ratio 
#  0.6152597 

# test for all sSNPs
eSNP=unique(c(as.character(eqtl$SNP), as.character(eqtl_gene$SNP))); length(eSNP)
table(eSNP %in% eqtl$SNP, eSNP %in% eqtl_gene$SNP)
#         FALSE TRUE
# FALSE     0 6724
# TRUE   5193    5
fisher.test(table(eSNP %in% eqtl$SNP, eSNP %in% eqtl_gene$SNP), alternative = 'great')
# Fisher's Exact Test for Count Data
# 
# data:  table(eSNP %in% eqtl$SNP, eSNP %in% eqtl_gene$SNP)
# p-value = 1
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
# 0 Inf
# sample estimates:
# odds ratio 
# 0 

# the 5 SNPs drive both circRNA and gene
inner_join(x=eqtl, y=eqtl_gene, by="SNP") %>% left_join(y=select(annotation_refilterd4eqtl,ID,chrom,start,end,exonCount,circType, geneName), by=c("gene.x"="ID"))  # all from ciRNA

# the top eQTL SNP
eqtl %>% filter(P<=1e-6, CHROM=="chr15") %>% inner_join(y=eqtl_gene, by="SNP") %>% left_join(y=select(annotation_refilterd4eqtl,ID,chrom,start,end,exonCount,circType, geneName), by=c("gene.x"="ID"))

## ================
## Q3: How many of eQTL circRNAs expressed in more than 20 samples?
## ================
Merge_circexp_norm= readRDS("~/projects/circRNA/data/Merge_circexplorer_BC106.normRPM.rds")
df=Merge_circexp_norm[as.character(annotation_refilterd4eqtl$ID),scan("~/neurogen/rnaseq_PD/results/merged/samplelist.HCILB_SNDA84",character())]; dim(df)
RPM_threshold = 0.001
eqtl %>% filter(gene %in% rownames(df)[rowSums(df>=RPM_threshold)>=45]) %>% arrange(P) %>% left_join(y=annotation_refilterd4eqtl,by = c("gene"="ID"))

# how many samples that eQTL circRNAs are expressed in?
d=data.frame(n_samples=rowSums(df[unique(as.character(eqtl$gene)),] >= RPM_threshold), type=annotation_refilterd4eqtl$circType[match(unique(as.character(eqtl$gene)), annotation_refilterd4eqtl$ID)])
library(ggplot2); 
require(scales)
mylog_trans <- function (base = exp(1), from = 0) 
{
  trans <- function(x) log(x, base) - from
  inv <- function(x) base^(x + from)
  trans_new("mylog", trans, inv, log_breaks(base = base), domain = c(base^from, Inf))
}
pdf("~/projects/circRNA/data/QTL_BC/eQTL.circularRNA.N84.expression.hist.pdf", width = 5, height = 3)
cols <- c("#ff0000", "#fdbb84"); names(cols) = paste0(names(table(d$type)), " (n=",table(d$type),")")
ggplot(d, aes(x=n_samples)) + geom_histogram(aes(fill="ciRNA (n=182)"), binwidth = 1,color='#ffffff') +  geom_histogram(data=subset(d,type=='circRNA'),binwidth = 1,color='#ffffff',aes(fill="circRNA (n=413)")) + scale_x_continuous(breaks=seq(0,85,5))+ scale_y_continuous(trans = mylog_trans(base=10, from=-1), breaks=c(1,10,100,1000),limits=c(0.1,500)) + geom_text(data=subset(d, n_samples > 40), aes(x=n_samples,y=1,label=rownames(subset(d, n_samples > 40))), angle=90, nudge_y=0.05, hjust=0) + theme_bw() + labs(x="Numbers of samples (out of 84 in total)", title='eQTL circular RNAs', y='Numbers of circular RNAs') + scale_fill_manual(name="circular RNA types",values=cols) + theme(legend.position = c(0.2, 0.9))
dev.off()

# vocano plot for significant eQTL
d= d %>% rownames_to_column('ID') %>% inner_join(y=eqtl, by=c("ID"="gene")) %>% mutate(geneName=annotation_refilterd4eqtl$geneName[match(ID, annotation_refilterd4eqtl$ID)])
p = ggplot(d) + geom_point(aes(x=beta, y=-log10(P), col=type, size=n_samples, text=paste0(ID, " (",geneName,")")), alpha=0.5) + theme_bw() + labs(x="Effect size (beta)", title='Significant cis-eQTL for circular RNAs', y='-log10(P value)') + geom_text(data=filter(d, beta >=1 | beta<= -1 | P<=1e-7), aes(x=beta,y=-log10(P),label=geneName), size=2,angle=0, nudge_x=0.01, hjust=0, check_overlap=T) 
ggsave("~/projects/circRNA/data/QTL_BC/eQTL.circularRNA.N84.volcano.pdf",p, width=8, height=7)

install.packages('plotly')
library(plotly)
p <- ggplotly(p, tooltip = c("text", "x",'y','col','size'))
p

# Create a shareable link to your chart
# Set up API credentials: https://plot.ly/r/getting-started
chart_link = plotly_POST(p, filename="geom_point/scatter")
chart_link

## FUTURE: combine with Shiny to hover boxplot image, see example below:
# https://plot.ly/r/shiny-coupled-hover-events/




