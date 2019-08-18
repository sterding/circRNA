## ===============================================
## get synaptic genes and test its enrichment
## ===============================================
library('dplyr')
library('reshape2')
setwd("~/projects/circRNA/data/") 

################################
## get synaptic genes
################################
## (1) 1461 genes in Synapse Proteomics Datasets in G2C (http://synsysdb.genes2cognition.org/db/GeneList/L00000069)
G2C_PSP=read.delim("SYNSYSDBdb_L00000069_BAYES-COLLINS-HUMAN-PSD-FULL.txt", header = T, stringsAsFactors = F)$Symbol %>% unique()
length(G2C_PSP)

## (2) genes with GO term including key words "synapse" or "synaptic";
library(biomaRt) #source("https://bioconductor.org/biocLite.R"); biocLite("biomaRt")
library(dplyr)
library(RCurl)
# ensembl_hs=useEnsembl(biomart="ensembl",GRCh=37,dataset='hsapiens_gene_ensembl')
# #listAttributes(ensembl_hs)
# listAttributes(ensembl_hs, page='feature_page')
# # get GO for all genes (if having GO annotation): N = 20475
# genesGO_all=getBM(attributes=c('ensembl_gene_id','hgnc_symbol', 'go_id','name_1006','namespace_1003'),
#                   filters =c('biotype','transcript_biotype','with_go'), 
#                   values =list("protein_coding","protein_coding",T),
#                   mart = ensembl_hs) %>% 
#   rename(go_name=name_1006, go_domain=namespace_1003) %>% distinct
# write.table(genesGO_all, file='genesGO_all.tab',sep = "\t", row.names = F)
genesGO_all=read.table(file='genesGO_all.tab',sep = "\t", header = T, stringsAsFactors = F)  # run above code ahead
genes_with_synapse_GO = genesGO_all %>% subset(grepl("synaptic|synapse", go_name), select=hgnc_symbol, drop=T) %>% unique()
length(genes_with_synapse_GO)

## (3) synaptomeDB
library(xlsx)
synaptomeDB=read.xlsx("SynaptomeDB.20190514.xlsx", sheetIndex = 1, header = T, stringsAsFactors = F)
head(synaptomeDB)

library(SuperExactTest); # install.packages('SuperExactTest')
plot(supertest(list(synaptomeDB=unique(synaptomeDB$symbol), GO=genes_with_synapse_GO, G2C=G2C_PSP), n=2000),degree=2:3, sort.by='size')

################################
## add SV annotation
################################
annotation_filtered_enriched=readRDS(file="Merge_circexplorer_BC106.filtered.enriched.annotation.bed14.rds")
dim(annotation_filtered_enriched)

# join "synapse" column
table(ifelse(annotation_filtered_enriched$geneName %in% c(G2C_PSP, genes_with_synapse_GO, synaptomeDB$symbol), "yes", "no"))
annotation_filtered_enriched$synapse=ifelse(annotation_filtered_enriched$geneName %in% c(G2C_PSP, genes_with_synapse_GO, synaptomeDB$symbol), "yes", "no")

total_circRNA_raw_reads_filtered = rowSums(Merge_circexp_raw[rowSums(Merge_circexp_raw)>=10, ])
total_circRNA_raw_reads_filtered_annotation=cbind(annotation[match(names(total_circRNA_raw_reads_filtered), annotation$ID),c('ID','geneName','synapse')], total_circRNA_raw_reads_filtered)

# all genes expressed in the 20 samples (RPKM>0.1)
fpkm=read.table("~/neurogen/rnaseq_PD/results/merged/genes.fpkm.cuffnorm.allSamples.uniq.xls", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
expressed_genes = gsub("\\..*","",rownames(fpkm)[rowMeans(fpkm[,sample20])>1])
expressed_genes_name = unique(genesGO_all$hgnc_symbol[genesGO_all$ensembl_gene_id %in% expressed_genes])
expressed_genes_ID = unique(genesGO_all$ensembl_gene_id[genesGO_all$ensembl_gene_id %in% expressed_genes])
write(expressed_genes_ID, file="pilot20_expressed_hostgenes.txt")
# write gene name to file for GSEA analysis
circRNA_hostgene = unique(as.character(total_circRNA_raw_reads_filtered_annotation$geneName))
write(circRNA_hostgene[circRNA_hostgene %in% expressed_genes_name], file="pilot20_circRNA_hostgenes.txt")

mytable=table(total_circRNA_raw_reads_filtered_annotation$synapse)
lbls=paste(names(mytable), "\n", mytable, sep="")
pie(mytable, labels = lbls, main="Pie chart of circRNAs of host genes with synaptic function")

mytable=table(unique(total_circRNA_raw_reads_filtered_annotation[,2:3])$synapse)
mytable

allgenes = unique(genesGO_all$hgnc_symbol)
synapse=ifelse(allgenes %in% c(G2C_PSP, genes_with_synapse_GO), "yes", "no")
circRNA=ifelse(allgenes %in% circRNA_hostgene, "yes", "no")
fisher.test(table(synapse, circRNA))


hist(total_circRNA_raw_reads_filtered_annotation$total_circRNA_raw_reads_filtered[total_circRNA_raw_reads_filtered_annotation$total_circRNA_raw_reads_filtered<=200], col="Red")
hist(total_circRNA_raw_reads_filtered_annotation$total_circRNA_raw_reads_filtered[total_circRNA_raw_reads_filtered_annotation$synapse=='yes'], breaks=seq(1,15000,1000), col="Blue", add=TRUE)