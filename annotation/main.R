## Main script to work on the output of circExplorer
library('dplyr')
library('reshape2')
library('tidyverse')
library('ggpubr') # install.packages("ggpubr")
setwd("~/projects/circRNA/data/") 

## See ~/projects/circRNA/src/main.sh for how to generate the input files: 
# Merge_circexplorer_BC221.annotation.bed14
# Merge_circexplorer_BC221.rawcount.long.txt

## sample size explaination
# total sample size in BRAINcode2 = 221 (see Google spreadsheet and FigS1a)
# after QC, outlier removal etc, n=197 (see FigS1a and Fig1b)
# Among them, neuronal samples n=190, healthy control n=109, 

## see pilot.R in ~/Dropbox/grant/2019R21/pilot.R to make figure for pilot study (for grant purpose)

## Log
# 20191028: change to use only neuronal samples for Fig1 (104 SNDA + 86 PY = 190) 

###########################################
################# load data  ##############
###########################################

## === annotation ===

annotation<- read.table("Merge_circexplorer_BC221.annotation.bed14",header=T, check.names = F) # include both ciRNA and circRNA
head(annotation); dim(annotation)

# add geneSymbol, geneType
GENCODEv19=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v19.annotation.genes.bed", col.names = c("chr","start","end","geneName","geneID","geneType","strand"), stringsAsFactors = F, header = F)
annotation$geneName=GENCODEv19$geneName[match(annotation$geneID, GENCODEv19$geneID)]
annotation$geneType=GENCODEv19$geneType[match(annotation$geneID, GENCODEv19$geneID)]

# add CDR1as
annotation$geneName[annotation$geneID=='CDR1as'] = 'CDR1as'
annotation$geneType[annotation$geneID=='CDR1as'] = 'antisense'
  
## === load raw reads count ===

Merge_circexp_raw <- read.table("Merge_circexplorer_BC221.rawcount.long.txt",sep="\t",check.names =F, col.names = c("ID","readsCount","sampleID")) %>% 
  group_by(ID, sampleID) %>% summarise(readsCount=sum(readsCount)) %>%
  spread(key = sampleID, value = readsCount, fill = 0) %>%
  column_to_rownames("ID")
head(Merge_circexp_raw); dim(Merge_circexp_raw)

###########################################
############## normalize to RPM ###########
###########################################
BRAINCODE2_Sequencing_Log.RNAseq_statistics="https://docs.google.com/spreadsheets/d/e/2PACX-1vQFQ4aQj0sD9oxIqaZ-cEgo7kWcCmNYGBH9emLw8iNu0f6TTjKE5Lte7IBfoMMy57cLjA4pXE0YlPY2/pub?gid=1049148747&single=true&output=tsv"
require(RCurl)
RNAseq_statistics<- read.delim(textConnection(getURL(BRAINCODE2_Sequencing_Log.RNAseq_statistics)), header=T, check.names = T, comment.char = "#", stringsAsFactors = F) 
readsNum_million <- as.numeric(gsub(",","",RNAseq_statistics$total_mapped_reads))/10^6; names(readsNum_million) = RNAseq_statistics$SOURCE_SAMPLE_ID
readsNum_million = readsNum_million[colnames(Merge_circexp_raw)]

Merge_circexp_norm<-sweep(Merge_circexp_raw,2,readsNum_million,"/")

# save
saveRDS(Merge_circexp_raw, file="Merge_circexplorer_BC221.rawcount.rds")
saveRDS(Merge_circexp_norm, file="Merge_circexplorer_BC221.normRPM.rds")
saveRDS(annotation, file="Merge_circexplorer_BC221.annotation.bed14.rds")

write.table(Merge_circexp_raw, file="Merge_circexplorer_BC221.rawcount.xls", quote=F, sep="\t", col.names = NA, row.names = TRUE); system("gzip Merge_circexplorer_BC221.rawcount.xls")
write.table(Merge_circexp_norm, file="Merge_circexplorer_BC221.normRPM.xls", quote=F, sep="\t", col.names = NA, row.names = TRUE); system("gzip Merge_circexplorer_BC221.normRPM.xls")

## Technical replication (Supplementary Fig. 1e)
# Rscript ~/pipeline/modules/_pairwise_compare.R HC_UWA616_SNDA_2_rep1 HC_UWA616_SN_6_rep1.amplified Merge_circexplorer_BC221.normRPM.rds
# Rscript ~/pipeline/modules/_pairwise_compare.R HC_BN12-44_TCPY_11_rep1 HC_BN12-44_TCPY_5_rep1 Merge_circexplorer_BC221.normRPM.rds 
###########################################
############ filter cirRNAs     ###########
###########################################

### (0) only  197 selected samples (104 SNDA + 86 PY + 7 NN = 197) 
BC197=filter(RNAseq_statistics, BRAINcode2.final.selection==1, CELL!='CSF') %>% pull(SOURCE_SAMPLE_ID); length(BC197)
Merge_circexp_raw_filtered <- Merge_circexp_raw[,BC197]
Merge_circexp_raw_filtered <- Merge_circexp_raw_filtered[rowSums(Merge_circexp_raw_filtered)>0,]; dim(Merge_circexp_raw_filtered)
saveRDS(Merge_circexp_raw_filtered, file="Merge_circexplorer_BC197.rawcount.rds")
saveRDS(annotation[annotation$ID %in% rownames(Merge_circexp_raw_filtered),], file="Merge_circexplorer_BC197.annotation.bed14.rds")

### (1) being expressed: with >=2 unique back-splicing reads in all samples (Zheng et al., doi:10.1038/ncomms11215)
### --------------------

Merge_circexp_raw_filtered <- Merge_circexp_raw_filtered[rowSums(Merge_circexp_raw_filtered)>=2,]
Merge_circexp_norm_filtered <- Merge_circexp_norm[rownames(Merge_circexp_raw_filtered), BC197]
annotation_filtered <- annotation[annotation$ID %in% rownames(Merge_circexp_norm_filtered),] 

dim(Merge_circexp_raw_filtered); dim(Merge_circexp_norm_filtered); dim(annotation_filtered)
# [1] 124771   197

# save
saveRDS(Merge_circexp_raw_filtered, file="Merge_circexplorer_BC197.filtered.rawcount.rds")
saveRDS(Merge_circexp_norm_filtered, file="Merge_circexplorer_BC197.filtered.normRPM.rds")
saveRDS(annotation_filtered, file="Merge_circexplorer_BC197.filtered.annotation.bed14.rds")
#write.table(annotation_filtered, file="Merge_circexplorer_BC197.filtered.annotation.bed14", sep = "\t", quote = F, row.names = F, col.names = F)

###########################################
############ being enriched     ###########
###########################################

# Definition of being enriched: at least 20 reads in RNase R treatment AND foldchange(RNase vs Mock) at least 2 in at least one sample. 
# Run main.RM.R first
Merge_circexp_raw_long_enriched.RM  = readRDS(file="Merge_circexp_raw_long_filterd.RM.rds")
length(unique(Merge_circexp_raw_long_enriched.RM$ID))
# [1] 50900
## circRNAs enriched in brain
length(unique(filter(Merge_circexp_raw_long_enriched.RM, sampleGroup3!="NN") %>% pull(ID)))
# [1] 39875  ## goes to Fig.1a
table(unique(Merge_circexp_raw_long_enriched.RM$ID) %in% rownames(Merge_circexp_raw_filtered))
# FALSE  TRUE 
# 37499  13401

Merge_circexp_norm_filtered_and_enriched=Merge_circexp_norm_filtered[intersect(rownames(Merge_circexp_raw_filtered), Merge_circexp_raw_long_enriched.RM$ID), ]
Merge_circexp_raw_filtered_and_enriched <- Merge_circexp_raw_filtered[rownames(Merge_circexp_norm_filtered_and_enriched), ]
annotation_filtered_enriched <- annotation[annotation$ID %in% rownames(Merge_circexp_raw_filtered_and_enriched),] 

dim(Merge_circexp_norm_filtered_and_enriched); dim(Merge_circexp_raw_filtered_and_enriched); dim(annotation_filtered_enriched)
# [1] 13401   197
# [1] 10431   125
# [1] 10017   106

# save
saveRDS(Merge_circexp_raw_filtered_and_enriched, file="Merge_circexplorer_BC197.filtered.enriched.rawcount.rds")
saveRDS(Merge_circexp_norm_filtered_and_enriched, file="Merge_circexplorer_BC197.filtered.enriched.normRPM.rds")
saveRDS(annotation_filtered_enriched, file="Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds")
annotation_filtered_enriched=readRDS(file="Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds")
write.table(annotation_filtered_enriched, file="Merge_circexplorer_BC197.filtered.enriched.annotation.bed14", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(Merge_circexp_norm_filtered_and_enriched, file="Merge_circexplorer_BC197.filtered.enriched.normRPM.tab", sep = "\t", quote = F, row.names = T, col.names = NA) # a blank column name is added for the topleft corner

# save annotation to different cell type
Merge_circexp_raw_filtered_and_enriched=readRDS("Merge_circexplorer_BC197.filtered.enriched.rawcount.rds")
annotation_filtered_enriched=readRDS(file="Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds")
Merge_circexp_raw_filtered_and_enriched %>% select(contains("_SNDA_")) %>% transmute(ID=rownames(.), rowsums=rowSums(.)) %>% filter(rowsums > 0) %>% 
  left_join(y=annotation_filtered_enriched,by = "ID") %>% mutate(itemRgb=ifelse(circType=="circRNA",ifelse(strand=="+", "255,0,0", "0,0,255"), ifelse(strand=="+", "255,100,100", "100,100,255")), score = floor(1000*rowsums/max(rowsums))) %>% 
  select(chrom, start, end, ID, score, strand, thickStart,  thickEnd, itemRgb, exonCount, exonSizes, exonOffsets) %>% arrange(chrom, start, end) %>% write.table(file="ucsc_tracks/Merge_circexplorer_BC197.filtered.enriched.annotation.SNDA.bed12", row.names = F, col.names = F, quote = F, sep = "\t")
Merge_circexp_raw_filtered_and_enriched %>% select(matches('_TCPY_|_MCPY_')) %>% transmute(ID=rownames(.), rowsums=rowSums(.)) %>% filter(rowsums > 0) %>% 
  left_join(y=annotation_filtered_enriched,by = "ID") %>% mutate(itemRgb=ifelse(circType=="circRNA",ifelse(strand=="+", "255,0,0", "0,0,255"), ifelse(strand=="+", "255,100,100", "100,100,255")), score = floor(1000*rowsums/max(rowsums))) %>% 
  select(chrom, start, end, ID, score, strand, thickStart,  thickEnd, itemRgb, exonCount, exonSizes, exonOffsets)%>% arrange(chrom, start, end) %>% write.table(file="ucsc_tracks/Merge_circexplorer_BC197.filtered.enriched.annotation.PY.bed12", row.names = F, col.names = F, quote = F, sep = "\t")
Merge_circexp_raw_filtered_and_enriched %>% select(matches('_PBMC_|_FB_')) %>% transmute(ID=rownames(.), rowsums=rowSums(.)) %>% filter(rowsums > 0) %>% 
  left_join(y=annotation_filtered_enriched,by = "ID") %>% mutate(itemRgb=ifelse(circType=="circRNA",ifelse(strand=="+", "255,0,0", "0,0,255"), ifelse(strand=="+", "255,100,100", "100,100,255")), score = floor(1000*rowsums/max(rowsums))) %>% 
  select(chrom, start, end, ID, score, strand, thickStart,  thickEnd, itemRgb, exonCount, exonSizes, exonOffsets) %>% arrange(chrom, start, end) %>% write.table(file="ucsc_tracks/Merge_circexplorer_BC197.filtered.enriched.annotation.NN.bed12", row.names = F, col.names = F, quote = F, sep = "\t")

# only the 109 HC
Merge_circexp_raw_filtered_and_enriched %>% select(matches("HC_.*_SNDA_")) %>% transmute(ID=rownames(.), rowsums=rowSums(.)) %>% filter(rowsums > 0) %>% 
  left_join(y=annotation_filtered_enriched,by = "ID") %>% mutate(itemRgb=ifelse(circType=="circRNA",ifelse(strand=="+", "255,0,0", "0,0,255"), ifelse(strand=="+", "255,100,100", "100,100,255")), score = floor(1000*rowsums/max(rowsums))) %>% 
  select(chrom, start, end, ID, score, strand, thickStart,  thickEnd, itemRgb, exonCount, exonSizes, exonOffsets) %>% arrange(chrom, start, end) %>% write.table(file="ucsc_tracks/Merge_circexplorer_BC109.filtered.enriched.annotation.SNDA.bed12", row.names = F, col.names = F, quote = F, sep = "\t")
Merge_circexp_raw_filtered_and_enriched %>% select(matches('HC_.*_[T,M]CPY_')) %>% transmute(ID=rownames(.), rowsums=rowSums(.)) %>% filter(rowsums > 0) %>% 
  left_join(y=annotation_filtered_enriched,by = "ID") %>% mutate(itemRgb=ifelse(circType=="circRNA",ifelse(strand=="+", "255,0,0", "0,0,255"), ifelse(strand=="+", "255,100,100", "100,100,255")), score = floor(1000*rowsums/max(rowsums))) %>% 
  select(chrom, start, end, ID, score, strand, thickStart,  thickEnd, itemRgb, exonCount, exonSizes, exonOffsets)%>% arrange(chrom, start, end) %>% write.table(file="ucsc_tracks/Merge_circexplorer_BC109.filtered.enriched.annotation.PY.bed12", row.names = F, col.names = F, quote = F, sep = "\t")
Merge_circexp_raw_filtered_and_enriched %>% select(matches('HC_.*_PBMC_|HC_.*_FB_')) %>% transmute(ID=rownames(.), rowsums=rowSums(.)) %>% filter(rowsums > 0) %>% 
  left_join(y=annotation_filtered_enriched,by = "ID") %>% mutate(itemRgb=ifelse(circType=="circRNA",ifelse(strand=="+", "255,0,0", "0,0,255"), ifelse(strand=="+", "255,100,100", "100,100,255")), score = floor(1000*rowsums/max(rowsums))) %>% 
  select(chrom, start, end, ID, score, strand, thickStart,  thickEnd, itemRgb, exonCount, exonSizes, exonOffsets) %>% arrange(chrom, start, end) %>% write.table(file="ucsc_tracks/Merge_circexplorer_BC109.filtered.enriched.annotation.NN.bed12", row.names = F, col.names = F, quote = F, sep = "\t")

## then go the ucsc_tracks/README.txt

#########################################
## impute the dropout zero using scImpute 
#########################################
# install 
# Sys.unsetenv("GITHUB_PAT"); # if Bad credentials error 
# install_github("Vivianstats/scImpute")
library('scImpute') 
scimpute(count_path = "Merge_circexplorer_BC197.filtered.enriched.rawcount.rds", # full path to raw count matrix
  infile = "rds",           # format of input file
  outfile = "rds",          # format of output file
  out_dir = "./scImpute/",   # full path to output directory
  labeled = FALSE,          # cell type labels not available
  drop_thre = 0.4,          # threshold set on dropout probability
  Kcluster = 2,             # 2 cell subpopulations
  ncores = 4)              # number of cores used in parallel computation

#########################################
## make the sebset for 190 neuronal samples: for Fig.1 
#########################################
BC190=filter(RNAseq_statistics, BRAINcode2.final.selection==1, CELL %in% c('SNDA',"MCPY","TCPY")) %>% pull(SOURCE_SAMPLE_ID); length(BC190)
Merge_circexp_raw_BC190 = Merge_circexp_raw[rowSums(Merge_circexp_raw[,BC190])>0,BC190]; 
saveRDS(Merge_circexp_raw_BC190, file="Merge_circexplorer_BC190.rawcount.rds")
saveRDS(annotation[annotation$ID %in% rownames(Merge_circexp_raw_BC190),], file="Merge_circexplorer_BC190.annotation.bed14.rds")
dim(Merge_circexp_raw_BC190) #  365145 (initial set of circRNAs called in brain neuron)

Merge_circexp_raw_filtered_BC190 = Merge_circexp_raw_filtered[rowSums(Merge_circexp_raw_filtered[,BC190])>=2,BC190]; 
saveRDS(Merge_circexp_raw_filtered_BC190, file="Merge_circexplorer_BC190.filtered.rawcount.rds")
Merge_circexp_norm_filtered_BC190 = Merge_circexp_norm_filtered[rownames(Merge_circexp_raw_filtered_BC190),BC190]; dim(Merge_circexp_norm_filtered_BC190)
saveRDS(Merge_circexp_norm_filtered_BC190, file="Merge_circexplorer_BC190.filtered.normRPM.rds")
annotation_filtered_BC190=filter(annotation_filtered, ID %in% rownames(Merge_circexp_raw_filtered_BC190)); dim(annotation_filtered_BC190)
saveRDS(annotation_filtered_BC190, file="Merge_circexplorer_BC190.filtered.annotation.bed14.rds")
dim(Merge_circexp_raw_filtered_BC190) # 111,419 (filtered set of circRNAs called in brain neuron)

Merge_circexp_raw_filtered_and_enriched_BC190 = Merge_circexp_raw_filtered_and_enriched[rowSums(Merge_circexp_raw_filtered_and_enriched[,BC190])>0,BC190];
saveRDS(Merge_circexp_raw_filtered_and_enriched_BC190, file="Merge_circexplorer_BC190.filtered.enriched.rawcount.rds")
Merge_circexp_norm_filtered_and_enriched_BC190 = Merge_circexp_norm_filtered_and_enriched[rownames(Merge_circexp_raw_filtered_and_enriched_BC190),BC190]
saveRDS(Merge_circexp_norm_filtered_and_enriched_BC190, file="Merge_circexplorer_BC190.filtered.enriched.normRPM.rds")
annotation_filtered_enriched_BC190=filter(annotation_filtered_enriched, ID %in% rownames(Merge_circexp_raw_filtered_and_enriched_BC190))
saveRDS(annotation_filtered_enriched_BC190, file="Merge_circexplorer_BC190.filtered.enriched.annotation.bed14.rds")
dim(Merge_circexp_norm_filtered_and_enriched_BC190); dim(Merge_circexp_raw_filtered_and_enriched_BC190); dim(annotation_filtered_enriched_BC190)
# [1] 11039  (validated set of circRNAs called in brain neuron)

#########################################
## Figure 1b: distribution of circRNAs supported by different number of reads
###########################################

# pie chart of circRNA among all circular RNAs
pdf("~/projects/circRNA/results/total_circRNA_raw_reads.BC190.pie.pdf", width=6, height = 2)
par(mfrow=c(1,3))
pie(table(readRDS("Merge_circexplorer_BC190.annotation.bed14.rds")$circType=='circRNA'),labels = NA, border='white', col=c('#fdbb84','red'),main="all")
pie(table(readRDS("Merge_circexplorer_BC190.filtered.annotation.bed14.rds")$circType=='circRNA'),labels = NA, border='white', col=c('#fdbb84','red'),main="filtered")
pie(table(readRDS("Merge_circexplorer_BC190.filtered.enriched.annotation.bed14.rds")$circType=='circRNA'),labels = NA, border='white', col=c('#fdbb84','red'),main="filtered+enriched")
dev.off()

table(readRDS("Merge_circexplorer_BC190.filtered.enriched.annotation.bed14.rds")[,c('geneType','circType')])
# geneType                   circRNA ciRNA
# 3prime_overlapping_ncrna       2     0
# antisense                     33     2
# lincRNA                       33     1
# polymorphic_pseudogene         1     0
# processed_transcript          18     3
# protein_coding             10708    188
# pseudogene                    48     0
# sense_intronic                 1     0
# sense_overlapping              1     0

Merge_circexp_raw=readRDS("Merge_circexplorer_BC190.rawcount.rds"); dim(Merge_circexp_raw)
Merge_circexp_raw_filtered_and_enriched=readRDS("Merge_circexplorer_BC190.filtered.enriched.rawcount.rds"); dim(Merge_circexp_raw_filtered_and_enriched)
annotation_filtered_enriched=readRDS(file="Merge_circexplorer_BC190.filtered.enriched.annotation.bed14.rds"); dim(annotation_filtered_enriched)

## -------------------
### overlap with PD-GWAS genes, AD-GWAS genes, and SynGO genes
## -------------------
ADgenes=unique(scan("~/neurogen/external_download/externalData/GWAS/AD/AD_gwas_associatedGene.txt", character())); length(ADgenes)
#PDgenes=unique(system("grep -v Candidate ~/neurogen/external_download/externalData/GWAS/PD.GWAS.Chang2017.table1n2.xls | cut -f3 | sed 's/,/\\\n/g' | sort -u", intern = TRUE)); length(PDgenes)  # Chang et al. 2017
PDgenes=unique(system("grep -v Nearest ~/neurogen/external_download/externalData/GWAS/PD/PD.GWAS.Nalls2019.TableS2.txt | cut -f4,5 | sed 's/\\\t$//;s/\\\t/\\\n/g' | sort -u", intern = TRUE)); length(PDgenes) # Nalls et al. 2019
synGO=unique(readxl::read_xlsx("syngo_genes.xlsx")$ensembl_id); length(synGO)  # N=1112; download from SynGO (https://www.syngoportal.org)
endocytosisGO=unique(read.delim("endocytosis.GO_0006897.Ensembl104.bioMart.tsv")$Gene.stable.ID); length(endocytosisGO)  # N=802; based on bioMart http://www.ensembl.org/biomart/martview/23b883ab8c5172aec6d515d3194108c6 
#endocytosisGO=unique(read.delim("endocytosis.GO_0006897.Ensembl75.bioMart.tsv")$Ensembl.Gene.ID); length(endocytosisGO)  # N=400; based on bioMart http://feb2014.archive.ensembl.org/biomart/martview/de88da23301b459dc48f92592696c222
## read the latest KEGG pathway (see note below)
#endocytosisKEGG=read.delim("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/KEGG/path2genes.tab") %>% filter(pathName=="Endocytosis") %>% pull(geneName) %>% unique(); length(endocytosisKEGG) # N=1082

annotation=readRDS(file="Merge_circexplorer_BC190.annotation.bed14.rds")
Merge_circexp_raw=readRDS("Merge_circexplorer_BC190.rawcount.rds")

# top PD genes
Merge_circexp_raw[as.character(annotation %>% filter(geneName %in% PDgenes) %>% pull(ID)),] %>% rownames_to_column() %>% mutate(rowsum=rowSums(select(.,contains("_")))) %>% select(rowsum, ID=rowname) %>% left_join(y=annotation, by="ID") %>% slice_max(rowsum, n = 20)
# top AD genes
Merge_circexp_raw[as.character(annotation %>% filter(geneName %in% ADgenes) %>% pull(ID)),] %>% rownames_to_column() %>% mutate(rowsum=rowSums(select(.,contains("_")))) %>% select(rowsum, ID=rowname) %>% left_join(y=annotation, by="ID") %>% slice_max(rowsum, n = 20) 
# top synGO genes
Merge_circexp_raw[as.character(annotation %>% mutate(geneID=gsub("\\..*","", geneID)) %>% filter(geneID %in% synGO) %>% pull(ID)),] %>% rownames_to_column() %>% mutate(rowsum=rowSums(select(.,contains("_")))) %>% select(rowsum, ID=rowname) %>% left_join(y=annotation, by="ID") %>% slice_max(rowsum, n = 20) 
# top EndocytosisGO genes
Merge_circexp_raw[as.character(annotation %>% mutate(geneID=gsub("\\..*","", geneID)) %>% filter(geneID %in% endocytosisGO) %>% pull(ID)),] %>% rownames_to_column() %>% mutate(rowsum=rowSums(select(.,contains("_")))) %>% select(rowsum, ID=rowname) %>% left_join(y=select(annotation,ID, strand,exonCount,circType,geneID, geneName, geneType), by="ID") %>% slice_max(rowsum, n = 20)
# top endocytosisKEGG genes
Merge_circexp_raw[as.character(annotation %>% filter(geneName %in% endocytosisKEGG) %>% pull(ID)),] %>% rownames_to_column() %>% mutate(rowsum=rowSums(select(.,contains("_")))) %>% select(rowsum, ID=rowname) %>% left_join(y=annotation, by="ID") %>% slice_max(rowsum, n = 20) 

# joint
rownames_to_column(Merge_circexp_raw) %>% mutate(rowsum=rowSums(select(.,contains("_")))) %>% select(rowsum, ID=rowname) %>% 
  left_join(y=select(annotation,ID, strand,exonCount,circType,geneID, geneName, geneType), by="ID") %>% 
  mutate(geneID=gsub("\\..*","", geneID)) %>%
  mutate(is.ADgene=geneName %in% ADgenes,
         is.PDgene=geneName %in% PDgenes,
         is.synGO =geneID %in% synGO,
         is.Endocy=geneID %in% endocytosisGO) %>%
  filter(is.PDgene, is.ADgene | is.synGO | is.Endocy) %>% 
  #filter(geneName=="NRXN3")
  slice_max(rowsum, n = 50) 


pdf("~/projects/circRNA/results/Merge_circexplorer_BC190.annotation.ADPDsynGOendoGO.pie.pdf", width=4, height = 4)
pie(table(PDgenes %in% annotation$geneName),labels = c(NA,paste0(table(PDgenes %in% annotation$geneName)[2],"/",sum(table(PDgenes %in% annotation$geneName)),"\n(",round(100*mean(PDgenes %in% annotation$geneName),2),"%)")), border='gray', col=c('white','gray'),main="PD-risk genes producing circRNAs")
pie(table(ADgenes %in% annotation$geneName),labels = c(NA,paste0(table(ADgenes %in% annotation$geneName)[2],"/",sum(table(ADgenes %in% annotation$geneName)),"\n(",round(100*mean(ADgenes %in% annotation$geneName),2),"%)")), border='gray', col=c('white','gray'),main="AD-risk genes producing circRNAs")
pie(table(synGO %in% gsub("\\..*","", annotation$geneID)),labels = c(NA,paste0(table(synGO %in% gsub("\\..*","", annotation$geneID))[2],"/",sum(table(synGO %in% gsub("\\..*","", annotation$geneID))),"\n(",round(100*mean(synGO %in% gsub("\\..*","", annotation$geneID)),2),"%)")), border='gray', col=c('white','gray'),main="Synaptic genes producing circRNAs")
pie(table(endocytosisGO %in% gsub("\\..*","", annotation$geneID)),labels = c(NA,paste0(table(endocytosisGO %in% gsub("\\..*","", annotation$geneID))[2],"/",sum(table(endocytosisGO %in% gsub("\\..*","", annotation$geneID))),"\n(",round(100*mean(endocytosisGO %in% gsub("\\..*","", annotation$geneID)),2),"%)")), border='gray', col=c('white','gray'),main="Endocytosis(GO) genes producing circRNAs")
#pie(table(endocytosisKEGG %in% annotation$geneName),labels = c(NA,paste0(table(endocytosisKEGG %in% annotation$geneName)[2],"/",sum(table(endocytosisKEGG %in% annotation$geneName)),"\n(",round(100*mean(endocytosisKEGG %in% annotation$geneName),2),"%)")), border='gray', col=c('white','gray'),main="Endocytosis(KEGG) genes producing circRNAs")
dev.off()

# Fisher's exact test
# GENCODE
# GENCODEv19=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v19.annotation.genes.bed", col.names = c("chr","start","end","geneName","geneID","geneType","strand"), stringsAsFactors = F, header = F)
dim(GENCODEv19)
fisher.test(table(is_synapse = gsub("\\..*","", GENCODEv19$geneID) %in% synGO, is_circRNA = GENCODEv19$geneID %in% annotation$geneID))

## Queston: Does the percentage of circRNA-producing synaptic genes decrease along the HC, ILB and PD?
#dim(annotation); dim(Merge_circexp_raw); head(Merge_circexp_raw)
#Merge_circexp_norm_filtered_and_enriched=readRDS("Merge_circexplorer_BC190.filtered.enriched.normRPM.rds"); dim(Merge_circexp_norm_filtered_and_enriched)
percentage_of_synGOgene_among_all_hostgenes = apply(Merge_circexp_raw, 2, function(x) {
  mean(unique(gsub("\\..*","", annotation$geneID[match(names(x[x>0]), annotation$ID)])) %in% synGO)
})
library(ggpubr)
as.data.frame(percentage_of_synGOgene_among_all_hostgenes) %>% rownames_to_column() %>% 
  separate(rowname, c("Dx",NA,"cellType"), sep="_", remove = F, extra='drop') %>% 
  filter(cellType!="MCPY") %>% 
  #filter(Dx!="ILB") %>% 
  mutate(Dx=ifelse(Dx=="ILB","HC",Dx)) %>%
  mutate(Dx=factor(Dx, levels = c("HC","ILB","PD","AD"))) %>%
  ggboxplot(x = "cellType", y = "percentage_of_synGOgene_among_all_hostgenes",
            ylab= "Percentage of synaptic genes among all host genes",
            color = "Dx", palette = "jco",
            add = "jitter",
            short.panel.labs = FALSE) + stat_compare_means(aes(group = Dx), method = "t.test") + 
  ggsave("../results/Merge_circexplorer_BC190.annotation.syn.vs.ADPD.boxplot.pdf", width = 4, height = 4)

## Question: How many AD/PD GWAS risk genes transcribed circRNAs in their corresponding neuron?
# PD >>
# all cirsRNAs
table(PDgenes %in% annotation$geneName)
# 96/109
# DA neuronal circRNAs
table(PDgenes %in% (annotation %>% filter(ID %in% (Merge_circexp_raw %>% select(contains("_SNDA_")) %>% rownames_to_column() %>% filter_at(vars(-rowname), any_vars(. != 0)) %>% pull(rowname))) %>% pull(geneName)))
# 86/109  

# AD >>
# all cirsRNAs
table(ADgenes %in% annotation$geneName)
# 159/217
# PY neuronal circRNAs
table(ADgenes %in% (annotation %>% filter(ID %in% (Merge_circexp_raw %>% select(contains("PY_")) %>% rownames_to_column() %>% filter_at(vars(-rowname), any_vars(. != 0)) %>% pull(rowname))) %>% pull(geneName)))
# 131/217  

## -------------------
### overlap with endocytosis genes etc. - related with new Fig3
## see comment from Clemens below:
# -	Try one alternative ways to visualize this in a more attractive manner. For example, we could show 4 pie graphs; each representing the totality of circRNAs. Then show the % of circRNAs in select key pathways 
# - A (eg endocytosis-related: endocytosis, axon guidance, adherence junction (check and correcct spelling in figure), long-term potentiation, long term depression => all related to endo/syn) as slices in the pies. 
# - B. pathways in cancers (eg pathways in cancer, thryodi cancer, prostate cncer,endometrial cncer etc)as another slice (…enriched in NN).
# -	C. Tight junction, gap junction
# -	D. Erb signalling
# -	E. Other sig pathways
## -------------------

# read "gene - pathway" relationship from MSigDB
library('fgsea')
genesets_list = fgsea::gmtPathways("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/msigdb_v7.2/msigdb_v7.2_GMTs/c2.cp.kegg.v7.2.symbols.gmt")
# genesets is a list, convert to a data frame
library(tidyverse)
genesets_df = data.frame(pathwayID=rep(names(genesets_list), lengths(genesets_list)), geneName=unlist(genesets_list, use.names = F))
# add pathway to circRNA annotation list
annotation_filtered_enriched=readRDS(file="Merge_circexplorer_BC190.filtered.enriched.annotation.bed14.rds"); dim(annotation_filtered_enriched)
annotation_filtered_enriched %>% left_join(y=genesets_df, by="geneName") %>% filter(!is.na(pathwayID)) %>% pull(ID) %>% unique() %>% length() ## only 2689 (out of 11039) circRNAs with host genes in KEGG pathway somehow (either due to the incomplete of KEGG or the gene name convertion)

# read the "gene - pathway" relationship from the latest KEGG download
kegg.enrichment <- function(circRNAset, output_filename="../results/circRNAs.hostgenes.kegg.enrichment.pie.pdf", path2genes="~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/KEGG/path2genes.tab"){
  # debug: circRNAset=annotation_filtered_enriched; output_filename="../results/Merge_circexplorer_BC190.filtered.enriched.annotation.hostgenes.kegg.enrichment.pie.pdf"; path2genes="~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/KEGG/path2genes.tab"
  # get KEGG path to gene relationship
  if(file.exists(path2genes)) path2genes = read.delim(path2genes, header = T) else {
    g=read.delim("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/KEGG/kegg.genes.list", sep = "\t", header = F, col.names = c("geneID", "geneNames","geneDesp"))
    p2g=read.delim("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/KEGG/pathway.genes.list", sep = "\t", header = F, col.names = c("pathID","geneID"))
    p=read.delim("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/KEGG/pathway.genes.tab", sep = "\t", header = F, col.names = c("pathID","pathName","geneNames"))
    path2genes = left_join(p2g,g,by="geneID") %>% group_by(pathID) %>% mutate(geneName = paste0(geneNames, collapse = ", ")) %>% select(pathID, geneName) %>% # tail()
      separate_rows(geneName, sep=", ") %>% mutate(pathID=sub("path:","",pathID)) %>% left_join(select(p,pathID, pathName), by="pathID") %>% distinct()    
    write.table(path2genes, file = "~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/KEGG/path2genes.tab", col.names = T, row.names = F, quote = F, sep = "\t")
  }
  head(path2genes)
  
  # annotation_filtered_enriched %>% left_join(y=path2genes, by="geneName") %>% filter(!is.na(pathID)) %>% pull(ID) %>% unique() %>% length() ## still, only 4318 (out of 11039) circRNAs with host genes found in KEGG pathway somehow (either due to the incomplete of KEGG or the gene name convertion)
  ## it seems that KEGG is not a good choice as it only covers a small subset of genes (total 5245 genes in mSigDbv7.2 KEGG set)
  ## we still use KEGG as temporary solution, but in long term, we might need to switch to GO 
  
  x=circRNAset %>% left_join(y=path2genes, by="geneName") %>% filter(!is.na(pathID)) %>% 
    mutate(group = ifelse(pathName %in%c("Endocytosis", "Axon guidance", "Adherens junction", "Long-term potentiation", "Long-term depression") | grepl("synapse", ignore.case =T, pathName), "A_endocytosis", 
                          ifelse(grepl("cancer",ignore.case =T, pathName), "B_cancer", 
                                 ifelse(grepl("junction",ignore.case =T, pathName), "C_junction", 
                                        ifelse(pathName %in%c("ErbB signaling pathway"), "D_Erbb", "E_others")))))
  ## prioritize in order of A->E
  df = select(x, ID, group) %>% distinct() %>% mutate(value=1) %>% pivot_wider(id_cols = ID, names_from=group, values_from=value)
  A=filter(df, A_endocytosis==1)
  B=filter(df, is.na(A_endocytosis)) %>% filter(B_cancer==1)
  C=filter(df, is.na(A_endocytosis), is.na(B_cancer)) %>% filter(C_junction==1)
  D=filter(df, is.na(A_endocytosis), is.na(B_cancer), is.na(C_junction)) %>% filter(D_Erbb==1)
  E=filter(df, is.na(A_endocytosis), is.na(B_cancer), is.na(C_junction), is.na(D_Erbb)) %>% filter(E_others==1)
  df=rbind(data.frame(ID=A$ID, group="A_endocytosis"),
          data.frame(ID=B$ID, group="B_cancer"),
          data.frame(ID=C$ID, group="C_junction"),
          data.frame(ID=D$ID, group="D_Erbb"),
          data.frame(ID=E$ID, group="E_others"))
  
  pdf(output_filename, width = 5, height = 5)
  par(mar=c(4,4,4,4))
  pie(table(df$group),radius = 1, cex = 0.4, labels = paste0(names(table(df$group))," (",round(100*table(df$group)/sum(table(df$group)),2),"%, ", table(df$group), "/", sum(table(df$group)),")"))
  par(mar=c(5,15,3,3))
  barplot(sort(table(x$pathName)), xlim=c(100,800),horiz = T, las=2, cex.names = 0.2)
  dev.off()
}

kegg.enrichment(annotation_filtered_enriched, output_filename="../results/Merge_circexplorer_BC190.filtered.enriched.annotation.hostgenes.kegg.enrichment.pie.pdf")
## go to cell-specific circRNAs
DF3=readRDS("Merge_circexplorer_BC109.cellspecific_heatmap.circRNA3.rds") ## see getTissueSpecificRNA.R
for(i in c("NN","PY","SNDA")) {
  message(paste("processing cell type",i,"..."));
  kegg.enrichment(annotation_filtered_enriched[as.character(annotation_filtered_enriched$ID) %in% unique(subset(DF3, celltype==i, select = 'gene',drop = T)),], output_filename=paste0("../results/Merge_circexplorer_BC109.cellspecific_heatmap5.genes3.",i,".kegg.enrichment.pie.pdf"))
}

## -------------------
## conserved to Drosophila
## -------------------
library('biomaRt') # BiocManager::install('biomaRt',ask=F); 
gmart = useEnsembl(biomart = "ensembl",version=75, dataset = "hsapiens_gene_ensembl") # Ensembl v75 is GRCh37.p13, which is same as GENCODE v19
listAttributes(gmart, page='homologs') %>% filter(grepl("dmelanogaster", name));
listFilters(gmart);
attributes = c("ensembl_gene_id",
               "dmelanogaster_homolog_ensembl_gene","dmelanogaster_homolog_associated_gene_name",
               "dmelanogaster_homolog_orthology_type",
               "dmelanogaster_homolog_perc_id", "dmelanogaster_homolog_perc_id_r1")
orth.human = getBM(attributes,mart = gmart,
                   filters="biotype",values='protein_coding', 
                   uniqueRows=TRUE)
dim(orth.human);head(orth.human); 
total_human_genes = length(unique(orth.human$ensembl_gene_id))
orth.human = mutate(orth.human, ave_perc_id=dmelanogaster_homolog_perc_id+dmelanogaster_homolog_perc_id_r1) %>% 
  group_by(ensembl_gene_id) %>% 
  filter(ave_perc_id==max(ave_perc_id)) %>% 
  filter(1:n() == 1) # make sure each gene only assigned one ortholog (the best recipocal one)
dim(orth.human);head(orth.human)
total_human_genes_with_fly_orth = length(unique(orth.human$ensembl_gene_id))

annotation_filtered_enriched$drosophila_orth_symbol=orth.human$dmelanogaster_homolog_associated_gene_name[match(gsub("\\..*","",annotation_filtered_enriched$geneID), orth.human$ensembl_gene_id)] # the unfound will be NA
head(annotation_filtered_enriched)
circRNA_fly_orth=table(select(annotation_filtered_enriched, geneID, geneName, geneType, drosophila_orth_symbol) %>% distinct() %>% filter(geneType=='protein_coding') %>% mutate(has_drosophila_orth=!is.na(drosophila_orth_symbol)) %>% pull(has_drosophila_orth))

x0= rbind(circRNA_hostgene = circRNA_fly_orth, all_genes=c(total_human_genes-total_human_genes_with_fly_orth, total_human_genes_with_fly_orth))
x=x0/rowSums(x0)

# those fly genes with circRNAs
# options(timeout=200); download.file("https://www.picb.ac.cn/rnomics/circpedia/static/download_cache/fly_dm6_All_circRNA.csv", destfile = "circpedia.fly_dm6_All_circRNA.csv")
CIRCpedia2_fly = read.csv("https://www.picb.ac.cn/rnomics/circpedia/static/download_cache/fly_dm6_All_circRNA.csv", header = F)
dim(CIRCpedia2_fly); head(CIRCpedia2_fly)

pdf("~/projects/circRNA/results/Merge_circexplorer_BC190.filtered.enriched.conserved2fly.pdf", width = 4, height = 5)
bp=barplot(t(x)[2:1,2:1], col=c('red','white'), ylab = "Proportion of genes conserved to Drosophila")
mtext(at=bp,rowSums(x0)[2:1])

pie(table(orth.human$dmelanogaster_homolog_associated_gene_name %in% CIRCpedia2_fly$V3), labels = NA, col=c('red','darkred'), border='red', lwd=4, main="among all genes with fly orthology")
pie(table(select(annotation_filtered_enriched, geneID, geneType, drosophila_orth_symbol) %>% distinct() %>% filter(geneType=='protein_coding', !is.na(drosophila_orth_symbol)) %>% pull(drosophila_orth_symbol) %in% CIRCpedia2_fly$V3), labels = NA, col=c('red','darkred'), border='red', lwd=4, main="among all circRNA host genes with fly orthology")
dev.off()

## -------------------
## draw histogram of supported reads per circRNAs
## -------------------

total_circRNA_raw_reads = rowSums(Merge_circexp_raw)

pdf("~/projects/circRNA/results/total_circRNA_raw_reads.BC190.pdf", width = 6, height = 5)
BREAK=300
layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(5,1), heights=c(1,1))
ht=hist(total_circRNA_raw_reads[total_circRNA_raw_reads<=BREAK], breaks = 0:BREAK, plot=F)
par(mar=c(4,4,0,1))
plot(log10(ht$counts+.1)+1, lwd=2,xlim=c(0,BREAK), type='h', ylim = c(0.9,6.6),
     col=c("#000000",rep('#aaaaaa',199)), 
     ylab='Number of circular RNAs', xlab="Number of backspliced reads",
     yaxt="n",xaxs="i", yaxs="i",bty="n")
y1=floor(range(log10(ht$counts+.1)+1))
pow <- seq(y1[1], y1[2])
ticksat <- as.vector(sapply(pow-1, function(p) log10((2:9)*10^p)))
axis(2, pow, labels = 10^(pow-1), las=1)
axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
# add the remained ones with the #2 filter: enriched
filtered_circRNA_raw_reads = rowSums(Merge_circexp_raw[rownames(Merge_circexp_raw_filtered_and_enriched),])
ht=hist(filtered_circRNA_raw_reads[filtered_circRNA_raw_reads<=BREAK], breaks = 0:BREAK, plot=F)
points(log10(ht$counts+.1)+1, lwd=2, type='h',col='#8B0000')

# add the remained ones with the #3 filter: circRNA only
circRNAonly=as.character(annotation_filtered_enriched$ID[annotation_filtered_enriched$circType=='circRNA'])
filtered_circRNA_raw_reads_circRNAonly = rowSums(Merge_circexp_raw[circRNAonly,])
ht=hist(filtered_circRNA_raw_reads_circRNAonly[filtered_circRNA_raw_reads_circRNAonly<=BREAK], breaks = 0:BREAK, plot=F)
points(log10(ht$counts+.1)+1, lwd=2, type='h',col='#ff0000')

legend("topleft",c(paste0("Distinct circular RNAs (n = ",format(length(total_circRNA_raw_reads),big.mark=","),")"),
                   paste0("- at least 2 reads in overall samples (n = ", format(sum(total_circRNA_raw_reads>=2),big.mark=","),")"),
                   paste0("-- enriched in RNase R (n = ", format(length(filtered_circRNA_raw_reads),big.mark=","),")"),
                   paste0("--- circRNA only (n = ", format(length(filtered_circRNA_raw_reads_circRNAonly),big.mark=","),")")),
       col=c("#000000",'#aaaaaa','#8B0000','#ff0000'),
       text.col = c("#000000",'#aaaaaa','#8B0000','#ff0000'),
       bty='n', xpd=T)

## for the BREAK- region
par(mar=c(4,1,0,1))
n200=total_circRNA_raw_reads[total_circRNA_raw_reads>BREAK]; 
MAX=max(total_circRNA_raw_reads)
ht=hist(n200, breaks = BREAK:MAX, plot=F)
plot(log10(ht$counts+.1)+1, lwd=2,xlim=range(total_circRNA_raw_reads), type='h', ylim = c(0.9,6.6),col='#aaaaaa', ylab="",xlab="",xaxt="n",yaxt='n',xaxs="i", yaxs="i",bty="n")
y1=floor(range(log10(ht$counts+.1)+1))
pow <- seq(y1[1], y1[2])
ticksat <- as.vector(sapply(pow-1, function(p) log10((2:9)*10^p)))
#axis(2, pow, labels = 10^(pow-1), las=1)
#axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
axis(1, c(0,seq(10000,MAX,20000)-BREAK), labels=c(BREAK,paste0(seq(10000,MAX,20000)/1000,"k")))

# add the remained ones with the #2 filter
n200=filtered_circRNA_raw_reads[filtered_circRNA_raw_reads>BREAK];
ht=hist(n200, breaks = BREAK:MAX, plot=F)
points(log10(ht$counts+.1)+1, lwd=2, type='h',col='#8B0000')

# add the remained ones with the #3 filter: circRNA only
ht=hist(filtered_circRNA_raw_reads_circRNAonly[filtered_circRNA_raw_reads_circRNAonly>BREAK], breaks = BREAK:MAX, plot=F)
points(log10(ht$counts+.1)+1, lwd=2, type='h',col='#ff0000')

dev.off()

## add gene highlight in the figure
annotation_filtered_enriched %>% dplyr::select(ID, circType, geneName) %>% 
  left_join(rownames_to_column(Merge_circexp_raw_filtered_and_enriched), by = c("ID"="rowname")) %>% 
  mutate(sumVar = rowSums(dplyr::select(., contains("_")))) %>% dplyr::select(-contains("_")) %>% 
  filter(circType=='circRNA') %>% arrange(-sumVar) %>% 
  filter(sumVar>10) %>%
  filter(geneName %in% c("SLC8A1", 'SNAP25', 'SV2A', 'RAB3A', 'SYNGR1','NRXN3', 'SYNPR','STXBP1', 'ERC2','CASK','PPFIA1','PPFIA4','UNC13B', "DNAJC6","RIMS1","RIMS2",
                         "SORL1","PSEN1","APP", "ABCA7","CLU","CR1","PICALM","PLD3","TREM2",'UCHL1','SRY',"SNCA","LRRK2","PINK1","ABCA7","CLU","CR1","PLD3","TREM2",'PARK7', 'PRKN', 'GBA'))

## BC190
# ID circType geneName sumVar
# 1     chr2_40655612_40657444  circRNA   SLC8A1   8556
# 2     chr3_55717821_56026278  circRNA     ERC2   1264
# 3     chr2_40655612_40657441  circRNA   SLC8A1   1188
# 4     chr6_73016960_73043538  circRNA    RIMS1    977
# 5     chr6_73005639_73043538  circRNA    RIMS1    974
# 6    chr11_85707868_85714494  circRNA   PICALM    769
# 7     chr1_65830317_65831879  circRNA   DNAJC6    676
# 8    chr11_85707868_85742653  circRNA   PICALM    553
# 9     chr3_55984452_56026278  circRNA     ERC2    316
# 10    chr3_55768798_56026278  circRNA     ERC2    239
# 11    chr2_40366540_40405633  circRNA   SLC8A1    226
# 12    chr6_72960032_72961071  circRNA    RIMS1    199
# 13  chr8_105080739_105161076  circRNA    RIMS2    196
# 14   chr14_80130120_80164280  circRNA    NRXN3    164
# 15    chr6_72960032_73043538  circRNA    RIMS1    159
# 16   chr11_85701292_85742653  circRNA   PICALM    130
# 17   chr11_85685750_85742653  circRNA   PICALM    123
# 18  chr8_105105698_105161076  circRNA    RIMS2    123
# 19    chr3_55922416_56026278  circRNA     ERC2    117
# 20   chr14_73614502_73614814  circRNA    PSEN1    103
# 21    chr6_73001636_73043538  circRNA    RIMS1    102
# 22  chr8_105053553_105161076  circRNA    RIMS2     93
# 23    chr3_55717821_55984588  circRNA     ERC2     85
# 24    chr1_65830317_65860715  circRNA   DNAJC6     84
# 25    chr3_55922416_56041349  circRNA     ERC2     71
# 26   chr21_27326903_27354790  circRNA      APP     70
# 27   chr14_73614502_73640415  circRNA    PSEN1     65
# 28    chr9_35295692_35313986  circRNA   UNC13B     65
# 29    chr3_55733405_55922577  circRNA     ERC2     52
# 30    chr3_55717821_55733540  circRNA     ERC2     48
# 31    chr3_55922416_55984588  circRNA     ERC2     48
# 32    chr3_55984452_56041349  circRNA     ERC2     48
# 33   chr14_73614502_73614802  circRNA    PSEN1     45
# 34    chr3_55717821_55768946  circRNA     ERC2     40
# 35   chr11_85718584_85742653  circRNA   PICALM     39
# 36    chr8_27462440_27464041  circRNA      CLU     33
# 37    chr3_55768798_55984588  circRNA     ERC2     31
# 38   chr11_85707868_85712201  circRNA   PICALM     25
# 39  chr8_105025669_105161076  circRNA    RIMS2     25
# 40    chrX_41428920_41469278  circRNA     CASK     25
# 41   chr14_73614502_73659572  circRNA    PSEN1     22
# 42   chr21_27347382_27372497  circRNA      APP     22
# 43    chr6_72952016_73043538  circRNA    RIMS1     20
# 44    chr6_73001636_73023375  circRNA    RIMS1     19
# 45   chr11_85707868_85718626  circRNA   PICALM     18
# 46   chr11_70179575_70184559  circRNA   PPFIA1     17
# 47    chr2_40387904_40405633  circRNA   SLC8A1     17
# 48   chr11_85685750_85695016  circRNA   PICALM     16
# 49   chr11_70176278_70185412  circRNA   PPFIA1     14
# 50   chr11_85692171_85742653  circRNA   PICALM     14
# 51   chr14_73614502_73664837  circRNA    PSEN1     14
# 52  chr8_104943490_105161076  circRNA    RIMS2     14
# 53    chr8_27462440_27462852  circRNA      CLU     14
# 54  chr1_155207131_155209868  circRNA      GBA     13
# 55   chr21_27326903_27372497  circRNA      APP     13
# 56    chr3_56173536_56183160  circRNA     ERC2     13
# 57 chr11_121348826_121367758  circRNA    SORL1     12
# 58   chr11_85722072_85742653  circRNA   PICALM     12
# 59    chr1_65830317_65845207  circRNA   DNAJC6     12
# 60    chr1_65830317_65871816  circRNA   DNAJC6     12
# 61  chr8_104922361_105026843  circRNA    RIMS2     12
# 62  chr8_105025669_105080840  circRNA    RIMS2     12

## ===============================================
## merge into gene  ## see Supplementary fig. ?
## ===============================================
setwd("~/projects/circRNA/data/") 
Merge_circexp_raw=readRDS("Merge_circexplorer_BC190.rawcount.rds")
Merge_circexp_raw_filtered_and_enriched=readRDS("Merge_circexplorer_BC190.filtered.enriched.rawcount.rds")
annotation=readRDS("Merge_circexplorer_BC190.annotation.bed14.rds")

genes_annotation = read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/annotation.genes.bed6+3", sep="\t", quote="", header = F, stringsAsFactors = F, 
                              col.names = c("chr","start","end","geneID","score","strand","geneName","geneType","geneDescription"));
annotation$geneDescription = genes_annotation$geneDescription[match(sub("\\..*","", annotation$geneID), genes_annotation$geneID)]
annotation$geneDescription = sub(" \\[Source:.*","", annotation$geneDescription)
annotation2gene = dplyr::select(annotation, ID, circType, geneID, geneName, geneType, geneDescription) %>% mutate(circID=paste0(sub("RNA","",as.character(circType)),as.character(geneName)))
dim(Merge_circexp_raw);  Merge_circexp_raw2gene = rownames_to_column(Merge_circexp_raw) %>% mutate(rowname=annotation2gene$circID[match(rowname, annotation2gene$ID)]) %>% group_by(rowname) %>% summarise_all(sum) %>% column_to_rownames(); dim(Merge_circexp_raw2gene)
dim(annotation2gene); annotation2gene = mutate(annotation2gene,ID=circID) %>% dplyr::select(-circID) %>% distinct(); dim(annotation2gene)
total_circRNA_raw_reads2gene = rowSums(Merge_circexp_raw2gene)

pdf("~/projects/circRNA/results/total_circRNA_raw_reads2genes.BC190.pdf", width = 6, height = 5)
BREAK=500
layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(5,1), heights=c(1,1))
ht=hist(total_circRNA_raw_reads2gene[total_circRNA_raw_reads2gene<=BREAK], breaks = 0:BREAK, plot=F)
MAXy=ceiling(max(log10(ht$counts+.1)+1))
par(mar=c(4,4,0,1))
plot(log10(ht$counts+.1)+1, lwd=2,xlim=c(0,BREAK), type='h', ylim = c(0.9,MAXy+0.5),col=c("#000000",rep('#aaaaaa',BREAK-1)), ylab='Number of circular RNAs', xlab="Number of backspliced reads",yaxt="n",xaxs="i", yaxs="i",bty="n")
y1=floor(range(log10(ht$counts+.1)+1))
pow <- seq(y1[1], y1[2])
ticksat <- as.vector(sapply(pow-1, function(p) log10((2:9)*10^p)))
axis(2, pow, labels = 10^(pow-1), las=1)
axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
# # add the remained ones with the #2 filter: enriched
# filtered_circRNA_raw_reads = rowSums(Merge_circexp_raw[rownames(Merge_circexp_raw_filtered_and_enriched),])
# ht=hist(filtered_circRNA_raw_reads[filtered_circRNA_raw_reads<=BREAK], breaks = 0:BREAK, plot=F)
# points(log10(ht$counts+.1)+1, lwd=2, type='h',col='#8B0000')
# 
# # add the remained ones with the #3 filter: circRNA only
# circRNAonly=as.character(annotation_filtered_enriched$ID[annotation_filtered_enriched$circType=='circRNA'])
# filtered_circRNA_raw_reads_circRNAonly = rowSums(Merge_circexp_raw[circRNAonly,])
# ht=hist(filtered_circRNA_raw_reads_circRNAonly[filtered_circRNA_raw_reads_circRNAonly<=BREAK], breaks = 0:BREAK, plot=F)
# points(log10(ht$counts+.1)+1, lwd=2, type='h',col='#ff0000')

legend("topleft",c(paste0("Distinct collapsed circular RNAs (n = ",format(length(total_circRNA_raw_reads2gene),big.mark=","),")"),
                   paste0("- at least 2 reads in overall samples (n = ", format(sum(total_circRNA_raw_reads2gene>1),big.mark=","),")")),
       col=c("#000000",'#aaaaaa'),
       text.col = c("#000000",'#aaaaaa'),
       bty='n', xpd=T)

## for the BREAK- region
par(mar=c(4,1,0,1))
n200=total_circRNA_raw_reads2gene[total_circRNA_raw_reads2gene>BREAK]; 
MAX=max(total_circRNA_raw_reads2gene)
ht=hist(n200, breaks = BREAK:MAX, plot=F)
plot(log10(ht$counts+.1)+1, lwd=2,xlim=range(total_circRNA_raw_reads2gene), type='h', ylim = c(0.9,MAXy+0.5),col='#aaaaaa', ylab="",xlab="",xaxt="n",yaxt='n',xaxs="i", yaxs="i",bty="n")
y1=floor(range(log10(ht$counts+.1)+1))
pow <- seq(y1[1], y1[2])
ticksat <- as.vector(sapply(pow-1, function(p) log10((2:9)*10^p)))
#axis(2, pow, labels = 10^(pow-1), las=1)
#axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
axis(1, c(0,seq(10000,MAX,20000)-BREAK), labels=c(BREAK,paste0(seq(10000,MAX,20000)/1000,"k")))

# # add the remained ones with the #2 filter
# n200=filtered_circRNA_raw_reads[filtered_circRNA_raw_reads>BREAK];
# ht=hist(n200, breaks = BREAK:MAX, plot=F)
# points(log10(ht$counts+.1)+1, lwd=2, type='h',col='#8B0000')
# 
# # add the remained ones with the #3 filter: circRNA only
# ht=hist(filtered_circRNA_raw_reads_circRNAonly[filtered_circRNA_raw_reads_circRNAonly>BREAK], breaks = BREAK:MAX, plot=F)
# points(log10(ht$counts+.1)+1, lwd=2, type='h',col='#ff0000')

dev.off()

## add gene highlight in the figure
annotation2gene %>% dplyr::select(ID, circType, geneName, geneDescription) %>% 
  left_join(rownames_to_column(Merge_circexp_raw2gene), by = c("ID"="rowname")) %>% 
  mutate(sumReads = rowSums(dplyr::select(., contains("_")))) %>% dplyr::select(-contains("_")) %>% 
  filter(circType=='circRNA') %>% arrange(-sumReads) %>% 
  #filter(sumReads>500) %>% 
  write.table("~/projects/circRNA/results/total_circRNA_raw_reads2genes.BC190.xls", col.names = T, row.names = F, sep = "\t", quote = F)
  filter(geneName %in% c("SLC8A1", 'SNAP25', 'SV2A', 'RAB3A', 'SYNGR1','NRXN3', 'SYNPR','STXBP1', 'ERC2','CASK','PPFIA1','PPFIA4','UNC13B', "DNAJC6","RIMS1","RIMS2",
                         "SORL1","PSEN1","APP", "ABCA7","CLU","CR1","PICALM","PLD3","TREM2",'UCHL1','SRY',"SNCA","LRRK2","PINK1","ABCA7","CLU","CR1","PLD3","TREM2",'PARK7', 'PRKN', 'GBA'))

## make sub-network of these circRNAs-hosting genes
hostGenes = read.table("~/projects/circRNA/results/total_circRNA_raw_reads2genes.BC190.xls", sep="\t", header=T, stringsAsFactors = F)
head(hostGenes); dim(hostGenes)
# read BioPlex
bioplex = read.table("~/projects/circRNA/data/BioPlex_interactionList_v4a.tsv", header = T, sep = "\t", stringsAsFactors = F)
head(bioplex); dim(bioplex)
length(unique(c(bioplex$SymbolA, bioplex$SymbolB)))
sum(as.character(hostGenes$geneName) %in% unique(c(bioplex$SymbolA, bioplex$SymbolB)))
inner_join(bioplex, dplyr::select(hostGenes, geneName, sumReadsA=sumReads), by=c("SymbolA"= "geneName")) %>% inner_join(y=dplyr::select(hostGenes, geneName, sumReadsB=sumReads), by=c("SymbolB"= "geneName")) %>% write.table("~/projects/circRNA/data/BioPlex_interactionList_v4a.circRNA_hostGene.tsv", quote = F, sep = "\t", row.names = F)
filter(bioplex, SymbolA %in% hostGenes$geneName, SymbolB %in% hostGenes$geneName) %>% write.table("~/projects/circRNA/data/BioPlex_interactionList_v4a.circRNA_hostGene.only2only.tsv", quote = F, sep = "\t", row.names = F)
filter(hostGenes, geneName %in% unique(c(bioplex$SymbolA, bioplex$SymbolB))) %>% dplyr::select(name=geneName, sumReads) %>% write.table("~/projects/circRNA/data/BioPlex_interactionList_v4a.circRNA_hostGene.only2only.annotation.tsv", quote = F, sep = "\t", row.names = F) 

## correlation between expression level: circRNA vs. hostgenes
# exon-circularized reads (no ciRNA)
df1 = read.table("~/projects/circRNA/results/total_circRNA_raw_reads2genes.BC190.xls", sep="\t", header=T, stringsAsFactors = F)
dim(df1); df1$geneName[df1$geneName=='CDR1as']='CDR1'
hostgenes = read.table("~/neurogen/rnaseq_PD/results/merged/genes.htseqcount.cufflinks.allSamples.BCv2.uniq.xls", sep="\t", header=T, stringsAsFactors = F, check.names = F)
df2 = mutate(hostgenes, tracking_id=gsub("\\..*","",tracking_id)) %>% mutate(tracking_id=genes_annotation$geneName[match(tracking_id, genes_annotation$geneID)]) %>% group_by(tracking_id) %>% summarise_all(sum) %>% filter(!is.na(tracking_id)) %>% column_to_rownames(var="tracking_id") %>% select(colnames(Merge_circexp_raw)) %>% rowSums()
# all circular reads (circRNA + ciRNA)
df11 = rownames_to_column(Merge_circexp_raw) %>% mutate(rowname=annotation$geneName[match(rowname, annotation$ID)]) %>% group_by(rowname) %>% summarise_all(sum) %>% filter(!is.na(rowname)) %>% column_to_rownames() %>% rowSums();
names(df11)[names(df11)=='CDR1as']='CDR1'

pdf("~/projects/circRNA/results/total_circRNA_raw_reads2genes.vs.mRNA.BC190.pdf", width = 4, height = 4)
cortest = cor.test(df1$sumReads, df2[df1$geneName],alternative = "two.sided", method = "spearman")
par(pty="s"); plot(df1$sumReads, df2[df1$geneName], pch=16, cex=0.8, col='#00000066', log = 'xy', asp=1,
                   xlab="exon-circularized reads (circRNA only) per gene", ylab="linear reads", 
                   main="Correlation btw circular and linear reads")
abline(a=0, b=1, col='red')
mtext(paste("Pearson's r =", round(cortest$estimate, 2), "(p-value:", format.pval(cortest$p.value),")"))

cortest = cor.test(df11, df2[names(df11)],alternative = "two.sided", method = "spearman")
par(pty="s"); plot(df11, df2[names(df11)], pch=16, cex=0.8, col='#00000066', log = 'xy', asp=1,
                   xlab="circular reads per gene", ylab="linear reads", 
                   main="Correlation btw circular and linear reads")
abline(a=0, b=1, col='red')
mtext(paste("Pearson's r =", round(cortest$estimate, 2), "(p-value:", format.pval(cortest$p.value),")"))
dev.off()

## ===============================================
## back-splicing reads vs. sample fraction
## similar to Fig. 3E in https://www-sciencedirect-com.ezp-prod1.hul.harvard.edu/science/article/pii/S0092867418316350?via%3Dihub#fig3
## ===============================================
Merge_circexp_raw_filtered_and_enriched=readRDS("Merge_circexplorer_BC190.filtered.enriched.rawcount.rds")
dim(Merge_circexp_raw_filtered_and_enriched)
df=Merge_circexp_raw_filtered_and_enriched
pdf("../results/Merge_circexplorer_BC190.filtered.enriched.nReads_vs_nSample.pdf", width = 3, height = 3)
plot(jitter(rowMeans(df>0))*100, rowSums(df), xlab="Sample fraction (%) ", ylab="Backspliced reads", pch=21, col='gray', bg='#00000033', log='y')
dev.off()
## only few circRNAs expressed in 90% samples in our case. So, it's applicable to divide into 'high' vs. 'low' like the paper above. 

## ===============================================
## Number of circRNAs per host gene
## ===============================================
require(scales)
mylog_trans <- function (base = exp(1), from = 0) 
  {
    trans <- function(x) log(x, base) - from
    inv <- function(x) base^(x + from)
    trans_new("mylog", trans, inv, log_breaks(base = base), domain = c(base^from, Inf))
}
annotation_filtered_enriched=readRDS("Merge_circexplorer_BC190.filtered.enriched.annotation.bed14.rds")
annotation_filtered_enriched_n = annotation_filtered_enriched %>% filter(circType=='circRNA') %>% group_by(geneName) %>% summarise(n=n()) %>% group_by(n) %>% summarise(N=n(), geneNames=paste(geneName, collapse="; "))
ggplot(annotation_filtered_enriched_n, aes(x=n, y=N)) + 
  geom_col() + 
  geom_text(data=subset(annotation_filtered_enriched_n, n > 17), aes(x=n,y=N,label=geneNames), size=2,angle=90, nudge_y=0.05, hjust=0) +
  xlab("Number of circRNAs in the host gene") + ylab("Count of host genes") +
  ggtitle("Histogram of number of circRNAs per host gene") +
  scale_y_continuous(trans = mylog_trans(base=10, from=-1), breaks=c(1,10,100,1000),limits=c(0.1,1500)) +
  scale_x_continuous(breaks=pretty_breaks())
ggsave("../results/Merge_circexplorer_BC190.filtered.enriched.circRNAs_per_hostgene.hist.pdf",width = 5, height = 3)

## vs. # of exon per host genes
longestTx_annotation=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v19.annotation.longestTx.bed12", header = F, col.names = c('chr','start','end','ID','score','strand','gstart','gend','rgb','nExon','lengths','starts'), check.names = F, stringsAsFactors=F)
longestTx_annotation=separate(longestTx_annotation, ID,c('geneName','geneID',NA,NA),sep='___')
annotation_filtered_enriched$nExon=longestTx_annotation$nExon[match(annotation_filtered_enriched$geneID, longestTx_annotation$geneID)]
annotation_filtered_enriched$nExon[annotation_filtered_enriched$geneID=='CDR1as']=1
head(annotation_filtered_enriched)
# remove NA case caused by inconsistance between RefSeq and GENCODE gene symbol
annotation_filtered_enriched %>% filter(!is.na(nExon), circType=='circRNA') %>% 
  group_by(geneID) %>% summarise(nExon=mean(nExon), nCircRNA=n()) %>% ungroup() %>% 
  group_by(nExon_interval=cut(nExon, breaks = c(0,10,20,30,1000))) %>% 
  ggplot(aes(nExon_interval,nCircRNA)) + geom_boxplot() + theme_classic() 
ggsave("../results/Merge_circexplorer_BC190.filtered.enriched.nCircRNA_vs_nExon.pdf", width = 3, height = 3)

## ===============================================
## exon and flanking intron length of circRNAs vs. all exons
## ===============================================
annotation_filtered_enriched=readRDS(file="Merge_circexplorer_BC190.filtered.enriched.annotation.bed14.rds") %>% filter(circType=='circRNA') %>% dplyr::select(1:12)
introns=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/introns.meta.bed", header = F, col.names = c('chr','start','end','ID','score','strand'), check.names = F, stringsAsFactors=F) %>% mutate(score=end-start) %>% unite("chrend",c("chr",'end'), remove =F) %>% unite("chrstart",c("chr",'start'), remove =F)
head(introns)
# run makeControl.sh beforehand. 
controls=read.table("Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.matched2", header = F, col.names = c('chrom','start','end','ID','score','strand','thickStart','thickEnd','itemRgb','exonCount','exonSizes','exonOffsets','exonindex','hostgene','matchedCircRNA'), stringsAsFactors = F) %>% filter(matchedCircRNA %in% annotation_filtered_enriched$ID) %>% dplyr::select(1:12) %>% mutate(itemRgb='0,0,0')
head(controls); dim(controls)
head(annotation_filtered_enriched)

# flanking intron lenth comparsion between called circRNAs and controls
df=bind_rows(mutate(annotation_filtered_enriched,type="circRNA"),mutate(controls,type="control")) %>% 
  unite('left_intron', c("chrom", "start"), remove =F) %>% 
  unite('right_intron', c("chrom", "end"), remove =F) %>% 
  mutate(left_intron_length=introns$score[match(left_intron, introns$chrend)], 
         right_intron_length=introns$score[match(right_intron, introns$chrstart)]) %>% 
  rowwise() %>% mutate(mean_intron_length=mean(c(left_intron_length,right_intron_length), na.rm = T))  
ggplot(df,aes(x=type, y=mean_intron_length, col=type)) + geom_boxplot() + scale_y_log10() + theme_classic() + 
  ggtitle(paste("Wilcox test P value =",signif(wilcox.test(mean_intron_length ~ as.factor(type), df)$p.value,2))) +
  scale_color_manual(guide=FALSE, name="type",values=c("circRNA"="red","control" = "grey"))  
ggsave("../results/Merge_circexplorer_BC190.filtered.enriched.length.flankingIntron.pdf", width = 2,height = 3)

## number of repetitive elements comparsion between called circRNAs and controls
# bash
# slopBed -g ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/hg19.chrom.size -i Merge_circexplorer_BC197.filtered.enriched.annotation.bed14 -b 1 | bedtools intersect -a ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/introns.meta.bed -b - -wo | awk '$23==1' | cut -f1-6,10 | bedtools intersect -a - -b ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Variation/rmsk.hg19.bed -wo | cut -f7 | sort | uniq -c | sed -e 's/^[ \t]*//' > Merge_circexplorer_BC197.filtered.enriched.nRepeats_in_flankingIntron.txt
# awk '{OFS="\t"; $4=$15; print $0,$15}' Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.matched2 | slopBed -g ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/hg19.chrom.size -b 1 | bedtools intersect -a ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/introns.meta.bed -b - -wo | awk '$23==1' | cut -f1-6,10 | bedtools intersect -a - -b ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Variation/rmsk.hg19.bed -wo | cut -f7 | sort | uniq -c  | sed -e 's/^[ \t]*//' > Merge_circexplorer_BC197.filtered.enriched.matched2.nRepeats_in_flankingIntron.txt
df=bind_rows(mutate(read.table('Merge_circexplorer_BC197.filtered.enriched.nRepeats_in_flankingIntron.txt', header = F, col.names = c('nRepeats','circRNAID')),type="circRNA"),
          mutate(read.table('Merge_circexplorer_BC197.filtered.enriched.matched2.nRepeats_in_flankingIntron.txt', header = F, col.names = c('nRepeats','circRNAID')),type="control"))
dim(df); df=filter(df, circRNAID %in% as.character(annotation_filtered_enriched$ID)); dim(df)
ggplot(df,aes(x=type, y=nRepeats, col=type)) + geom_boxplot() + scale_y_log10() + theme_classic() + 
  ggtitle(paste("Wilcox test P value =",signif(wilcox.test(nRepeats ~ as.factor(type), df)$p.value,2))) +
  scale_color_manual(guide=FALSE, name="type",values=c("circRNA"="red","control" = "grey"))
ggsave("../results/Merge_circexplorer_BC190.filtered.enriched.nRepeats_in_flankingIntron.pdf", width = 2,height = 3)

# single-exon circRNAs lenth vs. random exon controls
exons=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/exons.meta.bed", sep="\t", col.names = c("chrom","start","end","ID","score","strand"), stringsAsFactors = F, header = F) %>% 
  filter(grepl("protein_coding", ID)) %>% dplyr::select(1:3) %>% sample_n(5000) %>% mutate(type="background")
df=readRDS(file="Merge_circexplorer_BC190.filtered.enriched.annotation.bed14.rds") %>% 
  filter(exonCount==1, circType=='circRNA') %>% dplyr::select(1:3) %>% mutate(type="circRNA") %>% bind_rows(exons) %>%
  mutate(exon_length=end-start)
df$type <- factor(df$type, levels = c("circRNA","background"))
ggplot(df, aes(x=type, y=exon_length, col=type)) + geom_boxplot() + scale_y_log10() + theme_classic() + 
  ggtitle(paste("Wilcox test P value =",signif(wilcox.test(exon_length ~ as.factor(type), df)$p.value,2))) +
  scale_color_manual(guide=FALSE, name="type",values=c("circRNA"="red","background" = "grey")) 
ggsave("../results/Merge_circexplorer_BC190.filtered.enriched.length.single-exon-circRNAs.pdf", width = 2,height = 3)

# region length distribution for filtered.enriched circRNAs
readRDS(file="Merge_circexplorer_BC190.filtered.enriched.annotation.bed14.rds") %>% 
  ggplot(aes(x=end-start, fill=factor(circType, levels=c('ciRNA', 'circRNA')))) + 
  geom_histogram(position = 'stack',bins = 60)  + scale_color_brewer(palette="Dark2") + 
  theme_minimal()+theme_classic()+theme(legend.position="top")+ 
  scale_x_log10() + scale_fill_manual(name="circular RNA types",values=c("circRNA"="red","ciRNA" = "orange")) 
ggsave("../results/Merge_circexplorer_BC190.filtered.enriched.lengthDistribution.pdf", width = 4,height = 3)

# lenth distribution for all called circRNAs
readRDS(file="Merge_circexplorer_BC190.annotation.bed14.rds") %>% 
  ggplot(aes(x=end-start, fill=factor(circType, levels=c('ciRNA', 'circRNA')))) + 
  geom_histogram(position = 'stack',bins = 60)  + scale_color_brewer(palette="Dark2") + 
  theme_minimal()+theme_classic()+theme(legend.position="top")+ 
  scale_x_log10() + scale_fill_manual(name="circular RNA types",values=c("circRNA"="red","ciRNA" = "orange"))
ggsave("../results/Merge_circexplorer_BC190.unfiltered.lengthDistribution.pdf", width = 4,height = 3)

## ===============================================
## host genes function enrichment analysis
## ===============================================
### tryout 4:  --- USED
## GOseq of circRNA-hosting genes vs. all expressed genes (only in brain)
## ------------------------
library(goseq) # if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("goseq")
genes_fpkm=read.table("~/neurogen/rnaseq_PD/results/merged/genes.fpkm.cufflinks.allSamples.BCv2.uniq.xls", header=T,row.names = 1, check.names = F)
genes_annotation=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed", header = F, col.names = c('chr','start','end','EnsID','score','strand','symbol','type'), check.names = F, stringsAsFactors=F)
genes_fpkm=genes_fpkm[,BC190] # only the 190 brain samples
dim(genes_fpkm)
# filter "expressed genes" with at leat 30% samples with rpkm >1 (as https://www.biorxiv.org/content/biorxiv/early/2018/12/19/500991.full.pdf)
genes_expressed=genes_fpkm[rowMeans(genes_fpkm>1)>=0.3,]  
dim(genes_expressed)
# circRNA host genes
annotation_filtered_enriched=readRDS(file="Merge_circexplorer_BC190.filtered.enriched.annotation.bed14.rds")

# clip <- pipe("pbcopy", "w")                       
# dump(unique(annotation_filtered_enriched$geneName), file=clip)                               
# close(clip)

# limit to lincRNA and protein-coding
genes_expressed_symbol=subset(genes_annotation, EnsID %in% rownames(genes_expressed) & type %in% c('protein_coding', 'lincRNA'),select=symbol, drop =T)
length(genes_expressed_symbol)

source("~/neurogen/pipeline/RNAseq/bin/lib.R")
# # using the expressed genes as background
# ORA(inputGenes=unique(annotation_filtered_enriched$geneName), allGenes=genes_expressed_symbol, topN=10, output="Merge_circexplorer_BC190.filtered.enriched.circRNA.expressedinNeuron")
# topGOenrichment(unique(annotation_filtered_enriched$geneName), allGenes=genes_expressed_symbol, topN=10, pCutoff=0.01, type='all', output="Merge_circexplorer_BC190.filtered.enriched.circRNA.expressedinNeuron")

# using all genes as background  ## used in the final version
res=ORA(inputGenes=unique(annotation_filtered_enriched$geneName), allGenes=genes_annotation$symbol, topN=10, output="../results/Merge_circexplorer_BC190.filtered.enriched.circRNA.all")
res2=ORA_slim(res)

# Q: how many (N and %)  circRNAs annotated to KEGG endosytocis genes?
table(annotation_filtered_enriched$geneName %in% strsplit(res$genelist[res$V1=="KEGG_ENDOCYTOSIS"], ",")[[1]])
# FALSE  TRUE 
# 10822   217 

# Q: how many (N and %)  circRNAs annotated to synaptic vesicle endocytosis? (XD use both GO synapse definition here and secondly the ad hoc synapse gene set also used in Fig 3a)
# table(annotation_filtered_enriched$geneName %in% strsplit(filter(res, gene_set=='c5.go.cc', V1=='GO_SYNAPSE') %>% pull(genelist), split = ",")[[1]]) # synapse genes
# # FALSE  TRUE 
# # 9512  1527 
# ora0=read.table("../results/Merge_circexplorer_BC190.filtered.enriched.circRNA.all.ORA.all.xls", sep = "\t", header = T, row.names = 1)
# filter(ora0, gene_set=='c5.go.bp', V1=='GO_SYNAPTIC_VESICLE_ENDOCYTOSIS') ## Note that GO_SYNAPTIC_VESICLE_ENDOCYTOSIS is not included in the MySigDB v7.2. We changed to check v7.4
genelists=strsplit(system("grep GOBP_SYNAPTIC_VESICLE_EXOCYTOSIS ~/bioinformatics/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/msigdb_v7.4/msigdb_v7.4_GMTs/c5.go.bp.v7.4.symbols.gmt | cut -f3-", intern =T), split = "\t")[[1]]
mean(annotation_filtered_enriched$geneName %in% genelists)
# 0.01684935

## ===============================================
## Fig. 1f: add GDA (gene-disease assocaition) from DisGeNet (http://www.disgenet.org/)
## ===============================================
# library(disgenet2r) # library(devtools); install_bitbucket("ibi_group/disgenet2r")
# genes2diseases <- function(genes=genes) {
#   disgenet2r::extract(gene2disease(gene = genes,score = c(0.1, 1),warnings = F,database = "CURATED"))
# } 
# glist = unique(annotation_filtered_enriched$geneName)
# big.list.of.data.frames <- lapply(split(glist, ceiling(seq_along(glist)/100)), genes2diseases)
# big.data.frame <- do.call(rbind,big.list.of.data.frames)

# circRNA host genes
annotation_filtered_enriched=readRDS(file="Merge_circexplorer_BC190.filtered.enriched.annotation.bed14.rds")
glist = unique(annotation_filtered_enriched$geneName)

# enrichment
#disgenet_CURATED = read.table(gzcon(url("http://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz"),text=T), sep="\t", quote = "", header = T, stringsAsFactors = F)
#saveRDS(disgenet_CURATED,file = "~/neurogen/external_download/externalData/others/DisGeNET.curated_gene_disease_associations.RDS")
disgenet_CURATED = readRDS(file = "~/neurogen/external_download/externalData/others/DisGeNET.curated_gene_disease_associations.RDS")
colnames(disgenet_CURATED); head(disgenet_CURATED)
# disgenet_disease = read.table(gzcon(url("https://www.disgenet.org/static/disgenet_ap1/files/downloads/disease_mappings_to_attributes.tsv.gz"),text=T), sep="\t", quote = "", header = T, stringsAsFactors = F)
# saveRDS(disgenet_disease,file = "~/neurogen/external_download/externalData/others/DisGeNET.disease_mapping_to_attributes.RDS")
disgenet_disease = readRDS(file = "~/neurogen/external_download/externalData/others/DisGeNET.disease_mapping_to_attributes.RDS")
colnames(disgenet_disease); head(disgenet_disease)

source("../src/annotation/tools.R")  # rewrite some of DisGeNET functions
gsd=disgenet_CURATED %>% select(geneId, geneSymbol, diseaseId, diseaseName) %>% distinct()
dim(gsd); str(gsd)
res_enrich = disease_enrichment_v2(genes = glist, gdas = gsd)
write.table(cbind(celltype='all', res_enrich), paste0("../results/Merge_circexplorer_BC190.filtered.enriched.circRNA.expressedinNeuron.DisGeNet.ALL.xls"),sep="\t", na="", row.names=F) 
select(res_enrich, "ID","Description", "pvalue", "FDR", "GeneRatio",  "BgRatio","OR", "Count") %>% filter(Count>3, FDR<0.01)
select(res_enrich, "ID","Description", "pvalue", "FDR", "GeneRatio",  "BgRatio","OR", "Count") %>% filter(grepl("Parkinson|Alzheimer", Description))
pdf("../results/Merge_circexplorer_BC190.filtered.enriched.circRNA.expressedinNeuron.DisGeNet.pdf", width = 7, height = 7); 
plot_enrichment(select(res_enrich, "Description", "Count", "FDR", "OR"), cutoff = 0.01, count = 3, limit = 200); 
dev.off()

## add PD/AD gwas genes
ADgenes=unique(scan("~/neurogen/external_download/externalData/GWAS/AD/AD_gwas_associatedGene.txt", character())); length(ADgenes)
#PDgenes=unique(system("grep -v Candidate ~/neurogen/external_download/externalData/GWAS/PD.GWAS.Chang2017.table1n2.xls | cut -f3 | sed 's/,/\\\n/g' | sort -u", intern = TRUE)); length(PDgenes) # Chang et al. 2017
PDgenes=unique(system("grep -v Nearest ~/neurogen/external_download/externalData/GWAS/PD/PD.GWAS.Nalls2019.TableS2.txt | cut -f4,5 | sed 's/\\\t$//;s/\\\t/\\\n/g' | sort -u", intern = TRUE)); length(PDgenes) # Nalls et al. 2019
annotation=readRDS(file="Merge_circexplorer_BC190.annotation.bed14.rds")
GENCODEv19=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v19.annotation.genes.bed", 
                      col.names = c("chr","start","end","geneName","geneID","geneType","strand"), stringsAsFactors = F, header = F)

gwas_list = list(AD_risk_genes=ADgenes, PD_risk_genes=PDgenes)
genes=unique(annotation$geneName)
universe=unique(GENCODEv19$geneName)

data <- data.frame(ID = character(), Description = character(), 
                   GeneRatio = character(), BgRatio = character(), OR=numeric(), geneID = character(), 
                   pvalue = numeric(), Count = integer(), stringsAsFactors = FALSE)
for(i in 1:length(gwas_list)){
  aa <- gwas_list[[i]]
  inter <- length(intersect(aa, genes))
  t <- matrix(c(inter, length(aa) - inter, length(genes) - 
                  inter, length(universe) + inter - length(aa) - length(genes)), 
              nrow = 2, dimnames = list(module = c("in", "out"), 
                                        Pheno = c("phen", "nophen")))
  test <- fisher.test(t, alternative = "greater")
  pv <- test$p.value
  OR <- as.numeric(test$estimate) ## Odds ratio
  GeneRatio <- paste(inter, "/", length(genes), sep = "")
  BgRatio <- paste(length(aa), "/", length(universe), sep = "")
  geneID <- paste(intersect(aa, genes), collapse = "/") ## can be changed to geneSymbol 
  data[i, ] <- c(as.character(names(gwas_list)[i]), as.character(names(gwas_list)[i]), GeneRatio, BgRatio, OR, geneID, pv, inter)
}
data$pvalue <- as.numeric(data$pvalue)
data$FDR <- p.adjust(data$pvalue, method = "BH")
data$Count <- as.numeric(as.character(data$Count))
data$gg <- data$Count/length(genes)
data$OR <- as.numeric(data$OR)
data <- data[order(as.numeric(as.character(data$FDR))), ]

write.table(cbind(celltype='all', rbind(res_enrich, data)), paste0("../results/Merge_circexplorer_BC190.filtered.enriched.circRNA.expressedinNeuron.DisGeNet+gwas.ALL.xls"),sep="\t", na="", row.names=F) 

pdf("../results/Merge_circexplorer_BC190.filtered.enriched.circRNA.expressedinNeuron.DisGeNet+gwas.pdf", width = 8, height = 7); 
rbind(res_enrich, data) %>% select("Description", "Count", "FDR", "OR") %>% plot_enrichment(cutoff = 0.01, count = 3, title="", limit = 200); 
dev.off()

gsd2=disgenet_CURATED %>% select(geneId, geneSymbol, diseaseId=diseaseClass) %>% distinct() %>% 
  inner_join(y=select(disgenet_disease, diseaseId=diseaseClassMSH, diseaseName=diseaseClassNameMSH) %>% distinct(), by="diseaseId")
dim(gsd2); str(gsd2)
res_enrich = disease_enrichment_v2(genes = glist, gdas = gsd2)
filter(res_enrich, Count>3, FDR<0.01) %>% select("ID","Description", "pvalue", "FDR", "GeneRatio",  "BgRatio","OR")
plot_enrichment(res_enrich, cutoff = 0.01, count = 3, limit = 10)

gsd3=disgenet_CURATED %>% select(geneId, geneSymbol, diseaseId=diseaseClass) %>% distinct() %>% 
  mutate(diseaseId = strsplit(diseaseId, ";")) %>% unnest(diseaseId) %>% distinct() %>% 
  inner_join(y=select(disgenet_disease, diseaseId=diseaseClassMSH, diseaseName=diseaseClassNameMSH) %>% 
               mutate(diseaseId = strsplit(diseaseId, ";"), diseaseName = strsplit(diseaseName, "; ")) %>% 
               unnest(c(diseaseId, diseaseName)) %>% distinct(), 
             by="diseaseId")
dim(gsd3); head(gsd3)
res_enrich = disease_enrichment_v2(genes = glist, gdas = gsd3) 
filter(res_enrich, Count>=3, FDR<=0.05) %>% select("ID","Description", "pvalue", "FDR", "GeneRatio",  "BgRatio","OR")
pdf("../results/Merge_circexplorer_BC190.filtered.enriched.circRNA.expressedinNeuron.DisGeNet.class.pdf", width = 7, height = 2); 
plot_enrichment(res_enrich, cutoff = 0.01, count = 3, limit = 20); 
dev.off()


## how many ILB-DE circRNAs are from synGO genes?
de=read.table("../results/DE2gene_SNDA/DEresult.DE2gene_SNDA.CONDITION_ILB_vs_HC.xls.gz", sep="\t", header = T, row.names = 1)
head(de)
synGO=unique(readxl::read_xlsx("syngo_genes.xlsx")$hgnc_symbol); length(synGO)  # N=1112; download from SynGO (https://www.syngoportal.org)
filter(de[,1:12], pvalue<0.05, abs(log2FoldChange)>1, geneName %in% synGO)


## 99% of the AD-associated GWAS SNPs26 (2334/2357) and XX% (…) of the PD-associated GWAS SNPs were located in the proximity (i.e. within 1M bp) of a circRNA (data not shown)
# bash
# zcat ~/neurogen/external_download/externalData/GWAS/AD/AD_sumstats_Jansenetal.txt.gz | awk 'NR>1 && $8<5e-8{OFS="\t"; print "chr"$2, ($3<1000000)?0:($3-1000000), $3+1000000}' | intersectBed -a - -b ~/projects/circRNA/data/Merge_circexplorer_BC197.filtered.enriched.annotation.bed14 -u | wc -l
# cat ~/neurogen/external_download/externalData/GWAS/PD/nallsEtAl2019_allSamples_significantVariants5E-8.hg19.bed6 | awk '{OFS="\t"; print $1, ($3<1000000)?0:($3-1000000), $3+1000000}' | intersectBed -a - -b ~/projects/circRNA/data/Merge_circexplorer_BC197.filtered.enriched.annotation.bed14 -u | wc -l

## Q: XX% of all synaptic circRNAs were linked to brain diseases
annotation_filtered_enriched=readRDS(file="Merge_circexplorer_BC190.filtered.enriched.annotation.bed14.rds"); dim(annotation_filtered_enriched)
synGO=unique(readxl::read_xlsx("syngo_genes.xlsx")$hgnc_symbol); length(synGO)  # N=1112; download from SynGO (https://www.syngoportal.org)
disgenet_CURATED = readRDS(file = "~/neurogen/external_download/externalData/others/DisGeNET.curated_gene_disease_associations.RDS")
colnames(disgenet_CURATED); head(disgenet_CURATED)
disgenet_disease = readRDS(file = "~/neurogen/external_download/externalData/others/DisGeNET.disease_mapping_to_attributes.RDS")
gsd3=disgenet_CURATED %>% select(geneId, geneSymbol, diseaseId=diseaseClass) %>% distinct() %>% 
  mutate(diseaseId = strsplit(diseaseId, ";")) %>% unnest(diseaseId) %>% distinct() %>% 
  inner_join(y=select(disgenet_disease, diseaseId=diseaseClassMSH, diseaseName=diseaseClassNameMSH) %>% 
               mutate(diseaseId = strsplit(diseaseId, ";"), diseaseName = strsplit(diseaseName, "; ")) %>% 
               unnest(c(diseaseId, diseaseName)) %>% distinct(), 
             by="diseaseId")
dim(gsd3); head(gsd3)
filter(annotation_filtered_enriched, geneName %in% synGO) %>% dim() # 1362
filter(annotation_filtered_enriched, geneName %in% synGO, geneName %in% unique(filter(gsd3, diseaseName %in% c("Mental Disorders", "Nervous System Diseases")) %>% pull(geneSymbol))) %>% dim() # 837
837/1362
# 0.6145374

## ===============================================
## expression of linear RNA vs. circRNA
## ===============================================
library(tidyverse)
library(scales)
setwd("~/projects/circRNA/data/backup_circExplorer2_v20190819/") 
nS=read.table("Merge_circexplorer_BC.annotation.bed14.s3s5.gz", col.names=c("ID","type","sample","reads","count"), check.names = F, stringsAsFactors=F)
nU=read.table("Merge_circexplorer_BC.annotation.bed14.u3u5.gz", col.names=c("ID","type","sample","reads","count"), check.names = F, stringsAsFactors=F)
nC=read.table("Merge_circexplorer_BC.annotation.bed14.circReads.txt", header=T, check.names = F, row.names = 1, stringsAsFactors=F)
# nL0=read.table("Merge_circexplorer_BC.annotation.bed14.sum_u3u5s3s5", header=T, check.names = F, row.names = 1, stringsAsFactors=F) # from bedtools coverage (which might be buggy)
sample106=scan("~/neurogen/circRNA_seq_Rebeca_HC/BRAINCODE_circexp/BC.n106.samplelist",character())
filtered_enriched_annotation=readRDS("Merge_circexplorer_BC106.filtered.enriched.annotation.bed14.rds")$ID; length(filtered_enriched_annotation)

nLw=unite(rbind(nS,nU), temp, sample, type, sep = "__") %>% select(ID, temp, reads) %>% spread(key=temp, value=reads, fill=0) %>% column_to_rownames(var = "ID")
nS3=select(nLw, contains("__s3")) %>% rename_all(.funs = funs(sub("__.*", "", .))); dim(nS3)
nS5=select(nLw, contains("__s5")) %>% rename_all(.funs = funs(sub("__.*", "", .))); dim(nS5)
nU3=select(nLw, contains("__u3")) %>% rename_all(.funs = funs(sub("__.*", "", .))); dim(nU3)
nU5=select(nLw, contains("__u5")) %>% rename_all(.funs = funs(sub("__.*", "", .))); dim(nU5)
# nS3=filter(nS,type=='s3') %>% select(ID, sample, reads) %>% spread(key=sample, value=reads, fill = 0) %>% column_to_rownames(var = "ID")
# nS5=filter(nS,type=='s5') %>% select(ID, sample, reads) %>% spread(key=sample, value=reads, fill = 0) %>% column_to_rownames(var = "ID")
# nU3=filter(nU,type=='u3') %>% select(ID, sample, reads) %>% spread(key=sample, value=reads, fill = 0) %>% column_to_rownames(var = "ID")
# nU5=filter(nU,type=='u5') %>% select(ID, sample, reads) %>% spread(key=sample, value=reads, fill = 0) %>% column_to_rownames(var = "ID")
common_Rows=Reduce(intersect, list(filtered_enriched_annotation, nS$ID, nU$ID, rownames(nC))); length(common_Rows)
common_Cols=Reduce(intersect, list(colnames(nC), sample106, nS$sample, nU$sample)); length(common_Cols)
nC=nC[common_Rows, sample106]; nS3=nS3[common_Rows, sample106]; nS5=nS5[common_Rows, sample106]; nU3=nU3[common_Rows, sample106]; nU5=nU5[common_Rows, sample106]
save(nC,nS3,nS5,nU3,nU5, file = "Merge_circexplorer_BC.annotation.bed14.nS3nS5nU3nU5.RData")
load("Merge_circexplorer_BC.annotation.bed14.nS3nS5nU3nU5.RData")
nL=nS3+nS5+nU3+nU5  # not very fair to use the sum
## 3/21/2019: We could change to use max: nL=max(nS3,nS5,nU3,nU5)
# combine
nCL=inner_join(x=rownames_to_column(nC, var = 'ID') %>% gather(key="sample",value=nC, contains("_")),
               y=rownames_to_column(nL, var = 'ID') %>% gather(key="sample",value=nL, contains("_")),
               by=c("ID","sample")) %>% 
  filter(nC>0, nL>0) %>%
  mutate(celltype=gsub(".*_.*_(.*)_.*_.*","\\1",sample))

nCLrange=range(c(nCL$nL,nCL$nC))

# ratio boxplot
head(nCL); dim(nCL)
group_by(nCL, celltype) %>% summarise(rho=cor(nC,nL,method = 'spearman'))
# celltype    rho
# 1 FB       0.0654
# 2 MCPY     0.114 
# 3 PBMC     0.104 
# 4 SNDA     0.113 
# 5 TCPY     0.0895
with(nCL, cor(nC, nL, method = 'spearman')) # 0.1028416

# rCL difference between neuron vs. non-neuron
wilcox.test(rCL ~ is_neuron, data=mutate(nCL, rCL=nC/(nC+nL), is_neuron=as.factor(ifelse(celltype %in% c('SNDA','TCPY','MCPY'), 'yes', 'no'))))
# Wilcoxon rank sum test with continuity correction
# data:  rCL by is_neuron
# W = 38260758, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

mutate(nCL, rCL=nC/(nC+nL)) %>%
  ggplot(aes(factor(celltype, levels = c("SNDA","TCPY","MCPY","FB","PBMC")),rCL)) + 
  scale_y_log10(breaks = c(0,1e-5,1e-4,1e-3,0.01,.1,.5,1),labels = c(0,1e-5,1e-4,1e-3,0.01,.1,.5,1)) +
  geom_violin(aes(fill = factor(celltype, levels = c("SNDA","TCPY","MCPY","FB","PBMC"))), scale = "width", trim = FALSE) + 
  geom_boxplot(fill='white', width=0.1, outlier.shape = NA) + 
  scale_fill_manual(name="Cell types",
                     values=c("SNDA"="#F22A7B","TCPY" = "#3182bd","MCPY" = "#2659B2","FB" ="#BC9371","PBMC"="#D3CBCB")) +
  labs(x="Cell types", y="Circular-to-linear fraction (ratio = nC/total reads)") +
  theme_bw() + 
  theme(legend.position='none')
ggsave("../results/Merge_circexplorer_BC106.annotation.bed14.nC.vs.nL.boxplot.pdf", width=4.5, height = 4.5)


# nC vs. nL scatter plot
set.seed(1)
library(scales)
ggplot(nCL,aes(x = jitter(nL, amount=0.49), y = jitter(nC, amount=.49))) +
  geom_point(aes(colour = factor(celltype, levels = c("SNDA","TCPY","MCPY","FB","PBMC"))), shape=16, alpha=0.8) +
  geom_abline(intercept = 0, slope = 1, linetype=2) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(name="Cell types",
                     values=c("SNDA"="#F22A7B","TCPY" = "#3182bd","MCPY" = "#2659B2","FB" ="#BC9371","PBMC"="#D3CBCB")) +
  coord_fixed() +
  expand_limits(x=nCLrange,y=nCLrange) +
  labs(x="Number of linear reads at junction sites", y="Number of circular reads at junction sites") +
  theme_bw()
ggsave("../results/Merge_circexplorer_BC106.annotation.bed14.nC.vs.nL.pdf", width=6, height = 5, useDingbats=T)
ggsave("../results/Merge_circexplorer_BC106.annotation.bed14.nC.vs.nL.png", width=6, height = 5)

pdf("../results/Merge_circexplorer_BC106.annotation.bed14.nC.vs.nL.sub.pdf", width=6, height = 5)
p=ggplot(nCL,aes(x = nL, y = nC)) +
  geom_point(position=position_jitter(width = 0.1, height = 0.1), color='#aaaaaa', shape=19, size=1, alpha=0.8) + 
  geom_abline(intercept = 0, slope = 1, linetype=2) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  coord_fixed() +
  expand_limits(x=nCLrange,y=nCLrange) + 
  labs(x="Number of linear reads at junction sites", y="Number of circular reads at junction sites") +
  theme_bw()
print(p)
p1=ggplot(filter(nCL,celltype=="SNDA"),aes(x = nL, y = nC)) +
  geom_point(position=position_jitter(width = 0.1, height = 0.1),aes(colour = factor(celltype, levels = c("SNDA"))), shape=19, size=1, alpha=0.8) + 
  geom_abline(intercept = 0, slope = 1, linetype=2) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(name="Cell types",
                     values=c("SNDA" = "#F22A7B")) +
  coord_fixed() +
  expand_limits(x=nCLrange,y=nCLrange) + 
  labs(x="Number of linear reads at junction sites", y="Number of circular reads at junction sites") +
  theme_bw()
print(p1)
p2=ggplot(filter(nCL,celltype %in% c("TCPY","MCPY")),aes(x = nL, y = nC)) +
  geom_point(position=position_jitter(width = 0.1, height = 0.1),aes(colour = factor(celltype, levels = c("TCPY","MCPY"))), shape=19, size=1, alpha=0.8) + 
  geom_abline(intercept = 0, slope = 1, linetype=2) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(name="Cell types",
                     values=c("TCPY" = "#3182bd","MCPY" = "#2659B2")) +
  coord_fixed() +
  expand_limits(x=nCLrange,y=nCLrange) + 
  labs(x="Number of linear reads at junction sites", y="Number of circular reads at junction sites") +
  theme_bw()
print(p2)
p3=ggplot(filter(nCL,celltype %in% c("FB","PBMC")),aes(x = nL, y = nC)) +
  geom_point(position=position_jitter(width = 0.1, height = 0.1),aes(colour = factor(celltype, levels = c("FB","PBMC"))), shape=19, size=1, alpha=0.8) + 
  geom_abline(intercept = 0, slope = 1, linetype=2) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(name="Cell types",
                     values=c("FB" ="#BC9371","PBMC"="#D3CBCB")) +
  coord_fixed() +
  expand_limits(x=nCLrange,y=nCLrange) +
  labs(x="Number of linear reads at junction sites", y="Number of circular reads at junction sites") +
  theme_bw()
print(p3)
dev.off()

## Q: which circRNAs have significantly more circular reads than linear reads?
# use Wilcoxon Signed-Rank Test
length(grep("_SNDA_",names(nS3)))
readsNum_filtered<- read.table("~/neurogen/rnaseq_PD/run_output/linescounts.filtered.txt",row.names=1, header=F, check.names = F) 
readsNum_million<-(t(readsNum_filtered)[1,]/10^6)
# convert to RPM
SNDAsamples = names(nS3)[grep("_SNDA_",names(nS3))]
normalized_factor = readsNum_million[SNDAsamples]
rpmL<-sweep(nL[,SNDAsamples],2,normalized_factor,"/"); dim(rpmL)
rpmC<-sweep(nC[,SNDAsamples],2,normalized_factor,"/"); dim(rpmC)
# filter: at least 10% samples with RPM >0 in both circular and linear reads
expressed_rows = rowMeans(rpmL>0)>.01 & rowMeans(rpmC>0)>.01; sum(expressed_rows)
do.call(rbind, apply(rpmC[expressed_rows,] - rpmL[expressed_rows, ], 1, wilcox.test, alternative = "g", exact = F, correct=F)) %>% data.frame() %>% rownames_to_column() %>% filter(p.value<0.05)


# ## into 2D density (e.g. https://www.r-graph-gallery.com/2d-density-plot-with-ggplot2/)
# ggplot(nCL, aes(x=nL, y=nC) ) +
#   stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
#   scale_fill_distiller(palette="Spectral", direction=-1) +
#   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),expand = c(0, 0),
#                 labels = trans_format("log10", math_format(10^.x))) +
#   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),expand = c(0, 0),
#                 labels = trans_format("log10", math_format(10^.x))) +
#   coord_fixed() +
#   theme(legend.position='none')
