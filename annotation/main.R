## Main script to work on the output of circExplorer
library('dplyr')
library('reshape2')
setwd("~/projects/circRNA/data/") 

## See ~/projects/circRNA/README.Rmd for how to generate the input files: 
# Merge_circexplorer_BC.annotation.bed14
# Merge_circexplorer_BC.rawcount.txt

###########################################
################# load data  ##############
###########################################

## === annotation ===

annotation<- read.table("Merge_circexplorer_BC.annotation.bed14",header=T, check.names = F) # include both ciRNA and circRNA
head(annotation)

## === load raw reads count ===

#removed "#" on first line as well as ful name to just sample name
Merge_circexp_raw <- read.table("Merge_circexplorer_BC.rawcount.txt",sep="\t",check.names =F, header=TRUE, comment.char ='') 
colnames(Merge_circexp_raw)=gsub("#|_candidates.bed_circReads","",colnames(Merge_circexp_raw))

Merge_circexp_raw = mutate(Merge_circexp_raw, concated_column = paste(chrom,start,end, sep = '_'))
row.names(Merge_circexp_raw) = Merge_circexp_raw$concated_column
Merge_circexp_raw=select(Merge_circexp_raw,dplyr::contains('rep'))
dim(Merge_circexp_raw)
#[1] 189001    106
#[1] 417175    125

readsNum_filtered<- read.table("~/neurogen/rnaseq_PD/run_output/linescounts.filtered.txt",row.names=1, header=F, check.names = F) 
readsNum_million<-(t(readsNum_filtered)[1,]/10^6)
#only include the 125 (106 HCILB + 19 PD) samples
readsNum_million=readsNum_million[colnames(Merge_circexp_raw)]


###########################################
############## normalize to RPM ###########
###########################################

Merge_circexp_norm<-sweep(Merge_circexp_raw,2,readsNum_million,"/")

###########################################
############ filter cirRNAs     ###########
###########################################

# What can be a good filter?

### (1) being expressed
### --------------------

# Definition #1 of being expression: with RPM >= threshold in >=2 samples
pdf("~/projects/circRNA/results/filter_QC.pdf", width = 6, height = 5)
boxplot(1/readsNum_million, outcol="NA", ylab="1 RPM cutoff for each sample"); 
points(jitter(rep(1, length(readsNum_million))), 1/readsNum_million, pch=20, col=rgb(0,0,0,.6)) 
median(1/readsNum_million) # 0.006118284  --> 0.006
dev.off()
threshold<- signif(median(1/readsNum_million), 1)  # 0.006 
N = ncol(Merge_circexp_norm)
# with cutoff of 1/(min(readsNum_million)), 4030 remained; while RPM>=0.006, 9161 remained.
Merge_circexp_norm_filtered <- Merge_circexp_norm[rowSums(Merge_circexp_norm >= threshold) >= 2, ]
dim(Merge_circexp_norm_filtered)
# [1] 10279  125

# Definition #2 of being expression: with >=2 unique back-splicing reads in all samples (Zheng et al., doi:10.1038/ncomms11215)
Merge_circexp_raw_filtered <- Merge_circexp_raw[rowSums(Merge_circexp_raw)>=2, ]
dim(Merge_circexp_raw_filtered)
# [1] 65119   125

### (2) being enriched
### --------------------

# Definition of being enriched: at least 20 reads in RNase R treatment AND foldchange(RNase vs Mock) at least 2 in at least one sample. See main.RM.R first
Merge_circexp_raw_long_filterd.RM  = readRDS(file="Merge_circexp_raw_long_filterd.RM.rds")
length(unique(Merge_circexp_raw_long_filterd.RM$name))
# [1] 47994

### Final selection: being expressed (#2) and being enriched

Merge_circexp_norm_filtered=Merge_circexp_norm[intersect(rownames(Merge_circexp_raw_filtered), Merge_circexp_raw_long_filterd.RM$name), ]
dim(Merge_circexp_norm_filtered)
# [1] 10431  125

total_circRNA_raw_reads = rowSums(Merge_circexp_raw)
pdf("~/projects/circRNA/results/total_circRNA_raw_reads.pdf", width = 6, height = 5)
layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(5,1), heights=c(1,1))
ht=hist(total_circRNA_raw_reads[total_circRNA_raw_reads<=200], breaks = 0:200, plot=F)
par(mar=c(4,4,0,1))
plot(log10(ht$counts+.1)+1, lwd=2,xlim=c(0,200), type='h', ylim = c(0.9,6.6),col=c("#000000",rep('#aaaaaa',199)), ylab='Number of circular RNAs', xlab="Number of backspliced reads",yaxt="n",xaxs="i", yaxs="i",bty="n")
y1=floor(range(log10(ht$counts+.1)+1))
pow <- seq(y1[1], y1[2])
ticksat <- as.vector(sapply(pow-1, function(p) log10((2:9)*10^p)))
axis(2, pow, labels = 10^(pow-1), las=1)
axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
# add the remained ones with the #2 filter
filtered_circRNA_raw_reads = rowSums(Merge_circexp_raw[rownames(Merge_circexp_norm_filtered),])
ht=hist(filtered_circRNA_raw_reads[filtered_circRNA_raw_reads<=200], breaks = 0:200, plot=F)
points(log10(ht$counts+.1)+1, lwd=2, type='h',col='#ff0000')
legend("topleft",c(paste("Distinct circular RNAs (n =",format(length(total_circRNA_raw_reads),big.mark=","),")"),
                   paste("+ at least 2 reads in overall samples (n =", format(sum(total_circRNA_raw_reads>1),big.mark=","),")"),
                   paste("+ enriched in RNase R (n =", format(length(filtered_circRNA_raw_reads),big.mark=","),")")),
       col=c("#000000",'#aaaaaa','#ff0000'),
       text.col = c("#000000",'#aaaaaa','#ff0000'),
       bty='n')

## for the 200- region
par(mar=c(4,1,0,1))
ht=hist(total_circRNA_raw_reads[total_circRNA_raw_reads>200], breaks = 200:8200, plot=F)
plot(log10(ht$counts+.1)+1, lwd=2,xlim=c(0,8200), type='h', ylim = c(0.9,6.6),col='#aaaaaa', ylab="",xlab="",xaxt="n",yaxt='n',xaxs="i", yaxs="i",bty="n")
pow <- seq(y1[1], y1[2])
ticksat <- as.vector(sapply(pow-1, function(p) log10((2:9)*10^p)))
#axis(2, pow, labels = 10^(pow-1), las=1)
#axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
axis(1, c(0,seq(2000,8200,2000)-200), labels=c(200,seq(2000,8200,2000)))
# add the remained ones with the #2 filter
filtered_circRNA_raw_reads = rowSums(Merge_circexp_raw[rownames(Merge_circexp_norm_filtered),])
ht=hist(filtered_circRNA_raw_reads[filtered_circRNA_raw_reads>200], breaks = 200:8200, plot=F)
points(log10(ht$counts+.1)+1, lwd=2,xlim=c(0,8200), type='h',col='#ff0000')

dev.off()

# save
save(threshold, Merge_circexp_raw,Merge_circexp_norm, Merge_circexp_raw_filtered, Merge_circexp_norm_filtered, file="Merge_circexp.BC.Rdata")

write.table(annotation[annotation$ID %in% rownames(Merge_circexp_norm_filtered),], file="Merge_circexplorer_BC.annotation.Merge_circexp_norm_filtered.bed14", sep="\t", quote = F, row.names = F)

###########################################
## circRNA expressed in each cell #########
###########################################
#definition: “Expressed” circular RNAs: RPM > CUTOFF in at least 2 samples overall (see red above) and in at least 1 sample in a specific cell group 
load("Merge_circexp.BC.Rdata")
df = Merge_circexp_norm_filtered
df = df %>% mutate(gene=rownames(df)) %>% melt(id.vars="gene", variable.name = "sampleID", value.name='fpkm') %>%
  separate(sampleID, c("Dx", "subjectID","cellType5","batch","rep"), sep = "_", remove=FALSE) %>% 
  mutate(cellType3=ifelse(cellType5 %in% c("TCPY","MCPY"), "PY", ifelse(cellType5=="SNDA","SNDA","NN"))) %>% 
  filter(fpkm>threshold) %>% 
  select(gene, cellType3, cellType5) %>% arrange(gene, cellType3, cellType5)
df %>% group_by(gene) %>% summarize(celltype3 = paste(unique(cellType3), collapse = ','), celltype5 = paste(unique(cellType5), collapse = ',')) %>% ungroup() %>% merge(annotation, by.x='gene',by.y='ID', all.x=T, all.y=F) %>% write.table(file="../results/Merge_circexp_norm_filtered.celltype.annotation.xls", sep="\t", quote=F, row.names=F)

df %>% group_by(gene) %>% summarize(celltype3 = paste(unique(cellType3), collapse = ','), celltype5 = paste(unique(cellType5), collapse = ',')) %>% group_by(celltype3) %>% summarise(n=n())
# # A tibble: 7 x 2
# celltype3     n
# <chr> <int>
# 1         NN  1040
# 2      NN,PY   542
# 3 NN,PY,SNDA  1152
# 4    NN,SNDA  1664
# 5         PY   271
# 6    PY,SNDA  1843
# 7       SNDA  2649

# 1         NN  3021
# 2      NN,PY   352
# 3 NN,PY,SNDA  1071
# 4    NN,SNDA  1367
# 5         PY  1029
# 6    PY,SNDA   806
# 7       SNDA  2607

clip <- pipe("pbcopy", "w");
df %>% group_by(gene) %>% summarize(celltype3 = paste(unique(cellType3), collapse = ','), celltype5 = paste(unique(cellType5), collapse = ',')) %>% group_by(celltype5) %>% summarise(n=n()) %>% as.data.frame() %>% write.table(file=clip, sep="\t", row.names=F)
close(clip)

df %>% merge(annotation, by.x='gene',by.y='ID', all.x=T, all.y=F) %>% select(gene, cellType3, cellType5, circType) %>% group_by(cellType5, circType) %>% summarise(n=n_distinct(gene)) %>% dcast(cellType5 ~ circType, value.var="n") #%>% write.table(file=clip, sep="\t", row.names=F)
# cellType5 circRNA ciRNA
# 1        FB    1476   298
# 2      MCPY    1074   385
# 3      PBMC    3421   409
# 4      SNDA    4540  2768
# 5      TCPY    2292   727
df %>% merge(annotation, by.x='gene',by.y='ID', all.x=T, all.y=F) %>% select(gene, cellType3, cellType5, circType) %>% group_by(cellType3, circType) %>% summarise(n=n_distinct(gene)) %>% dcast(cellType3 ~ circType, value.var="n") #%>% write.table(file=clip, sep="\t", row.names=F)
# cellType3 circRNA ciRNA
# 1        NN    3792   606
# 2        PY    2777  1031
# 3      SNDA    4540  2768

## call eulerAPE to generate Merge_circexp_norm_filtered.celltype3.pdf

# 5-ways venn diagram (not area-proportional)
library(gplots)
pdf("../results/Merge_circexp_norm_filtered.celltype5.pdf")
venn(list(SNDA = unique(df$gene[df$cellType5=='SNDA']),
          TCPY = unique(df$gene[df$cellType5=='TCPY']),
          MCPY = unique(df$gene[df$cellType5=='MCPY']),
          PBMC = unique(df$gene[df$cellType5=='PBMC']),
          FB   = unique(df$gene[df$cellType5=='FB'])
          ))
dev.off()


###########################################
############ cell specificity   ###########
###########################################
source("../src/getTissueSpecificRNA.R")
# output: ../results/cellspecific_heatmap.pdf, and ../results/cellspecific_heatmap+.pdf

######################################################
####### validated in RNase R treated RNAseq   ########
######################################################



######################################################
#######     QTL of  circRNA in SNDA          ########
######################################################
# for the filterd circRNAs in SNDA
load("Merge_circexp.BC.Rdata")
df = Merge_circexp_norm_filtered

# get circRNA annotatio
library(readr)
annotation %>% select(ID, chrom, start, end) %>% rename(geneid=`ID`, chr=`chrom`, s1=`start`, s2=`end`) %>% write_tsv(path="../results/Merge_circexplorer_BC.annotation.loci.txt")

######################################################
####get circula ratio (CR) of circRNAs in SNDA  ######
######################################################

######################################################
#######   crQTL of filtered circRNA in SNDA   ########
######################################################


















#head(Merge_circexp_norm_filtered[order(-Merge_circexp_norm_filtered$columnGreat),][ 85:91])
head(Merge_circexp_norm_filtered[order(-Merge_circexp_norm_filtered$columnGreat),][ (samplenum_BC+1):(samplenum_BC+7)])

#start       end        name score strand columnGreat
#134333  40655612  40657444 circ_134333 11647      -          64
#117125  20107645  20109225 circ_117125  4208      +          53
#16800  117944807 117963271  circ_16800  7360      +          48
#95547   65941524  65972074  circ_95547  2711      +          44
#97736    9182379   9221997  circ_97736  7390      +          41
#11802   27056141  27059283  circ_11802  1910      +          33

circ_count<-colSums(Merge_circexp_norm_filtered[1:samplenum_BC]>threshold) #add how many samples have values greater than 0, so at least one circRNA

circ_count_BC<-circ_count
op <- par(mar=c(12,6,4,2))
barplot(circ_count,width=0.8,las=2,cex.names=0.7,
        cex.axis=0.8, ylab="Number of circRNAs or ciRNAs", main="circRNAs and ciRNAs count Braincode with CircExplorer Normalized and Filtered")
rm(op)

#Subset by cell type
SNDA_samples<-(Merge_circexp_norm_filtered[grepl( "SNDA|chrom|start|end|name|score|strand" , names(Merge_circexp_norm_filtered))])
TCPY_samples<-(Merge_circexp_norm_filtered[grepl( "TCPY|chrom|start|end|name|score|strand" , names(Merge_circexp_norm_filtered))])
PBMC_samples<-(Merge_circexp_norm_filtered[grepl( "PBMC|chrom|start|end|name|score|strand" , names(Merge_circexp_norm_filtered))])
MCPY_samples<-(Merge_circexp_norm_filtered[grepl( "MCPY|chrom|start|end|name|score|strand" , names(Merge_circexp_norm_filtered))])
FB_samples<-(Merge_circexp_norm_filtered[grepl( "FB|chrom|start|end|name|score|strand" , names(Merge_circexp_norm_filtered))])
neuronal_samples<-(Merge_circexp_norm_filtered[grepl( "SNDA|TCPY|MCPY|chrom|start|end|name|score|strand" , names(Merge_circexp_norm_filtered))])
NonNeuronal_samples<-(Merge_circexp_norm_filtered[grepl( "PBMC|FB|chrom|start|end|name|score|strand" , names(Merge_circexp_norm_filtered))])

#count sample number and make variable
num_SNDA= 62 
num_TCPY= 9 
num_PBMC=4 
num_MCPY =3 
num_FB = 3 
num_neuronal=74
num_nonNeuronal=7

#How many in each cell type 
#At least 1 sample greater than threshold

SNDA_samples$columnGreat<-rowSums(SNDA_samples[1:num_SNDA]>threshold) #add how many rows have cols greater than threshold
SNDA_samples_filtered<-subset(SNDA_samples, columnGreat>0) #keep those rows with 1 sample or more greater than threshold
nrow(SNDA_samples_filtered)
#[1] 2426 # how many have at least one greater than threshold. 

TCPY_samples$columnGreat<-rowSums(TCPY_samples[1:num_TCPY]>threshold) #add how many rows have cols greater than threshold
TCPY_samples_filtered<-subset(TCPY_samples, columnGreat>0) #keep those rows with 1 sample or more greater than threshold
nrow(TCPY_samples_filtered)
#[1] 1377

PBMC_samples$columnGreat<-rowSums(PBMC_samples[1:num_PBMC]>threshold) #add how many rows have cols greater than threshold
PBMC_samples_filtered<-subset(PBMC_samples, columnGreat>0)#keep those rows with 1 sample or more greater than threshold
nrow(PBMC_samples_filtered)
#[1]  2236

MCPY_samples$columnGreat<-rowSums(MCPY_samples[1:num_MCPY]>threshold) #add how many rows have cols greater than threshold
MCPY_samples_filtered<-subset(MCPY_samples, columnGreat>0) #keep those rows with 1 sample or more greater than threshold
nrow(MCPY_samples_filtered)
#[1] 708

FB_samples$columnGreat<-rowSums(FB_samples[1:num_FB]>threshold) #add how many rows have cols greater than threshold
FB_samples_filtered<-subset(FB_samples, columnGreat>0)#keep those rows with 1 sample or more greater than threshold
nrow(FB_samples_filtered)
#[1] 967

neuronal_samples$columnGreat<-rowSums(neuronal_samples[1:num_neuronal]>threshold) #add how many rows have cols greater than threshold
neuronal_samples_filtered<-subset(neuronal_samples, columnGreat>0) #keep those rows with 1 sample or more greater than threshold
nrow(neuronal_samples_filtered)
#[1] 2917 # how many have at least one greater than threshold. 

NonNeuronal_samples$columnGreat<-rowSums(NonNeuronal_samples[1:num_nonNeuronal]>threshold) #add how many rows have cols greater than threshold
NonNeuronal_samples_filtered<-subset(NonNeuronal_samples, columnGreat>0) #keep those rows with 1 sample or more greater than threshold
nrow(NonNeuronal_samples_filtered)
#[1] 2528
# how many have at least one greater than threshold. 


#find unique to these subsets
#It'll be easier to just use the location information to view gene replicates

SNDA_samples_filtered_pos<-subset(SNDA_samples_filtered, select = c(chrom, start,end,name,strand))
TCPY_samples_filtered_pos<-subset(TCPY_samples_filtered, select = c(chrom, start,end,name,strand))
PBMC_samples_filtered_pos<-subset(PBMC_samples_filtered, select = c(chrom, start,end,name,strand))
MCPY_samples_filtered_pos<-subset(MCPY_samples_filtered, select = c(chrom, start,end,name,strand))
FB_samples_filtered_pos<-subset(FB_samples_filtered, select = c(chrom, start,end,name,strand))
neuronal_samples_filtered_pos<-subset(neuronal_samples_filtered, select = c(chrom, start,end,name,strand))
NonNeuronal_samples_filtered_pos<-subset(NonNeuronal_samples_filtered, select = c(chrom, start,end,name,strand))


SNDA_samples_filtered_pos$Coder <- "S"
TCPY_samples_filtered_pos$Coder <- "T"
PBMC_samples_filtered_pos$Coder <- "P"
MCPY_samples_filtered_pos$Coder<- "M"
FB_samples_filtered_pos$Coder<- "F"


df_unique <- rbind(SNDA_samples_filtered_pos, TCPY_samples_filtered_pos, PBMC_samples_filtered_pos,MCPY_samples_filtered_pos,FB_samples_filtered_pos)    


#check size
#nrow(df_unique)
#7714
#2426 +1377+2236+708+967
#7714

# Find the rows which have duplicates in a different group.
dupRows <- dupsBetweenGroups(df_unique, "Coder")
# Print it alongside the data frame
cbind(df_unique, dup=dupRows)
#find unique rows
cbind(df_unique, unique=!dupRows)

# Store the results in df
dfUniq <- cbind(df_unique, unique=!dupRows)

#split with results
SNDA_samples_filtered_pos <- subset(dfUniq, Coder=="S", select=-Coder)
TCPY_samples_filtered_pos<- subset(dfUniq, Coder=="T", select=-Coder)
PBMC_samples_filtered_pos<- subset(dfUniq, Coder=="P", select=-Coder)
MCPY_samples_filtered_pos<- subset(dfUniq, Coder=="M", select=-Coder)
FB_samples_filtered_pos<- subset(dfUniq, Coder=="F", select=-Coder)

#count unique circRNAs per group

sum(SNDA_samples_filtered_pos$unique)
#[1] 394
sum(TCPY_samples_filtered_pos$unique)
#[1] 59
sum(PBMC_samples_filtered_pos$unique)
#[1] 335
sum(MCPY_samples_filtered_pos$unique)
#[1] 3
sum(FB_samples_filtered_pos$unique)
#[1] 100

#count neuronal specific

neuronal_samples_filtered_pos$Coder<- "N"
NonNeuronal_samples_filtered_pos$Coder<- "O"

df_unique_neuro <- rbind(neuronal_samples_filtered_pos,NonNeuronal_samples_filtered_pos)    

# Find the rows which have duplicates in a different group.
dupRows_neuro <- dupsBetweenGroups(df_unique_neuro, "Coder")
# Print it alongside the data frame
cbind(df_unique_neuro, dup=dupRows_neuro)
#find unique rows
cbind(df_unique_neuro, unique=!dupRows_neuro)

# Store the results in df
dfUniq_neuro <- cbind(df_unique_neuro, unique=!dupRows_neuro)

#split with results
neuronal_samples_filtered_pos<- subset(dfUniq_neuro, Coder=="N", select=-Coder)
NonNeuronal_samples_filtered_pos<- subset(dfUniq_neuro, Coder=="O", select=-Coder)

#count unique circRNAs per group

sum(neuronal_samples_filtered_pos$unique)
#1058
sum(NonNeuronal_samples_filtered_pos$unique)
#669

#subset
SNDA_unique_circRNAs<-SNDA_samples_filtered_pos[ which(SNDA_samples_filtered_pos$unique=='TRUE'), ]
TCPY_unique_circRNAs<-TCPY_samples_filtered_pos[ which(TCPY_samples_filtered_pos$unique=='TRUE'), ]
PBMC_unique_circRNAs<-PBMC_samples_filtered_pos[ which(PBMC_samples_filtered_pos$unique=='TRUE'), ]
MCPY_unique_circRNAs<-MCPY_samples_filtered_pos[ which(MCPY_samples_filtered_pos$unique=='TRUE'), ]
FB_unique_circRNAs<-FB_samples_filtered_pos[ which(FB_samples_filtered_pos$unique=='TRUE'), ]
neuronal_unique_circRNAs<-neuronal_samples_filtered_pos[ which(neuronal_samples_filtered_pos$unique=='TRUE'), ]
NonNeuronal_unique_circRNAs<-NonNeuronal_samples_filtered_pos[ which(NonNeuronal_samples_filtered_pos$unique=='TRUE'), ]

#example to test unique
#df_unique[grep("circ_499", df_unique$name), ] #confirm unique
#Merge_circexp_norm_filtered[grep("circ_499",Merge_circexp_norm_filtered$name), ] #present in FB and PBMC but at lower levels PBMC

#Write tables

#write.table(SNDA_unique_circRNAs, "SNDA_unique_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)
#write.table(TCPY_unique_circRNAs, "TCPY_unique_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)
#write.table(PBMC_unique_circRNAs, "PBMC_unique_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)
#write.table(MCPY_unique_circRNAs, "MCPY_unique_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)
#write.table(FB_unique_circRNAs, "FB_unique_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)
#write.table(neuronal_unique_circRNAs, "neuronal_unique_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)
#write.table(NonNeuronal_unique_circRNAs, "NonNeuronal_unique_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)


#Shared circRNAs per sample SNDA

hist(SNDA_samples_filtered$columnGreat, breaks=40, xlim=c(0,60),ylim=c(0,400), col = "blue",main="Shared expression of circRNAs in SNDA Samples",xlab="Number of SNDA Samples", ylab="CircRNAS")
axis(1, at=seq(0,60,5))
#shared_SNDA<-merge(circexp, mapsplice, by=c("V1","V3")) #end is same, mapsplice isnt' in 0 based so not in bed

head(SNDA_samples_filtered[order(-SNDA_samples_filtered$columnGreat),][ (num_SNDA+1):(num_SNDA+7)])
#chrom     start       end        name score strand columnGreat
#134333  chr2  40655612  40657444 circ_134333 11647      -          49
#117125 chr17  20107645  20109225 circ_117125  4208      +          46
#95547  chr17  65941524  65972074  circ_95547  2711      +          37
#16800   chr1 117944807 117963271  circ_16800  7360      +          36
#97736  chr18   9182379   9221997  circ_97736  7390      +          36
#11802   chr1  27056141  27059283  circ_11802  1910      +          35

#obtain row mean (not just unique ones but ones expressed in cell type)

SNDA_samples_filtered$meanRPM<-rowMeans(SNDA_samples_filtered[,1:num_SNDA])
TCPY_samples_filtered$meanRPM<-rowMeans(TCPY_samples_filtered[,1:num_TCPY])
PBMC_samples_filtered$meanRPM<-rowMeans(PBMC_samples_filtered[,1:num_PBMC])
MCPY_samples_filtered$meanRPM<-rowMeans(MCPY_samples_filtered[,1:num_MCPY])
FB_samples_filtered$meanRPM<-rowMeans(FB_samples_filtered[,1:num_FB])
neuronal_samples_filtered$meanRPM<-rowMeans(neuronal_samples_filtered[,1:num_neuronal])
NonNeuronal_samples_filtered$meanRPM<-rowMeans(NonNeuronal_samples_filtered[,1:num_nonNeuronal])


head(SNDA_samples_filtered[order(-SNDA_samples_filtered$meanRPM),][(num_SNDA+1):(num_SNDA+8)])

 head(SNDA_samples_filtered[order(-SNDA_samples_filtered$meanRPM),][ 65:72])
#chrom     start       end        name score strand columnGreat    meanRPM
#134333  chr2  40655612  40657444 circ_134333 11647      -          49 0.38800363
#117125 chr17  20107645  20109225 circ_117125  4208      +          46 0.15094636
#95547  chr17  65941524  65972074  circ_95547  2711      +          37 0.09107188
#16800   chr1 117944807 117963271  circ_16800  7360      +          36 0.07969346
#99398   chr1  67356836  67371058  circ_99398  1733      -          19 0.07647262
#11802   chr1  27056141  27059283  circ_11802  1910      +          35 0.07447913


 head(TCPY_samples_filtered[order(-TCPY_samples_filtered$meanRPM),][(num_TCPY+1):(num_TCPY+8)])
#chrom     start       end        name score strand columnGreat   meanRPM
#16800   chr1 117944807 117963271  circ_16800  7360      +           9 0.4902523
#134333  chr2  40655612  40657444 circ_134333 11647      -           7 0.4760345
#147714 chr16  14328002  14334341 circ_147714   675      +           1 0.3849846
#28242   chr4  54280781  54310270  circ_28242   682      +           2 0.2917123
#102372  chr1 243708811 244006584 circ_102372  1594      -           8 0.2875935
#136878  chr3  55717821  56026278 circ_136878   872      -           6 0.2706504

head(PBMC_samples_filtered[order(-PBMC_samples_filtered$meanRPM),][(num_PBMC+1):(num_PBMC+8)])
#chrom     start       end        name score strand columnGreat   meanRPM
#97736  chr18   9182379   9221997  circ_97736  7390      +           4 3.4252354
#16800   chr1 117944807 117963271  circ_16800  7360      +           4 2.0697196
#117125 chr17  20107645  20109225 circ_117125  4208      +           4 1.3099024
#134333  chr2  40655612  40657444 circ_134333 11647      -           4 1.1318929
#58370   chr4  48371865  48385801  circ_58370  1808      +           4 1.0060072
#95547  chr17  65941524  65972074  circ_95547  2711      +           4 0.8931477

head(PBMC_samples_filtered[order(-PBMC_samples_filtered$columnGreat),][ (num_PBMC+1):(num_PBMC+8)])
chrom     start       end       name score strand columnGreat    meanRPM
#79     chr1 171492359 171502100    circ_79    70      +           4 0.06902097
#2332   chr3 169831147 169854453  circ_2332   140      -           4 0.14961726
#2786   chr3 169694733 169706147  circ_2786   444      +           4 0.28212613
#11802  chr1  27056141  27059283 circ_11802  1910      +           4 0.51532506
#13589  chr9   5968018   5988545 circ_13589   687      -           4 0.19036014
#16800  chr1 117944807 117963271 circ_16800  7360      +           4 2.06971960


head(MCPY_samples_filtered[order(-MCPY_samples_filtered$meanRPM),][(num_MCPY+1):(num_MCPY+8)])
#chrom     start       end        name score strand columnGreat   meanRPM
#16800   chr1 117944807 117963271  circ_16800  7360      +           3 0.6834573
#66682   chr2 214174782 214239843  circ_66682   749      +           1 0.4038502
#179409  chr1  21220009  21268823 circ_179409   918      -           2 0.3534628
#106882  chr1 240370098 240421332 circ_106882   466      +           1 0.3439924
#126458  chr1 211192205 211264032 circ_126458   419      -           2 0.2967044
#117125 chr17  20107645  20109225 circ_117125  4208      +           3 0.2847736

 head(MCPY_samples_filtered[order(-MCPY_samples_filtered$columnGreat),][ (num_MCPY+1):(num_MCPY+8)])
#chrom     start       end       name score strand columnGreat    meanRPM
#2421   chr1  32381495  32385259  circ_2421  2186      -           3 0.15636047
#13029  chr1  54059780  54065951 circ_13029   109      -           3 0.03468401
#16373  chr2 200173482 200246543 circ_16373   218      -           3 0.05659342
#16800  chr1 117944807 117963271 circ_16800  7360      +           3 0.68345733
#32053  chr7   8043537   8110761 circ_32053   615      +           3 0.23198566
#36767 chr22  50810448  50832564 circ_36767   102      +           3 0.02319532

head(FB_samples_filtered[order(-FB_samples_filtered$meanRPM),][(num_FB+1):(num_FB+8)])
#chrom     start       end        name score strand columnGreat   meanRPM
#97736  chr18   9182379   9221997  circ_97736  7390      +           3 0.7618473
#134333  chr2  40655612  40657444 circ_134333 11647      -           3 0.6766341
#60104   chr5  64747301  64769779  circ_60104   400      -           3 0.5015491
#16800   chr1 117944807 117963271  circ_16800  7360      +           3 0.4355838
#11802   chr1  27056141  27059283  circ_11802  1910      +           3 0.3292680
#11093   chr3 171965322 171969331  circ_11093  1235      +           3 0.2525699

head(FB_samples_filtered[order(-FB_samples_filtered$columnGreat),][(num_FB+1):(num_FB+8)])
chrom     start       end      name score strand columnGreat    meanRPM
#1232  chr6 131276244 131277639 circ_1232   502      -           3 0.14701978
#2421  chr1  32381495  32385259 circ_2421  2186      -           3 0.10054845
#3200 chr15  39877764  39878158 circ_3200    33      +           3 0.04615300
#3480 chr15  39877764  39878167 circ_3480    35      +           3 0.04949432
#4928 chr15  39877764  39878176 circ_4928    92      +           3 0.11337987
#5912 chr15  39877764  39878140 circ_5912    22      +           3 0.02886838

head(neuronal_samples_filtered[order(-neuronal_samples_filtered$meanRPM),][(num_neuronal+1):(num_neuronal+8)])
chrom     start       end        name score strand columnGreat    meanRPM
#134333  chr2  40655612  40657444 circ_134333 11647      -          59 0.38948571
#117125 chr17  20107645  20109225 circ_117125  4208      +          54 0.15451849
#16800   chr1 117944807 117963271  circ_16800  7360      +          48 0.15410321
#95547  chr17  65941524  65972074  circ_95547  2711      +          47 0.10489177
#11802   chr1  27056141  27059283  circ_11802  1910      +          41 0.08680448
#97736  chr18   9182379   9221997  circ_97736  7390      +          43 0.08578242


head(neuronal_samples_filtered[order(-neuronal_samples_filtered$columnGreat),][ (num_neuronal+1):(num_neuronal+8)])
chrom     start       end        name score strand columnGreat    meanRPM
#134333  chr2  40655612  40657444 circ_134333 11647      -          59 0.38948571
#117125 chr17  20107645  20109225 circ_117125  4208      +          54 0.15451849
#16800   chr1 117944807 117963271  circ_16800  7360      +          48 0.15410321
#95547  chr17  65941524  65972074  circ_95547  2711      +          47 0.10489177
#97736  chr18   9182379   9221997  circ_97736  7390      +          43 0.08578242
#2421    chr1  32381495  32385259   circ_2421  2186      -          42 0.05366963

head(NonNeuronal_samples_filtered[order(-NonNeuronal_samples_filtered$meanRPM),][(num_nonNeuronal+1):(num_nonNeuronal+8)])
chrom     start       end        name score strand columnGreat   meanRPM
#97736  chr18   9182379   9221997  circ_97736  7390      +           7 2.2837833
#16800   chr1 117944807 117963271  circ_16800  7360      +           7 1.3693757
#134333  chr2  40655612  40657444 circ_134333 11647      -           7 0.9367820
#117125 chr17  20107645  20109225 circ_117125  4208      +           7 0.8087355
#58370   chr4  48371865  48385801  circ_58370  1808      +           7 0.6262296
#177406  chr4 178274461 178281831 circ_177406   879      +           7 0.6003218

head(NonNeuronal_samples_filtered[order(-NonNeuronal_samples_filtered$columnGreat),][(num_nonNeuronal+1):(num_nonNeuronal+8)])
#chrom     start       end       name score strand columnGreat   meanRPM
#11802  chr1  27056141  27059283 circ_11802  1910      +           7 0.4355863
#16800  chr1 117944807 117963271 circ_16800  7360      +           7 1.3693757
#24347 chr16  68155889  68157024 circ_24347   588      +           7 0.2023212
#44355 chr16  68155889  68160513 circ_44355   771      +           7 0.3812181
#58351  chr5 122881110 122893258 circ_58351   630      +           7 0.3290225
#58370  chr4  48371865  48385801 circ_58370  1808      +           7 0.6262296


#Write tables with circRNAs non unique

write.table((SNDA_samples_filtered[(num_SNDA+1):(num_SNDA+8)]), "SNDA_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table((TCPY_samples_filtered[(num_TCPY+1):(num_TCPY+8)]), "TCPY_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table((PBMC_samples_filtered[(num_PBMC+1):(num_PBMC+8)]), "PBMC_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table((MCPY_samples_filtered[(num_MCPY+1):(num_MCPY+8)]), "MCPY_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table((FB_samples_filtered[(num_FB+1):(num_FB+8)]), "FB_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table((neuronal_samples_filtered[(num_neuronal+1):(num_neuronal+8)]), "neuronal_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table((NonNeuronal_samples_filtered[(num_nonNeuronal+1):(num_nonNeuronal+8)]), "NonNeuronal_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)

#save unique with mean RPM

SNDA_merge_unique_filtered<-merge(SNDA_unique_circRNAs, (SNDA_samples_filtered[(num_SNDA+1):(num_SNDA+8)]),by=c("chrom","start","end","name","strand") )
TCPY_merge_unique_filtered<-merge(TCPY_unique_circRNAs, (TCPY_samples_filtered[(num_TCPY+1):(num_TCPY+8)]),by=c("chrom","start","end","name","strand") )
PBMC_merge_unique_filtered<-merge(PBMC_unique_circRNAs, (PBMC_samples_filtered[(num_PBMC+1):(num_PBMC+8)]),by=c("chrom","start","end","name","strand") )
MCPY_merge_unique_filtered<-merge(MCPY_unique_circRNAs, (MCPY_samples_filtered[(num_MCPY+1):(num_MCPY+8)]),by=c("chrom","start","end","name","strand") )
FB_merge_unique_filtered<-merge(FB_unique_circRNAs, (FB_samples_filtered[(num_FB+1):(num_FB+8)]),by=c("chrom","start","end","name","strand") )
neuronal_merge_unique_filtered<-merge(neuronal_unique_circRNAs, (neuronal_samples_filtered[(num_neuronal+1):(num_neuronal+8)]),by=c("chrom","start","end","name","strand") )
NonNeuronal_merge_unique_filtered<-merge(NonNeuronal_unique_circRNAs, (NonNeuronal_samples_filtered[(num_nonNeuronal+1):(num_nonNeuronal+8)]),by=c("chrom","start","end","name","strand") )


write.table(SNDA_merge_unique_filtered, "SNDA_unique_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(TCPY_merge_unique_filtered, "TCPY_unique_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(PBMC_merge_unique_filtered, "PBMC_unique_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(MCPY_merge_unique_filtered, "MCPY_unique_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(FB_merge_unique_filtered, "FB_unique_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(neuronal_merge_unique_filtered, "neuronal_unique_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(NonNeuronal_merge_unique_filtered, "NonNeuronal_unique_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt", sep="\t", quote=FALSE, row.names=FALSE)

###END BRAINCODE ###

#################################################################
#################################################################
###### RNASE R samples ananysis No Fibroblast No Outliers########
#################################################################
#################################################################

setwd("~/projects/circRNA/")
#removed "#" on first line as well as ful name to just sample name
Merge_circexp_RNAse <- read.table("Merge_circexplorer_RNAse.bed",sep="\t",header=TRUE, comment.char ='') 
colnames(Merge_circexp_RNAse)=gsub("X.|_candidates.bed_circReads","",colnames(Merge_circexp_RNAse))
readsNum_RNAse<- read.table("ReadsNum_RNAseR_filtered.txt",sep="\t",header=TRUE) 
#keep rows we need
readsNum_RNAse$HC_WGC082362_SN_R_b1 <- NULL
readsNum_RNAse$HC_WGC082363_SN_M_b1 <- NULL
readsNum_RNAse$HC_WGC082364_SN_R_b2. <- NULL
readsNum_RNAse$HC_WGC082365_SN_M_b2 <- NULL

nrow(Merge_circexp_RNAse)
#[1] 226090

readsNum_RNAse_million<-(readsNum_RNAse/10^6)
samplenum_RNase<-ncol(readsNum_RNAse_million) #sample number


#Note: still use the minimal sequencing depth in Braincode samples as cutoff threshold. 
readsNum_filtered<- read.table("ReadsNum_Braincode_filtered.txt",sep="\t",header=TRUE) 
readsNum_million<-(readsNum_filtered/10^6)
#remove file that didn't run
readsNum_million$HC_BN12.44_TCPY_5_rep1<- NULL
#remove outlier samples
readsNum_million$HC_MD5028_SNDA_1_rep1<- NULL
readsNum_million$HC_MGH1026_SNDA_1_rep1<- NULL
#keep amp PBMC but not SN since it is whole tissue
readsNum_million$HC_UWA616_SN_6_rep1.amplified<- NULL
min(readsNum_million)
# 64.42181
threshold<-1/(min(readsNum_million))
#filtered
# 0.01552269

Merge_circexp_RNAse_norm<-sweep(as.matrix(Merge_circexp_RNAse[7:(6+samplenum_RNase)]),2,as.matrix(readsNum_RNAse_million),"/")
Merge_circexp_RNAse_norm<-as.data.frame(Merge_circexp_RNAse_norm) #back to df
#add back other columns 
Merge_circexp_RNAse_norm$chrom=Merge_circexp_RNAse$chrom
Merge_circexp_RNAse_norm$start=Merge_circexp_RNAse$start
Merge_circexp_RNAse_norm$end=Merge_circexp_RNAse$end
Merge_circexp_RNAse_norm$name=Merge_circexp_RNAse$name
Merge_circexp_RNAse_norm$score=Merge_circexp_RNAse$score
Merge_circexp_RNAse_norm$strand=Merge_circexp_RNAse$strand

#count how many rows have samples beyond threshold

Merge_circexp_RNAse_norm$columnGreat<-rowSums(Merge_circexp_RNAse_norm[1:samplenum_RNase]>threshold) #add how many rows have cols greater than threshold
#choose value for pseudocounts, this will be smallest power to detect, so 1/largest library
pseudocount<-1/max(readsNum_RNAse_million)
# 0.001823273
Merge_circexp_RNAse_norm_filtered<-subset(Merge_circexp_RNAse_norm, columnGreat>0) #keep those rows with at least one sample at threshold

# save dataset
save(Merge_circexp_RNAse,Merge_circexp_RNAse_norm, Merge_circexp_RNAse_norm_filtered, file="Merge_circexp.RNAse.Rdata")

nrow(Merge_circexp_RNAse_norm)
#[1] 226090

nrow(Merge_circexp_RNAse_norm_filtered) #count at least one sample above threshold
#[1] 73728

head(Merge_circexp_RNAse_norm_filtered[order(-Merge_circexp_RNAse_norm_filtered$columnGreat),][ 12:18],352)

#chrom     start       end        name  score strand columnGreat
#76     chr18   9195548   9221997     circ_76  26893      +          11
#321     chr1  28362054  28374937    circ_321  34018      -          11
#771     chr2  99802639  99812219    circ_771  68640      +          11
#1098   chr11  34111725  34113603   circ_1098   7142      +          11
#1252   chr15  91030185  91035943   circ_1252   4427      +          11
#1934   chr15  68434283  68457142   circ_1934  67825      +          11
#2154    chr6 131276244 131277639   circ_2154  17056      -          11
#2413   chr10 104445572 104465246   circ_2413  11027      -          11


circ_count_RNAse<-colSums(Merge_circexp_RNAse_norm_filtered[1:11]>threshold) #add how many samples have values greater than 0, so at least one circRNA

#HC_WGC082362_SN_R_b1_r1 HC_WGC082363_SN_M_b1_r1 HC_WGC082364_SN_R_b2_r1 HC_WGC082365_SN_M_b2_r1 
#5026                    5370                    5009                    4817 
#HC_WGC082366_TC_R_b1    HC_WGC082367_TC_M_b1    HC_WGC082368_TC_R_b2    HC_WGC082369_TC_M_b2 
#19536                    7828                   42052                   13308 
#HC_WGC082371_FB_M_b1  HC_WGC082372_PBMC_R_b1  HC_WGC082373_PBMC_M_b1 
#2940                   16228                    7765 

op <- par(mar=c(12,6,4,2))
barplot(circ_count_RNAse,width=0.8,las=2,cex.names=0.9,
        cex.axis=0.8, ylab="Number of circRNAs and ciRNAs", main="circRNAs and ciRNAs count RNAse and Mock \n with CircExplorer Normalized and Filtered",col=c("red","darkblue","red","darkblue","red","darkblue","red","darkblue","darkblue","red","darkblue"))
legend( "topright", inset = c(-0.2,0),cex = 1.5, bty = "n", legend = c("RNase R", "Mock"), col = c("red", "darkblue"), pt.bg = c("red","blue"),pch = c(22,22),x.intersp = .3 ,y.intersp = .9)
rm(op)

#PLot circ count by depth

#all three
plot(as.numeric(readsNum_RNAse)/(10^6),as.numeric(circ_count_RNAse),col=c("red","cyan","red","cyan","red","cyan","red","cyan","cyan","red","cyan"),type = 'p', pch=16,cex = 1, main="Number of circRNAs by Sequencing Depth",xlim=c(0,600),ylim=c(0,45000),xlab="Number of Raw Reads (Million)",ylab="Number of circRNAS")
points(as.numeric(readsNum_filtered[1:81])/(10^6),as.numeric(circ_count_BC),col="black",type = 'p', pch=16,cex = 1)
legend( "topleft", cex = 1, bty = "n", legend = c("RNase R", "Mock","BRAINCODE"), col = c("red", "cyan","Black"), pt.bg = c("red","cyan","black"),pch = c(22,22,22),x.intersp = .5 ,y.intersp = 1)

#BC
plot(as.numeric(readsNum_filtered[1:81])/(10^6),as.numeric(circ_count_BC),col="black",type = 'p', pch=16,cex = .9, main="Number of circRNAs by Sequencing Depth \n BRAINCODE",xlab="Number of Raw Reads (Million)",ylab="Number of circRNAS")


#Paired T testwithout FB
RNAse_group = c(circ_count_RNAse[1], circ_count_RNAse[3], circ_count_RNAse[5], circ_count_RNAse[7], circ_count_RNAse[10])
mean(RNAse_group )
#17570.2
Mock_group= c(circ_count_RNAse[2], circ_count_RNAse[4], circ_count_RNAse[6], circ_count_RNAse[8], circ_count_RNAse[11])
mean(Mock_group)
#7817.6

shapiro.test(RNAse_group)
#Shapiro-Wilk normality test

#data:  RNAse_group
#W = 0.85948, p-value = 0.2264

shapiro.test(Mock_group)
#Shapiro-Wilk normality test

#data:  Mock_group
#W = 0.86166, p-value = 0.2343


t.test(RNAse_group,Mock_group, paired=TRUE,alternative="greater")
#alternative = "greater" is the alternative that x has a larger mean than y.

#Paired t-test

#data:  RNAse_group and Mock_group
#t = 1.8439, df = 4, p-value = 0.06949
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -1523.1     Inf
#sample estimates:
#  mean of the differences 
#9752.6 

#violin plot of paired values
library(vioplot)
plot(1, 1, xlim = c(0, 4), ylim = c(0,45000), type = 'n', xlab = '', ylab = '', xaxt = 'n', main= "circRNAs and ciRNAs count \n RNAse and Mock Distributions")
vioplot(RNAse_group, at = 1, add = T, col = 'red')
vioplot(Mock_group, at = 3, add = T, col = 'darkblue')
axis(1, at = c(1,3), labels = c('RNase R', 'Mock'))


#confirm they're all unique
coords<-(Merge_circexp_RNAse_norm_filtered[grepl( "chrom|start|end|name|score|strand" , names(Merge_circexp_RNAse_norm_filtered))]) #select necessary columns

anyDuplicated(coords) 
#[1] 0

#Are the pair circRNAs overlapping?

Merge_circexp_RNAse_norm_filtered_noFB=Merge_circexp_RNAse_norm_filtered
Merge_circexp_RNAse_norm_filtered_noFB$HC_WGC082371_FB_M_b1 <-NULL

nrow(Merge_circexp_RNAse_norm_filtered_noFB)
#[1] 66023

Merge_circexp_RNAse_norm_filtered_noFB$columnGreat<-rowSums(Merge_circexp_RNAse_norm_filtered_noFB[1:10]>threshold) #add how many rows have cols greater than threshold
Merge_circexp_RNAse_norm_filtered_noFB<-subset(Merge_circexp_RNAse_norm_filtered_noFB, columnGreat>0) #keep those rows with at least one sample at threshold

nrow(Merge_circexp_RNAse_norm_filtered_noFB)
#[1] 65681

Merge_circexp_RNAse_norm_filtered_noFB$MeanRNase=rowMeans(Merge_circexp_RNAse_norm_filtered_noFB[,c(1,3,5,7,9)])
Merge_circexp_RNAse_norm_filtered_noFB$MeanMock=rowMeans(Merge_circexp_RNAse_norm_filtered_noFB[,c(2,4,6,8,10)])
Merge_circexp_RNAse_norm_filtered_noFB$Log2FC_R_M<-((log2(Merge_circexp_RNAse_norm_filtered_noFB$MeanRNase+pseudocount))-(log2(Merge_circexp_RNAse_norm_filtered_noFB$MeanMock+pseudocount)))
Merge_circexp_RNAse_norm_filtered_noFB$MeanAll<-rowMeans(Merge_circexp_RNAse_norm_filtered_noFB[,c(1:10)])
plot(Merge_circexp_RNAse_norm_filtered_noFB$MeanAll,Merge_circexp_RNAse_norm_filtered_noFB$Log2FC_R_M,cex=.2, main="MA Plot SN, TC, and PBMC \n RNAse vs. Mock", xlab="Mean Normalized RPM", ylab="Mean Log2 FC (RNase R / Mock)")
# paired t test

Merge_circexp_RNAse_norm_filtered_noFB$p_value<- apply(Merge_circexp_RNAse_norm_filtered_noFB[,1:10], 1, function (x) {t.test(x[c(1,3,5,7,9)],x[c(2,4,6,8,10)],paired=TRUE,,alternative="greater")$p.value})
Merge_circexp_RNAse_norm_filtered_noFB$FDR<-p.adjust(Merge_circexp_RNAse_norm_filtered_noFB$p_value, method="fdr")  

low_pval <- Merge_circexp_RNAse_norm_filtered_noFB[Merge_circexp_RNAse_norm_filtered_noFB$p_value <= 0.05, ]
hist(Merge_circexp_RNAse_norm_filtered_noFB$FDR, main="FDR Values Distribution in \n 65681 circRNAs", xlab="FDR Values", ylab="Frequency", breaks=50)
hist(Merge_circexp_RNAse_norm_filtered_noFB$p_value,main="P-values Distribution in \n 65681 circRNAs", xlab="P-value", ylab="Frequency", breaks=50,ylim=c(0,40000),xaxt='n')
axis(side=1, at=seq(0, 1, .05))


####select Tissue regions###

SN_pairs<-(Merge_circexp_RNAse_norm_filtered_noFB[grepl( "SN|chrom|start|end|name|score|strand" , names(Merge_circexp_RNAse_norm_filtered_noFB))])

#select only circRNAs that are expressed above threshold in at least one of the 2

SN_pairs$columnGreat<-rowSums(SN_pairs[1:4]>threshold) #add how many rows have cols greater than threshold
SN_pairs_filtered<-subset(SN_pairs, columnGreat>0) #keep those rows with at least one sample at threshold

nrow(SN_pairs_filtered)
#[1] 13340 #at least one SN sample (either mock or rnaseR above threshold)

#log2FC
#add small value to prevent INF

SN_pairs_filtered$MeanRNase_SN=rowMeans(SN_pairs_filtered[,c(1,3)])
SN_pairs_filtered$MeanMock_SN=rowMeans(SN_pairs_filtered[,c(2,4)])
SN_pairs_filtered$Log2FC_R_M_SN<-((log2(SN_pairs_filtered$MeanRNase_SN+pseudocount))-(log2(SN_pairs_filtered$MeanMock_SN+pseudocount)))
SN_pairs_filtered$MeanAll_SN<-rowMeans(SN_pairs_filtered[,c(1:4)])
plot(SN_pairs_filtered$MeanAll_SN,SN_pairs_filtered$Log2FC_R_M_SN,cex=.2, main="MA Plot Substantia Nigra \n RNAse vs. Mock", xlab="Mean Normalized RPM", ylab="Mean Log2 Fold Change (RNase R / Mock)")
abline(h=0,col="red",lwd = 3)
hist(SN_pairs_filtered$Log2FC_R_M_SN, breaks=30, main="SN Both Pairs Fold Change \n RNaseR treated vs. Mock frequency", xlab="Mean Log2 Fold Change (RNase R / Mock)")

#view highest average ones

head(SN_pairs_filtered[order(-SN_pairs_filtered$MeanAll_SN),])

#HC_WGC082362_SN_R_b1_r1 HC_WGC082363_SN_M_b1_r1 HC_WGC082364_SN_R_b2_r1 HC_WGC082365_SN_M_b2_r1 chrom     start
#155085               0.7084684               0.7565177              53.4436175                2.909665  chr1  65830317
#11549                2.2828427               3.1126948              48.0378168                3.853621 chr13  78293666
#4290                 2.9575745               2.5871749              32.8677629                3.417949  chr1  32381495
#148087               7.2983493               8.7894648              15.5246678                8.049555 chr15  30053341
#30149               15.4119995              16.5047900               0.9360158                5.160636  chr1 117944807
#149586               0.0000000               0.0000000              35.3459370                1.338135 chr11  46098304
#end        name  score strand columnGreat MeanRNase_SN MeanMock_SN Log2FC_R_M_SN MeanAll_SN
#155085  65831879 circ_155085  70426      +           4    27.076043   1.8330914     3.8833292  14.454567
#11549   78327493  circ_11549  47622      +           4    25.160330   3.4831578     2.8520327  14.321744
#4290    32385259   circ_4290 142054      -           4    17.912669   3.0025619     2.5759854  10.457615
#148087  30065560 circ_148087  50445      -           4    11.411509   8.4195099     0.4385995   9.915509
#30149  117963271  circ_30149  90673      +           4     8.174008  10.8327130    -0.4062001   9.503360
#149586  46113774 circ_149586 446993      -           2    17.672968   0.6690674     4.7194717   9.171018


#view lowest average ones

head(SN_pairs_filtered[order(SN_pairs_filtered$MeanAll_SN),])

#HC_WGC082362_SN_R_b1_r1 HC_WGC082363_SN_M_b1_r1 HC_WGC082364_SN_R_b2_r1 HC_WGC082365_SN_M_b2_r1 chrom     start
#793                        0                       0                       0              0.01555971 chr12 124158164
#903                        0                       0                       0              0.01555971  chr3 197707196
#952                        0                       0                       0              0.01555971 chr15  90972898
#1168                       0                       0                       0              0.01555971 chr17  66303083
#1523                       0                       0                       0              0.01555971 chr13  42249350
#1632                       0                       0                       0              0.01555971 chr14  71428942
#end      name score strand columnGreat MeanRNase_SN MeanMock_SN Log2FC_R_M_SN  MeanAll_SN
#793  124172724  circ_793    51      +           1            0 0.007779853     -2.396974 0.003889927
#903  197726213  circ_903  1964      +           1            0 0.007779853     -2.396974 0.003889927
#952   90973118  circ_952     3      +           1            0 0.007779853     -2.396974 0.003889927
#1168  66397591 circ_1168    72      +           1            0 0.007779853     -2.396974 0.003889927
#1523  42267108 circ_1523    60      -           1            0 0.007779853     -2.396974 0.003889927
#1632  71479919 circ_1632    12      +           1            0 0.007779853     -2.396974 0.003889927

#TC pair:

TC_pairs<-(Merge_circexp_RNAse_norm_filtered_noFB[grepl( "TC|chrom|start|end|name|score|strand" , names(Merge_circexp_RNAse_norm_filtered_noFB))])

#select only circRNAs that are expressed above threshold in at least one of the 2

TC_pairs$columnGreat<-rowSums(TC_pairs[1:4]>threshold) #add how many rows have cols greater than threshold
TC_pairs_filtered<-subset(TC_pairs, columnGreat>0) #keep those rows with at least one sample at threshold

nrow(TC_pairs_filtered)
#50407 #at least one TC sample (either mock or rnaseR above threshold)

#log2FC
#add small value to prevent INF
TC_pairs_filtered$MeanRNase_TC=rowMeans(TC_pairs_filtered[,c(1,3)])
TC_pairs_filtered$MeanMock_TC=rowMeans(TC_pairs_filtered[,c(2,4)])
TC_pairs_filtered$Log2FC_R_M_TC<-((log2(TC_pairs_filtered$MeanRNase_TC+pseudocount))-(log2(TC_pairs_filtered$MeanMock_TC+pseudocount)))
TC_pairs_filtered$MeanAll_TC<-rowMeans(TC_pairs_filtered[,c(1:4)])
plot(TC_pairs_filtered$MeanAll_TC,TC_pairs_filtered$Log2FC_R_M_TC,cex=.2, main="MA Plot Temporal Cortex \n RNAse vs. Mock", xlab="Mean Normalized RPM", ylab="Mean Log2 Fold Change (RNase R / Mock)")
abline(h=0,col="red",lwd = 3)
hist(TC_pairs_filtered$Log2FC_R_M_TC,breaks=30, main="TC Both Pairs Fold Change \n RNaseR treated vs. Mock frequency", xlab="Mean Log2 Fold Change (RNase R / Mock)")



#view highest average ones

head(TC_pairs_filtered[order(-TC_pairs_filtered$MeanAll_TC),])
#HC_WGC082366_TC_R_b1 HC_WGC082367_TC_M_b1 HC_WGC082368_TC_R_b2 HC_WGC082369_TC_M_b2 chrom    start      end
#149586            419.96312             8.048388             453.3160            12.630941 chr11 46098304 46113774
#31267             199.66032             3.711795             372.1208            11.650319  chr6 54013853 54095715
#132921            106.48195             1.985379             244.1783             6.387065  chr2 72945231 72960247
#50482             125.69931             1.644206             198.2953             4.629754 chr15 59204761 59209198
#127463            117.94111             1.981268             186.4187             3.627436  chr3 47079155 47088111
#158814             92.92388             2.606066             190.5148             6.695136 chr17 26449628 26499644
#name  score strand columnGreat MeanRNase_TC MeanMock_TC Log2FC_R_M_TC MeanAll_TC
#149586 circ_149586 446993      -           4     436.6396   10.339665      5.399933  223.48961
#31267   circ_31267 287545      +           4     285.8906    7.681057      5.217681  146.78582
#132921 circ_132921 197971      -           4     175.3301    4.186222      5.387668   89.75818
#50482   circ_50482 359595      -           4     161.9973    3.136980      5.689628   82.56714
#127463 circ_127463 175688      -           4     152.1799    2.804352      5.761046   77.49212
#158814 circ_158814 187461      +           4     141.7194    4.650601      4.928929   73.18498

#low expression

head(TC_pairs_filtered[order(TC_pairs_filtered$MeanAll_TC),])
#HC_WGC082366_TC_R_b1 HC_WGC082367_TC_M_b1 HC_WGC082368_TC_R_b2 HC_WGC082369_TC_M_b2 chrom     start       end
#43                      0           0.01644206                    0                    0  chr1  97235258  97272514
#1193                    0           0.01644206                    0                    0  chr3 179082905 179093129
#1948                    0           0.01644206                    0                    0  chr6 152734486 152757236
#2505                    0           0.01644206                    0                    0 chr12  15656841  15704605
#3930                    0           0.01644206                    0                    0  chr8 100729398 100796704
#5026                    0           0.01644206                    0                    0 chr11 120328788 120330029
#name score strand columnGreat MeanRNase_TC MeanMock_TC Log2FC_R_M_TC  MeanAll_TC
#43     circ_43     4      +           1            0 0.008221029     -2.461775 0.004110515
#1193 circ_1193     4      +           1            0 0.008221029     -2.461775 0.004110515
#1948 circ_1948     4      -           1            0 0.008221029     -2.461775 0.004110515
#2505 circ_2505     4      +           1            0 0.008221029     -2.461775 0.004110515
#3930 circ_3930     4      +           1            0 0.008221029     -2.461775 0.004110515
#5026 circ_5026     4      +           1            0 0.008221029     -2.461775 0.004110515

#PBMC select

PBMC_pairs<-(Merge_circexp_RNAse_norm_filtered_noFB[grepl( "PBMC|chrom|start|end|name|score|strand" , names(Merge_circexp_RNAse_norm_filtered_noFB))])

#select only circRNAs that are expressed above threshold in at least one of the 2

PBMC_pairs$columnGreat<-rowSums(PBMC_pairs[1:2]>threshold) #add how many rows have cols greater than threshold
PBMC_pairs_filtered<-subset(PBMC_pairs, columnGreat>0) #keep those rows with at least one sample at threshold

nrow(PBMC_pairs_filtered)
#19885 #at least one SN sample (either mock or rnaseR above threshold)

#log2FC
#add small value to prevent INF

PBMC_pairs_filtered$MeanRNase_PBMC=PBMC_pairs_filtered$HC_WGC082372_PBMC_R_b1
PBMC_pairs_filtered$MeanMock_PBMC=PBMC_pairs_filtered$HC_WGC082373_PBMC_M_b1
PBMC_pairs_filtered$Log2FC_R_M_PBMC<-((log2(PBMC_pairs_filtered$MeanRNase_PBMC+pseudocount))-(log2(PBMC_pairs_filtered$MeanMock_PBMC+pseudocount)))
PBMC_pairs_filtered$MeanAll_PBMC<-rowMeans(PBMC_pairs_filtered[,c(1:2)])
plot(PBMC_pairs_filtered$MeanAll_PBMC,PBMC_pairs_filtered$Log2FC_R_M_PBMC,cex=.2, main="MA Plot Peripheral Blood Mononuclear Cell \n RNAse vs. Mock", xlab="Mean Normalized RPM", ylab="Log2 Fold Change (RNase R / Mock)")
abline(h=0,col="red",lwd = 3)
hist(PBMC_pairs_filtered$Log2FC_R_M_PBMC, breaks=30, main="PBMC Both Pairs Fold Change \n RNaseR treated vs. Mock frequency", xlab="Log2 Fold Change (RNase R / Mock)")

#view highest average ones

head(PBMC_pairs_filtered[order(-PBMC_pairs_filtered$MeanAll_PBMC),])

#HC_WGC082372_PBMC_R_b1 HC_WGC082373_PBMC_M_b1 chrom    start      end        name  score strand columnGreat
#50482                354.6904              12.391063 chr15 59204761 59209198  circ_50482 359595      -           2
#43926                267.7768               7.889758 chr20  2944917  2945848  circ_43926 206598      +           2
#103554               178.0791               3.892480 chr20 32207322 32211102 circ_103554 207228      +           2
#102142               156.3092               4.361574 chr11 85707868 85714494 circ_102142 136091      -           2
#4839                 137.2596               3.907451 chr10 98703869 98711953   circ_4839 102014      +           2
#3622                 132.6358               3.054100 chr13 28748408 28752072   circ_3622  84135      +           2
#MeanRNase_PBMC MeanMock_PBMC Log2FC_R_M_PBMC MeanAll_PBMC
#50482        354.6904     12.391063        4.838984    183.54073
#43926        267.7768      7.889758        5.084582    137.83328
#103554       178.0791      3.892480        5.515023     90.98577
#102142       156.3092      4.361574        5.162824     80.33538
#4839         137.2596      3.907451        5.133882     70.58354
#3622         132.6358      3.054100        5.439738     67.84495

#lowest

head(PBMC_pairs_filtered[order(PBMC_pairs_filtered$MeanAll_PBMC),])
#HC_WGC082372_PBMC_R_b1 HC_WGC082373_PBMC_M_b1 chrom     start       end      name score strand columnGreat
#491              0.01640946                      0 chr12 108929134 108936928  circ_491     9      -           1
#595              0.01640946                      0 chr15  39874953  39875269  circ_595     9      +           1
#769              0.01640946                      0  chrX  48462778  48463010  circ_769    17      +           1
#954              0.01640946                      0  chr1 153516942 153517129  circ_954     9      -           1
#1062             0.01640946                      0 chr12  55356020  55356214 circ_1062     9      -           1
#1438             0.01640946                      0 chr10  25140316  25160992 circ_1438     9      -           1
#MeanRNase_PBMC MeanMock_PBMC Log2FC_R_M_PBMC MeanAll_PBMC
#491      0.01640946             0        3.321928  0.008204728
#595      0.01640946             0        3.321928  0.008204728
#769      0.01640946             0        3.321928  0.008204728
#954      0.01640946             0        3.321928  0.008204728
#1062     0.01640946             0        3.321928  0.008204728
#1438     0.01640946             0        3.321928  0.008204728

#merge Braincode with RNase R and Mock samples

TCPY_BC_samples_filtered<- read.table("TCPY_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt",sep="\t",header=TRUE) 
TCPY_BC_merge_unique_filtered<- read.table("TCPY_unique_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt" ,sep="\t",header=TRUE) 
SN_BC_samples_filtered<- read.table("SNDA_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt"   ,sep="\t",header=TRUE) 
SN_BC_merge_unique_filtered<- read.table("SNDA_unique_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt"  ,sep="\t",header=TRUE) 
PBMC_BC_samples_filtered<- read.table("PBMC_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt" ,sep="\t",header=TRUE) 
PBMC_BC_merge_unique_filtered<- read.table("PBMC_unique_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt"   ,sep="\t",header=TRUE) 
NonNeuronal_BC_samples_filtered<- read.table("NonNeuronal_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt" ,sep="\t",header=TRUE) 
NonNeuronal_BC_merge_unique_filtered<- read.table("NonNeuronal_unique_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt"   ,sep="\t",header=TRUE) 
neuronal_BC_samples_filtered<- read.table("neuronal_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt" ,sep="\t",header=TRUE) 
neuronal_BC_merge_unique_filtered<- read.table("neuronal_unique_circRNAs_Braincode_circexcp_rpmfilt_noOut.txt",sep="\t",header=TRUE) 


nrow(TCPY_BC_samples_filtered)
#[1] 1377          
nrow(TCPY_BC_merge_unique_filtered)
#[1] 59
nrow(SN_BC_samples_filtered)
#[1] 2426
nrow(SN_BC_merge_unique_filtered)
#[1] 394
nrow(PBMC_BC_samples_filtered)
#[1]2236
nrow(PBMC_BC_merge_unique_filtered)
#[1] 335
nrow(NonNeuronal_BC_samples_filtered)
#[1] 2528
nrow(NonNeuronal_BC_merge_unique_filtered)
#669
nrow(neuronal_BC_samples_filtered)
#[1] 2917
nrow(neuronal_BC_merge_unique_filtered)
#[1] 1058


nrow(TC_pairs_filtered )# Mock union RNaseR
#50407


############# Merge
#BC intersect  RNAse union Mock 
#TCPY_braincode_rnase_mock_merge<-merge(TCPY_BC_samples_filtered, TC_pairs_filtered ,by=c("chrom","start","end","strand") )
#nrow(TCPY_braincode_rnase_mock_merge)
#[1] 601

#TC RNase only:
TC_Rnase<-(Merge_circexp_RNAse_norm_filtered_noFB[grepl( "TC_R|chrom|start|end|name|score|strand" , names(Merge_circexp_RNAse_norm_filtered_noFB))])
#select only circRNAs that are expressed above threshold in at least one of the 2
TC_Rnase$columnGreat<-rowSums(TC_Rnase[1:2]>threshold) #add how many rows have cols greater than threshold
TC_Rnase_filtered<-subset(TC_Rnase, columnGreat>0) #keep those rows with at least one sample at threshold
nrow(TC_Rnase_filtered)
#[1] 47204 #TC RNAse Only

#TC Mock only:
TC_Mock<-(Merge_circexp_RNAse_norm_filtered_noFB[grepl( "TC_M|chrom|start|end|name|score|strand" , names(Merge_circexp_RNAse_norm_filtered_noFB))])
#select only circRNAs that are expressed above threshold in at least one of the 2
TC_Mock$columnGreat<-rowSums(TC_Mock[1:2]>threshold) #add how many rows have cols greater than threshold
TC_Mock_filtered<-subset(TC_Mock, columnGreat>0) #keep those rows with at least one sample at threshold
nrow(TC_Mock_filtered)
#15863 #TC RNAse Only

#Mock and RNase R intersection
TC_Mock_RnaseR<-merge(TC_Mock_filtered, TC_Rnase_filtered ,by=c("chrom","start","end","strand") )
nrow(TC_Mock_RnaseR)
#[1] 12660

#BC intersect  RNAse intersect Mock 
TCPY_braincode_rnase_mock_merge_intersect<-merge(TCPY_BC_samples_filtered, TC_Mock_RnaseR ,by=c("chrom","start","end","strand") )
nrow(TCPY_braincode_rnase_mock_merge_intersect)
#1163

#BC and mock intersection
TCPY_braincode_mock_merge<-merge(TCPY_BC_samples_filtered, TC_Mock_filtered ,by=c("chrom","start","end","strand") )
nrow(TCPY_braincode_mock_merge)
#1221

#BC and RNaseR intersection
TCPY_braincode_rnase_merge<-merge(TCPY_BC_samples_filtered, TC_Rnase_filtered,by=c("chrom","start","end","strand") )
nrow(TCPY_braincode_rnase_merge)
#1242

############# TC Venn diagram

library(VennDiagram)

grid.newpage()
draw.triple.venn(area1 =1377, area2 = 47204, area3 = 15863, n12 = 1242, n23 = 12660, n13 = 1221, 
                 n123 = 1163, category = c("BRAINCODE\n (n=1377)", "RNase R \n(n=47,204)", "Mock (n=15,863)"),  lty = rep("solid", 3), col = rep("white", 3), 
                 fill = c("yellow", "red", "blue"),cat.cex=1.5,cex=2.5)


##BRAINCODE AND SN OVERLAP

#SN RNase only:
SN_Rnase<-(Merge_circexp_RNAse_norm_filtered_noFB[grepl( "SN_R|chrom|start|end|name|score|strand" , names(Merge_circexp_RNAse_norm_filtered_noFB))])
#select only circRNAs that are expressed above threshold in at least one of the 2
SN_Rnase$columnGreat<-rowSums(SN_Rnase[1:2]>threshold) #add how many rows have cols greater than threshold
SN_Rnase_filtered<-subset(SN_Rnase, columnGreat>0) #keep those rows with at least one sample at threshold
nrow(SN_Rnase_filtered)
#8580 #SN RNAse Only

#SN Mock only:
SN_Mock<-(Merge_circexp_RNAse_norm_filtered_noFB[grepl( "SN_M|chrom|start|end|name|score|strand" , names(Merge_circexp_RNAse_norm_filtered_noFB))])
#select only circRNAs that are expressed above threshold in at least one of the 2
SN_Mock$columnGreat<-rowSums(SN_Mock[1:2]>threshold) #add how many rows have cols greater than threshold
SN_Mock_filtered<-subset(SN_Mock, columnGreat>0) #keep those rows with at least one sample at threshold
nrow(SN_Mock_filtered)
#8238 #RN mock Only

#Mock and RNase R intersection
SN_Mock_RnaseR<-merge(SN_Mock_filtered, SN_Rnase_filtered ,by=c("chrom","start","end","strand") )
nrow(SN_Mock_RnaseR)
#3478

#BC intersect  RNAse intersect Mock 
SN_braincode_rnase_mock_merge_intersect<-merge(SN_BC_samples_filtered, SN_Mock_RnaseR ,by=c("chrom","start","end","strand") )
nrow(SN_braincode_rnase_mock_merge_intersect)
#1283

#BC and mock intersection
SN_braincode_mock_merge<-merge(SN_BC_samples_filtered, SN_Mock_filtered ,by=c("chrom","start","end","strand") )
nrow(SN_braincode_mock_merge)
#1529

#BC and RNaseR intersection
SN_braincode_rnase_merge<-merge(SN_BC_samples_filtered, SN_Rnase_filtered,by=c("chrom","start","end","strand") )
nrow(SN_braincode_rnase_merge)
#1525

############# TC Venn diagram

library(VennDiagram)

grid.newpage()
draw.triple.venn(area1 =2426, area2 = 8580, area3 =8238, n12 = 1525, n23 = 3478, n13 = 1529, 
                 n123 = 1283, category = c("BRAINCODE\n (n=2,426)", "RNase R \n(n=8,580)", "Mock (n=8,238)"),,  lty = rep("solid", 3), col = rep("white", 3),  
                 fill = c("yellow", "red", "blue"),cat.cex=1.5,cex=2.5)

# PMCB Only

PBMC_Rnase<-(Merge_circexp_RNAse_norm_filtered_noFB[grepl( "PBMC_R|chrom|start|end|name|score|strand" , names(Merge_circexp_RNAse_norm_filtered_noFB))])
#select only circRNAs that are expressed above threshold in at least one of the 2
PBMC_Rnase$columnGreat<-rowSums(PBMC_Rnase[1]>threshold) #add how many rows have cols greater than threshold
PBMC_Rnase_filtered<-subset(PBMC_Rnase, columnGreat>0) #keep those rows with at least one sample at threshold
nrow(PBMC_Rnase_filtered)
#16228 #PBMC RNAse Only

#PBMC Mock only:
PBMC_Mock<-(Merge_circexp_RNAse_norm_filtered_noFB[grepl( "PBMC_M|chrom|start|end|name|score|strand" , names(Merge_circexp_RNAse_norm_filtered_noFB))])
#select only circRNAs that are expressed above threshold in at least one of the 2
PBMC_Mock$columnGreat<-rowSums(PBMC_Mock[1]>threshold) #add how many rows have cols greater than threshold
PBMC_Mock_filtered<-subset(PBMC_Mock, columnGreat>0) #keep those rows with at least one sample at threshold
nrow(PBMC_Mock_filtered)
#7765 #PMBC mock Only

#Mock and RNase R intersection
PBMC_Mock_RnaseR<-merge(PBMC_Mock_filtered, PBMC_Rnase_filtered ,by=c("chrom","start","end","strand") )
nrow(PBMC_Mock_RnaseR)
#4108

#BC intersect  RNAse intersect Mock 
PBMC_braincode_rnase_mock_merge_intersect<-merge(PBMC_BC_samples_filtered, PBMC_Mock_RnaseR ,by=c("chrom","start","end","strand") )
nrow(PBMC_braincode_rnase_mock_merge_intersect)
#1092

#BC and mock intersection
PBMC_braincode_mock_merge<-merge(PBMC_BC_samples_filtered, PBMC_Mock_filtered ,by=c("chrom","start","end","strand") )
nrow(PBMC_braincode_mock_merge)
#1671

#BC and RNaseR intersection
PBMC_braincode_rnase_merge<-merge(PBMC_BC_samples_filtered, PBMC_Rnase_filtered,by=c("chrom","start","end","strand") )
nrow(PBMC_braincode_rnase_merge)
#1203

############# PBMC Venn diagram

library(VennDiagram)

grid.newpage()
draw.triple.venn(area1 =2236, area2 = 16228, area3 =7765, n12 = 1203, n23 = 4108, n13 = 1671, 
                 n123 = 1092, category = c("BRAINCODE\n (n=2236)", "RNase R \n(n=16,228)", "Mock (n=7,765)"),,  lty = rep("solid", 3), col = rep("white", 3), 
                 fill = c("yellow", "red", "blue"),cat.cex=1.5,cex=2.5)


# Stacked Plots

#obtain log2 FC for intersectin of BC,mock, rnaser

TCPY_braincode_rnase_mock_merge_intersect$MeanRNase_TC=rowMeans(TCPY_braincode_rnase_mock_merge_intersect[,c(14,15)])
TCPY_braincode_rnase_mock_merge_intersect$MeanMock_TC=rowMeans(TCPY_braincode_rnase_mock_merge_intersect[,c(9,10)])
TCPY_braincode_rnase_mock_merge_intersect$Log2FC_R_M_TC<-((log2(TCPY_braincode_rnase_mock_merge_intersect$MeanRNase_TC+pseudocount))-(log2(TCPY_braincode_rnase_mock_merge_intersect$MeanMock_TC+pseudocount)))

TC_higher2log2FC <- TCPY_braincode_rnase_mock_merge_intersect[ TCPY_braincode_rnase_mock_merge_intersect$Log2FC_R_M_TC>2 , ]
nrow(TC_higher2log2FC)
#910

TC_lowermin2log2FC <- TCPY_braincode_rnase_mock_merge_intersect[ TCPY_braincode_rnase_mock_merge_intersect$Log2FC_R_M_TC < (-2) , ]
nrow(TC_lowermin2log2FC)
#24

TC_eqlog2FC <- TCPY_braincode_rnase_mock_merge_intersect[  TCPY_braincode_rnase_mock_merge_intersect$Log2FC_R_M_TC <= 2 & TCPY_braincode_rnase_mock_merge_intersect$Log2FC_R_M_TC >= (-2) , ]
nrow(TC_eqlog2FC )
#229

##SN

SN_braincode_rnase_mock_merge_intersect$MeanRNase_SN=rowMeans(SN_braincode_rnase_mock_merge_intersect[,c(14,15)])
SN_braincode_rnase_mock_merge_intersect$MeanMock_SN=rowMeans(SN_braincode_rnase_mock_merge_intersect[,c(9,10)])
SN_braincode_rnase_mock_merge_intersect$Log2FC_R_M_SN<-((log2(SN_braincode_rnase_mock_merge_intersect$MeanRNase_SN+pseudocount))-(log2(SN_braincode_rnase_mock_merge_intersect$MeanMock_SN+pseudocount)))

SN_higher2log2FC <- SN_braincode_rnase_mock_merge_intersect[ SN_braincode_rnase_mock_merge_intersect$Log2FC_R_M_SN>2 , ]
nrow(SN_higher2log2FC)
#419

SN_lowermin2log2FC <- SN_braincode_rnase_mock_merge_intersect[ SN_braincode_rnase_mock_merge_intersect$Log2FC_R_M_SN < (-2) , ]
nrow(SN_lowermin2log2FC)
#75

SN_eqlog2FC <- SN_braincode_rnase_mock_merge_intersect[  SN_braincode_rnase_mock_merge_intersect$Log2FC_R_M_SN <= 2 & SN_braincode_rnase_mock_merge_intersect$Log2FC_R_M_SN >= (-2) , ]
nrow(SN_eqlog2FC )
#789


##PBMC

PBMC_braincode_rnase_mock_merge_intersect$Log2FC_R_M_PBMC<-((log2(PBMC_braincode_rnase_mock_merge_intersect$HC_WGC082372_PBMC_R_b1+0.01))-(log2(PBMC_braincode_rnase_mock_merge_intersect$HC_WGC082373_PBMC_M_b1+0.01)))

PBMC_higher2log2FC <- PBMC_braincode_rnase_mock_merge_intersect[ PBMC_braincode_rnase_mock_merge_intersect$Log2FC_R_M_PBMC >2 , ]
nrow(PBMC_higher2log2FC)
#898

PBMC_lowermin2log2FC <- PBMC_braincode_rnase_mock_merge_intersect[ PBMC_braincode_rnase_mock_merge_intersect$Log2FC_R_M_PBMC < (-2) , ]
nrow(PBMC_lowermin2log2FC)
#21


PBMC_eqlog2FC <- PBMC_braincode_rnase_mock_merge_intersect[  PBMC_braincode_rnase_mock_merge_intersect$Log2FC_R_M_PBMC <= 2 & PBMC_braincode_rnase_mock_merge_intersect$Log2FC_R_M_PBMC >= (-2) , ]
nrow(PBMC_eqlog2FC)
#173

dat <- read.table(text = "SN   TC   PBMC
                  A 419 910 898
                  B 789 229 173
                  C 75 24 21", header = TRUE)

library(reshape2)

#dat$row <- seq_len(nrow(dat))
#dat$row <- c("> 2 log2 FC",">= -2 and <2 log2 FC","< -2 log2 FC")
dat$row <- c("Enriched by RNAse R","Unaffected by RNAse R","Depleted by RNAse R")
dat2 <- melt(dat, id.vars = "row")

library(ggplot2)

ggplot(dat2, aes(x=variable, y=value, fill=row)) + 
  geom_bar(stat="identity") +
  xlab("\nRegion") +
  ylab("\nNumber of CircRNAs\n") +
  #guides(fill=FALSE) +
  guides(color=guide_legend) +
  labs(fill="")+
  ggtitle("RNAse R Effects in \n BRAINCODE Confirmed circRNAS ") +
  #scale_fill_manual(values=c("red","darkgreen", "gray")) +
  coord_flip()+
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=27),legend.text=element_text(size=15),plot.title = element_text(size=27))
#theme_bw()


#binomial test for INTERSECTION

SN_success=419
SN_total=(419+789+75)

binom.test(SN_success,SN_total, p = 0.5,
           alternative = "greater",
           conf.level = 0.95)

#p-value = 1
#
TN_success=910
TN_total=(910+229+24)

binom.test(TN_success,TN_total, p = 0.5,
           alternative = "greater",
           conf.level = 0.95)

#p-value < 2.2e-16

PBMC_success=898
PBMC_total=(898+173+21)

binom.test(PBMC_success,PBMC_total, p = 0.5,
           alternative = "greater",
           conf.level = 0.95)

#p-value < 2.2e-16

#Do a summary of results

#merge present in RNAse union Mock and BC TCPY

TCPY_braincode_mockUNIONrnase_merge<-merge(TCPY_BC_samples_filtered, TC_pairs_filtered,by=c("chrom","start","end","strand") )
nrow(TCPY_braincode_mockUNIONrnase_merge)
#1300 #other 77 are the ones not present in RNAse or Mock samples

#See which ones enriched, depleted, same, not present
#all intersected with BC
TCPY_higher2log2FC_pairs <-TCPY_braincode_mockUNIONrnase_merge[ TCPY_braincode_mockUNIONrnase_merge$Log2FC_R_M_TC>2 , ]
nrow(TCPY_higher2log2FC_pairs)
#985

TCPY_eqlog2FC <- TCPY_braincode_mockUNIONrnase_merge[  TCPY_braincode_mockUNIONrnase_merge$Log2FC_R_M_TC <= 2 & TCPY_braincode_mockUNIONrnase_merge$Log2FC_R_M_TC >= (-2) , ]
nrow(TCPY_eqlog2FC )
#233

TCPY_lowermin2log2FC_pairs <-TCPY_braincode_mockUNIONrnase_merge[ TCPY_braincode_mockUNIONrnase_merge$Log2FC_R_M_TC < (-2) , ]
nrow(TCPY_lowermin2log2FC_pairs)
#82


#SN

SN_braincode_mockUNIONrnase_merge<-merge(SN_BC_samples_filtered, SN_pairs_filtered,by=c("chrom","start","end","strand") )
nrow(SN_braincode_mockUNIONrnase_merge)
#1771 #other 655 are the ones not present in RNAse or Mock samples

#See which ones enriched, depleted, same, not present
#all intersected with BC
SN_higher2log2FC_pairs <-SN_braincode_mockUNIONrnase_merge[ SN_braincode_mockUNIONrnase_merge$Log2FC_R_M_SN>2 , ]
nrow(SN_higher2log2FC_pairs)
#639


SN_eqlog2FC <- SN_braincode_mockUNIONrnase_merge[  SN_braincode_mockUNIONrnase_merge$Log2FC_R_M_SN <= 2 & SN_braincode_mockUNIONrnase_merge$Log2FC_R_M_SN >= (-2) , ]
nrow(SN_eqlog2FC )
#824

SN_lowermin2log2FC_pairs <-SN_braincode_mockUNIONrnase_merge[ SN_braincode_mockUNIONrnase_merge$Log2FC_R_M_SN < (-2) , ]
nrow(SN_lowermin2log2FC_pairs)
#308



#PBMC

PBMC_braincode_mockUNIONrnase_merge<-merge(PBMC_BC_samples_filtered, PBMC_pairs_filtered,by=c("chrom","start","end","strand") )
nrow(PBMC_braincode_mockUNIONrnase_merge)
#1782 #other 454 are the ones not present in RNAse or Mock samples

#See which ones enriched, depleted, same, not present
#all intersected with BC
PBMC_higher2log2FC_pairs <-PBMC_braincode_mockUNIONrnase_merge[ PBMC_braincode_mockUNIONrnase_merge$Log2FC_R_M_PBMC>2 , ]
nrow(PBMC_higher2log2FC_pairs)
#1013

PBMC_eqlog2FC <- PBMC_braincode_mockUNIONrnase_merge[  PBMC_braincode_mockUNIONrnase_merge$Log2FC_R_M_PBMC <= 2 & PBMC_braincode_mockUNIONrnase_merge$Log2FC_R_M_PBMC >= (-2) , ]
nrow(PBMC_eqlog2FC)
#169

PBMC_lowermin2log2FC_pairs <-PBMC_braincode_mockUNIONrnase_merge[ PBMC_braincode_mockUNIONrnase_merge$Log2FC_R_M_PBMC < (-2) , ]
nrow(PBMC_lowermin2log2FC_pairs)
#600


#Plot summary of these results

sum <- read.table(text = "SN   TC   PBMC
                  A 639 985 1013
                  B 824 233 169
                  C 308 82 600
                  D 655 77 454", header = TRUE)


library(reshape2)

#dat$row <- seq_len(nrow(dat))
#dat$row <- c("> 2 log2 FC",">= -2 and <2 log2 FC","< -2 log2 FC")
sum$row <- c("Enriched by RNAse R","Unaffected by RNAse R","Depleted by RNAse R","Not Expressed in RNAse R or Mock")
sum2 <- melt(sum, id.vars = "row")

library(ggplot2)

ggplot(sum2, aes(x=variable, y=value, fill=row)) + 
  geom_bar(stat="identity") +
  xlab("\nRegion") +
  ylab("\nNumber of CircRNAs\n") +
  #guides(fill=FALSE) +
  guides(color=guide_legend) +
  labs(fill="")+
  ggtitle("RNAse R Effects in all BRAINCODE  \n Cell Type Expressed circRNAS ") +
  scale_fill_manual(values=c("coral1","springgreen3", "orange","steelblue1")) +
  coord_flip()+
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=27),legend.text=element_text(size=15),plot.title = element_text(size=27))

####Here
#binomial test for all BRAINCODE

SN_success=639
SN_total=(639+824+308+655)

binom.test(SN_success,SN_total, p = 0.5,
           alternative = "greater",
           conf.level = 0.95)

#p-value = 1


TC_success=985
TC_total=(985+233+82+77)

binom.test(TC_success,TC_total, p = 0.5,
           alternative = "greater",
           conf.level = 0.95)

#p-value < 2.2e-16 

PBMC_success=1013
PBMC_total=(1013+169+600+454)

binom.test(PBMC_success,PBMC_total, p = 0.5,
           alternative = "greater",
           conf.level = 0.95)

#p-value = 1

#theme_bw()

#Fischer's Test for Enrichment

#BC total= 1169

##Finished repetition here ##

####Visualization

#Select Most highly expressed of enriched, and view in genome browser plus view a fiew reads. Also see what 
#else has been published on this

head(TCPY_higher2log2FC_pairs[order(-TCPY_higher2log2FC_pairs$MeanRNase_TC),])

#chrom    start      end strand      name.x score.x columnGreat.x     meanRPM HC_WGC082366_TC_R_b1 HC_WGC082367_TC_M_b1
#254  chr11 46098304 46113774      - circ_110515    1350             7 0.176583496            419.96312             8.048388
#1064  chr6 54013853 54095715      + circ_124241     353             5 0.226805792            199.66032             3.711795
#712   chr2 72945231 72960247      - circ_165727     458             4 0.019080940            106.48195             1.985379
#458  chr15 59204761 59209198      -  circ_90756     123             1 0.007774251            125.69931             1.644206
#835   chr3 47079155 47088111      -  circ_26755     563             1 0.002329606            117.94111             1.981268
#520  chr17 26449628 26499644      +  circ_42796     646             2 0.024706750             92.92388             2.606066
#HC_WGC082368_TC_R_b2 HC_WGC082369_TC_M_b2      name.y score.y columnGreat.y MeanRNase_TC MeanMock_TC Log2FC_R_M_TC
#254              453.3160            12.630941 circ_149586  446993             4     436.6396   10.339665      5.399933
#1064             372.1208            11.650319  circ_31267  287545             4     285.8906    7.681057      5.217681
#712              244.1783             6.387065 circ_132921  197971             4     175.3301    4.186222      5.387668
#458              198.2953             4.629754  circ_50482  359595             4     161.9973    3.136980      5.689628
#835              186.4187             3.627436 circ_127463  175688             4     152.1799    2.804352      5.761046
#520              190.5148             6.695136 circ_158814  187461             4     141.7194    4.650601      4.928929
#MeanAll_TC
#254   223.48961
#1064  146.78582
#712    89.75818
#458    82.56714
#835    77.49212
#520    73.18498

head(SN_higher2log2FC_pairs[order(-SN_higher2log2FC_pairs$MeanRNase_SN),])
#chrom     start       end strand      name.x score.x columnGreat.x     meanRPM HC_WGC082362_SN_R_b1_r1
#157  chr1  65830317  65831879      +  circ_86835     540            20 0.020721624               0.7084684
#473 chr13  78293666  78327493      +  circ_97514    2534            29 0.052399581               2.2828427
#123  chr1  32381495  32385259      -   circ_2421    2186            32 0.038430763               2.9575745
#310 chr11  46098304  46113774      - circ_110515    1350            20 0.029576386               0.0000000
#888  chr2 230723487 230744844      - circ_165224     723             8 0.011575307               1.9792134
#653 chr16  56385295  56388993      +  circ_30021     292             7 0.008985813               1.6306019
#HC_WGC082363_SN_M_b1_r1 HC_WGC082364_SN_R_b2_r1 HC_WGC082365_SN_M_b2_r1      name.y score.y columnGreat.y MeanRNase_SN
#157               0.7565177                53.44362                2.909665 circ_155085   70426             4     27.07604
#473               3.1126948                48.03782                3.853621  circ_11549   47622             4     25.16033
#123               2.5871749                32.86776                3.417949   circ_4290  142054             4     17.91267
#310               0.0000000                35.34594                1.338135 circ_149586  446993             2     17.67297
#888               1.8191074                30.76482                1.768620 circ_132032   68025             4     16.37202
#653               0.9008913                27.53618                1.664889  circ_54025   33685             4     14.58339
#MeanMock_SN Log2FC_R_M_SN MeanAll_SN
#157   1.8330914      3.883329  14.454567
#473   3.4831578      2.852033  14.321744
#123   3.0025619      2.575985  10.457615
#310   0.6690674      4.719472   9.171018
#888   1.7938637      3.188785   9.082940
#653   1.2828899      3.504989   7.933141

head(PBMC_higher2log2FC_pairs[order(-PBMC_higher2log2FC_pairs$MeanRNase_PBMC),])
#chrom    start      end strand      name.x score.x columnGreat.x    meanRPM HC_WGC082372_PBMC_R_b1
#616  chr15 59204761 59209198      -  circ_90756     123             1 0.07877683               354.6904
#1007 chr20  2944917  2945848      +  circ_24400     776             4 0.55155845               267.7768
#1008 chr20 32207322 32211102      + circ_149041     563             2 0.04023203               178.0791
#364  chr11 85707868 85714494      -  circ_56926    1126             4 0.22268973               156.3092
#281  chr10 98703869 98711953      +  circ_93724     447             3 0.28524826               137.2596
#474  chr13 28748408 28752072      +   circ_2074     100             3 0.06985378               132.6358
#HC_WGC082373_PBMC_M_b1      name.y score.y columnGreat.y MeanRNase_PBMC MeanMock_PBMC Log2FC_R_M_PBMC MeanAll_PBMC
#616               12.391063  circ_50482  359595             2       354.6904     12.391063        4.838984    183.54073
#1007               7.889758  circ_43926  206598             2       267.7768      7.889758        5.084582    137.83328
#1008               3.892480 circ_103554  207228             2       178.0791      3.892480        5.515023     90.98577
#364                4.361574 circ_102142  136091             2       156.3092      4.361574        5.162824     80.33538
#281                3.907451   circ_4839  102014             2       137.2596      3.907451        5.133882     70.58354
#474                3.054100   circ_3622   84135             2       132.6358      3.054100        5.439738     67.84495


PBMC_higher2log2FC_pairs[grep("circ_110515", PBMC_higher2log2FC_pairs$name.x), ] #confirm unique

#chrom    start      end strand      name.x score.x columnGreat.x    meanRPM HC_WGC082372_PBMC_R_b1
#334 chr11 46098304 46113774      - circ_110515    1350             1 0.05079652               26.37364
#HC_WGC082373_PBMC_M_b1      name.y score.y columnGreat.y MeanRNase_PBMC MeanMock_PBMC Log2FC_R_M_PBMC MeanAll_PBMC
#334                      0 circ_149586  446993             1       26.37364             0        13.82038     13.18682



#Plot RPM and FC for first circRNA

barplot(c(4.719,5.399,13.82),col=c("orange","dodgerblue","purple"),names.arg=c("SN", "TC", "PBMC"), ylab="Mean Log2 FC (RNase R / Mock)", xlab="Brain Region",main="PHF21A circRNA RNase R \n Expression Changes")

barplot(c(17.67,0.669,436.6396,10.339,26.37,pseudocount),col=c("orange","sandybrown","dodgerblue","paleturquoise1","purple","plum"),names.arg=c("SN\n RNase R","SN\n Mock", "TC \nRNase R","TC\n Mock", "PBMC \nRNase R", "PBMC \n Mock"), ylab="Mean Reads Per Million Raw Reads (RPM)", xlab="Brain Region",main="PHF21A circRNA Expression")

##########################################################################################
##########################################################################################
################################   Obtain tables:    #####################################
##########################################################################################
##########################################################################################

# 5 tables
#Table 1. Must have: chr,start,end,id(chr_start_end_minus/plus),strand,host_gene for union of BC/RNASE/MC

#Table 2. Raw reads and RPM for all circRNAs in Union of circexplorer output (not filtered)

#Table 3,4, and 5. for each tissue- id, mean_m, mean_r, log2 FC

#Step 1. Obtain the union of BC/R/M circRNAs 

setwd("~/projects/circRNA")
load("Merge_circexp.BC.Rdata")
load("Merge_circexp.RNAse.Rdata")

All_circexpOut_BC_RM<- read.table("concatenated_circOut_unique_BC_RM.txt",sep="\t",header=FALSE) 
nrow(All_circexpOut_BC_RM)
#[1] 380201
names(All_circexpOut_BC_RM)<- c("chrom","start","end","strand","exonNum","ExonSize","ciRNA","GeneName")

#check for duplicates
sum(duplicated(All_circexpOut_BC_RM))
#[1] 0

#create column that says m (minus) or p (plus) for each row to do ID
All_circexpOut_BC_RM$strand_letter <- ifelse(All_circexpOut_BC_RM$strand=="+", "p", "m")

#create IDs
All_circexpOut_BC_RM <- within(All_circexpOut_BC_RM,  id <- paste(chrom, start, end, strand_letter, sep="_"))

#remove extra column 
All_circexpOut_BC_RM$strand_letter<- NULL

#ciRNA or circRNA

All_circexpOut_BC_RM$type <- ifelse(All_circexpOut_BC_RM$ciRNA=="Yes", "ciRNA", "circRNA")

#remove extra column 
All_circexpOut_BC_RM$ciRNA<- NULL


#Table 2. 

#create column that says m (minus) or p (plus) for each row to do ID
Merge_circexp$strand_letter <- ifelse(Merge_circexp$strand=="+", "p", "m")
Merge_circexp_RNAse$strand_letter <- ifelse(Merge_circexp_RNAse$strand=="+", "p", "m")
Merge_circexp_norm$strand_letter <- ifelse(Merge_circexp_norm$strand=="+", "p", "m")
Merge_circexp_RNAse_norm$strand_letter <- ifelse(Merge_circexp_RNAse_norm$strand=="+", "p", "m")

#create IDs
Merge_circexp <- within(Merge_circexp,  id <- paste(chrom, start, end, strand_letter, sep="_"))
Merge_circexp_RNAse <- within(Merge_circexp_RNAse,  id <- paste(chrom, start, end, strand_letter, sep="_"))
Merge_circexp_norm <- within(Merge_circexp_norm,  id <- paste(chrom, start, end, strand_letter, sep="_"))
Merge_circexp_RNAse_norm <- within(Merge_circexp_RNAse_norm,  id <- paste(chrom, start, end, strand_letter, sep="_"))

#remove extra column 
Merge_circexp$strand_letter<- NULL
Merge_circexp_RNAse$strand_letter<- NULL
Merge_circexp_norm$strand_letter<- NULL
Merge_circexp_RNAse_norm$strand_letter<- NULL

#Change column names to Raw for raw counts
colnames(Merge_circexp)[7:(ncol(Merge_circexp)-1)] <- paste(colnames(Merge_circexp[7:(ncol(Merge_circexp)-1)]) , "Raw", sep = "_") 
colnames(Merge_circexp_RNAse)[7:(ncol(Merge_circexp_RNAse)-1)] <- paste(colnames(Merge_circexp_RNAse[7:(ncol(Merge_circexp_RNAse)-1)]) , "Raw", sep = "_") 


#Change column names to RPM for normalized RPM counts

colnames(Merge_circexp_norm)[1:81] <- paste(colnames(Merge_circexp_norm[1:81]) , "RPM", sep = "_") 
colnames(Merge_circexp_RNAse_norm)[1:12] <- paste(colnames(Merge_circexp_RNAse_norm[1:12]) , "RPM", sep = "_") 



Merge_circexp_Raw_RPM<-merge(Merge_circexp, Merge_circexp_norm,by=c("chrom","start","end","name","strand","id") )

Merge_circexp_Raw_RPM$score.x <- NULL
Merge_circexp_Raw_RPM$score.y <- NULL
Merge_circexp_Raw_RPM$columnGreat <- NULL

Merge_circexp_RNaser_Raw_RPM<-merge(Merge_circexp_RNAse, Merge_circexp_RNAse_norm,by=c("chrom","start","end","name","strand","id") )

Merge_circexp_RNaser_Raw_RPM$score.x <- NULL
Merge_circexp_RNaser_Raw_RPM$score.y <- NULL
Merge_circexp_RNaser_Raw_RPM$columnGreat <- NULL


Merge_circexp_BCMR_Raw_RPM<- merge(Merge_circexp_Raw_RPM, Merge_circexp_RNaser_Raw_RPM,by=c("chrom","start","end","strand","id") )
#intersect is  26632

Merge_circexp_BCMR_Raw_RPM<- merge(Merge_circexp_Raw_RPM, Merge_circexp_RNaser_Raw_RPM,by=c("chrom","start","end","strand","id"), all =TRUE)
nrow(Merge_circexp_BCMR_Raw_RPM)
[1] 380201
#Union should be 380201 circRNAs 

Merge_circexp_BCMR_Raw_RPM$name.x <- NULL
Merge_circexp_BCMR_Raw_RPM$name.y <- NULL

#replace NA with zero
Merge_circexp_BCMR_Raw_RPM[is.na(Merge_circexp_BCMR_Raw_RPM)] <- 0

#Merge Table 1 and 2 to make sure they're the same and are in order

Table_anno_expr<- merge(All_circexpOut_BC_RM, Merge_circexp_BCMR_Raw_RPM,by=c("chrom","start","end","strand","id"))
nrow(Table_anno_expr)
[1] 380201

Table_annotation<-Table_anno_expr[,c(5,1,2,3,4,8,9,6,7)]

Table_expression<- Table_anno_expr[,c(5,10:ncol(Table_anno_expr))]

#write annotation and expression table

write.table(Table_annotation, "Table_annotation.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(Table_expression, "Table_expression.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

#Select Tissue from expression table
#SN
SN_pairs_nofilt<-(Table_expression[grepl( "SN_R|SN_M|id" , names(Table_expression))])
SN_pairs_nofilt<-(SN_pairs_nofilt[grepl( "RPM|id" , names(SN_pairs_nofilt))])

#Obtain log2 FC
SN_pairs_nofilt$MeanRNase=rowMeans(SN_pairs_nofilt[,c(2,4)])
SN_pairs_nofilt$MeanMock=rowMeans(SN_pairs_nofilt[,c(3,5)])
SN_pairs_nofilt$Log2FC_R_M<-((log2(SN_pairs_nofilt$MeanRNase+pseudocount))-(log2(SN_pairs_nofilt$MeanMock+pseudocount)))

#order by RNAse RPM and then FC
SN_pairs_nofilt<-(SN_pairs_nofilt[order(-SN_pairs_nofilt$MeanRNase,-SN_pairs_nofilt$Log2FC_R_M),])
write.table(SN_pairs_nofilt, "SN_BC_R_M_union_FoldChange.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

#TC
TC_pairs_nofilt<-(Table_expression[grepl( "TC_R|TC_M|id" , names(Table_expression))])
TC_pairs_nofilt<-(TC_pairs_nofilt[grepl( "RPM|id" , names(TC_pairs_nofilt))])

#Obtain log2 FC
TC_pairs_nofilt$MeanRNase=rowMeans(TC_pairs_nofilt[,c(2,4)])
TC_pairs_nofilt$MeanMock=rowMeans(TC_pairs_nofilt[,c(3,5)])
TC_pairs_nofilt$Log2FC_R_M<-((log2(TC_pairs_nofilt$MeanRNase+pseudocount))-(log2(TC_pairs_nofilt$MeanMock+pseudocount)))

#order by RNAse RPM and then FC
TC_pairs_nofilt<-(TC_pairs_nofilt[order(-TC_pairs_nofilt$MeanRNase,-TC_pairs_nofilt$Log2FC_R_M),])
write.table(TC_pairs_nofilt, "TC_BC_R_M_union_FoldChange.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

#PBMC
PBMC_pairs_nofilt<-(Table_expression[grepl( "PBMC_R|PBMC_M|id" , names(Table_expression))])
PBMC_pairs_nofilt<-(PBMC_pairs_nofilt[grepl( "RPM|id" , names(PBMC_pairs_nofilt))])

#Obtain log2 FC
PBMC_pairs_nofilt$MeanRNase=PBMC_pairs_nofilt[,2] # PBMC_pairs_nofilt$HC_WGC082372_PBMC_R_b1_RPM
PBMC_pairs_nofilt$MeanMock=PBMC_pairs_nofilt[,3] # PBMC_pairs_nofilt$HC_WGC082373_PBMC_M_b1_RPM
PBMC_pairs_nofilt$Log2FC_R_M<-((log2(PBMC_pairs_nofilt$MeanRNase+pseudocount))-(log2(PBMC_pairs_nofilt$MeanMock+pseudocount)))

#order by RNAse RPM and then FC
PBMC_pairs_nofilt<-(PBMC_pairs_nofilt[order(-PBMC_pairs_nofilt$MeanRNase,-PBMC_pairs_nofilt$Log2FC_R_M),])
write.table(PBMC_pairs_nofilt, "PBMC_BC_R_M_union_FoldChange.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

#FB
FB_pairs_nofilt<-(Table_expression[grepl( "FB_R|FB_M|id" , names(Table_expression))])
FB_pairs_nofilt<-(FB_pairs_nofilt[grepl( "RPM|id" , names(FB_pairs_nofilt))])

#Obtain log2 FC  
FB_pairs_nofilt$MeanRNase=FB_pairs_nofilt[,2]
FB_pairs_nofilt$MeanMock=FB_pairs_nofilt[,3]
FB_pairs_nofilt$Log2FC_R_M<-((log2(FB_pairs_nofilt$MeanRNase+pseudocount))-(log2(FB_pairs_nofilt$MeanMock+pseudocount)))

#order by RNAse RPM and then FC
FB_pairs_nofilt<-(FB_pairs_nofilt[order(-FB_pairs_nofilt$MeanRNase,-FB_pairs_nofilt$Log2FC_R_M),])
write.table(FB_pairs_nofilt, "FB_BC_R_M_union_FoldChange.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

