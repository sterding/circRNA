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

setwd("~/projects/circRNA/data")
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
#intersect is  20305

Merge_circexp_BCMR_Raw_RPM<- merge(Merge_circexp_Raw_RPM, Merge_circexp_RNaser_Raw_RPM,by=c("chrom","start","end","strand","id"), all =TRUE)
nrow(Merge_circexp_BCMR_Raw_RPM)
#[1] 371620
#Union should be 371620 circRNAs 

Merge_circexp_BCMR_Raw_RPM$name.x <- NULL
Merge_circexp_BCMR_Raw_RPM$name.y <- NULL

#replace NA with zero
Merge_circexp_BCMR_Raw_RPM[is.na(Merge_circexp_BCMR_Raw_RPM)] <- 0

#Merge Table 1 and 2 to make sure they're the same and are in order

Table_anno_expr<- merge(All_circexpOut_BC_RM, Merge_circexp_BCMR_Raw_RPM,by=c("chrom","start","end","strand","id"))
nrow(Table_anno_expr)
#[1] 371620

Table_annotation<-Table_anno_expr[,c(5,1,2,3,4,8,9,6,7)]
Table_expression<- Table_anno_expr[,c(5,10:ncol(Table_anno_expr))]

table(Table_annotation$type)
# circRNA   ciRNA 
#  97031  274589

dim(Table_expression[rowSums(Table_expression[,grep("_Raw",colnames(Table_expression))])==0,])

#write annotation and expression table

write.table(Table_annotation, "Table_annotation.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(Table_expression, "Table_expression.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

#Select Tissue from expression table
# =====================================
pseudocount=0.001823273
#SN
SN_pairs_nofilt<-(Table_expression[grepl( "SN_R|SN_M|id" , names(Table_expression))])
SN_pairs_nofilt<-(SN_pairs_nofilt[grepl( "RPM|id" , names(SN_pairs_nofilt))])

# the SN sample1 RNase_R is not working well, as there is no enrichment on circRNA. So we just remove it when counting the fold change
colSums(Table_expression[,grep("HC_.*_SN_.*_Raw", colnames(Table_expression))])
# HC_WGC082362_SN_R_b1_r1_Raw HC_WGC082363_SN_M_b1_r1_Raw HC_WGC082364_SN_R_b2_r1_Raw HC_WGC082365_SN_M_b2_r1_Raw 
# 101893                      100030                      896267                      122185

#Obtain log2 FC
SN_pairs_nofilt$MeanRNase=SN_pairs_nofilt[,4] # rowMeans(SN_pairs_nofilt[,c(2,4)])
SN_pairs_nofilt$MeanMock=SN_pairs_nofilt[,5] # rowMeans(SN_pairs_nofilt[,c(3,5)])
SN_pairs_nofilt$Log2FC_R_M<-((log2(SN_pairs_nofilt$MeanRNase+pseudocount))-(log2(SN_pairs_nofilt$MeanMock+pseudocount)))

#order by RNAse RPM and then FC
SN_pairs_nofilt<-(SN_pairs_nofilt[order(-SN_pairs_nofilt$MeanRNase,-SN_pairs_nofilt$Log2FC_R_M),])
write.table(SN_pairs_nofilt, "SN_BC_R_M_union_FoldChange.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

#TC
TC_pairs_nofilt<-(Table_expression[grepl( "TC_R|TC_M|id" , names(Table_expression))])
TC_pairs_nofilt<-(TC_pairs_nofilt[grepl( "RPM|id" , names(TC_pairs_nofilt))])

colSums(Table_expression[,grep("HC_.*_TC_.*_Raw", colnames(Table_expression))])  # not the case for TC samples

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

# sort Annotation according to mean RNase values in each cell type [For Yunfei]
write.table(cbind(Table_annotation[match(SN_pairs_nofilt$id,Table_annotation$id),],SN_pairs_nofilt[,c('MeanRNase', 'MeanMock', 'Log2FC_R_M')]), "~/Google\ Drive/circRNA/Table_annotation.sorted.by.meanRNase.in.SN.xls", sep="\t", quote=FALSE, row.names=1:nrow(Table_annotation), col.names=NA)
write.table(cbind(Table_annotation[match(TC_pairs_nofilt$id,Table_annotation$id),],TC_pairs_nofilt[,c('MeanRNase', 'MeanMock', 'Log2FC_R_M')]), "~/Google\ Drive/circRNA/Table_annotation.sorted.by.meanRNase.in.TC.xls", sep="\t", quote=FALSE, row.names=1:nrow(Table_annotation), col.names=NA)
write.table(cbind(Table_annotation[match(PBMC_pairs_nofilt$id,Table_annotation$id),],PBMC_pairs_nofilt[,c('MeanRNase', 'MeanMock', 'Log2FC_R_M')]), "~/Google\ Drive/circRNA/Table_annotation.sorted.by.meanRNase.in.PBMC.xls", sep="\t", quote=FALSE, row.names=1:nrow(Table_annotation), col.names=NA)
write.table(cbind(Table_annotation[match(FB_pairs_nofilt$id,Table_annotation$id),],PBMC_pairs_nofilt[,c('MeanRNase', 'MeanMock', 'Log2FC_R_M')]), "~/Google\ Drive/circRNA/Table_annotation.sorted.by.meanRNase.in.FB.xls", sep="\t", quote=FALSE, row.names=1:nrow(Table_annotation), col.names=NA)
