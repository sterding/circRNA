##========================================================================
# 4. get circulization ratio == back-spliced reads / (back-spliced reads + linear reads)
##========================================================================
args<-commandArgs(TRUE)

circRNA_annotation=args[1] # circRNA_annotation="Merge_circexplorer_BC.annotation.bed14"
sample_list_subset=args[2] # "~/neurogen/rnaseq_PD/results/merged/samplelist.HCILB_SNDA84"

#message(paste0("test|", file_test("-f",sample_list_subset),"|"))
#quit('no')

setwd("~/projects/circRNA/data/")

# read data
nC=read.table(paste0(circRNA_annotation, ".sum_u3u5s3s5"), header=T, check.names = F, row.names = 1, stringsAsFactors=F)
dim(nC)
nL=read.table(paste0(circRNA_annotation, ".circReads.txt"), header=T, check.names = F, row.names = 1, stringsAsFactors=F)
dim(nL)
common_Cols=intersect(colnames(nC), colnames(nL))
if(file_test("-f",sample_list_subset)) common_Cols=intersect(common_Cols, scan(sample_list_subset,character()))
common_Rows=intersect(rownames(nC), rownames(nL))
nC=nC[common_Rows, common_Cols]
nL=nL[common_Rows, common_Cols]
nSum=nC+nL
dim(nSum)
## filter and transformation learned from LeafCutter (https://github.com/davidaknowles/leafcutter/blob/master/scripts/prepare_phenotype_table.py)
# 1. If ratio is missing for over 50% of the samples, skip the row
nC=nC[rowMeans(nSum==0)<=.5,]; nSum=nSum[rowMeans(nSum==0)<=.5,]
dim(nSum)
# 2. ratio
rC = (nC+0.5) / (nSum+0.5)  # add a 0.5 pseudocount
# 3. Set missing values as the mean of all values
nSum_1NA=ifelse(nSum==0, NA, 1)
rC_NA=nSum_1NA * rC  # now rC_NA has NA for the cell where nSum is 0
impute.mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))  # https://stackoverflow.com/a/9322975
rC = t(apply(rC_NA, 1, impute.mean))
dim(rC)
# 4. If there is too little variation, skip (it might be a bug of fastqtl which doesn't handle cases with no variation)
rC = rC[apply(rC, 1, sd, na.rm = TRUE) >= 0.005, ]
dim(rC)
# 5. scale normalize on the rows
rC = t(scale(t(rC)))
# 6. qqnorms on the columns
rC = apply(rC, 2, rank, ties.method = "average")
a=ifelse(nrow(rC)<=10, 3/8, 0.5) # from LeafCutter. No idea why
rC = qnorm((rC -a) / (nrow(rC)+1-2*a));  # to standard normalization

# save for eQTL
saveRDS(rC, file=paste0(circRNA_annotation, ".cRatio.rds"))
saveRDS(nC[rownames(rC),], file=paste0(circRNA_annotation, ".nC.rds"))
saveRDS(nL[rownames(rC),], file=paste0(circRNA_annotation, ".nL.rds"))

message("The following three files are generated:")
message(paste0(circRNA_annotation,c(".cRatio.rds",".nC.rds",".nL.rds"), collapse="\n"))