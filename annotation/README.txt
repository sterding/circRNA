# cross-sample normalize the RPM, for RNase and Mock separately
# R
setwd("~/Google\ Drive/circRNA")
df=list()
for(i in c('SN', 'PBMC', 'TC', 'FB')){
    print(i) # i='SN'
    df[[i]]=read.table(paste('Table_annotation.sorted.by.meanRNase.in',i,'xls',sep="."))
}

# normalize to the minimal one
normalized_factor_RNase = sapply(c('SN', 'PBMC', 'TC', 'FB'), function(x) mean(df[[x]]$MeanRNase)) / min(sapply(c('SN', 'PBMC', 'TC', 'FB'), function(x) mean(df[[x]]$MeanRNase)))
normalized_factor_Mock  = sapply(c('SN', 'PBMC', 'TC', 'FB'), function(x) mean(df[[x]]$MeanMock)) / min(sapply(c('SN', 'PBMC', 'TC', 'FB'), function(x) mean(df[[x]]$MeanMock)))

PSEU=0.0158 #?
for(i in c('SN', 'PBMC', 'TC', 'FB')){
    print(i)
    df[[i]]$MeanRNase = df[[i]]$MeanRNase / normalized_factor_RNase[i]
    df[[i]]$MeanMock = df[[i]]$MeanMock / normalized_factor_Mock[i]
    df[[i]]$Log2FC_R_M = log2((df[[i]]$MeanRNase + PSEU) / (df[[i]]$MeanMock + PSEU))
    # write
    write.table(df[[i]],file=paste('Table_annotation.sorted.by.meanRNase.in',i,'norm.xls',sep="."), quote=F,sep="\t",na="",col.names = NA,row.names = TRUE)
}


## intersect PD associated genes
cd ~/Google\ Drive/circRNA/
dos2unix PD_associated_genes.txt; for i in SN PBMC TC FB; do echo $i; awk '$11>=0.015' Table_annotation.sorted.by.meanRNase.in.$i.norm.xls | grep -wf <(sed 's/,/\n/g;s/ //g;s/(/\n/g;s/)//g;/^$/d' PD_associated_genes.txt | sort -u) > PD_associated_genes.circRNA.$i.txt; done