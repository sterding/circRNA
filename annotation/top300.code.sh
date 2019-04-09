# sample1 is not properly done, so we changed to use the sample2 RNase_R value for sorting, instead of MeanRNase
#sort -k6,6gr SN_BC_R_M_union_FoldChange.txt | head -n300 | awk '{OFS="\t"; print NR,$1}' > top300.by.meanRNase.in.SN_BC_R_M_union_FoldChange.txt

./rowsToCols Table_expression.txt stdout | grep _SN_ | grep Raw | ./rowsToCols stdin stdout  | grep -v HC | datamash sum 1 sum 2 sum 3 sum 4
101893	100030	896267	122185  # the sample1 RNase_R is not working well, as no enrichment on circRNA

mv SN_BC_R_M_union_FoldChange.txt SN_BC_R_M_union_FoldChange.old.txt
awk '$7==0' SN_BC_R_M_union_FoldChange.old.txt | head # calculate the pseudocount Rebeca added
awk '{OFS="\t"; if(NR>1) {seu=0.001823273; $6=$4;$7=$5;$8=log((seu+$4)/(seu+$5))/log(2);} print $1,$4,$5,$6,$7,$8}' SN_BC_R_M_union_FoldChange.old.txt > SN_BC_R_M_union_FoldChange.txt

sort -k4,4gr SN_BC_R_M_union_FoldChange.txt | head -n300 | awk '{OFS="\t"; print NR,$1}' > top300.by.meanRNase.in.SN_BC_R_M_union_FoldChange.txt
fgrep -wf <(cut -f2 top300.by.meanRNase.in.SN_BC_R_M_union_FoldChange.txt) Table_annotation.txt | sort -k1,1 | paste <(sort -k2,2 top300.by.meanRNase.in.SN_BC_R_M_union_FoldChange.txt) - | sort -k1,1n | cut -f1,3- > top300.by.meanRNase.in.SN_BC_R_M_union_FoldChange.annotation.txt

sort -k6,6gr TC_BC_R_M_union_FoldChange.txt | head -n300 | awk '{OFS="\t"; print NR,$1}' > top300.by.meanRNase.in.TC_BC_R_M_union_FoldChange.txt
fgrep -wf <(cut -f2 top300.by.meanRNase.in.TC_BC_R_M_union_FoldChange.txt) Table_annotation.txt | sort -k1,1 | paste <(sort -k2,2 top300.by.meanRNase.in.TC_BC_R_M_union_FoldChange.txt) - | sort -k1,1n | cut -f1,3- > top300.by.meanRNase.in.TC_BC_R_M_union_FoldChange.annotation.txt

sort -k4,4gr PBMC_BC_R_M_union_FoldChange.txt | head -n300 | awk '{OFS="\t"; print NR,$1}' > top300.by.meanRNase.in.PBMC_BC_R_M_union_FoldChange.txt
fgrep -wf <(cut -f2 top300.by.meanRNase.in.PBMC_BC_R_M_union_FoldChange.txt) Table_annotation.txt | sort -k1,1 | paste <(sort -k2,2 top300.by.meanRNase.in.PBMC_BC_R_M_union_FoldChange.txt) - | sort -k1,1n | cut -f1,3- > top300.by.meanRNase.in.PBMC_BC_R_M_union_FoldChange.annotation.txt


# TO RUN top 1000
for i in SN TC PBMC FB; do 
    echo $i;
    sort -k4,4gr ${i}_BC_R_M_union_FoldChange.txt | awk '{OFS="\t"; print NR,$1}' > sorted.by.meanRNase.in.${i}_BC_R_M_union_FoldChange.txt    
    fgrep -wf <(cut -f2 sorted.by.meanRNase.in.${i}_BC_R_M_union_FoldChange.txt) Table_annotation.txt | sort -k1,1 | paste <(sort -k2,2 sorted.by.meanRNase.in.${i}_BC_R_M_union_FoldChange.txt) - | sort -k1,1n | cut -f1,3- > sorted.by.meanRNase.in.${i}_BC_R_M_union_FoldChange.annotation.txt
done

