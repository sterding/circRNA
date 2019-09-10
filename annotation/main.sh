################################
## get all circular RNAs from circExlorer3 for all noon-CSF BRAINcode2 samples (N=221)
################################
# $ wc -l samplelist.braincode2.*
#    87 samplelist.braincode2.CSF.N87.txt
#   221 samplelist.braincode2.woCSF.allN221.txt
#   197 samplelist.braincode2.woCSF.selectedN197.txt

cd ~/projects/circRNA/data

## for CSF
> Merge_circexplorer_CSF87.rawcount.long.txt
cat samplelist.braincode2.CSF.N87.txt | while read i; do echo $i; awk -vi=$i '{OFS="\t"; print $1"_"$2"_"$3, $13, i}' ~/neurogen/rnaseq_CSF/run_output/$i/circularRNA_known3.txt >> Merge_circexplorer_CSF87.rawcount.long.txt; done

echo -e "chrom\tstart\tend\tID\tscore\tstrand\tthickStart\tthickEnd\titemRgb\texonCount\texonSizes\texonOffsets\tcircType\tgeneID" > Merge_circexplorer_CSF87.annotation.bed14
cat `ls ~/neurogen/rnaseq_CSF/run_output/*/circularRNA_known3.txt | grep -f samplelist.braincode2.CSF.N87.txt -` | awk '{OFS="\t"; $5=0;$4=$1"_"$2"_"$3; print}' | cut -f1-12,14-15 | sort -u >> Merge_circexplorer_CSF87.annotation.bed14
wc -l Merge_circexplorer_CSF87.annotation.bed14
# 42,424

## for BC
> Merge_circexplorer_BC221.rawcount.long.txt
cat samplelist.braincode2.woCSF.allN221.txt | while read i; do echo $i; awk -vi=$i '{OFS="\t"; print $1"_"$2"_"$3, $13, i}' ~/neurogen/rnaseq_PD/run_output/$i/circularRNA_known3.txt >> Merge_circexplorer_BC221.rawcount.long.txt; done

echo -e "chrom\tstart\tend\tID\tscore\tstrand\tthickStart\tthickEnd\titemRgb\texonCount\texonSizes\texonOffsets\tcircType\tgeneID" > Merge_circexplorer_BC221.annotation.bed14
cat `ls ~/neurogen/rnaseq_PD/run_output/*/circularRNA_known3.txt | grep -f samplelist.braincode2.woCSF.allN221.txt -` | awk '{OFS="\t"; $5=0;$4=$1"_"$2"_"$3; print}' | cut -f1-12,14-15 | sort -u >> Merge_circexplorer_BC221.annotation.bed14
wc -l Merge_circexplorer_BC221.annotation.bed14
# 503,255

## for RM
echo "HC_WGC082364_SN_R_b2_r1
HC_WGC082365_SN_M_b2_r1
HC_MD6326_MB_M_b3
HC_MD6326_MB_R_b3
HC_WGC082366_TC_R_b1
HC_WGC082367_TC_M_b1
HC_WGC082368_TC_R_b2
HC_WGC082369_TC_M_b2
HC_WGC082370_FB_R_b1
HC_WGC082371_FB_M_b1
HC_WGC082372_PBMC_R_b1
HC_WGC082373_PBMC_M_b1" > samplelist.RM.passedN12.txt # excluded the failed ones: HC_WGC082362_SN_R_b1_r1 and HC_WGC082363_SN_M_b1_r1
> Merge_circexplorer_RM12.rawcount.long.txt
cat samplelist.RM.passedN12.txt | while read i; do echo $i; awk -vi=$i '{OFS="\t"; print $1"_"$2"_"$3, $13, i}' ~/neurogen/circRNA_seq_Rebeca_HC/run_output/$i/circularRNA_known3.txt | sed 's/WGC082364/MC3290/;s/WGC082365/MC3290/;s/WGC082366/TCKY1247/;s/WGC082367/TCKY1247/;s/WGC082368/TCKY1217/;s/WGC082369/TCKY1217/;s/WGC082370/ND34770/;s/WGC082371/ND34770/;s/WGC082372/H02018/;s/WGC082373/H02018/' >> Merge_circexplorer_RM12.rawcount.long.txt; done

echo -e "chrom\tstart\tend\tID\tscore\tstrand\tthickStart\tthickEnd\titemRgb\texonCount\texonSizes\texonOffsets\tcircType\tgeneID" > Merge_circexplorer_RM12.annotation.bed14
cat `ls ~/neurogen/circRNA_seq_Rebeca_HC/run_output/*/circularRNA_known3.txt | grep -f samplelist.RM.passedN12.txt -` | awk '{OFS="\t"; $5=0;$4=$1"_"$2"_"$3; print}' | cut -f1-12,14-15 | sort -u >> Merge_circexplorer_RM12.annotation.bed14
wc -l Merge_circexplorer_RM12.annotation.bed14
# 267,398

## BC+RM
cat Merge_circexplorer_BC221.annotation.bed14 Merge_circexplorer_RM12.annotation.bed14 | grep -v chrom | sort -k1,1 -k2,2n -u  >  Merge_circexplorer_BC_RM.annotation.bed14


## track for UCSC 
awk '{OFS="\t";$9=($13=="circRNA")?(($6=="+")?"255,0,0":"0,0,255"):(($6=="+")?"255,100,100":"100,100,255");print}' Merge_circexplorer_BC_RM.annotation.bed14 | cut -f1-12 >  Merge_circexplorer_BC_RM.annotation.bed12
bedToBigBed -type=bed12 Merge_circexplorer_BC_RM.annotation.bed12 ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/hg19.chrom.size Merge_circexplorer_BC_RM.annotation.bb
chmod 644 Merge_circexplorer_BC_RM.annotation.bb;
scp Merge_circexplorer_BC_RM.annotation.bb xd010@panda.dipr.partners.org:~/public_html/tracks
# add the following description to http://panda.partners.org/~xd010/myHub/hg19/trackDb.circRNA.txt
    track BC_BM_circRNAs
    bigDataUrl http://pd:brain@panda.partners.org/~xd010/tracks/Merge_circexplorer_BC_RM.annotation.bb
    shortLabel BC_RM_circularRNA
    longLabel Circular RNAs from all RM (n=12) and BC (n=221) samples (see main.sh in src)
    visibility pack
    itemRgb on
    type bigBed 12
    parent circRNA_braincode

