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
cat Merge_circexplorer_BC221.annotation.bed14 Merge_circexplorer_RM12.annotation.bed14 | grep -v chrom | sort -k1,1 -k2,2n -k3,3n -u  >  Merge_circexplorer_BC_RM.annotation.bed14

#########################################
## tracks for UCSC Genome Browser
#########################################
## track for UCSC 
mkdir ~/projects/circRNA/data/ucsc_tracks;
cd ~/projects/circRNA/data/ucsc_tracks;
awk '{OFS="\t";$9=($13=="circRNA")?(($6=="+")?"255,0,0":"0,0,255"):(($6=="+")?"255,100,100":"100,100,255");print}' ../Merge_circexplorer_BC_RM.annotation.bed14 | cut -f1-12 | sort -k1,1 -k2,2n >  Merge_circexplorer_BC_RM.annotation.bed12
bedToBigBed -type=bed12 Merge_circexplorer_BC_RM.annotation.bed12 ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/hg19.chrom.size Merge_circexplorer_BC_RM.annotation.bb
rsync -a --chmod=u+rw,go+r Merge_circexplorer_BC_RM.annotation.bb xd010@panda.dipr.partners.org:~/public_html/tracks

# separate by cell type
for i in SNDA PY NN; do echo $i; bedToBigBed -type=bed12 Merge_circexplorer_BC197.filtered.enriched.annotation.$i.bed12 ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/hg19.chrom.size Merge_circexplorer_BC197.filtered.enriched.annotation.$i.bb; done
rsync -a --chmod=u+rw,go+r Merge_circexplorer_BC197.filtered.enriched.annotation.*.bb xd010@panda.dipr.partners.org:~/public_html/tracks


## bigBarChart for 197 samples
# require to run main.R to get Merge_circexplorer_BC197.filtered.enriched.normRPM.tab and Merge_circexplorer_BC197.filtered.enriched.annotation.bed14
# - matrixFile: 
ln -fs ../Merge_circexplorer_BC197.filtered.enriched.normRPM.tab Merge_circexplorer_BC197.filtered.enriched.normRPM.tab
# - sampleFile: 
head -n1 Merge_circexplorer_BC197.filtered.enriched.normRPM.tab | rowsToCols stdin stdout | awk '{OFS="\t"; split($1,a,"_"); print $1,a[1]"_"a[3];}' > Merge_circexplorer_BC197.filtered.enriched.sampleFile
# - bedFile (bed6+1): 
cut -f1-6,15 ../Merge_circexplorer_BC197.filtered.enriched.annotation.bed14 > Merge_circexplorer_BC197.filtered.enriched.annotation.bed6+1
# - GROUPORDERFILE
echo HC_SNDA ILB_SNDA PD_SNDA HC_TCPY AD_TCPY HC_MCPY HC_FB HC_PBMC | rowsToCols stdin GROUPORDERFILE
expMatrixToBarchartBed  --useMean --groupOrderFile GROUPORDERFILE Merge_circexplorer_BC197.filtered.enriched.sampleFile Merge_circexplorer_BC197.filtered.enriched.normRPM.tab Merge_circexplorer_BC197.filtered.enriched.annotation.bed6+1 Merge_circexplorer_BC197.filtered.enriched.BarChart.bed
wget https://genome.ucsc.edu/goldenPath/help/examples/barChart/barChartBed.as
sort -k1,1 -k2,2n Merge_circexplorer_BC197.filtered.enriched.BarChart.bed > Merge_circexplorer_BC197.filtered.enriched.BarChart.bed.sorted;
bedToBigBed -as=barChartBed.as -type=bed6+5 Merge_circexplorer_BC197.filtered.enriched.BarChart.bed.sorted ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/hg19.chrom.size Merge_circexplorer_BC197.filtered.enriched.BarChart.bb

## bigbed for 197 samples, separately by celltype
Merge_circexplorer_BC197.filtered.enriched.normRPM.tab
# add the following description to http://panda.partners.org/~xd010/myHub/hg19/trackDb.circRNA.txt

echo "
track circRNA
shortLabel circRNA
longLabel BRAINcode circRNA datasets
dataVersion Version 1 (June 2020)
type bed 3
visibility full
boxedCfg on
priority 24
superTrack on show

    track published_circRNA
    compositeTrack on
    parent circRNA
    shortLabel published_circRNA
    type bigBed 12
    allButtonPair on
    
        track jeck_circRNAs
        bigDataUrl http://bimsbstatic.mdc-berlin.de/hubs/rajewsky/circBase/hg19/Jeck2013_sorted_coloured.bb
        shortLabel jeck2013
        longLabel Jeck et al. 2013 circular RNA candidates
        itemRgb off
        type bigBed 12
        parent published_circRNA
        
        track salzman_circRNAs
        bigDataUrl http://bimsbstatic.mdc-berlin.de/hubs/rajewsky/circBase/hg19/Salzman2013_sorted.bb
        shortLabel salzman2013
        longLabel Salzman et al. 2013 circular RNA candidates
        itemRgb off
        type bigBed 12
        parent published_circRNA
        
        track memczak_circRNAs
        bigDataUrl http://bimsbstatic.mdc-berlin.de/hubs/rajewsky/circBase/hg19/Memczak2013_circRNAs.bb
        shortLabel memczak2013
        longLabel Memczak et al. 2013 circular RNA candidates
        itemRgb off
        type bigBed 12
        parent published_circRNA
        
        track zhang_circRNAs
        bigDataUrl http://bimsbstatic.mdc-berlin.de/hubs/rajewsky/circBase/hg19/Zhang2013_circRNAs.bb
        shortLabel zhang2013
        longLabel Zhang et al. 2013 circular RNA candidates
        itemRgb off
        type bigBed 12
        parent published_circRNA
        
        track rybak_circRNAs
        bigDataUrl http://bimsbstatic.mdc-berlin.de/hubs/rajewsky/circBase/hg19/Rybak2015_circRNAs.bb
        shortLabel rybak2015
        longLabel Rybak et al. 2015 circular RNA candidates
        itemRgb off
        type bigBed 12
        parent published_circRNA

    track circRNA_BC2_filtered
    compositeTrack on
    parent circRNA
    shortLabel circRNA_BC2_filtered
    type bigBed 12
    allButtonPair on
    
        track circRNA_BC2_filtered_SNDA
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/tracks/Merge_circexplorer_BC197.filtered.enriched.annotation.SNDA.bb
        shortLabel circRNA_BC2_filtered_SNDA
        longLabel Filtered and enriched circRNAs in BRAINcode2 (dopamine neuron, N=104)
        itemRgb on
        type bigBed 12
        parent circRNA_BC2_filtered

        track circRNA_BC2_filtered_PY
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/tracks/Merge_circexplorer_BC197.filtered.enriched.annotation.PY.bb
        shortLabel circRNA_BC2_filtered_PY
        longLabel Filtered and enriched circRNAs in BRAINcode2 (pyramidal neuron, N=86)
        itemRgb on
        type bigBed 12
        parent circRNA_BC2_filtered
        
        track circRNA_BC2_filtered_NN
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/tracks/Merge_circexplorer_BC197.filtered.enriched.annotation.NN.bb
        shortLabel circRNA_BC2_filtered_NN
        longLabel Filtered and enriched circRNAs in BRAINcode2 (non-neuronal cells, N=7)
        itemRgb on
        type bigBed 12
        parent circRNA_BC2_filtered

        track circRNA_BC2_filtered_barchart
        type bigBarChart
        visibility full
        shortLabel circRNA_BC2_filtered_barchart
        longLabel Filtered and enriched circRNAs in BRAINcode2 (N=197)
        barChartBars HC_SNDA ILB_SNDA PD_SNDA HC_TCPY AD_TCPY HC_MCPY HC_FB HC_PBMC
        barChartColors #F22A7B #E44892 #A24B9C #3182BD #ED6120 #2659B2 #BC9371 #D3CBCB
        barChartLabel Celltype
        barChartMetric mean
        barChartUnit RPM
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/tracks/Merge_circexplorer_BC197.filtered.enriched.BarChart.bb
        barChartMatrixUrl http://pd:brain@panda.partners.org/~xd010/tracks/Merge_circexplorer_BC197.filtered.enriched.normRPM.tab
        barChartSampleUrl http://pd:brain@panda.partners.org/~xd010/tracks/Merge_circexplorer_BC197.filtered.enriched.sampleFile
        parent circRNA_BC2_filtered
        
    track circularRNA_BC2+RM_all
    bigDataUrl http://pd:brain@panda.partners.org/~xd010/tracks/Merge_circexplorer_BC_RM.annotation.bb
    shortLabel circularRNA_BC2+RM_all
    longLabel Circular RNAs from all RM (n=12) and BC2 (n=221) samples (see main.sh in src)
    visibility pack
    itemRgb on
    type bigBed 12
    parent circRNA

" > trackDb.circRNA.txt
rsync -a --chmod=u+rw,go+r --copy-links trackDb.circRNA.txt xd010@panda.dipr.partners.org:~/public_html/myHub/hg19/
rsync -a --chmod=u+rw,go+r --copy-links *.bb xd010@panda.dipr.partners.org:~/public_html/tracks/
rsync -a --chmod=u+rw,go+r --copy-links *.tab xd010@panda.dipr.partners.org:~/public_html/tracks/
rsync -a --chmod=u+rw,go+r --copy-links *.sampleFile xd010@panda.dipr.partners.org:~/public_html/tracks/



## for Bennett
> Merge_circexplorer_Bennett.rawcount.long.txt
cat samplelist.Bennett.RNAseq.allN42.txt | while read i; do echo $i; awk -vi=$i '{OFS="\t"; print $1"_"$2"_"$3, $13, i}' ~/neurogen/rnaseq_Bennett/run_output/$i/circularRNA_known3.txt >> Merge_circexplorer_Bennett.rawcount.long.txt; done

echo -e "chrom\tstart\tend\tID\tscore\tstrand\tthickStart\tthickEnd\titemRgb\texonCount\texonSizes\texonOffsets\tcircType\tgeneID" > Merge_circexplorer_Bennett.annotation.bed14
cat `ls ~/neurogen/rnaseq_Bennett/run_output/*/circularRNA_known3.txt | grep -f samplelist.Bennett.RNAseq.allN42.txt -` | awk '{OFS="\t"; $5=0;$4=$1"_"$2"_"$3; print}' | cut -f1-12,14-15 | sort -u >> Merge_circexplorer_Bennett.annotation.bed14
wc -l Merge_circexplorer_Bennett.annotation.bed14
# 88699
