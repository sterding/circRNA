# R script
# check RNA binding protein expression and binding

## check cell-specificity of RBP genes

library(tidyverse)
library(readxl)

setwd("~/projects/circRNA/data")

# RBP source: https://rbpbase.shiny.embl.de
# wget https://rbpbase.shiny.embl.de/data/RBPbase_Hs_DescriptiveID.xlsx 
rbpbase = read_excel("RBPbase_Hs_DescriptiveID.xlsx",sheet=1)

# Superset of 3,470 known and identified RNA binding proteins (Gebauer et al. 2021)
RBP = rbpbase %>% filter(`humanRBPs-2021\nRBPANNO000000078.1`=='YES') %>% select(1:3) 
# for those with multiple IDs (e.g. "ENSG00000275700|ENSG00000276072"), just take the first one, as from example review, the rest are (Human Alternative sequence Gene) on scaffold in Ensembl 
RBP = mutate(RBP, ensID=str_extract(ID, "\\w+"))

#Q1: any RBP genes are cell-specific expressed?
geneSpecificity_group3mean=readRDS("~/projects/circRNA/data/geneSpecificity.gene_group3mean.rds")
geneSpecificity_group3mean %>% mutate(gene=str_extract(gene, "\\w+")) %>% 
 inner_join(RBP, by=c("gene" = "ensID")) %>% filter(Private_or_not==1) %>% 
 select(celltype, UNIQUE) %>%
 group_by(celltype) %>% mutate(n=length(UNIQUE), UNIQUE=paste0(UNIQUE, collapse = "; ")) %>% distinct()
# celltype UNIQUE                                                              n
# 1 SNDA     RPL10L; RBMXL3; BOLA2; ZCCHC13; AKAP17A; RBMY1A1; RBMY1B; CALR3     8
# 2 NN       SLC25A6; HIST1H1B                                                   2
# 3 PY       RBMY1F; CALML5; DDX53; RNASE13; RBMY1E                              5

#Q2: any predicted RBP binding sites on the cell-specific circRNAs?
#bash to get fasta sequencing of circRNAs
cut -f1-12 Merge_circexplorer_BC197.filtered.enriched.annotation.bed14 | bedtools getfasta -fi ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -bed - -nameOnly -split -s -tab -fo Merge_circexplorer_BC197.filtered.enriched.annotation.fasta.tab
# cell-specific ones
awk 'NR>1{print $4}' Merge_circexplorer_BC109.cellspecific_heatmap.circRNA3.txt | fgrep -f - Merge_circexplorer_BC197.filtered.enriched.annotation.fasta.tab | awk 'NR<5000{split($1,a,"("); print ">"a[1]"\n"$2;}' > Merge_circexplorer_BC109.cellspecific_heatmap.circRNA3.part1.fa
awk 'NR>1{print $4}' Merge_circexplorer_BC109.cellspecific_heatmap.circRNA3.txt | fgrep -f - Merge_circexplorer_BC197.filtered.enriched.annotation.fasta.tab | awk 'NR>=5000{split($1,a,"("); print ">"a[1]"\n"$2;}' > Merge_circexplorer_BC109.cellspecific_heatmap.circRNA3.part2.fa
# go to http://www.csbio.sjtu.edu.cn/bioinf/RBPsuite/ http://rbpmap.technion.ac.il to predict
# download reports for split files
cat cellspecificcircRNAs*.RBPmap_Predictions.txt | awk -v OFS="\t" 'BEGIN{print "ID","RBP","SeqPosition","Motif","Kmer","Zscore","Pvalue";}{OFS="\t"; if($1 ~ /^chr/) id=$1; if($1=="Protein:") {match($2,/(.*)\(/,arr);gene=arr[1];} if($1 ~ /^[0-9]+$/) print id,gene,$0;}' | sed 's/ //g' > cellspecificcircRNAs_RBPmap_Predictions.txt
grep -E "Motif|SLC25A6|HIST1H1B|RBMY1F|CALML5|DDX53|RNASE13|RBMY1E|RPL10L|RBMXL3|BOLA2|ZCCHC13|AKAP17A|RBMY1A1|RBMY1B|CALR3" cellspecificcircRNAs_RBPmap_Predictions.txt

#Q2b: any predicted RBP binding sites on the flanking intron of cell-specific circRNAs?
#bash to get fasta sequencing of circRNAs flanking intron
cut -f1-12 Merge_circexplorer_BC197.filtered.enriched.annotation.bed14 | bedtools flank -b 300 -i - -g ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/hg19.chrom.size | awk '{OFS="\t"; i=(NR%2==1)?"up200":"dw200"; $4=$4"."i; print}' | bedtools getfasta -fi ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -bed - -nameOnly -s -tab -fo Merge_circexplorer_BC197.filtered.enriched.annotation.flankingintron.fasta.tab
# cell-specific ones
awk 'NR>1{print $4}' Merge_circexplorer_BC109.cellspecific_heatmap.circRNA3.txt | fgrep -f - Merge_circexplorer_BC197.filtered.enriched.annotation.flankingintron.fasta.tab | awk '{split($1,a,"("); print ">"a[1]"\n"$2;}' > Merge_circexplorer_BC109.cellspecific_heatmap.flankingintron.circRNA3.fa
# go to http://rbpmap.technion.ac.il to predict
# download reports for split files
cat part*_res/All_Predictions.txt | awk -v OFS="\t" 'BEGIN{print "ID","RBP","SeqPosition","Motif","Kmer","Zscore","Pvalue";}{OFS="\t"; if($1 ~ /^chr/) id=$1; if($1=="Protein:") {match($2,/(.*)\(/,arr);gene=arr[1];} if($1 ~ /^[0-9]+$/) print id,gene,$0;}' | sed 's/ //g' > cellspecificcircRNAs_RBPmap_Predictions.flankingintron.txt
grep -E "Motif|SLC25A6|HIST1H1B|RBMY1F|CALML5|DDX53|RNASE13|RBMY1E|RPL10L|RBMXL3|BOLA2|ZCCHC13|AKAP17A|RBMY1A1|RBMY1B|CALR3" cellspecificcircRNAs_RBPmap_Predictions.flankingintron.txt

#Q3: any predicted RBP binding sites on the flanking introns of KANSL1 circRNAs?
grep KANSL1 Merge_circexplorer_BC197.filtered.enriched.annotation.bed14 | cut -f1-12 | bedtools flank -b 300 -i - -g ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/hg19.chrom.size | awk '{OFS="\t"; i=(NR%2==1)?"up200":"dw200"; $4=$4"."i; print}' | bedtools getfasta -fi ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -bed - -nameOnly -fo circKANSL1.flanking200bp.fasta
# go to http://rbpmap.technion.ac.il to predict
# download reports for split files, eg. wget http://rbpmap.technion.ac.il/1675731465/All_Predictions.txt -O circKANSL1.flanking200bp.RBPmap_Predictions.txt 
cat circKANSL1.flanking200bp.RBPmap_Predictions.txt | awk -v OFS="\t" 'BEGIN{print "ID","RBP","SeqPosition","Motif","Kmer","Zscore","Pvalue";}{OFS="\t"; if($1 ~ /^chr/) id=$1; if($1=="Protein:") {match($2,/(.*)\(/,arr);gene=arr[1];} if($1 ~ /^[0-9]+$/) print id,gene,$0;}' | sed 's/ //g' > circKANSL1.flanking200bp_RBPmap_Predictions.txt

awk 'NR>1{OFS="\t"; split($1,a,"."); split(a[1],b,"_"); if(a[2]=="up200" && $3>100) print b[1],b[2]-300+$3,b[2]-300+$3+length($4),$2,$6; if(a[2]=="dw200" && $3<200) print b[1],b[3]+$3,b[3]+$3+length($4),$2,$6}' circKANSL1.flanking200bp_RBPmap_Predictions.txt | sort -u | wc -l # 10076
awk 'NR>1{OFS="\t"; split($1,a,"."); split(a[1],b,"_"); if(a[2]=="up200" && $3>100) print b[1],b[2]-300+$3,b[2]-300+$3+length($4),$2,$6; if(a[2]=="dw200" && $3<200) print b[1],b[3]+$3,b[3]+$3+length($4),$2,$6}' circKANSL1.flanking200bp_RBPmap_Predictions.txt | sort -u | cut -f4 | sort -u | wc -l # 131

grep KANSL1 Merge_circexplorer_BC197.filtered.enriched.annotation.bed14 | cut -f1-12 | bedtools flank -b 200 -i - -g ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/hg19.chrom.size | awk '{OFS="\t"; i=(NR%2==1)?"up200":"dw200"; $4=$4"."i; print}' | cut -f1-6 > circKANSL1.flanking200bp.hg19.bed
liftOver <(cut -f1-6 RBP_eCLIP_ENCODE/ENCODE.eCLIP.significant.midpoint.bed2) ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/hg38ToHg19.over.chain.gz RBP_eCLIP_ENCODE/ENCODE.eCLIP.significant.midpoint.hg19.bed unmapped -bedPlus=10  # signal and pvalue column are not integer
#liftOver circKANSL1.flanking200bp.hg19.bed ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/hg19ToHg38.over.chain.gz circKANSL1.flanking200bp.hg38.bed unmapped
intersectBed -a RBP_eCLIP_ENCODE/ENCODE.eCLIP.significant.midpoint.hg19.bed -b circKANSL1.flanking200bp.hg19.bed -u | sort -u > circKANSL1.flanking200bp.hg19.ENCODE.eCLIP.bed

# cell-specificity score of these RBPs
# R
GENCODEv19=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v19.annotation.genes.bed", col.names = c("chr","start","end","geneName","geneID","geneType","strand"), stringsAsFactors = F, header = F)
geneSpecificity.gene_group3mean=readRDS(file="~/projects/circRNA/data/geneSpecificity.gene_group3mean.rds")
circKANSL1.flanking200bp.hg19.ENCODE.eCLIP = read.table("circKANSL1.flanking200bp.hg19.ENCODE.eCLIP.bed", header = F) %>% 
  separate(V4, c("RBP_cellline","RBP_target","RBP_rep")) 

inner_join(geneSpecificity.gene_group3mean, GENCODEv19, by = c("gene" = "geneID")) %>% 
  filter(geneName %in% circKANSL1.flanking200bp.hg19.ENCODE.eCLIP$RBP_target) %>%
  ggplot(aes(y=S, x=geneName, size=mean)) + 
  geom_point(shape=21, fill="white", color="black") +
  geom_hline(yintercept = 0.5, color='blue', linetype = "dashed") +
  ylim(0,1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename = "circKANSL1.flanking200bp.hg19.ENCODE.eCLIP.RBP.cellspecificity.dotplot.pdf", width = 6, height = 3)
  


## RBP binding sites from ENCODE eCLIP

```{bash}

cd ~/projects/circRNA/data/RBP_eCLIP_ENCODE

# download narrowpeak files for 209 ENCODE eCLIP experiments:
# https://www.encodeproject.org/search/?type=Experiment&control_type!=*&status=released&perturbed=false&target.investigated_as=RNA+binding+protein&assay_title=eCLIP&files.file_type=bed+narrowPeak&assembly=GRCh38
# IMPORTANT NOTE: Choose "Download default analysis files", instead of "Download default files", as some experiments have wrongly assigned rep1_vs_rep2 into default file. So to get rep1 and rep2 orignial files, you need to get all analysis files, then filter out them in the metatable. 
head -n1 files.txt | xargs -L 1 curl -O -J -L # save metadata.tsv
awk -v FS="\t" 'NR>1 && $35 !~ /,/' metadata.tsv | wc -l # n = 418 (209*2) --> correct
awk -v FS="\t" 'NR>1 && $35 !~ /,/{print $48}' metadata.tsv > files.reps.txt # only the rep1 and rep2
xargs -L 1 curl -O -J -L < files.reps.txt

# filter: We generally recommend that for a high stringency set of high-confidence peaks, cutoffs should be set at fold-enrichment ≥ 8 and p-value ≤ 10^-3
# in the narrowpeak file, Signal value (e.g., log2(fold-enrichment)) for RBP IP compared to size- matched input for the listed peak region
# P -value (typically reported as –log10(p-value)), determined by Fisher’s Exact test (or Chi-Square approximation where appropriate) comparing RBP IP versus size-matched input for the listed peak region
# see https://www.encodeproject.org/documents/9422b49f-7105-4d90-86f8-1201e008eb58/@@download/attachment/eCLIP_experimental_guidelines_doc_v1.0.pdf

# still some files (rep1 or rep2) have missing (.) value in the 4th column
# zcat *.bed.gz | awk '$7>=3 && $8>=3' > ENCODE.eCLIP.significant.bed  # or awk '$5==1000'
# # mid point
# awk '{OFS="\t"; s=int(($3+$2)/2); $2=s; $3=s+1; print;}' ENCODE.eCLIP.significant.bed > ENCODE.eCLIP.significant.midpoint.bed

ENCODE.eCLIP.significant.midpoint.bed2; 
awk -v FS="\t" -v OFS="\t" 'NR>1 && $35 !~ /,/{print $1, $11"_"$23"_rep0"$35}' metadata.tsv | sed 's/-human//g;s/ //g' | while read id target; do echo $id $target;  zcat $id.bed.gz | awk -v target=$target '$5==1000{$4=target; OFS="\t"; s=int(($3+$2)/2); $2=s; $3=s+1; print}' >> ENCODE.eCLIP.significant.midpoint.bed2; done

# liftover from hg38 to hg19
#liftOver ENCODE.eCLIP.significant.midpoint.bed ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/hg38ToHg19.over.chain.gz ENCODE.eCLIP.significant.midpoint.hg19.bed unmapped -bedPlus=10  # signal and pvalue column are not integer

liftOver ../Merge_circexplorer_BC197.filtered.enriched.annotation.bed14 ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/hg19ToHg38.over.chain.gz ../Merge_circexplorer_BC197.filtered.enriched.annotation.hg38.bed14 unmapped -bedPlus=10 


# circRNA upstream 500bp
bedtools flank -i ../Merge_circexplorer_BC197.filtered.enriched.annotation.hg38.bed14 -g ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/hg38.chrom.size -l 500 -r 0 -s | intersectBed -a - -b ENCODE.eCLIP.significant.midpoint.bed2 -s -wo | awk '{OFS="\t"; d=($6=="+")?($18-$3):($2-$19); print $4,$13,$14,$15,$16,d,$20,$23,$24}' > Merge_circexplorer_BC197.filtered.enriched.annotation.up500bp.RBP_ENCODEeCLIPmidpoint.tab
# circRNA downstream 500bp
bedtools flank -i ../Merge_circexplorer_BC197.filtered.enriched.annotation.hg38.bed14 -g ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/hg38.chrom.size -r 500 -l 0 -s | intersectBed -a - -b ENCODE.eCLIP.significant.midpoint.bed2 -s -wo | awk '{OFS="\t"; d=($6=="+")?($19-$2):($3-$18); print $4,$13,$14,$15,$16,d,$20,$23,$24}' > Merge_circexplorer_BC197.filtered.enriched.annotation.dw500bp.RBP_ENCODEeCLIPmidpoint.tab
# circRNA itself (exon only) -- rrS: relative position in 1-100 scale, relative to the 5SS of exon 
awk '{OFS="\t"; $4=$4"."$13; print}' ../Merge_circexplorer_BC197.filtered.enriched.annotation.hg38.bed14 | cut -f1-12 | intersectBed -a - -b ENCODE.eCLIP.significant.midpoint.bed2 -s -split -wo | awk '{OFS="\t"; split($11,l,","); split($12,s,","); L=0; rS=0; for(i=1;i<=$10;i++) {if($14>=($2+s[i]) && $14<($2+s[i]+l[i])) rS=L+($14-$2-s[i]); L=L+l[i];} rrS=($18=="+")?(rS/L):(1-rS/L);split($4,id,".");print id[1],id[2],rS,L,$18,int(100*rrS),$16,$19,$20;}'  > Merge_circexplorer_BC197.filtered.enriched.annotation.exon.RBP_ENCODEeCLIPmidpoint.tab

```

## plot relative position 
library(tidyverse)
setwd("~/projects/circRNA/data/RBP_eCLIP_ENCODE")
RBP_ENCODEeCLIPmidpoint=rbind(cbind(read.table("Merge_circexplorer_BC197.filtered.enriched.annotation.up500bp.RBP_ENCODEeCLIPmidpoint.tab", stringsAsFactors=F),segment="up500bp"),
                              cbind(read.table("Merge_circexplorer_BC197.filtered.enriched.annotation.exon.RBP_ENCODEeCLIPmidpoint.tab", stringsAsFactors=F),segment="exon100"),
                              cbind(read.table("Merge_circexplorer_BC197.filtered.enriched.annotation.dw500bp.RBP_ENCODEeCLIPmidpoint.tab", stringsAsFactors=F),segment="dw500bp"))
names(RBP_ENCODEeCLIPmidpoint) = c("ID","circType","hostID","hostSymbol","hostType","RBP_midpoint","RBP","log2FC","log10pValue","segment")

p=RBP_ENCODEeCLIPmidpoint %>% separate(RBP, c("RBP_cellline","RBP_target","RBP_rep")) %>%
 # filter(RBP_rep=="rep1") %>% # optional
 filter(RBP_cellline!="adrenalgland") %>% # optional
 filter(abs(RBP_midpoint)<200) %>%
 mutate(RBP_midpoint=ifelse(segment=="dw500bp",RBP_midpoint+100,RBP_midpoint)) %>%
 ggplot(aes(RBP_midpoint, fill = RBP_target)) + 
 geom_histogram(position = "stack", binwidth=1) +
 scale_x_continuous(breaks=seq(-200,300,100),labels=c(-200, -100, "start","end",100,  200)) +
 facet_grid(vars(RBP_cellline), vars(RBP_rep), scales = "fixed") +
guides(fill = guide_legend(ncol= 19, title = NULL))+ theme(legend.key.size = unit(5, 'points'),legend.position="bottom")
ggsave(p,filename = "RBP_ENCODEeCLIPmidpoint.rep.pdf", width = 14, height = 6)

p=RBP_ENCODEeCLIPmidpoint %>% separate(RBP, c("RBP_cellline","RBP_target","RBP_rep")) %>%
  # filter(RBP_rep=="rep1") %>% # optional
  filter(RBP_cellline!="adrenalgland") %>% # optional
  filter(abs(RBP_midpoint)<200) %>%
  mutate(RBP_midpoint=ifelse(segment=="dw500bp",RBP_midpoint+100,RBP_midpoint)) %>%
  ggplot(aes(RBP_midpoint, fill = RBP_target)) + 
  geom_histogram(position = "stack", binwidth=1) +
  scale_x_continuous(breaks=seq(-200,300,100),labels=c(-200, -100, "start","end",100,  200)) +
  facet_grid(vars(RBP_cellline), vars(circType), scales = "fixed") +
  guides(fill = guide_legend(ncol= 19, title = NULL))+ theme(legend.key.size = unit(5, 'points'),legend.position="bottom")
ggsave(p,filename = "RBP_ENCODEeCLIPmidpoint.circType.pdf", width = 14, height = 6)

## Q: If more RBP binding site on circRNAs vs. intron
RBP_ENCODEeCLIPmidpoint %>% separate(RBP, c("RBP_cellline","RBP_target","RBP_rep")) %>%
  filter(RBP_cellline!="adrenalgland") %>% # optional
  filter(abs(RBP_midpoint)<200) %>%
  mutate(len=as.numeric(ifelse(segment=="exon100",hostSymbol, 200)),
         segment2=ifelse(segment=="exon100","exon","intron")) %>%
  group_by(ID, segment2, len) %>% summarise(n=n()) %>%
  mutate(m=n/len) %>% wilcox.test(data =., m ~ segment2)
# p-value < 2.2e-16

RBP_ENCODEeCLIPmidpoint %>% separate(RBP, c("RBP_cellline","RBP_target","RBP_rep")) %>%
  filter(RBP_cellline!="adrenalgland") %>% # optional
  filter(abs(RBP_midpoint)<200) %>%
  mutate(len=as.numeric(ifelse(segment=="exon100",hostSymbol, 200)),
         segment2=ifelse(segment=="exon100","exon","intron")) %>%
  group_by(ID, segment2, len) %>% summarise(n=n()) %>%
  mutate(m=n/len) %>% ggplot(aes(x=segment2, y=m)) +
  geom_boxplot() + scale_y_continuous(trans='log10')
ggsave(filename = "RBP_ENCODEeCLIPmidpoint.ExonvsIntron.boxplot.pdf", width = 4, height = 4)

## RBP binding at random non-circularized exons (using Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.matched2 from makeControl.sh)
# bash
bedtools intersect -a ../Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.background.gz -b ../Merge_circexplorer_BC197.filtered.enriched.annotation.bed14 -s -v | shuf -n 13401 --random-source=<(yes 42) > ../Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.background.random
liftOver ../Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.background.random ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/hg19ToHg38.over.chain.gz ../Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.background.random.hg38 unmapped -bedPlus=12 
bedtools flank -i <(cut -f1-12 ../Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.background.random.hg38) -g ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/hg38.chrom.size -l 500 -r 0 -s | intersectBed -a - -b ENCODE.eCLIP.significant.midpoint.bed2 -s -wo | awk '{OFS="\t"; d=($6=="+")?($14-$3):($2-$15); print $4,d,$16,$19,$20}' > Merge_circexplorer_BC197.filtered.enriched.annotation.matched2.up500bp.RBP_ENCODEeCLIPmidpoint.tab
bedtools flank -i <(cut -f1-12 ../Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.background.random.hg38) -g ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/hg38.chrom.size -r 500 -l 0 -s | intersectBed -a - -b ENCODE.eCLIP.significant.midpoint.bed2 -s -wo | awk '{OFS="\t"; d=($6=="+")?($15-$2):($3-$14); print $4,d,$16,$19,$20}' > Merge_circexplorer_BC197.filtered.enriched.annotation.matched2.dw500bp.RBP_ENCODEeCLIPmidpoint.tab
# circRNA itself (exon only) -- rrS: relative position in 1-100 scale, relative to the 5SS of exon 
cut -f1-12 ../Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.background.random.hg38 | intersectBed -a - -b ENCODE.eCLIP.significant.midpoint.bed2 -s -split -wo | awk '{OFS="\t"; split($11,l,","); split($12,s,","); L=0; rS=0; for(i=1;i<=$10;i++) {if($14>=($2+s[i]) && $14<($2+s[i]+l[i])) rS=L+($14-$2-s[i]); L=L+l[i];} rrS=($18=="+")?(rS/L):(1-rS/L);print $4,int(100*rrS),$16,$19,$20;}'  > Merge_circexplorer_BC197.filtered.enriched.annotation.matched2.exon.RBP_ENCODEeCLIPmidpoint.tab

# R
setwd("~/projects/circRNA/data/RBP_eCLIP_ENCODE")
randomRBP_ENCODEeCLIPmidpoint=rbind(cbind(read.table("Merge_circexplorer_BC197.filtered.enriched.annotation.matched2.up500bp.RBP_ENCODEeCLIPmidpoint.tab", stringsAsFactors=F),segment="up500bp"),
                              cbind(read.table("Merge_circexplorer_BC197.filtered.enriched.annotation.matched2.exon.RBP_ENCODEeCLIPmidpoint.tab", stringsAsFactors=F),segment="exon100"),
                              cbind(read.table("Merge_circexplorer_BC197.filtered.enriched.annotation.matched2.dw500bp.RBP_ENCODEeCLIPmidpoint.tab", stringsAsFactors=F),segment="dw500bp"))
names(randomRBP_ENCODEeCLIPmidpoint) = c("ID","RBP_midpoint","RBP","log2FC","log10pValue","segment")
  
df=rbind(RBP_ENCODEeCLIPmidpoint %>% 
          separate(RBP, c("RBP_cellline","RBP_target","RBP_rep")) %>%
          filter(circType=="circRNA") %>% 
          select("ID","RBP_midpoint","RBP_cellline","RBP_target","RBP_rep","segment") %>%
          mutate(type="circRNA-forming exon"),
        randomRBP_ENCODEeCLIPmidpoint %>% 
          separate(RBP, c("RBP_cellline","RBP_target","RBP_rep")) %>%       
          select("ID","RBP_midpoint","RBP_cellline","RBP_target","RBP_rep","segment") %>%
          mutate(type="random exon")
        )

saveRDS(df, file="df.rds")

p=filter(df, abs(RBP_midpoint)<200) %>%
  filter(RBP_cellline!="adrenalgland") %>% # optional
  mutate(RBP_midpoint=ifelse(segment=="dw500bp",RBP_midpoint+100,RBP_midpoint)) %>%
 ggplot(aes(RBP_midpoint, fill = RBP_target)) + 
 geom_histogram(position = "stack", binwidth=2) +
 scale_x_continuous(breaks=seq(-200,300,100),labels=c(-200, -100, "start","end",100,  200)) +
 facet_grid(vars(RBP_cellline), vars(type), scales = "free") +
 guides(fill = guide_legend(ncol= 19, title = NULL))+ theme(legend.key.size = unit(5, 'points'),legend.position="bottom")

ggsave(p,filename = "ENCODEeCLIPmidpoint.AllvsRandom.pdf", width = 14, height = 6)

## RBP binding site on circRNAs vs. control exons
filter(df, segment=="exon100", RBP_cellline!="adrenalgland") %>%
  group_by(type) %>% summarize(n=n(),m=n_distinct(ID),k=n()/n_distinct(ID))
# type                     n     m     k
# circRNA-forming exon 267106  9898  27.0
# random exon          500852  7300  68.6
filter(df, segment=="exon100", RBP_cellline!="adrenalgland") %>%
  group_by(type, ID) %>% summarize(n=n()) %>% wilcox.test(data =., n ~ type)
# p-value < 2.2e-16

filter(df, segment=="exon100", RBP_cellline!="adrenalgland") %>%
  group_by(type, ID) %>% 
  summarize(n=n()) %>% 
  #group_by(type) %>% summarise(mean=mean(n))
  ggplot(aes(x=type, y=n)) +
  geom_boxplot() + scale_y_continuous(trans='log10')
ggsave(filename = "RBP_ENCODEeCLIPmidpoint.AllvsRandom.exon.boxplot.pdf", width = 4, height = 4)
# flanking intron
filter(df, segment!="exon100", RBP_cellline!="adrenalgland") %>%
  group_by(type, ID) %>% 
  summarize(n=n()) %>% 
  #group_by(type) %>% summarise(mean=mean(n))
  ggplot(aes(x=type, y=n)) +
  geom_boxplot() + scale_y_continuous(trans='log10')
ggsave(filename = "RBP_ENCODEeCLIPmidpoint.AllvsRandom.flankingIntron.boxplot.pdf", width = 4, height = 4)
# conclution: circRNA are bound by less RBPs than non-circRNA exons, both on the exon and the flanking introns

## Q1: whether cell-specific circRNAs are less likely harbored with RBP binding?
# 2x2 table for all circRNAs: cell-specific vs. RBP binding 
all_circRNAs = readRDS("../Merge_circexplorer_BC109.filtered.enriched.groupmean_s3.rds")
circRNAs_withRBP = unique(as.character(filter(RBP_ENCODEeCLIPmidpoint, segment=="exon100") %>% pull(ID)))
withRBP=ifelse(all_circRNAs$gene %in% circRNAs_withRBP, 1, 0)
(Fsh <- fisher.test(table(all_circRNAs$Private_or_not, withRBP), alternative = "less"))
#                 withRBP
# Private_or_not  0    1
#             0  329 1613
#             1 2249 7445
# p = 3.168e-10
# OR = 0.6752279
# Conclusion: cell-specific circRNAs are less likely bound with RBP binding than non-private circRNAs

# test for each RBP
all_circRNAs = readRDS("../Merge_circexplorer_BC109.filtered.enriched.groupmean_s3.rds")
Fsh.result = c()
for(i in unique(df$RBP_target)){
  circRNAs_withRBP = unique(as.character(filter(df, type=="circRNA-forming exon", segment=="exon100", RBP_target==i, RBP_cellline!="adrenalgland") %>% pull(ID)))
  if(length(circRNAs_withRBP)<20) next;
  withRBP=ifelse(all_circRNAs$gene %in% circRNAs_withRBP, 1, 0)
  Fsh <- fisher.test(table(all_circRNAs$Private_or_not, withRBP))
  Fsh.result=rbind(Fsh.result, c(i, sum(withRBP),Fsh$p.value, Fsh$estimate))
  message(paste(i, sum(withRBP), Fsh$p.value, Fsh$estimate))
}
colnames(Fsh.result) = c("RBP","n","p","OR")
Fsh.result=as.data.frame(Fsh.result, stringsAsFactors=F)

mutate(Fsh.result, n=as.numeric(n),p=as.numeric(p),OR=as.numeric(OR)) %>% 
  mutate(padj=p.adjust(p, method="fdr")) %>%
  filter(padj<0.05) %>%  
  arrange(desc(OR)) %>% mutate(RBP = factor(RBP, unique(as.character(RBP)))) %>% 
  ggplot(aes(y=-OR, x=RBP, size=n,  color=-log10(padj))) + 
  geom_point(shape=16) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename = "RBP_ENCODEeCLIPmidpoint.private_circRNAs.dotplot.pdf", width = 6, height = 8)

# cell-specificity score of these RBPs
GENCODEv19=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v19.annotation.genes.bed", col.names = c("chr","start","end","geneName","geneID","geneType","strand"), stringsAsFactors = F, header = F)
geneSpecificity.gene_group3mean=readRDS(file="~/projects/circRNA/data/geneSpecificity.gene_group3mean.rds")
geneSpecificity.gene_group3mean = inner_join(geneSpecificity.gene_group3mean, GENCODEv19, by = c("gene" = "geneID"))
mutate(Fsh.result, n=as.numeric(n),p=as.numeric(p),OR=as.numeric(OR)) %>% 
  mutate(padj=p.adjust(p, method="fdr")) %>%
  filter(padj<0.05) %>%  
  inner_join(y=geneSpecificity.gene_group3mean, by = c("RBP" = "geneName")) %>%
  arrange(desc(OR)) %>% mutate(RBP = factor(RBP, unique(as.character(RBP)))) %>% 
  #ggplot(aes(y=S, x=RBP, size=mean,  color=celltype)) + 
  ggplot(aes(y=S, x=RBP, size=mean)) + 
  geom_point(shape=21, fill="white", color="black") +
  geom_hline(yintercept = 0.5, color='blue', linetype = "dashed") +
  ylim(0,1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename = "RBP_ENCODEeCLIPmidpoint.private_circRNAs.RBPcellspecificity.dotplot.pdf", width = 6, height = 8)

mutate(Fsh.result, n=as.numeric(n),p=as.numeric(p),OR=as.numeric(OR)) %>% 
  filter(p<0.05) %>% arrange(-n)
#       RBP    n            p        OR
# 1    BUD13 4850 2.475595e-03 0.8582739
# 2     PPIG 4437 2.124078e-02 0.8885877
# 3   ZNF622 4319 8.785242e-07 0.7788141
# 4    UCHL5 4019 2.263018e-04 0.8259146
# 5    DDX24 4017 4.992595e-02 0.9032809
# 6    DDX43 3657 5.343631e-03 0.8620730
# 7    U2AF2 3499 3.121846e-03 0.8533879
# 8     FXR2 2810 2.768981e-03 0.8431312
# 9    PRPF8 2804 8.153994e-03 0.8593075
# 10     AQR 2731 1.166450e-02 0.8650046
# 11    YBX3 2349 1.707411e-02 0.8647603
# 12    UTP3 2158 6.557401e-03 0.8441850
# 13 IGF2BP1 2115 2.311377e-08 0.7078650
# 14   SRSF1 1810 4.414125e-03 0.8276771
# 15   GPKOW 1742 2.575532e-02 0.8595087
# 16   RBM15 1718 2.643555e-04 0.7807791
# 17  LIN28B 1630 3.453143e-02 0.8620025
# 18    FMR1 1622 6.464231e-04 0.7900743
# 19    PUM1 1290 3.587644e-02 0.8514854
# 20 IGF2BP2 1210 1.473014e-09 0.6332934
# 21 EXOSC10 1166 2.039303e-02 0.8306773
# 22   TRA2A 1063 2.517150e-03 0.7801815
# 23    HLTF  997 1.297335e-04 0.7238864
# 24   DDX3X  764 1.968474e-04 0.7033971
# 25   NOLC1  552 7.037441e-03 0.7413523
# 26   PCBP2  351 1.459863e-04 0.6065723
# 27  SMNDC1  268 9.929945e-03 0.6728494
# 28 IGF2BP3  235 9.965139e-04 0.5901866
# 29  ZRANB2  226 2.847560e-03 0.6161038
# 30   NCBP2  189 1.796012e-02 0.6550471
# 31   RBM22  128 3.090577e-02 0.6230994
# 32   SRSF9  117 2.417795e-02 0.6043386
# 33  DROSHA   87 4.181365e-02 0.5891624
# 34 FAM120A   47 4.784579e-03 0.3862188
# 35    SLTM   38 4.997586e-02 0.4903232
# 36   DDX51   36 8.651945e-05 0.2489044
# 37    SAFB   20 2.844204e-03 0.2440312

# Q2: Do cell-specific circRNAs (n=9694) have different RBP binding pattern than the non-cell-specific ones?

p=RBP_ENCODEeCLIPmidpoint %>% filter(ID %in% all_circRNAs$gene) %>% 
  mutate(Private_or_not = all_circRNAs$Private_or_not[match(ID, all_circRNAs$gene)]) %>%
  separate(RBP, c("RBP_cellline","RBP_target","RBP_rep")) %>%
  filter(RBP_cellline!="adrenalgland") %>% # optional
  filter(circType=="circRNA") %>% # optional
  filter(abs(RBP_midpoint)<200) %>%
  mutate(RBP_midpoint=ifelse(segment=="dw500bp",RBP_midpoint+100,RBP_midpoint)) %>%
  ggplot(aes(RBP_midpoint, fill = RBP_target)) + 
  geom_histogram(position = "stack", binwidth=1) +
  scale_x_continuous(breaks=seq(-200,300,100),labels=c(-200, -100, "start","end",100,  200)) +
  facet_grid(vars(Private_or_not), vars(RBP_cellline), scales = "free") +
  guides(fill = guide_legend(ncol= 19, title = NULL))+ theme(legend.key.size = unit(5, 'points'),legend.position="bottom")
ggsave(p,filename = "RBP_ENCODEeCLIPmidpoint.private_circRNAs.pdf", width = 14, height = 6)

# rows are cellline and columns are circRNA type
p=RBP_ENCODEeCLIPmidpoint %>% filter(ID %in% all_circRNAs$gene) %>% 
  mutate(Private_or_not = all_circRNAs$Private_or_not[match(ID, all_circRNAs$gene)]) %>%
  separate(RBP, c("RBP_cellline","RBP_target","RBP_rep")) %>%
  filter(RBP_cellline!="adrenalgland") %>% # optional
  filter(circType=="circRNA") %>% # optional
  filter(abs(RBP_midpoint)<200) %>%
  mutate(RBP_midpoint=ifelse(segment=="dw500bp",RBP_midpoint+100,RBP_midpoint)) %>%
  ggplot(aes(RBP_midpoint, fill = RBP_target)) + 
  geom_histogram(position = "stack", binwidth=1) +
  scale_x_continuous(breaks=seq(-200,300,100),labels=c(-200, -100, "start","end",100,  200)) +
  facet_grid(vars(RBP_cellline), vars(Private_or_not), scales = "free") +
  guides(fill = guide_legend(ncol= 19, title = NULL))+ theme(legend.key.size = unit(5, 'points'),legend.position="bottom")
ggsave(p,filename = "RBP_ENCODEeCLIPmidpoint.private_circRNAs.flipped.pdf", width = 14, height = 6)


# Q: Would be great to also go on step deeper according to the Reviewer’s suggestion and examine the association between RBP activity, predicted binding sites, and expression in cell type specific circRNA isoforms separately (e.g. the alternatively backspliced class of circRNAs: one host gene -> multiple distinct circRNAs in distinct cell types) vs non-isoform cell type specific circRNAs (e.g. the on/off class of circRNAs: one host gene -> expression of one circRNA in a specific cell type) vs ALL non-cell type specific circRNAs.
# The expectation would be that a few RBP have differential activity in circRNA isoforms; another few in non-isoform circRNAs; compared to non-cell type specific circRNAs.

# replicate the number in Fig2A
# filter(all_circRNAs, Private_or_not==1) %>% group_by(hostgene, celltype) %>% summarise(n = n()) %>% arrange(hostgene, celltype) %>% group_by(hostgene) %>% mutate(n=n_distinct(celltype), celltypes=paste0(celltype, collapse = ";")) %>% dplyr::select(-celltype) %>% distinct() %>% group_by(celltypes, n) %>% summarise(N=n())

all_circRNAs2= filter(all_circRNAs, Private_or_not==1) %>% group_by(hostgene, celltype) %>% summarise(n = n()) %>% arrange(hostgene, celltype) %>% group_by(hostgene) %>% mutate(n=n_distinct(celltype), celltypes=paste0(celltype, collapse = ";")) %>% dplyr::select(-celltype) %>% distinct() %>% 
  right_join(y=all_circRNAs, by='hostgene') %>% 
  mutate(hostgenetype = factor(ifelse(Private_or_not==0,"non-cell type specific circRNAs", ifelse(n==1,"non-isoform cell type specific circRNAs","isoform cell type specific circRNAs")),
                               levels = c("isoform cell type specific circRNAs", "non-isoform cell type specific circRNAs", "non-cell type specific circRNAs"))) 

p=RBP_ENCODEeCLIPmidpoint %>% filter(ID %in% all_circRNAs2$gene) %>% 
  mutate(hostgenetype = all_circRNAs2$hostgenetype[match(ID, all_circRNAs2$gene)]) %>%
  separate(RBP, c("RBP_cellline","RBP_target","RBP_rep")) %>%
  filter(RBP_cellline!="adrenalgland") %>% # optional
  filter(circType=="circRNA") %>% # optional
  filter(abs(RBP_midpoint)<200) %>%
  mutate(RBP_midpoint=ifelse(segment=="dw500bp",RBP_midpoint+100,RBP_midpoint)) %>%
  ggplot(aes(RBP_midpoint, fill = RBP_target)) + 
  geom_histogram(position = "stack", binwidth=1) +
  scale_x_continuous(breaks=seq(-200,300,100),labels=c(-200, -100, "start","end",100,  200)) +
  facet_grid(vars(RBP_cellline), vars(hostgenetype), scales = "free") +
  guides(fill = guide_legend(ncol= 19, title = NULL))+ theme(legend.key.size = unit(5, 'points'),legend.position="bottom")
ggsave(p,filename = "RBP_ENCODEeCLIPmidpoint.private_circRNAs.v2.flipped.pdf", width = 14, height = 6)

# test for each RBP
plot_list = list()
plot1_list = list()
plot0_list = list()
j=1
for(hostgenetypes in c("isoform cell type specific circRNAs", "non-isoform cell type specific circRNAs")){
  Fsh.result=data.frame(matrix(ncol = 5, nrow = 0))
  
  for(i in unique(df$RBP_target)){
    circRNAs_withRBP = unique(as.character(filter(df, type=="circRNA-forming exon", segment=="exon100", RBP_target==i, RBP_cellline!="adrenalgland") %>% pull(ID)))
    if(length(circRNAs_withRBP)<20) next;
    
    # isoform cell-specific circRNAs vs. non-cellspecific circRNAs # N = 7970
    # non-isoform cell-specific circRNAs vs. non-cellspecific circRNAs # N = 5608
    all_circRNAs2_subset = filter(all_circRNAs2, hostgenetype %in% c(hostgenetypes, "non-cell type specific circRNAs")) 
    if_isoform = ifelse(all_circRNAs2_subset$hostgenetype==hostgenetypes,1,0)
    if_withRBP = ifelse(all_circRNAs2_subset$gene %in% circRNAs_withRBP, 1, 0)
    Fsh <- fisher.test(table(if_isoform, if_withRBP))
    Fsh.result=rbind(Fsh.result, c(i, hostgenetypes, sum(if_withRBP),Fsh$p.value, Fsh$estimate))
    message(paste(i, hostgenetypes, sum(if_withRBP), Fsh$p.value, Fsh$estimate))
  }
  colnames(Fsh.result) = c("RBP","hostgenetype","n","p","OR")
  
  Fsh.result = mutate(Fsh.result, n=as.numeric(n),p=as.numeric(p),OR=as.numeric(OR)) %>% 
    mutate(padj=p.adjust(p, method="fdr")) %>%
    filter(padj<0.05) %>%  
    arrange(desc(OR)) %>% mutate(RBP = factor(RBP, unique(as.character(RBP)))) 
  
  write.table(Fsh.result, paste0("RBP_ENCODEeCLIPmidpoint.",gsub(" ","-",hostgenetypes),".xls"), quote = F, sep = "\t", col.names = F)

  p = ggplot(Fsh.result, aes(y=OR, x=RBP, size=n,  color=-log10(padj))) + 
    geom_point(shape=16) +
    scale_y_reverse() +
    geom_hline(yintercept = 1, color='blue', linetype = "dashed") +
    ggtitle(paste(hostgenetypes, "vs. non-cell type specific circRNAs")) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  plot_list[[j]] = p
  #ggsave(filename = paste0("RBP_ENCODEeCLIPmidpoint.",gsub(" ","-",hostgenetypes),".dotplot.pdf"), width = 6, height = 8)
  
  # y axis from 1 (top) to ~0.1 (bottom)
  p = ggplot(Fsh.result, aes(y=OR, x=RBP, size=n,  color=-log10(padj))) + 
    geom_point(shape=16) +
    scale_y_log10(limits=c(0.1,10)) +
    geom_hline(yintercept = 1, color='blue', linetype = "dashed") +
    ggtitle(paste(hostgenetypes, "vs. non-cell type specific circRNAs")) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  plot1_list[[j]] = p
  
  p0=inner_join(x=Fsh.result, y=geneSpecificity.gene_group3mean, by = c("RBP" = "geneName")) %>%
    arrange(desc(OR)) %>% mutate(RBP = factor(RBP, unique(as.character(RBP)))) %>%
    ggplot(aes(y=S, x=RBP, size=mean)) + 
    geom_point(shape=21, fill="white", color="black") +
    geom_hline(yintercept = 0.5, color='blue', linetype = "dashed") +
    ylim(0,1) +
    labs(size='Mean expression') + ylab("Cell specificity score") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  plot0_list[[j]] = p0
  
  j=j+1;
}
library('cowplot')
plot_grid(plotlist=plot1_list, labels = "AUTO", rel_widths = c(1, 3))
ggsave(filename = "RBP_ENCODEeCLIPmidpoint.private_circRNAs.v4.dotplot.pdf", width = 15, height = 5)

plot_grid(plotlist=plot_list, labels = "AUTO", rel_widths = c(1, 3))
ggsave(filename = "RBP_ENCODEeCLIPmidpoint.private_circRNAs.v2.dotplot.pdf", width = 15, height = 5)

plot_grid(plotlist=plot0_list, labels = "AUTO", rel_widths = c(1, 3))
ggsave(filename = "RBP_ENCODEeCLIPmidpoint.private_circRNAs.RBPcellspecificity.v2.dotplot.pdf", width = 15, height = 5)

## comparison: "non-isoform cell type specific circRNAs" vs. "isoform cell type specific circRNAs"
# test for each RBP
Fsh.result=data.frame(matrix(ncol = 5, nrow = 0))
for(i in unique(df$RBP_target)){
  circRNAs_withRBP = unique(as.character(filter(df, type=="circRNA-forming exon", segment=="exon100", RBP_target==i, RBP_cellline!="adrenalgland") %>% pull(ID)))
  if(length(circRNAs_withRBP)<20) next;
  
  all_circRNAs2_subset = filter(all_circRNAs2, hostgenetype %in% c("isoform cell type specific circRNAs", "non-isoform cell type specific circRNAs")) 
  if_isoform = ifelse(all_circRNAs2_subset$hostgenetype=="non-isoform cell type specific circRNAs",1,0)
  if_withRBP = ifelse(all_circRNAs2_subset$gene %in% circRNAs_withRBP, 1, 0)
  Fsh <- fisher.test(table(if_isoform, if_withRBP))
  Fsh.result=rbind(Fsh.result, c(i, "non-isoform vs. isoform", sum(if_withRBP),Fsh$p.value, Fsh$estimate))
  message(paste(i, sum(if_withRBP), Fsh$p.value, Fsh$estimate))
}
colnames(Fsh.result) = c("RBP","hostgenetype","n","p","OR")

Fsh.result = mutate(Fsh.result, n=as.numeric(n),p=as.numeric(p),OR=as.numeric(OR)) %>% 
  mutate(padj=p.adjust(p, method="fdr")) %>%
  filter(padj<0.05) %>%  
  arrange(desc(OR)) %>% mutate(RBP = factor(RBP, unique(as.character(RBP)))) 

write.table(Fsh.result, "RBP_ENCODEeCLIPmidpoint.nonisoform-vs-isoform.xls", quote = F, sep = "\t", col.names = F)

# y axis from 1 (top) to ~0.1 (bottom)
p = ggplot(Fsh.result, aes(y=OR, x=RBP, size=n,  color=-log10(padj))) + 
  geom_point(shape=16) +
  scale_y_log10(limits=c(0.1,10)) +
  geom_hline(yintercept = 1, color='blue', linetype = "dashed") +
  ggtitle("non-isoform vs. isoform cell type specific circRNAs") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(plot=p, filename = "RBP_ENCODEeCLIPmidpoint.nonisoform-vs-isoform.dotplot.pdf", width = 15, height = 5)

p0=inner_join(x=Fsh.result, y=geneSpecificity.gene_group3mean, by = c("RBP" = "geneName")) %>%
  arrange(desc(OR)) %>% mutate(RBP = factor(RBP, unique(as.character(RBP)))) %>%
  ggplot(aes(y=S, x=RBP, size=mean)) + 
  geom_point(shape=21, fill="white", color="black") +
  geom_hline(yintercept = 0.5, color='blue', linetype = "dashed") +
  ylim(0,1) +
  labs(size='Mean expression') + ylab("Cell specificity score") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(plot=p0, filename = "RBP_ENCODEeCLIPmidpoint.nonisoform-vs-isoform.RBPcellspecificity.dotplot.pdf", width = 12, height = 5)


## RBP binding site on three types of circRNAs vs. control exons

left_join(x=filter(df, segment=="exon100", RBP_cellline!="adrenalgland"),
          y=mutate(all_circRNAs2, hostgenetype=as.character(hostgenetype)), 
          by = c("ID" = "gene")) %>% 
  mutate(hostgenetype=ifelse(is.na(hostgenetype), "control",hostgenetype)) %>% 
  #filter(type=="random exon") 
  group_by(hostgenetype) %>% summarize(n=n(),m=n_distinct(ID),k=n()/n_distinct(ID)) 
# hostgenetype                                 n     m     k
# 1 control                                 531117  8507  62.4
# 2 isoform cell type specific circRNAs     134828  4681  28.8
# 3 non-cell type specific circRNAs          44022  1553  28.3
# 4 non-isoform cell type specific circRNAs  57991  2457  23.6

my_comparisons <- list( c("non-cell type specific circRNAs", "isoform cell type specific circRNAs"),
                        c("non-cell type specific circRNAs", "non-isoform cell type specific circRNAs"),
                        c("non-isoform cell type specific circRNAs", "isoform cell type specific circRNAs"),
                        c("control", "non-cell type specific circRNAs"), 
                        c("control", "isoform cell type specific circRNAs"), 
                        c("control", "non-isoform cell type specific circRNAs")
                        )

library(ggpubr)
left_join(x=filter(df, segment=="exon100", RBP_cellline!="adrenalgland"),
          y=mutate(all_circRNAs2, hostgenetype=as.character(hostgenetype)), 
          by = c("ID" = "gene")) %>% 
  mutate(hostgenetype=ifelse(is.na(hostgenetype), "control",hostgenetype)) %>%
  mutate(hostgenetype=factor(hostgenetype, levels = c("control","non-cell type specific circRNAs", "isoform cell type specific circRNAs", "non-isoform cell type specific circRNAs"))) %>%
  group_by(hostgenetype, ID) %>% summarize(n=n()) %>%
  #ggboxplot(x = "hostgenetype", y = "n", add = "jitter", add.params = list(color = "grey"))+ 
  ggboxplot(x = "hostgenetype", y = "n")+ 
  scale_y_continuous(trans='log10') +
  #stat_compare_means(label = "p.signif", comparisons = my_comparisons) +
  stat_compare_means(label = "p.format", comparisons = my_comparisons) + # default: wilcox test
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

ggsave(filename = "RBP_ENCODEeCLIPmidpoint.AllvsRandom.boxplot.v2.pdf", width = 4, height = 8)

# flanking intron
left_join(x=filter(df, segment!="exon100", RBP_cellline!="adrenalgland"),
          y=mutate(all_circRNAs2, hostgenetype=as.character(hostgenetype)), 
          by = c("ID" = "gene")) %>% 
  mutate(hostgenetype=ifelse(is.na(hostgenetype), "control",hostgenetype)) %>%
  mutate(hostgenetype=factor(hostgenetype, levels = c("control","non-cell type specific circRNAs", "isoform cell type specific circRNAs", "non-isoform cell type specific circRNAs"))) %>%
  group_by(hostgenetype, ID) %>% summarize(n=n()) %>%
  #ggboxplot(x = "hostgenetype", y = "n", add = "jitter", add.params = list(color = "grey"))+ 
  ggboxplot(x = "hostgenetype", y = "n")+ 
  scale_y_continuous(trans='log10') +
  #stat_compare_means(label = "p.signif", comparisons = my_comparisons) +
  stat_compare_means(label = "p.format", comparisons = my_comparisons) + # default: wilcox test
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

ggsave(filename = "RBP_ENCODEeCLIPmidpoint.AllvsRandom.flankingIntron.boxplot.v2.pdf", width = 4, height = 8)
