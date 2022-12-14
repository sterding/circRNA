# bash script to generate input file for test the SNP co-localization enrichment analysis

type='SNAP'

cd ~/projects/circRNA/data/
#circRNA_annotation=Merge_circexplorer_BC106.filtered.enriched.annotation.bed14  # BRAINCODE filted only
circRNA_annotation=Merge_circexplorer_BC109.filtered.enriched.annotation.bed14  # BRAINCODEv2 filted only


## pre-steps to get LD data for GWAS SNPs: 
## Ref ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/GWASCatalog/README.txt

[ "$type" == "SNAP" ] && snps_in_LD=$GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.SNAP.LD_w250.r2_0.8.bed  # in the final figure we used SNAP
[ "$type" == "PLINK" ] && snps_in_LD=$GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.PLINK.LD_w250.r2_0.8.bed

## extract all autosomal.associations
#[ -e $snps_in_LD.autosomal.associations.bed ] || awk 'BEGIN{FS="\t"; OFS="\t";}{split($8,a,"|");  n=split(a[2],b,";"); print $1,$2,$3,$7; for(i=1;i<n;i++) print $1,b[i]-1,b[i],$7;}' $snps_in_LD | grep -v chrX | grep -v chrY | sort -u > $snps_in_LD.autosomal.associations.bed
[ -e $snps_in_LD.autosomal.associations.bed ] || awk 'BEGIN{FS="\t"; OFS="\t";}$5>=-log(5e-8)/log(10){split($8,a,"|");  n=split(a[2],b,";"); print $1,$2,$3,$7; for(i=1;i<n;i++) print $1,b[i]-1,b[i],$7;}' $snps_in_LD | grep -v chrX | grep -v chrY | sort -u > $snps_in_LD.autosomal.associations.bed

# number of gwas SNPs
wc -l $snps_in_LD
# number of diseases/traits
cut -f7 $snps_in_LD |sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort -u | wc -l
# number of associations
wc -l $snps_in_LD.autosomal.associations.bed

## circRNA-matched exons
# see makeControl.sh --> control2

## circRNA-matched exons from the same gene
# see makeControl.sh --> control3

# circRNA-flanking intron
bedtools closest -io -s -t all -a <(sortBed -i $circRNA_annotation) -b <(sortBed -i $GENOME/Annotation/Genes/introns.meta.bed) | awk '($6=="-" && $16==$3) || ($6=="+" && $17==$2){OFS="\t"; print $15,$16,$17,$18,$19,$20,$4}' > $circRNA_annotation.flankingIntron.5prime.bed
bedtools closest -io -s -t all -a <(sortBed -i $circRNA_annotation) -b <(sortBed -i $GENOME/Annotation/Genes/introns.meta.bed) | awk '($6=="+" && $16==$3) || ($6=="-" && $17==$2){OFS="\t"; print $15,$16,$17,$18,$19,$20,$4}' > $circRNA_annotation.flankingIntron.3prime.bed
# circRNA-internal intron
bedtools intersect -s -F 1.0 -a <(sortBed -i $circRNA_annotation) -b <(sortBed -i $GENOME/Annotation/Genes/introns.meta.bed) -wo | awk '{OFS="\t"; print $15,$16,$17,$18,$19,$20,$4}' > $circRNA_annotation.internalIntron.bed

echo "## overlapped SNPs with each dataset"
### ##################
echo "# all"
cat $snps_in_LD.autosomal.associations.bed | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.all
echo "# circRNA"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b $circRNA_annotation -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.circRNA
echo "# circRNA-SNDA"
awk '$2~/SNDA/' ../results/Merge_circexplorer_BC.annotation_per_cell.xls | cut -f4- | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.circRNA-SNDA
echo "# circRNA-PY"
awk '$2~/PY/' ../results/Merge_circexplorer_BC.annotation_per_cell.xls | cut -f4- | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.circRNA-PY
echo "# circRNA-NN"
awk '$2~/NN/' ../results/Merge_circexplorer_BC.annotation_per_cell.xls | cut -f4- | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.circRNA-NN
echo "# circRNA-Neuron"
awk 'NR>1 && $2!~/NN/' ../results/Merge_circexplorer_BC.annotation_per_cell.xls | cut -f4- | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.circRNA-Neuron
echo "# circRNA-SNDA-private"
awk '$2=="SNDA"' ../results/Merge_circexplorer_BC.annotation_per_cell.xls | cut -f4- | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.circRNA-SNDA-private
echo "# circRNA-PY-private"
awk '$2=="PY"' ../results/Merge_circexplorer_BC.annotation_per_cell.xls | cut -f4- | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.circRNA-PY-private
echo "# circRNA-NN-private"
awk '$2=="NN"' ../results/Merge_circexplorer_BC.annotation_per_cell.xls | cut -f4- | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.circRNA-NN-private
echo "# circRNA-matched exons"
cut -f4 $circRNA_annotation | fgrep -w -f - Merge_circexplorer_BC.annotation.bed14.matched2 | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.circRNA-matched-exons
echo "# circRNA-matched exons from the same gene"
cut -f4 $circRNA_annotation | fgrep -w -f - Merge_circexplorer_BC.annotation.bed14.matched3 | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.circRNA-matched-exons-same-host
echo "# circRNA-flanking intron"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b $circRNA_annotation.flankingIntron.3prime.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.circRNA-flanking-intron3
intersectBed -a $snps_in_LD.autosomal.associations.bed -b $circRNA_annotation.flankingIntron.5prime.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.circRNA-flanking-intron5
echo "# circRNA-internal intron"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b $circRNA_annotation.internalIntron.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.circRNA-internal-introns
echo "# circRNA-matched internal intron"
bedtools intersect -s -F 1.0 -a <(sortBed -i Merge_circexplorer_BC.annotation.bed14.matched2) -b <(sortBed -i $GENOME/Annotation/Genes/introns.meta.bed) -wo | awk '{OFS="\t"; print $16,$17,$18,$19,$20,$21,$4}'| intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.circRNA-matched-introns

echo "## total overlapped SNPs count with each dataset"
### ##################
echo "all" `wc -l $GENOME/Annotation/Variation/snp137.bed.groupped.SNP | cut -f1 -d' '` > SNP.$type.counts.summary
echo "circRNA" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b $circRNA_annotation -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "circRNA-SNDA" `awk '$2~/SNDA/' ../results/Merge_circexplorer_BC.annotation_per_cell.xls | cut -f4- | intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b - -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "circRNA-PY" `awk '$2~/PY/' ../results/Merge_circexplorer_BC.annotation_per_cell.xls | cut -f4- | intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b - -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "circRNA-NN" `awk '$2~/NN/' ../results/Merge_circexplorer_BC.annotation_per_cell.xls | cut -f4- | intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b - -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "circRNA-Neuron" `awk 'NR>1 && $2!~/NN/' ../results/Merge_circexplorer_BC.annotation_per_cell.xls | cut -f4- | intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b - -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "circRNA-SNDA-private" `awk '$2=="SNDA"' ../results/Merge_circexplorer_BC.annotation_per_cell.xls | cut -f4- | intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b - -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "circRNA-PY-private" `awk '$2=="PY"' ../results/Merge_circexplorer_BC.annotation_per_cell.xls | cut -f4- | intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b - -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "circRNA-NN-private" `awk '$2=="NN"' ../results/Merge_circexplorer_BC.annotation_per_cell.xls | cut -f4- | intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b - -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "circRNA-matched-exons" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b Merge_circexplorer_BC.annotation.bed14.matched2 -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "circRNA-matched-exons-same-host" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b Merge_circexplorer_BC.annotation.bed14.matched3 -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "circRNA-flanking-intron5" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b $circRNA_annotation.flankingIntron.5prime.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "circRNA-flanking-intron3" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b $circRNA_annotation.flankingIntron.3prime.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "circRNA-internal-introns" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b $circRNA_annotation.internalIntron.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "circRNA-matched-introns" `bedtools intersect -s -F 1.0 -a <(sortBed -i Merge_circexplorer_BC.annotation.bed14.matched2) -b <(sortBed -i $GENOME/Annotation/Genes/introns.meta.bed) -wo | awk '{OFS="\t"; print $16,$17,$18,$19,$20,$21,$4}' | intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b - -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary

#echo "## Fisher test and make plot"  # move out of the script now
### ##################
#Rscript $pipeline_path/src/circRNA.SNP.enrichment.R