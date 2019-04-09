# main script for QTL

### prepare genotype
## BC 
# # Tao's imputation + QC (remove imputation R2 <0.3, MAF < 0.05, multiple allelic SNPs) --> subset 84 sample --> --maf 0.05 --remove-filtered-all
# # on hms:
  # cd hms:/home/tw83/twang/AMP/BrainCode/Imputation/QC/
## combine chrX.male and chrX.female, see README
# for i in chr{[0-9]*,X}.dos.postQC.vcf.gz; do echo $i; bsub -q short -W 12:00 -n 1 -J ${i/.*/} bash subset.sh $i; done
# # on hpc:
# cd ~/projects/circRNA/data/QTL_BC/genotype
# rsync -azv xd20@transfer.med.harvard.edu:/home/tw83/twang/AMP/BrainCode/Imputation/QC/84samples/chr* .

###
### prepare phenotype
###
cd ~/projects/circRNA/data/QTL_BC/phenotype
# circRNA expression, gene expression, PCI
Rscript ~/projects/circRNA/src/qtl/QTL_prepare_phenotype.R
# PSI
cd ~/neurogen/rnaseq_PD/results/sQTL
mkdir bamfiles142 
cd bamfiles142; for i in /PHShome/xd010/neurogen/rnaseq_PD/run_output/*/accepted_hits.bam; do jj=${i/\/accepted_/.accepted_}; ln -fs $i `basename $jj`; done; cd ..
ls bamfiles142/* | grep -f ~/neurogen/rnaseq_PD/results/merged/samplelist.HCILB_SNDA - > bamfiles86.txt
ls bamfiles142/* > bamfiles142.txt
bsub -q long -n 1 -J leafcutter python $HOME/bin/leafcutter/example_data/run_sQTL.py -b bamfiles142.txt -d $HOME/bin/leafcutter -a $HOME/bin/leafcutter/clustering/gencode.v19.annotation.gtf.gz -t sqtl_bamfiles142
bsub -q long -n 1 -J leafcutter python $HOME/bin/leafcutter/example_data/run_sQTL.py -b bamfiles86.txt -d $HOME/bin/leafcutter -a $HOME/bin/leafcutter/clustering/gencode.v19.annotation.gtf.gz -t sqtl_bamfiles86
cd ~/projects/circRNA/data/QTL_BC/phenotype/
for i in `seq 1 22`; do echo $i;  cat ~/neurogen/rnaseq_PD/results/sQTL/sqtl_bamfiles86/leafcutter_perind.counts.gz.qqnorm_chr$i | sed 's/.accepted_hits.bam//g;s/HC_//g;s/ILB_//g;s/_SNDA_[1-9]_rep[12]//g' | awk '{OFS="\t"; s=0; for(i=5;i<=NF;i++) if($i>0) s++; if(NR==1 || s>=1) {if(NR>1 && $1 !~/^chr/) {$1="chr"$1; $4="chr"$4;} print}}' > circRNA.PSI.chr$i.bed & done
# gzip and index *.bed 
for i in *bed; do echo $i; bgzip -f $i && tabix -p bed $i.gz; done

###
### QTL analysis
###
cd ~/projects/circRNA/data/QTL_BC

# eQTL for all circular RNA
for i in `seq 1 22` X; do
  bsub -q short -n 1 -J eQTL.nominal$i fastQTL --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/circularRNA.expression.bed.gz -L _eQTL.nominal$i --include-samples sample.incl --region chr$i --out eQTL.nominal.chr$i.txt --normal --seed 123 --cov covs.fastqtl.txt
  bsub -q short -n 1 -J eQTL.permute$i fastQTL --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/circularRNA.expression.bed.gz -L _eQTL.permute$i --include-samples sample.incl --region chr$i --permute 100 100000 --out eQTL.permutations.chr$i.txt --normal --seed 123 --cov covs.fastqtl.txt
done

# eQTL for all circular RNA (with ciRNAs/intron lariats from the same donor site merged into one group)
for i in `seq 1 22` X; do
  [ -e eQTLcircularRNAmerged.nominal.chr$i.txt ] || bsub -q short -n 1 -J eQTL.nominal$i fastQTL --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/circularRNAmerged.expression.bed.gz -L _eQTLcircularRNAmerged.nominal$i --include-samples sample.incl --region chr$i --out eQTLcircularRNAmerged.nominal.chr$i.txt --normal --seed 123 --cov covs.fastqtl.txt
  [ -e eQTLcircularRNAmerged.permutations.chr$i.txt ] || bsub -q short -n 1 -J eQTL.permute$i fastQTL --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/circularRNAmerged.expression.bed.gz -L _eQTLcircularRNAmerged.permute$i --include-samples sample.incl --region chr$i --permute 100 100000 --out eQTLcircularRNAmerged.permutations.chr$i.txt --normal --seed 123 --cov covs.fastqtl.txt
  [ -e eQTLcircularRNAmerged.permutations2.chr$i.txt ] || bsub -q short -n 1 -J eQTL.permute$i QTLtools cis --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/circularRNAmerged.expression.bed.gz --include-samples sample.incl --region chr$i --permute 1000 100000 --out eQTLcircularRNAmerged.permutations2.chr$i.txt --normal --seed 123 --cov covs.fastqtl.txt
done

# # eQTL for circRNA only (not ciRNA)
# for i in `seq 1 22` X; do
#   bsub -q short -n 1 -J eQTL.nominal$i fastQTL --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/circRNA.expression.bed.gz -L _eQTL.nominal$i --include-samples sample.incl --region chr$i --out eQTLcircRNA.nominal.chr$i.txt --normal --seed 123 --cov covs.fastqtl.txt
#   bsub -q short -n 1 -J eQTL.permute$i fastQTL --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/circRNA.expression.bed.gz -L _eQTL.permute$i --include-samples sample.incl --region chr$i --permute 100 100000 --out eQTLcircRNA.permutations.chr$i.txt --normal --seed 123 --cov covs.fastqtl.txt
# done

# eQTL for gene
for i in `seq 1 22` X; do
  bsub -q short -n 1 -J eQTLgene.nominal$i fastQTL --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/genes.expression.bed.gz -L _eQTLgene.nominal$i --include-samples sample.incl --region chr$i --out eQTLgene.nominal.chr$i.txt --normal --seed 123 --cov covs.fastqtl.txt
  bsub -q short -n 1 -J eQTLgene.permute$i fastQTL --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/genes.expression.bed.gz -L _eQTLgene.permute$i --include-samples sample.incl --region chr$i --permute 100 100000 --out eQTLgene.permutations.chr$i.txt --normal --seed 123 --cov covs.fastqtl.txt
done

## TODO: add chrX for phenotype/circRNA.PSI.chrX.bed.gz
# sQTL
for i in `seq 1 22`; do
  bsub -q short -n 1 -J sQTL.nominal$i fastQTL --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/circRNA.PSI.chr$i.bed.gz -L _sQTL.nominal$i --include-samples sample.incl --exclude-phenotypes PSI.excl --region chr$i --out sQTL.nominal.chr$i.txt --seed 123 --window 1e5 --cov covs.fastqtl.txt
  bsub -q short -n 1 -J sQTL.permute$i fastQTL --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/circRNA.PSI.chr$i.bed.gz -L _sQTL.permute$i --include-samples sample.incl --exclude-phenotypes PSI.excl --region chr$i --permute 100 100000 --out sQTL.permutations.chr$i.txt --seed 123 --window 1e5 --cov covs.fastqtl.txt
done

# cQTL
for i in `seq 1 22`; do
  bsub -q short -n 1 -J cQTL.nominal$i fastQTL --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/circRNA.PCI.bed.gz -L _cQTL.nominal$i --include-samples sample.incl --exclude-phenotypes PCI.excl --region chr$i --out cQTL.nominal.chr$i.txt --seed 123 --window 1e5 --cov covs.fastqtl.txt
  bsub -q short -n 1 -J cQTL.permute$i fastQTL --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/circRNA.PCI.bed.gz -L _cQTL.permute$i --include-samples sample.incl --exclude-phenotypes PCI.excl --region chr$i --permute 100 100000 --out cQTL.permutations.chr$i.txt --seed 123 --window 1e5 --cov covs.fastqtl.txt
done

ll *QTL*.permutations.chr* | wc -l
ll *QTL*.nominal.chr* | wc -l

# remove file with 0 size
find . -size  0 -delete

cat eQTL.permutations.chr*.txt | gzip > eQTL.permutations.txt.gz
cat eQTL.nominal.chr*.txt | gzip > eQTL.nominal.txt.gz
cat eQTLcircularRNAmerged.permutations.chr*.txt | gzip > eQTLcircularRNAmerged.permutations.txt.gz
cat eQTLcircularRNAmerged.nominal.chr*.txt | gzip > eQTLcircularRNAmerged.nominal.txt.gz
cat eQTLcircularRNAmerged.permutations2.chr*.txt | gzip > eQTLcircularRNAmerged.permutations2.txt.gz
cat eQTLgene.permutations.chr*.txt | gzip > eQTLgene.permutations.txt.gz
cat eQTLgene.nominal.chr*.txt | gzip > eQTLgene.nominal.txt.gz
cat sQTL.permutations.chr*.txt | gzip > sQTL.permutations.txt.gz
cat sQTL.nominal.chr*.txt | gzip > sQTL.nominal.txt.gz
cat cQTL.permutations.chr*.txt | gzip > cQTL.permutations.txt.gz
cat cQTL.nominal.chr*.txt | gzip > cQTL.nominal.txt.gz

rm *QTL*.permutations.chr* *QTL*.nominal.chr*
rm _*.log

# Rscript to eGene and QTL cutoff
for i in eQTL eQTLcircularRNAmerged eQTLcircRNA eQTLgene sQTL cQTL; do 
  echo $i; 
  Rscript ~/projects/circRNA/src/qtl/post_fastQTL.R $i > $i.post_fastQTL.log
  # get all SNP-gene pairs
  awk 'FNR==NR { if(NR>1) array[$1]=$15; next;} $1 in array { if($4<=array[$1]) print; }' <(zcat $i.egenes.txt.gz) <(zcat $i.nominal.txt.gz) | gzip -c > $i.allsignificantpairs.txt.gz
  
  # make manhattan plot
  p_cutoff=.01; [[ $i == "cQTL" ]] && p_cutoff=1
  Psignifiance=`grep threshold $i.post_fastQTL.log | cut -f2 -d':'`
  Rscript ~/projects/circRNA/src/qtl/_makeManhattanPlot.R -i $i.nominal.txt.gz -p $p_cutoff -P $Psignifiance -f fastQTL

done

## boxplot for top eQTL

zcat eQTL.nominal.txt.gz | sort -k4,4g | head
zcat eQTLcircRNA.nominal.txt.gz | sort -k4,4g | head
zcat eQTLcircularRNAmerged.nominal.txt.gz | awk '$4<1e-6' > eQTLcircularRNAmerged.nominal.txt.p1e-6.list

Rscript ~/projects/circRNA/src/qtl/_eQTL_boxplot.R -i eQTLcircularRNAmerged.nominal.txt.p1e-6.list -p ~/projects/circRNA/data/QTL_BC -g phenotype/circularRNAmerged.expression.bed.gz -e eQTLcircularRNAmerged.nominal.txt.gz -f fastQTL


# mark the genes on each tower
# focus on circRNA example
awk '$6<.01' final.cis.eQTL.xls | cut -f2 | fgrep -wf - ~/projects/circRNA/data/Merge_circexplorer_BC.annotation.bed14 | grep circRNA

zcat eQTL.nominal.txt.gz | awk '$4<1e-6' | join -1 1 -2 1 - <(cut -f4,13,14 ~/projects/circRNA/data/Merge_circexplorer_BC.annotation.bed14)


## RTC using QTLtools
# awk '{FS="\t"; OFS="\t"; print $1":"$3, $7}' ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.SNAP.LD_w250.r2_0.8.bed > GWAS.b37.txt
# wget http://jungle.unige.ch/QTLtools_examples/hotspots_b37_hg19.bed
QTLtools rtc --vcf genotype/chr15.dos.postQC.vcf.gz --bed phenotype/circularRNAmerged.expression.bed.gz --cov covs.fastqtl.txt --hotspot hotspots_b37_hg19.bed --gwas-cis GWAS.b37.txt eQTLcircularRNAmerged.permutations2.significant.txt --normal --out rtc_results.txt


####################################
## eQTL using matrix-eQTL
####################################
cd ~/projects/circRNA/results/eQTL/HCILB_SNDA/
Rscript ~/projects/circRNA/src/qtl/_SVA.eQTL.R HCILB_SNDA  # --> ~/projects/circRNA/results/eQTL/HCILB_SNDA/final.cis.eQTL.xls
Psignifiance=`sort -k5,5gr final.cis.eQTL.FDR.05.xls | head -n1 | cut -f5`
Rscript ~/projects/circRNA/src/qtl/_makeManhattanPlot.R -i final.cis.eQTL.xls -p 0.01 -P $Psignifiance -f matrixeQTL
# RTC
Rscript ~/projects/circRNA/src/qtl/_RTC.R ~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.txt expression.postSVA.xls final.cis.eQTL.FDR.05.xls circRNA  # --> final.cis.eQTL.FDR.05.xls.RTC.xls
Rscript ~/projects/circRNA/src/qtl/_RTC.R ~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.txt expression.postSVA.xls final.cis.eQTL.d1e6.p1e-2.xls circRNA  # --> final.cis.eQTL.d1e6.p1e-2.xls.RTC.xls

## TODO
# mark the genes on each tower
# focus on circRNA example
awk '$6<.01' final.cis.eQTL.xls | cut -f2 | fgrep -wf - ~/projects/circRNA/data/Merge_circexplorer_BC.annotation.bed14 | grep circRNA

####################################
### post-eQTL analysis
####################################
## Q1: Are the eQTL circRNAs more likely hosted in an eQTL gene?
## Q2: Are circRNA and their hostgene likely associated with the same set of SNPs?

####################################
### SNPs at the splicing sites
####################################
awk 'NR>1{OFS="\t"; print $2,$3-1,$3,$1,$4}' ~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID | intersectBed -a - -b <(awk '$13=="circRNA"{OFS="\t"; print $1,$2-2,$2,$4,($6=="+")?5:3,$6; print $1,$3,$3+2,$4,($6=="+")?3:5,$6;}' Merge_circexplorer_BC.annotation.bed14) -wo
# chr10	32807433	32807434	exm818678_G:T	0	chr10	32807433	32807435	chr10_32761408_32807433	3	+	1
# chr13	21735927	21735928	exm1056293_C:T	0	chr13	21735926	21735928	chr13_21735928_21742538	3	-	1
# chr13	21735927	21735928	exm1056293_C:T	0	chr13	21735926	21735928	chr13_21735928_21746820	3	-	1
# chr20	18142461	18142462	exm1527244_A:G	0	chr20	18142461	18142463	chr20_18142463_18162490	5	+	1
# chr8	124154696	124154697	exm718677_G:T	0	chr8	124154696	124154698	chr8_124089350_124154696	3	+	1
# chr8	124154696	124154697	exm718677_G:T	0	chr8	124154696	124154698	chr8_124141305_124154696	3	+	1
# chr8	124154696	124154697	exm718677_G:T	0	chr8	124154696	124154698	chr8_124153000_124154696	3	+	1
# chr5	68513571	68513572	rs2242351:68513572:C:G_C:G	1	chr5	68513570	68513572	chr5_68513572_68513698	5	+	1
awk 'NR>1{OFS="\t"; print $2,$3-1,$3,$1,$4}' ~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID | intersectBed -b - -a <(awk '$13=="circRNA"{OFS="\t"; print $1,$2-2,$2,$4,($6=="+")?5:3,$6; print $1,$3,$3+2,$4,($6=="+")?3:5,$6;}' Merge_circexplorer_BC.annotation.bed14) -wo | cut -f4,10 > SNPs.on.circRNA.splicingSites.list

