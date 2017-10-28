# script for QTL

### prepare genotype
## BC 
# # Tao's imputation + QC (remove imputation R2 <0.3, MAF < 0.05, multiple allelic SNPs) --> subset 84 sample --> --maf 0.05 --remove-filtered-all
# # on hms:
# cd hms:/home/tw83/twang/AMP/BrainCode/Imputation/QC/
# for i in chr*dos.postQC.vcf.gz; do echo $i; bsub -q short -W 12:00 -n 1 -J ${i/.*/} bash subset.sh $i; done
# # on hpc:
# cd ~/projects/circRNA/data/QTL/genotype
# rsync -azv xd20@transfer.orchestra.med.harvard.edu:/home/tw83/twang/AMP/BrainCode/Imputation/QC/84samples/chr* .

### prepare phenotype
cd ~/projects/circRNA/data/QTL_BC/phenotype
Rscript ~/projects/circRNA/src/QTL_prepare_phenotype.R
# sQTL
# BC
cd ~/neurogen/rnaseq_PD/results/sQTL
mkdir bamfiles142; cd bamfiles142
for i in /PHShome/xd010/neurogen/rnaseq_PD/run_output/*/accepted_hits.bam; do jj=${i/\/accepted_/.accepted_}; ln -fs $i `basename $jj`; done
cd ..
ls bamfiles142/* | grep -f ~/neurogen/rnaseq_PD/results/merged/samplelist.HCILB_SNDA - > bamfiles86.txt
ls bamfiles142/* > bamfiles142.txt
bsub -q long -n 1 -J leafcutter python $HOME/bin/leafcutter/example_data/run_sQTL.py -b bamfiles142.txt -d $HOME/bin/leafcutter -a $HOME/bin/leafcutter/clustering/gencode.v19.annotation.gtf.gz -t sqtl_bamfiles142

### QTL
cd ~/projects/circRNA/data/QTL_BC

# eQTL
for i in `seq 1 22`; do
  bsub -q short -n 1 -J eQTL.nominal$i fastQTL --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/circularRNA104.expression.bed.gz -L _eQTL.nominal$i --include-samples sample.incl --region $i --out eQTL.nominal.chr$i.txt.gz --normal --seed 123 --cov covs.fastqtl.txt
  bsub -q short -n 1 -J eQTL.permute$i fastQTL --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/circularRNA104.expression.bed.gz -L _eQTL.permute$i --include-samples sample.incl --region $i --permute 100 100000 --out eQTL.permutations.chr$i.txt.gz --normal --seed 123 --cov covs.fastqtl.txt
done
# eQTL for gene
for i in `seq 1 22`; do
  bsub -q short -n 1 -J eQTLgene.nominal$i fastQTL --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/genes.expression.bed.gz -L _eQTLgene.nominal$i --include-samples sample.incl --region $i --out eQTLgene.nominal.chr$i.txt.gz --normal --seed 123 --cov covs.fastqtl.txt
  bsub -q short -n 1 -J eQTLgene.permute$i fastQTL --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/genes.expression.bed.gz -L _eQTLgene.permute$i --include-samples sample.incl --region $i --permute 100 100000 --out eQTLgene.permutations.chr$i.txt.gz --normal --seed 123 --cov covs.fastqtl.txt
done

# sQTL
for i in `seq 1 22`; do
  bsub -q short -n 1 -J sQTL.nominal$i fastQTL --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/circRNA.PSI.chr$i.bed.gz -L _sQTL.nominal$i --include-samples sample.incl --exclude-phenotypes PSI.excl --region $i --out sQTL.nominal.chr$i.txt.gz --seed 123 --window 1e5 --cov covs.fastqtl.txt
  bsub -q short -n 1 -J sQTL.permute$i fastQTL --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/circRNA.PSI.chr$i.bed.gz -L _sQTL.permute$i --include-samples sample.incl --exclude-phenotypes PSI.excl --region $i --permute 100 100000 --out sQTL.permutations.chr$i.txt.gz --seed 123 --window 1e5 --cov covs.fastqtl.txt
done

# cQTL
for i in `seq 1 22`; do
  bsub -q short -n 1 -J cQTL.nominal$i fastQTL --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/circRNA.PCI.bed.gz -L _cQTL.nominal$i --include-samples sample.incl --exclude-phenotypes PCI.excl --region $i --out cQTL.nominal.chr$i.txt.gz --seed 123 --window 1e5 --cov covs.fastqtl.txt
  bsub -q short -n 1 -J cQTL.permute$i fastQTL --vcf genotype/chr$i.dos.postQC.vcf.gz --bed phenotype/circRNA.PCI.bed.gz -L _cQTL.permute$i --include-samples sample.incl --exclude-phenotypes PCI.excl --region $i --permute 100 100000 --out cQTL.permutations.chr$i.txt.gz --seed 123 --window 1e5 --cov covs.fastqtl.txt
done

ll *QTL*.permutations.chr* | wc -l
ll *QTL*.nominal.chr* | wc -l

cat eQTL.permutations.chr*.txt.gz > eQTL.permutations.txt.gz
cat eQTL.nominal.chr*.txt.gz > eQTL.nominal.txt.gz
cat eQTLgene.permutations.chr*.txt.gz > eQTLgene.permutations.txt.gz
cat eQTLgene.nominal.chr*.txt.gz > eQTLgene.nominal.txt.gz
cat sQTL.permutations.chr*.txt.gz > sQTL.permutations.txt.gz
cat sQTL.nominal.chr*.txt.gz > sQTL.nominal.txt.gz
cat cQTL.permutations.chr*.txt.gz > cQTL.permutations.txt.gz
cat cQTL.nominal.chr*.txt.gz > cQTL.nominal.txt.gz

rm *QTL*.permutations.chr* *QTL*.nominal.chr*
rm _*.log

# Rscript to eGene and QTL cutoff
for i in eQTL eQTLgene sQTL cQTL; do 
  echo $i; 
  Rscript ~/projects/circRNA/src/post_fastQTL.R $i > $i.post_fastQTL.log
  # get all SNP-gene pairs
  awk 'FNR==NR { if(NR>1) array[$1]=$15; next;} $1 in array { if($4<=array[$1]) print; }' <(zcat $i.egenes.txt.gz) <(zcat $i.nominal.txt.gz) | gzip -c > $i.allsignificantpairs.txt.gz
  
  # make manhattan plot
  Rscript ~/projects/circRNA/src/_makeManhattanPlot.R $i.nominal.txt.gz

done


zcat eQTL.nominal.txt.gz | sort -k4,4g | head