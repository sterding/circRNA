## # script for QTL for MSBB dataset @ hms

cd ~/projects/circRNA/MSBB/
echo individualIdentifier DNA_SampleName RNA_SampleName | sed 's/ /\t/g' > Samples_hasGenoExpr.ID.DNA.RNA.tab 
join -1 2 -2 1 <(cut -f2,3 Genotype/Plink/Imputed/Download/QC/Samples_hasGenoExpr.txt | sort -k2,2) <(cat Job_RIN5_FrontalPole_Aug1717.tsv| sed 's/ /_/g;s/"//g' | cut -f8,14 | sort -k1,1) | sed 's/ /\t/g' >> Samples_hasGenoExpr.ID.DNA.RNA.tab


### prepare genotype
### phenotype
cd ~/projects/circRNA/MSBB/genotype
for i in ../Genotype/Plink/Imputed/Download/chr*dose.vcf.gz*; do echo $i; ln -fs $i `basename $i`; done
## PSI
cd ~/projects/circRNA/MSBB/phenotypes/leafcutter
ls ~/projects/circRNA/MSBB/run_output/*/*Aligned.out.bam > bamfiles.txt
bsub -q medium -n 2 -W 120:00 -J leafcutter python $HOME/leafcutter/example_data/run_sQTL.py -b bamfiles.txt -d $HOME/leafcutter -a $HOME/leafcutter/clustering/gencode.v19.annotation.gtf.gz -t sqtl_tmpdir
cd ~/projects/circRNA/MSBB/phenotypes
for i in `seq 1 22`; do echo $i;  cat leafcutter/sqtl_tmpdir/leafcutter_perind.counts.gz.qqnorm_chr$i | sed 's/.accepted_hits.bam//g;s/.Aligned.out.bam//g;s/HC_//g;s/ILB_//g;s/_SNDA_[1-9]_rep[12]//g' | awk '{OFS="\t"; s=0; for(i=5;i<=NF;i++) if($i>0) s++; if(NR==1 || s>=1) {if(NR>1 && $1 !~/^chr/) {$1="chr"$1; $4="chr"$4;} print}}' > circRNA.PSI.chr$i.bed & done
for i in *bed; do echo $i; bgzip -f $i && tabix -f -p bed $i.gz; done
## PCI
cd ~/projects/circRNA/MSBB/
#1. run STAR --> circExplorer --> circ.txt in each run_output folder (done)
for i in Fastq/*/*gz; do sample=`basename $i .accepted_hits.sort.coord.combined.fastq.gz`; [ ! -s outputs/$sample/circ.txt_circ.txt ] && bsub -q mcore -n 4 -W 12:00 -J $sample -R "rusage[mem=10000]" bash run-star.sh $i; done
#2. merge all circ.txt into one big matrix of circ expression, only for samples with Genotype+phrenotyp
mkdir phenotypes/circexplorer; cd phenotypes/circexplorer;
ls ~/projects/circRNA/MSBB/run_output/*/*circ.txt | grep -f <(cut -f3 ~/projects/circRNA/MSBB/Samples_hasGenoExpr.ID.DNA.RNA.tab) - | while read line; do i=${line/*\//}; ii=${i/circ.txt/candidates.bed}; echo $ii; ln -fs $line $ii; done
bsub -q medium -n 1 -R "rusage[mem=20000]" -W 3:00 -J circ_merge "python $HOME/projects/circRNA/MSBB/merge_circexplorer.py *candidates.bed > ~/projects/circRNA/MSBB/Merge_circexplorer_MSBB.rawcount.txt"
echo -e "chrom\tstart\tend\tID\tscore\tstrand\tthickStart\tthickEnd\titemRgb\texonCount\texonSizes\texonOffsets\tcircType\tgeneName" > ~/projects/circRNA/MSBB/Merge_circexplorer_MSBB.annotation.bed14
awk -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3,1,$6,$7,$8,$9,$10,$11,$12,($14=="Yes")?"ciRNA":"circRNA",$15}' *_candidates.bed | sort -u >>  ~/projects/circRNA/MSBB/Merge_circexplorer_MSBB.annotation.bed14
#3. in hpc: run main.MSBB.R, then QTL_prepare_phenotype.R to get circRNA.PCI.bed, circRNA.expression.bed etc. 
ssh to hpc
rsync -azv xd20@transfer.orchestra.med.harvard.edu:~/projects/circRNA/MSBB/Merge_circexp* ~/projects/circRNA/data
rsync -azv xd20@transfer.orchestra.med.harvard.edu:~/projects/circRNA/MSBB/Metadata ~/neurogen/ROSMAP/MSBB
rsync -azv xd20@transfer.orchestra.med.harvard.edu:~/projects/circRNA/MSBB/Samples_hasGenoExpr.ID.DNA.RNA.tab ~/neurogen/ROSMAP/MSBB


### QTL
cd ~/projects/circRNA/MSBB/

# eQTL
for i in `seq 1 22`; do
  bsub -q short -n 1 -J eQTL.nominal$i fastQTL --vcf genotype/chr$i.dose.vcf.gz --bed phenotype/circRNA.expression.bed.gz -L _eQTL.nominal$i --include-samples sample.incl --region $i --out eQTL.nominal.chr$i.txt.gz --normal --seed 123 --cov covs.fastqtl.txt
  bsub -q short -n 1 -J eQTL.permute$i fastQTL --vcf genotype/chr$i.dose.vcf.gz --bed phenotype/circRNA.expression.bed.gz -L _eQTL.permute$i --include-samples sample.incl --region $i --permute 100 100000 --out eQTL.permutations.chr$i.txt.gz --normal --seed 123 --cov covs.fastqtl.txt
done

# sQTL
for i in `seq 1 22`; do
  bsub -q short -n 1 -J sQTL.nominal$i fastQTL --vcf genotype/chr$i.dose.vcf.gz --bed phenotype/circRNA.PSI.chr$i.bed.gz -L _sQTL.nominal$i --include-samples sample.incl --exclude-phenotypes PSI.excl --region $i --out sQTL.nominal.chr$i.txt.gz --seed 123 --window 1e5 --cov covs.fastqtl.txt
  bsub -q short -n 1 -J sQTL.permute$i fastQTL --vcf genotype/chr$i.dose.vcf.gz --bed phenotype/circRNA.PSI.chr$i.bed.gz -L _sQTL.permute$i --include-samples sample.incl --exclude-phenotypes PSI.excl --region $i --permute 100 100000 --out sQTL.permutations.chr$i.txt.gz --seed 123 --window 1e5 --cov covs.fastqtl.txt
done

# cQTL
for i in `seq 1 22`; do
  bsub -q short -n 1 -J cQTL.nominal$i fastQTL --vcf genotype/chr$i.dose.vcf.gz --bed phenotype/circRNA.PCI.bed.gz -L _cQTL.nominal$i --include-samples sample.incl --exclude-phenotypes PCI.excl --region $i --out cQTL.nominal.chr$i.txt.gz --seed 123 --window 1e5 --cov covs.fastqtl.txt
  bsub -q short -n 1 -J cQTL.permute$i fastQTL --vcf genotype/chr$i.dose.vcf.gz --bed phenotype/circRNA.PCI.bed.gz -L _cQTL.permute$i --include-samples sample.incl --exclude-phenotypes PCI.excl --region $i --permute 100 100000 --out cQTL.permutations.chr$i.txt.gz --seed 123 --window 1e5 --cov covs.fastqtl.txt
done

ll *QTL.permutations.chr* | wc -l
ll *QTL.nominal.chr* | wc -l

cat eQTL.permutations.chr*.txt.gz > eQTL.permutations.txt.gz
cat eQTL.nominal.chr*.txt.gz > eQTL.nominal.txt.gz
cat sQTL.permutations.chr*.txt.gz > sQTL.permutations.txt.gz
cat sQTL.nominal.chr*.txt.gz > sQTL.nominal.txt.gz
cat cQTL.permutations.chr*.txt.gz > cQTL.permutations.txt.gz
cat cQTL.nominal.chr*.txt.gz > cQTL.nominal.txt.gz

rm *QTL.permutations.chr* *QTL.nominal.chr*
rm _*.log

# Rscript to eGene and QTL cutoff
for i in eQTL sQTL cQTL; do 
  echo $i; 
  Rscript ~/projects/circRNA/src/post_fastQTL.R $i > $i.post_fastQTL.log
  # get all SNP-gene pairs
  awk 'FNR==NR { if(NR>1) array[$1]=$15; next;} $1 in array { if($4<=array[$1]) print; }' <(zcat $i.egenes.txt.gz) <(zcat $i.nominal.txt.gz) | gzip -c > $i.allsignificantpairs.txt.gz
done