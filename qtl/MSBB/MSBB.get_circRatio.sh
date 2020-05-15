# script to compute the circularization ratio of circRNA
cd ~/projects/circRNA/data/

circRNA_annotation=Merge_circexplorer_MSBB.annotation.bed14  # MSBB

##========================================================================
# 1. get the 1nt--1nt around the two ends of circRNA
##========================================================================
awk '{OFS="\t"; s=($6=="+")?$2:$3; e=($6=="+")?$3:$2; if(NR>1) {print $1,s-1,s+1,$4; print $1,e-1,e+1,$4;}}' $circRNA_annotation | sortBed -faidx hg19.hms.genome > $circRNA_annotation.2ntAtEnds.bed

# ##========================================================================
# # 1. get linear spliced reads on the two ends of circRNA: s3+s5
# ##========================================================================
# > $circRNA_annotation.sum_s3s5
# for i in ~/neurogen/rnaseq_PD/run_output/*/junctions.bed; do
#   sampleName=`echo $i | sed 's/.*output\/\(.*\)\/junction.*/\1/g'`
#   echo "Processing $sampleName ..."
#   # convert junctions.bed from Tophat to 1nt--1nt pattern
#   cat $i | \
#   awk '{OFS="\t"; split($11,a,","); split($12,b,","); $2=$2+a[1]-1; $3=$3-a[2]+1; $7=$2;$8=$3; $11="1,1"; $12="0,"(b[2]-a[1]+1); if($1!="track") print}' | \
#   intersectBed -b - -a $circRNA_annotation.2ntAtEnds.bed -split -wao | \
#   cut -f4,9,17 | awk -v sampleName=$sampleName '{OFS="\t"; if($2==-1) $2=0;print $1,sampleName,$2,$3}' | \
#   /apps/lib-osver/bedtools/2.23.0/bin/bedtools groupby -g 1,2 -c 3,4 -o sum,sum >> $circRNA_annotation.sum_s3s5
# done

##========================================================================
# 2. get the linear spliced and unspliced reads on the two ends of circRNA: s3+s5+u3+u5
##========================================================================
## Note: *sorted* input bam files are required. So, either change STAR option "--outSAMtype BAM SortedByCoordinate" or sort bam with samtools sort
## cd ~/neurogen/AMPAD/MSBB/rnaseq_runoutput/; for i in `find *.Aligned.out.bam -mtime -365 -mtime +2`;  do echo $i; bsub -q short -n 3 -J `basename $i .Aligned.out.bam` "samtools sort -m 2G -@ 3 -o $i.sorted $i; mv $i.sorted $i; samtools index $i;"; done
for i in ~/neurogen/AMPAD/MSBB/rnaseq_runoutput/*out.bam; do
  sampleName=`echo $i | sed 's/.*output\/\(.*\)\..*Aligned.*/\1/g'`
  [ -s $circRNA_annotation.sum_u3u5s3s5.$sampleName ] || bsub -q big-multi -n 4 -M 8000 -J $sampleName "module load bedtools/2.26.0; bedtools coverage -split -b $i -a $circRNA_annotation.2ntAtEnds.bed -sorted -counts -g hg19.hms.genome | LC_ALL=C sort -k4,4 --parallel=4 --buffer-size=1G | ~/bin/bedtools223 groupby -g 4 -c 5 -o sum > $circRNA_annotation.sum_u3u5s3s5.$sampleName"
done
# merge all together, for individual tissue
cd MSBB.frontal_pole
echo "geneID" `ls $circRNA_annotation.sum_u3u5s3s5.* | sed 's/.*3u5s3s5.//g' | rowsToCols stdin stdout` > header
tmp=$(mktemp)
tmp2=$(mktemp)
files=($circRNA_annotation.sum_u3u5s3s5.*)
cp "${files[0]}" $tmp2
for file in "${files[@]:1}"; do
    echo $file
    join --nocheck-order $tmp2 $file > $tmp && mv $tmp $tmp2
done
cat header $tmp2 | tr ' ' '\t' > $circRNA_annotation.sum_u3u5s3s5
rm $circRNA_annotation.sum_u3u5s3s5.*

##========================================================================
# 3. get back-spliced reads on the circRNA
##========================================================================
cat ${circRNA_annotation/annotation.bed14/}rawcount.txt | awk '{OFS="\t"; if(NR>1) $4=$1"_"$2"_"$3; print}' | cut -f4,7- | sed 's/_candidates.bed_circReads//g' > $circRNA_annotation.circReads.txt

##========================================================================
# 4. get circulization ratio == back-spliced reads / (back-spliced reads + linear reads)
##========================================================================
Rscript ../src/get_circRatio.R $circRNA_annotation
