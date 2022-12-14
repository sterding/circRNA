# script to compute the circularization ratio of circRNA
cd ~/projects/circRNA/data/

#circRNA_annotation=Merge_circexplorer_BC.annotation.bed14  # BRAINCODE only
#circRNA_annotation=Merge_circexplorer_BC_RM.annotation.bed14  # BRAINCODE+RM
#circRNA_annotation=Merge_circexplorer_BC.annotation.bed14.matched  # BC exon number- and length-matched controls
circRNA_annotation=exons.internal.meta.pc.bed  # internal meta exons in the genome background

##========================================================================
# 1. get the 1nt--1nt around the two ends of circRNA
##========================================================================
#awk '{OFS="\t"; s=($6=="+")?$2:$3; e=($6=="+")?$3:$2; if(NR>1) {print $1,s-1,s+1,$4; print $1,e-1,e+1,$4;}}' $circRNA_annotation | sortBed -faidx hg19.hpc.genome > $circRNA_annotation.2ntAtEnds.bed
# use 5 to indicate 5' and 3 for 3' 
awk '{OFS="\t"; s=($6=="+")?$2:$3; e=($6=="+")?$3:$2; if(NR>1) {print $1,s-1,s+1,$4,5,$6; print $1,e-1,e+1,$4,3,$6;}}' $circRNA_annotation | sortBed -faidx hg19.hpc.genome > $circRNA_annotation.2ntAtEnds.bed6  

# ##========================================================================
# # 1. get linear spliced reads on the two ends of circRNA: s3+s5
# ##========================================================================
# ## Note: junctions.bed has strand information
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

## record s5,s3,u5,s3 separately
# s3, s5
## Note: below code only consider the exon-circularized circRNA, not the intron lariat
> $circRNA_annotation.s3s5
for i in ~/neurogen/rnaseq_PD/run_output/*/junctions.bed; do
  sampleName=`echo $i | sed 's/.*output\/\(.*\)\/junction.*/\1/g'`
  echo "Processing $sampleName ..."
  # convert junctions.bed from Tophat to 1nt--1nt pattern
  cat $i | \
  awk '{OFS="\t"; split($11,a,","); split($12,b,","); $2=$2+a[1]-1; $3=$3-a[2]+1; $7=$2;$8=$3; $11="1,1"; $12="0,"(b[2]-a[1]+1); if($1!="track") print}' | \
  intersectBed -b - -a $circRNA_annotation.2ntAtEnds.bed6 -split -s -wo | \
  awk -v sampleName=$sampleName '{OFS="\t"; if(($6=="+" && $5==5 && $3==$9) || ($6=="+" && $5==3 && $2==$8) || ($6=="-" && $5==5 && $2==$8) || ($6=="-" && $5==3 && $3==$9)) print $4,"s"$5,sampleName,$11,$19}' | sort -k1,1 -k2,2 | \
  /apps/lib-osver/bedtools/2.23.0/bin/bedtools groupby -g 1,2,3 -c 4,5 -o sum,sum >> $circRNA_annotation.s3s5
done
gzip $circRNA_annotation.s3s5

# u3, u5
## Note a bug in "bedtools coverage -split -f -counts" (see https://github.com/arq5x/bedtools2/issues/673)
# for i in ~/neurogen/rnaseq_PD/run_output/*/accepted_hits.bam; do
#   sampleName=`echo $i | sed 's/.*output\/\(.*\)\/accepted_hits.*/\1/g'`
#   [ -s $TMPDIR/$circRNA_annotation.u3u5.$sampleName ] || bsub -q big -n 1 -J $sampleName "module load bedtools/2.26.0; bedtools coverage -split -b $i -a $circRNA_annotation.2ntAtEnds.bed6 -sorted -counts -g hg19.hpc.genome -f 1.0 > $TMPDIR/$circRNA_annotation.u3u5.$sampleName"
# done

# alternative way
awk -vi=$circRNA_annotation.2ntAtEnds.bed6 '{print > i"."$1}' $circRNA_annotation.2ntAtEnds.bed6
for i in ~/neurogen/rnaseq_PD/run_output/*/accepted_hits.bam; do
  sampleName=`echo $i | sed 's/.*output\/\(.*\)\/accepted_hits.*/\1/g'`
  for y in `seq 1 22` X Y; do
    chr="chr"$y;
    [ -s $TMPDIR/$circRNA_annotation.2ntAtEnds.bed6.$sampleName.$chr ] || bsub -q short -n 1 -J $sampleName.$chr bash ~/neurogen/pipeline/RNAseq/bin/_bam_over_bed.sh $i $circRNA_annotation.2ntAtEnds.bed6 hg19.hpc.genome $chr    
  done
done
# combine
for i in ~/neurogen/rnaseq_PD/run_output/*/accepted_hits.bam; do
  sampleName=`echo $i | sed 's/.*output\/\(.*\)\/accepted_hits.*/\1/g'`; echo $sampleName
  cat $TMPDIR/$circRNA_annotation.2ntAtEnds.bed6.$sampleName.chr* > $TMPDIR/$circRNA_annotation.2ntAtEnds.bed6.$sampleName && \
  rm $TMPDIR/$circRNA_annotation.2ntAtEnds.bed6.$sampleName.chr*
done

> $circRNA_annotation.u3u5
find $TMPDIR/$circRNA_annotation.2ntAtEnds.bed6.* -cmin -160 | while read i; do 
  sampleName=`echo $i | sed 's/.*bed6.\(.*\)/\1/g'`; echo $sampleName
  awk -vsample=$sampleName '{OFS="\t"; if($7>0) print $4,"u"$5,sample,$7,1}' $i >> $circRNA_annotation.u3u5
done
gzip $circRNA_annotation.u3u5

##========================================================================
# 2. get the linear spliced and unspliced reads on the two ends of circRNA: s3+s5+u3+u5
##========================================================================
# for i in ~/neurogen/rnaseq_PD/run_output/[HI]*/accepted_hits.bam; do
#   sampleName=`echo $i | sed 's/.*output\/\(.*\)\/accepted_.*/\1/g'`
#   #[ -s $circRNA_annotation.sum_u3u5s3s5.$sampleName ] || bsub -q normal -n 1 -J $sampleName "bedtools multicov -split -bams $i -bed $circRNA_annotation.2ntAtEnds.bed | sort -k4,4 | bedtools223 groupby -g 4 -c 5 -o sum > $circRNA_annotation.sum_u3u5s3s5.$sampleName"
#   [ -s $circRNA_annotation.sum_u3u5s3s5.$sampleName ] || bsub -q normal -n 1 -R "rusage[mem=1000]" -J $sampleName "module load bedtools/2.26.0; bedtools coverage -split -b $i -a $circRNA_annotation.2ntAtEnds.bed -sorted -counts -g hg19.hpc.genome | sort -k4,4 | ~/bin/bedtools223 groupby -g 4 -c 5 -o sum > $circRNA_annotation.sum_u3u5s3s5.$sampleName"
# done

ls ~/neurogen/rnaseq_PD/run_output/*/accepted_hits.bam | grep -f ~/neurogen/circRNA_seq_Rebeca_HC/BRAINCODE_circexp/BC.n125.samplelist - | while read i; do
  sampleName=`echo $i | sed 's/.*output\/\(.*\)\/accepted_.*/\1/g'`
  [ -s $circRNA_annotation.sum_u3u5s3s5.$sampleName ] || bsub -q big -n 1 -R "rusage[mem=16000:swp=16000]" -J $sampleName "module load bedtools/2.26.0; bedtools coverage -split -b $i -a $circRNA_annotation.2ntAtEnds.bed -sorted -counts -g hg19.hpc.genome | sort -k4,4 | ~/bin/bedtools223 groupby -g 4 -c 5 -o sum > $circRNA_annotation.sum_u3u5s3s5.$sampleName"
done

# merge all together
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
Rscript ../src/qtl/get_circRatio.R $circRNA_annotation ~/neurogen/rnaseq_PD/results/merged/samplelist.HCILB_SNDA84
#Rscript ../src/get_circRatio.R $circRNA_annotation ~/neurogen/circRNA_seq_Rebeca_HC/BRAINCODE_circexp/BC.n125.samplelist
