# script to make exon length-matched control set of circRNAs
cd ~/projects/circRNA/data/

circRNA_annotation=Merge_circexplorer_BC197.filtered.enriched.annotation.bed14  # BRAINCODE only
#circRNA_annotation=Merge_circexplorer_BC_RM.annotation.bed14  # BRAINCODE+RM

genes_annotation=~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v19.annotation.bed12

#fgrep protein_coding___protein_coding $genes_annotation | awk '{OFS="\t"; split($11,L,",");split($12,S,","); N=$10; for(n=1;n<=N;n++) for(i=1;(i+n-1)<=N;i++) {l=0; s=S[i];e=0;for(j=1;j<=n;j++) {l=l+L[j+i-1];e=L[j+i-1]+S[i+j-1]} print n,l,$1"_"(s+$2)"_"($2+e),$4;}}' > $genes_annotation.protein_coding.allExonCombinations.txt

[ -e $genes_annotation.protein_coding.allExonCombinations.txt ] || fgrep protein_coding___protein_coding $genes_annotation | awk '{OFS="\t"; split($11,L,",");split($12,S,","); N=$10; split($4,G,"___"); for(n=1;n<=N;n++) for(i=1;(i+n-1)<=N;i++) {ll=""; ss=""; s=S[i];e=0;for(j=1;j<=n;j++) {ll=ll""L[j+i-1]","; ss=ss""(S[i+j-1]-s)","; e=L[j+i-1]+S[i+j-1]} print $1,s+$2,$2+e,$1"_"(s+$2)"_"($2+e),$5,$6,s+$2,$2+e,$9,n,ll,ss,i"/"N,G[2];}}' > $genes_annotation.protein_coding.allExonCombinations.txt

cut -f4 $circRNA_annotation | fgrep -v -f - $genes_annotation.protein_coding.allExonCombinations.txt | sort -u > $circRNA_annotation.background

############################
# control1: controls with matched numbers of exon and total lenght of exons
############################
# make a fake BED file using exonCount as chr and exonLength as start, then find the cloest feature in the background for each circRNA. 
# In that way, we match circRNAs with the closest combination of exons with the same number of exons and similar length of exons 
bedtools closest -a <(awk '{OFS="\t"; if(NR>1) {split($11,a,",");s=0;for(i in a) s+=a[i]; print "chr"$10,s-1,s,$0;}}' $circRNA_annotation | cut -f1-17 | sortBed) -b <(awk '{OFS="\t"; if(NR>1) {split($11,a,",");s=0;for(i in a) s+=a[i]; print "chr"$10,s-1,s,$0;}}' $circRNA_annotation.background | sortBed) -t all > /scratch/x/xd010/$circRNA_annotation.matchedall

## pick up one randomly from the ties
## HOW: add a random number at the end and then sort by  ID, then the rand number ...
awk 'BEGIN{srand(1);} {printf "%06d\t%s\n", rand()*1000000, $0;}' /scratch/x/xd010/$circRNA_annotation.matchedall | LC_ALL=C sort --parallel=6 --buffer-size=1G -k8,8 -k1,1n | awk '{OFS="\t";if($8!=id) {print $0,$8;id=$8;}}' | cut -f22- > $circRNA_annotation.matched

## Next
# run get_circRatio.sh for $circRNA_annotation.matched to get $circRNA_annotation.matched.u3u5 and $circRNA_annotation.matched.s3s5

############################
# control2: controls with matched numbers of exon and matched lenght of region
############################
# make a fake BED file using exonCount as chr and exonLength as start, then find the cloest feature in the background for each circRNA. 
# choose the first hit from the ties
# In that way, we match circRNAs with the closest combination of exons with the same number of exons and similar length of region 
bedtools closest -a <(awk 'NR>1 {OFS="\t"; print "chr"$10,$3-$2,$3-$2+1,$0;}' $circRNA_annotation | cut -f1-17 | sortBed) -b <(awk 'NR>1{OFS="\t"; print "chr"$10,$3-$2,$3-$2+1,$0;}' $circRNA_annotation.background | sortBed) -t all > /scratch/x/xd010/$circRNA_annotation.matched2all
awk 'BEGIN{srand(1);} {printf "%06d\t%s\n", rand()*1000000, $0;}' /scratch/x/xd010/$circRNA_annotation.matched2all | LC_ALL=C sort --parallel=6 --buffer-size=1G -k8,8 -k1,1n | awk '{OFS="\t";if($8!=id) {print $0,$8;id=$8;}}' | cut -f22- > $circRNA_annotation.matched2

# alternative: exact length match: many circRNAs cannot find exact lenght-mathced controls
# bedtools closest -a <(awk 'NR>1 {OFS="\t"; len=$3-$2; print "chr"len,$10,$10+1,$0;}' $circRNA_annotation | sortBed) -b <(awk 'NR>1{OFS="\t"; len=$3-$2; print "chr"len,$10,$10+1,$0;}' $circRNA_annotation.background | sortBed) -t first | awk '{OFS="\t";print $0,$7;}' | cut -f21- > $circRNA_annotation.matched3

############################
# control3: exons from the same gene
############################
# make a fake BED file using genename as chr and total length as start, then find the cloest feature in the background for each circRNA. 
# choose the first hit from the ties
bedtools closest -a <(awk 'NR>1 {OFS="\t"; print "chr"$14,$3-$2,$3-$2+1,$0;}' $circRNA_annotation | cut -f1-17 | sortBed) -b <(awk '{OFS="\t"; print "chr"$14,$3-$2,$3-$2+1,$0;}' $circRNA_annotation.background | sortBed) -t all > /scratch/x/xd010/$circRNA_annotation.matched3all
awk 'BEGIN{srand(1);} $18!="." {printf "%06d\t%s\n", rand()*1000000, $0;}' /scratch/x/xd010/$circRNA_annotation.matched3all | LC_ALL=C sort --parallel=6 --buffer-size=1G -k8,8 -k1,1n | awk '{OFS="\t";if($8!=id) {print $0,$8;id=$8;}}' | cut -f22- > $circRNA_annotation.matched3
