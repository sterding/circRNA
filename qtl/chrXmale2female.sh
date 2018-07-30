# script for convert manoploit male chrX to a pseudo diploit chrX (like female) in order to combine the chrX male and chrX female.
# usage: bash chrXmale2female.sh chrX.male.dos.postQC.vcf.gz chrX.female.dos.postQC.vcf.gz | bgzip -c > chrX.dos.postQC.vcf.gz

chrXmale_vcf_gz=$1
chrXfemale_vcf_gz=$2

zcat $chrXmale_vcf_gz | awk '{ OFS="\t";
if($0 ~ /^#/) {
 if($0 ~ /ID=DS/) $0="##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]\">";
 if($0 ~ /ID=GP/) $0="##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1\">";
}else{
  for(i=10;i<=NF;i++){
    split($i, a, ":");
    GT=a[1];DS=a[2];GP=a[3]; # Our data have GT:DS:GP format. This should be changed if a different format
    GT=GT"/"GT;
    DS=DS*(1-DS) + 2*DS*DS;
    split(GP,b,","); GP=b[1]",0,"b[2];
    $i=GT":"DS":"GP;
  }
}
  print $0;
}' |  bgzip -c > /tmp/chrXmale.vcf.gz
vcf-merge /tmp/chrXmale.vcf.gz $chrXfemale_vcf_gz