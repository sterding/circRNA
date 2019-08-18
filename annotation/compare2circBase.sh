#### comparsion between circBase, MiOncoCirc, and circRRAIN
cd ~/projects/circRNA/data/

## circBase

curl -s http://www.circbase.org/download/hsa_hg19_circRNA.txt | awk 'NR>1{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3}' | sort -u > circBase.hg19.bed 

## MiOncoCirc

# download 
# go to https://mioncocirc.github.io/download/ --> MiOncoCirc.v0.1.release.txt.gz  (coordinates in bed format already, reads >= 2, 2000+ libraries)
# liftover to hg19
zcat MiOncoCirc.v0.1.release.txt.gz | awk 'NR>1{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3"_hg38";}' | sort | uniq -c | awk '{OFS="\t"; print $2,$3,$4,$5,$1;}' | liftOver stdin ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/GWASCatalog/hg38ToHg19.over.chain.gz MiOncoCirc.v0.1.release.txt.hg19.bed unmapped

BC=Merge_circexplorer_BC_RM.annotation.bed14