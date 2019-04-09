#!/usr/bin/env bash
# align FASTQ using STAR and then call circRNA using CIRCExplore
# Note: this works for datasets: MSBB
# Aug 18, 2017
# Modified by XD
# Usage: 
# for i in Fastq/*gz; do sample=`basename $i .accepted_hits.sort.coord.combined.fastq.gz`; [ ! -s outputs/$sample/$sample.circ.txt ] && bsub -q mcore -n 4 -W 12:00 -J $sample -R "rusage[mem=10000]" bash run-star.sh $i; done

# $1 = gzipped FASTQ file with unmapped reads (full path optional; containing
#      folder hardcoded below)

module load seq/STAR/2.5.3a

sample=`basename $1 .accepted_hits.sort.coord.combined.fastq.gz`

# Define paths
rootdir="/home/xd20/neurogen/circRNA/MSBB"
fastqdir="${rootdir}/Fastq"

# Reference files
index="$GENOME/Homo_sapiens/UCSC/hg19/Sequence/starIndex"

# Specify output folder
outdir="${rootdir}/run_output/${sample}"
if [[ ! -e "$outdir" ]]; then
    mkdir -p $outdir
fi
cd $outdir

# Align reads with STAR
[ ! -s Chimeric.out.junction ] && \
STAR \
    --chimSegmentMin 10 \
    --chimJunctionOverhangMin 10 \
    --runMode alignReads \
    --runThreadN 4 \
    --genomeDir $index \
    --sjdbGTFfile "$GENOME/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf" \
    --readFilesIn "${fastqdir}/${sample}.accepted_hits.sort.coord.combined.fastq.gz" \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped Within \
    --quantMode GeneCounts \
    --twopassMode Basic \
    --readFilesCommand zcat

# convert Chimeric.out.junction to fusion_junction.txt
[ ! -f fusion_junction.txt ] && \
star_parse.py Chimeric.out.junction fusion_junction.txt

# parse fusion_junction.txt
# update with new reference file (add CDR1as) - 20171026
[ ! -f ${sample}_circ.txt ] && \
CIRCexplorer.py \
  -j fusion_junction.txt \
  -g $index/genome.fa \
  -r /home/xd20/neurogen/circRNA/MSBB/refFlat_plus_CDR1as.txt \
  -o $sample
  
# ## After STAR is done, copy the output *bam and *circ.txt to HPC and sort/index there
# cd ~/neurogen/ROSMAP/MSBB/rnaseq_runoutput
# rsync -azv xd20@transfer.orchestra.med.harvard.edu:~/projects/circRNA/MSBB/run_output/*/*circ.txt .
# rsync -azv xd20@transfer.orchestra.med.harvard.edu:~/projects/circRNA/MSBB/run_output/*/*.out.bam .
# # if not sorted/indexed
# for i in *.out.bam; do echo $i; bsub -q short -n 3 -J `basename $i .Aligned.out.bam` "samtools sort -m 2G -@ 3 -o $i.sorted $i; mv $i.sorted $i; samtools index $i;"; done