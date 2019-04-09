#!/usr/bin/env bash
# take input of unmapped.bam and then call circRNA using CIRCExplore
# Note: this works for datasets: CommonMind
# Nov 13, 2018
# Modified by XD
# Usage: 
# cd /home/xd20/neurogen/circRNA/CMC
# for i in RNAseq/*.unmapped.bam; do sample=`basename $i .unmapped.bam`; [ ! -s /n/scratch2/xd20/$sample/${sample}_circ.txt ] && sbatch -p short -c 4 -t 12:00:00 -J $sample -o /n/scratch2/xd20/$sample.log -e /n/scratch2/xd20/$sample.log --open-mode=append --mem-per-cpu=8G --wrap="bash ~/pipeline/circRNA/annotation/unmappedBAM2circRNA.sh $i"; done
# cp /n/scratch2/xd20/*/*circ.txt /home/xd20/neurogen/circRNA/CMC/circRNA

# $1 = gzipped FASTQ file with unmapped reads (full path optional; containing
#      folder hardcoded below)

#module load tophat/2.1.1  # use tophat v2.0.10 instead
#module load bowtie/1.2.1.1  # only bowtie v1.0.1 works
module load samtools/0.1.19

sample=`basename $1 .unmapped.bam`

# Define paths
rootdir="/home/xd20/neurogen/circRNA/CMC"
fastqdir="${rootdir}/RNAseq"

# Specify output folder
outdir="/n/scratch2/xd20/${sample}"
if [[ ! -e "$outdir" ]]; then
    mkdir -p $outdir
fi
cd $outdir

# Convert unmapp.bam to fastq
[ ! -s ${sample}.unmapped.fastq ] && \
bamToFastq -i ${fastqdir}/${sample}.unmapped.bam -fq ${sample}.unmapped.fastq

echo "# Align unmapp reads with tophat_fusion"

[ ! -s tophat_fusion/accepted_hits.bam ] && \
tophat \
    -o tophat_fusion \
    -p 4 \
    --fusion-search \
    --fusion-min-dist 200 \
    --fusion-ignore-chromosomes chrM \
    --keep-fasta-order \
    --bowtie1 \
    --no-coverage-search \
    /n/groups/shared_databases/igenome/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome \
    ${sample}.unmapped.fastq

echo "# calling CIRCexplorer"
# update with new reference file (add CDR1as) - 20171026
[ ! -f ${sample}_circ.txt ] && \
CIRCexplorer.py \
  -f tophat_fusion/accepted_hits.bam \
  -g /n/groups/shared_databases/igenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \
  -r /home/xd20/neurogen/circRNA/MSBB/refFlat_plus_CDR1as.txt \
  -o $sample
  
# ## After STAR is done, copy the output *bam and *circ.txt to HPC and sort/index there
# cd ~/neurogen/ROSMAP/MSBB/rnaseq_runoutput
# rsync -azv xd20@transfer.orchestra.med.harvard.edu:~/projects/circRNA/MSBB/run_output/*/*circ.txt .
# rsync -azv xd20@transfer.orchestra.med.harvard.edu:~/projects/circRNA/MSBB/run_output/*/*.out.bam .
# # if not sorted/indexed
# for i in *.out.bam; do echo $i; bsub -q short -n 3 -J `basename $i .Aligned.out.bam` "samtools sort -m 2G -@ 3 -o $i.sorted $i; mv $i.sorted $i; samtools index $i;"; done