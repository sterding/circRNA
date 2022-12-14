###########################################
## bash script for running DE analysis using count data and covariance table
###########################################
#!/bin/sh

cd ~/projects/circRNA/results/
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C CONDITION:PD:HC -s
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C CONDITION:ILB:HC -s
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C CONDITION:PD:ILB -s
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O Mmi -C CONDITION2:PD:HC -s
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C PDpathologygroup:pd:hc -s
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C PD.pathology.group:early:no -s
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C PD.pathology.group:late:early -s
#Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O mi -C PD.pathology.group:late:no -s
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.AD.TCPY.pathology.covariates.xls -o DE_TCPY -O Mmi -C CONDITION:AD:HC -s
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.AD.TCPY.pathology.covariates.xls -o DE_TCPY -O i -C Braak_Braak_stage -s
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA -O i -C MUSS -s

## PD: collapse to gene
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE2gene_SNDA -O Mmi -C CONDITION:PD:HC -k -s
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE2gene_SNDA -O Mmi -C CONDITION:ILB:HC -k -s
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE2gene_SNDA -O Mmi -C CONDITION2:PD:HC -k -s
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE2gene_SNDA -O i -C MUSS -k -s
## AD: collapse to gene
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.AD.TCPY.pathology.covariates.xls -o DE2gene_TCPY -O Mmi -C CONDITION:AD:HC -k -s
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.AD.TCPY.pathology.covariates.xls -o DE2gene_TCPY -O i -C Braak_Braak_stage -k -s

# for i in DEresult.DE2gene_SNDA.*xls; do awk '$6<=0.05' $i | cut -f1-6; done | cut -f1 | sort -u > DEgenes.list
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE2gene_SNDA -O Mmi -C CONDITION2:PD:HC -k -l DEgenes.list -s
# awk '$2<=0.1 && $4<=0.05 && $6<=0.05{print $1}' summaryTable.VMB.meta.p0.05.xls > summaryTable.VMB.meta.p0.05.list
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_BC197.filtered.enriched.rawcount.rds -c Table.PD.SNDA.pathology.covariates.xls -o DE2gene_SNDA -O Mmi -C CONDITION:ILB:HC -k -l summaryTable.VMB.meta.p0.05.list -s

# CSF
Rscript ~/projects/circRNA/src/DE/_DE2gene.R -i ../data/Merge_circexplorer_CSF87.rawcount.rds -c Table.PD.CSF.pathology.covariates.xls -o DE2gene_CSF -O Mmi -C CONDITION2:PD:HC -a ~/projects/circRNA/data/Merge_circexplorer_CSF87.annotation.bed14.rds -k -s

## Bennett's dataset
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_Bennett_VMB.filtered.rawcount.rds -c Table.Bennett.PD.VMP.pathology.covariates.xls -o DE2gene_VMB -O Mmi -C CONDITION:PD:HC -a ~/projects/circRNA/data/Merge_circexplorer_Bennett_VMB.filtered.annotation.bed14.rds -k
Rscript ~/projects/circRNA/src/DE/_DE.R -i ../data/Merge_circexplorer_Bennett_FCX.filtered.rawcount.rds -c Table.Bennett.AD.FCX.pathology.covariates.xls -o DE2gene_FCX -O Mmi -C CONDITION:AD:HC -a ~/projects/circRNA/data/Merge_circexplorer_Bennett_FCX.filtered.annotation.bed14.rds -k

## pathway analysis
cd ~/projects/circRNA/results/DE2gene_SNDA; Rscript ~/projects/circRNA/src/DE/_pathway.R -i DEresult.DE2gene_SNDA.CONDITION2_PD_vs_HC.xls -F loose
cd ~/projects/circRNA/results/DE2gene_TCPY; Rscript ~/projects/circRNA/src/DE/_pathway.R -i DEresult.DE2gene_TCPY.CONDITION_AD_vs_HC.xls -F loose

## WGCNA analysis
cd ~/projects/circRNA/results/DE_SNDA; Rscript ~/projects/circRNA/src/DE/_WGCNA.R -i ~/projects/circRNA/data/Merge_circexplorer_BC197.filtered.enriched.normRPM.rds -c 

## AD vs. PD DE comparison
# see comparison.R


##======================
## gene
##======================
cd ~/neurogen/rnaseq_PD/results/; Rscript ~/projects/circRNA/src/DE/_DE.R -i merged/genes.htseqcount.cufflinks.allSamples.BCv2.uniq.xls -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA.gene -O Mmi -C CONDITION2:PD:HC -t gene
cd ~/neurogen/rnaseq_PD/results/; Rscript ~/projects/circRNA/src/DE/_DE.R -i merged/genes.htseqcount.cufflinks.allSamples.BCv2.uniq.xls -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA.gene -O Mmi -C CONDITION:ILB:HC -t gene
cd ~/neurogen/rnaseq_PD/results/; Rscript ~/projects/circRNA/src/DE/_DE.R -i merged/genes.htseqcount.cufflinks.allSamples.BCv2.uniq.xls -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA.gene -O Mmi -C MUSS -t gene
cd ~/neurogen/rnaseq_PD/results/; Rscript ~/projects/circRNA/src/DE/_DE.R -i merged/genes.htseqcount.cufflinks.allSamples.BCv2.uniq.xls -c Table.AD.TCPY.pathology.covariates.xls -o DE_TCPY.gene -O Mmi -C CONDITION:AD:HC -t gene
cd ~/neurogen/rnaseq_PD/results/; Rscript ~/projects/circRNA/src/DE/_DE.R -i merged/genes.htseqcount.cufflinks.allSamples.BCv2.uniq.xls -c Table.AD.TCPY.pathology.covariates.xls -o DE_TCPY.gene -O Mmi -C Braak_Braak_stage -t gene
cd ~/neurogen/rnaseq_CSF/results; Rscript ~/projects/circRNA/src/DE/_DE.R -i merged/genes.htseqcount.cufflinks.CSF.uniq.xls -c Table.PD.CSF.pathology.covariates.xls -o DE_CSF.gene -O Mmi -C CONDITION:PD:HC -t gene

# make plot for example genes ERIC1, DNAJC6, DYM
cd ~/neurogen/rnaseq_PD/results/; 
grep -Ew "ERC1|DNAJC6|DYM" DE_SNDA.gene/DEresult.DE_SNDA.gene.CONDITION2_PD_vs_HC.f1-17.xls | sort -k1,1 | cut -f1 | sort -u > ERC1.DNAJC6.DYM.txt
Rscript ~/projects/circRNA/src/DE/_DE.R -i merged/genes.htseqcount.cufflinks.allSamples.BCv2.uniq.xls -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA.gene -O Mmi -C CONDITION2:PD:HC -t gene -l ERC1.DNAJC6.DYM.txt
Rscript ~/projects/circRNA/src/DE/_DE.R -i merged/genes.htseqcount.cufflinks.allSamples.BCv2.uniq.xls -c Table.PD.SNDA.pathology.covariates.xls -o DE_SNDA.gene -O Mmi -C CONDITION:ILB:HC -t gene -l ERC1.DNAJC6.DYM.txt


## pathway analysis
cd ~/neurogen/rnaseq_PD/results/DE_SNDA.gene; Rscript ~/projects/circRNA/src/DE/_pathway.R -i DEresult.DE_SNDA.gene.CONDITION2_PD_vs_HC.xls -F medium
cd ~/neurogen/rnaseq_PD/results/DE_TCPY.gene; Rscript ~/projects/circRNA/src/DE/_pathway.R -i DEresult.DE_TCPY.gene.CONDITION_AD_vs_HC.xls -F medium

## co-expression heatmap
cd ~/neurogen/rnaseq_PD/results/DE_SNDA.gene; Rscript ~/pipeline/modules/_expression2heatmap.R -i genes.htseqcount.cufflinks.allSamples.BCv2.uniq.xls.filtered.vsd_adjusted_log2.xls -l DEresult.DE_SNDA.gene.MUSS.padj0.05.xls