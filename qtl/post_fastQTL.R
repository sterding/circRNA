# Rscript to check fastQTL output files
# Ref: https://github.com/francois-a/fastqtl/blob/master/R/calculateSignificanceFastQTL.R
# Usage: Rscript post_fastQTL.R eQTL
library(qvalue) 

args<-commandArgs(TRUE)
prefix=args[1] 
# prefix='eQTLcircRNA'

q_cutoff=0.05

permutation_file=paste0(prefix,".permutations.txt.gz")
D = read.table(permutation_file, head=F, stringsAsFactors=F)
colnames(D) = c("gene_id", "num_var", "beta_shape1", "beta_shape2", "dummy", "variant_id", "dist", "pval_nominal", "slope","pval_perm", "pval_beta")

# remove genes w/o variation
nanrows <- is.na(D[, 'pval_beta'])
D <- D[!nanrows, ]
cat("  * Number of genes tested: ", nrow(D), " (excluding ", sum(nanrows), " genes w/o variation)\n", sep="")
cat("  * Correlation between Beta-approximated and empirical p-values: ", round(cor(D[, 'pval_perm'], D[, 'pval_beta']), 4), "\n", sep="")

pdf(paste0(prefix,".QC.correlation.pdf"), width=5, height = 5)
plot(D$pval_perm, D$pval_beta, xlab="Empirical p-values", ylab="Beta approximated p-values", main=paste0("Number of genes tested: ", nrow(D), "\n(excluding ", sum(nanrows), " genes w/o variants)"))
abline(0, 1, col="red")
legend("topleft", paste("R =",round(cor(D[, 'pval_perm'], D[, 'pval_beta']), 4)), bty='n')
dev.off()

# only for those eGene with nominal pvalue <= 0.05
# D=D[D$pval_beta<=0.05,]

# multiple correction
D$bonferroni = p.adjust(D$pval_beta, method="bonferroni")
D$fdr = p.adjust(D$pval_beta, method="fdr")
D$qval = qvalue(D$pval_beta)$qvalues
#arrange(D, pval_nominal, qval) %>% head()

cat("  * Proportion of significant phenotypes (1-pi0): " , round((1 - qvalue(D$pval_beta)$pi0), 2), "\n", sep="")
cat("  * eGenes @ FDR ", q_cutoff, ":   ", sum(D$qval<q_cutoff), "\n", sep="")

# determine global min(p) significance threshold and calculate nominal p-value threshold for each gene
ub <- sort(D[D$qval >= q_cutoff, 'pval_beta'])[1]  # smallest p-value above FDR
lb <- -sort(-D[D$qval <= q_cutoff, 'pval_beta'])[1]  # largest p-value below FDR
pthreshold <- mean(c(ub, lb), na.rm = T)
cat("  * min p-value threshold @ FDR ", q_cutoff, ": ", pthreshold, "\n", sep="")
D$pval_nominal_threshold <- signif(qbeta(pthreshold, D$beta_shape1, D$beta_shape2, ncp=0, lower.tail=TRUE, log.p=FALSE), 6)

write.table(D[D$qval<=q_cutoff,], gzfile(paste0(prefix,".egenes.txt.gz")), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
