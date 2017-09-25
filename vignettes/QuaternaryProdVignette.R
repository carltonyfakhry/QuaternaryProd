## ---- message=FALSE------------------------------------------------------
library(QuaternaryProd)

# Compute the probability mass function
pmf <- QP_Pmf(q_p = 20, q_m = 20, q_z = 20, q_r = 0, n_p = 20, n_m = 20, n_z = 20)

# Plot the mass function
plot(names(pmf), pmf, col="blue", xlab = "scores", ylab = "probabilities")
lines(names(pmf), pmf, col = "blue")

## ---- message=FALSE------------------------------------------------------
# Get the p-value of score 5
pval <- QP_Pvalue(score = 5, q_p = 20, q_m = 20, q_z = 20, q_r = 0, 
                                                     n_p = 20, n_m = 20, n_z = 20)
pval

# Compue the p-value only if it is statistically significant otherwise
# return -1
pval <- QP_SigPvalue(score = 5, q_p = 20, q_m = 20, q_z = 20, q_r = 0, 
                                                     n_p = 20, n_m = 20, n_z = 20)
pval

## ---- message=FALSE------------------------------------------------------
library(QuaternaryProd)
library(fdrtool)

# Gene expression data
evidence1 <- system.file("extdata", "e2f3_sig.txt", package = "QuaternaryProd")
evidence1 <- read.table(evidence1, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
evidence2 <- system.file("extdata", "myc_sig.txt", package = "QuaternaryProd")
evidence2 <- read.table(evidence2, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
evidence3 <- system.file("extdata", "ras_sig.txt", package = "QuaternaryProd")
evidence3 <- read.table(evidence3, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Remove duplicated entrez ids in evidence and rename column names appropriately
names(evidence1) <- c("entrez", "pvalue", "fc")
evidence1 <- evidence1[!duplicated(evidence1$entrez),]

names(evidence2) <- c("entrez", "pvalue", "fc")
evidence2 <- evidence2[!duplicated(evidence2$entrez),]

names(evidence3) <- c("entrez", "pvalue", "fc")
evidence3 <- evidence3[!duplicated(evidence3$entrez),]

# Run Quaternary CRE for entire Knowledge base on new evidence
# which computes the statistic for each of the source proteins

CRE_results <- RunCRE_StringDB(evidence1, is.Logfc = TRUE)
# Get FDR corrected p-values
CRE_results$pvalue <- fdrtool(CRE_results$pvalue, "pvalue", FALSE,
                                                      FALSE, FALSE, "fndr")$q
head(CRE_results[order(CRE_results$pvalue), c("uid","name","pvalue")])

