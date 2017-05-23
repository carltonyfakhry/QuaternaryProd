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
library(readr)
library(org.Hs.eg.db)
library(dplyr)
library(stringr)
library(fdrtool)

# Get the full file name containing the STRINGdb relations
ff <- system.file("extdata", "9606.protein.actions.v10.txt.gz", package="QuaternaryProd") 
all_rels <- read_tsv(gzfile(ff), col_names = TRUE)


## ---- message=FALSE------------------------------------------------------
# Set new names for columns
names(all_rels) <- c("srcuid", "trguid", "mode", "action", "direction","score")
Rels <- all_rels[, c("srcuid", "trguid", "mode", "direction")]

# Get all rows with causal relations
Rels <- Rels[Rels$mode %in% c("activation", "inhibition","expression"),]

# Get causal relations where direction is not specified, and consider reversed 
# direction of causality as a valid causal relation
Bidirectional <- Rels[Rels$direction == 0 , c("trguid", "srcuid", "mode", "direction")]
names(Bidirectional) <- c("srcuid", "trguid", "mode", "direction")
Rels <- unique(bind_rows(Rels, Bidirectional))
Rels$direction <- NULL

# Rename activation as increases, inhibition as decreases, expression
# as regulates
Rels$mode <- sub("activation", "increases", Rels$mode)
Rels$mode <- sub("inhibition", "decreases", Rels$mode)
Rels$mode <- sub("expression", "regulates", Rels$mode)
Rels <- unique(Rels)

# Get a subset of the network: Skip this step if you want the p-values 
# of the scores corresponding to the source nodes computed over the 
# entire network. 
Rels <- Rels[sample(1:nrow(Rels), 40000, replace=FALSE),]

## ---- message=FALSE------------------------------------------------------
# Get all unique protein ensemble ids in the causal network
allEns <- unique(c(Rels$srcuid, Rels$trguid))

# Map ensemble protein ids to entrez gene ids 
map <- org.Hs.egENSEMBLPROT2EG
id <- unlist(mget(sub("9606.","",allEns), map, ifnotfound=NA))
id[is.na(id)] <- "-1"
uid <- paste("9606.", names(id), sep="")

# Function to map entrez ids to gene symbols
map <- org.Hs.egSYMBOL
symbol <- unlist(mget(id, map, ifnotfound=NA))
symbol[is.na(symbol)] <- "-1"

# Create data frame of STRINGdb protein Id, entrez id and gene symbol and type of entity
Ents <- data_frame(uid, id, symbol, type="protein")
Ents <- Ents[Ents$uid %in% allEns,]

# Remove ensemble ids in entities with duplicated entrez id
Ents <- Ents[!duplicated(Ents$id),]

# Add mRNAs to entities
uid <- paste("mRNA_", Ents$uid, sep = "")
mRNAs <- data_frame(uid=uid, id=Ents$id, symbol=Ents$symbol, type="mRNA")
Ents <- bind_rows(Ents, mRNAs)

## ---- message=FALSE------------------------------------------------------
# Get all unique relations
Rels$trguid <- paste("mRNA_", Rels$trguid, sep="")
Rels <- Rels[Rels$srcuid %in% Ents$uid & Rels$trguid %in% Ents$uid,]
Rels <- unique(Rels)

# Leave source proteins which contain at least 10 edges
sufficientRels <- group_by(Rels, srcuid) %>% summarise(count=n()) 
sufficientRels <- sufficientRels %>% filter(count > 10)
Rels <- Rels %>% filter(srcuid %in% sufficientRels$srcuid)

## ---- message=FALSE------------------------------------------------------
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

CRE_results <- BioQCREtoNet(Rels, evidence1, Ents, is.Logfc = TRUE)
# Get FDR corrected p-values
CRE_results$pvalue <- fdrtool(CRE_results$pvalue, "pvalue", FALSE, 
                                                      FALSE, FALSE, "fndr")$q
head(CRE_results[order(CRE_results$pvalue), c("uid","name","pvalue")])

