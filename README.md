# QuaternaryProd
An R package for computing the Quaternary Product Scoring Statistic for signed and unsigned causal graphs.

## Installation
Before installing this package, make sure you have the latest version of *Rstudio*, *R* and the *devtools* package. You can install this R pacakge using the following:
```{R}
library(devtools)
install_github("carltonyfakhry/QuaternaryProd", build_vignettes = TRUE, local = FALSE)
```
## Usage
For an introduction to the Quaternary Product Scoring Statistic and for an example on how to compute it over the publicly available network *Stringdb*, please see 
the *Vignette* for this package using the following:
```{R}
browseVignettes("QuaternaryProd")
```
---
title: "QuaternaryProd"
author: "Carl Tony Fakhry, Ping Chen and Kourosh Zarringhalam"
date: "`r Sys.Date()`"
output: rmarkdown::pdf_document
        

vignette: >
  %\VignetteIndexEntry{<span style="color:red">QuaternaryProdVignette</span>}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

\(\indent\) \(\indent\) \(\indent\) \(\indent\) A signed causal graph is a directed graph where the edges are signed and the signs indicate the direction of effect of the source node on the target node (the signs are either + or -). **QuaternaryProd** is a package for computing the Quaternary Product Scoring Statistic for signed causal graphs. The Quaternary Product Scoring Statistic is a generalization of the Ternary Product statistic which allows for ambiguities to arise in a signed causal graph. Ambiguities arise when a source node can affect a target node in two different ways or if the direction of causality is unknown. We will first provide some background, and then we will apply the statistic to Stringdb which is a publicly available biological network.

## Introduction
\(\indent\) \(\indent\) \(\indent\) \(\indent\)  The Quaternary Product Scoring Statistic [1] is a goodness of fit test for examining how well the predictions of a signed and directed causal graph predict on newly realized experimental data. Given a source node \(s\) in a signed causal graph, let \(q_p\), \(q_m\) and \(q_r\) denote the number of target nodes which are increased, decreased and regulated by the source node respectively. Similarly, let \(q_z\) denote the set of target nodes in the causal network which do not share a relation with \(s\) i.e which are not affected by \(s\). Regulated relations occur when a source node regulates a target node without knowing the direction of causality or if an ambiguity in direction of causality occurs. An ambiguity can occur if a source node, according to a given network, shares both increase and decrease relations with the same target node. 
Next, Suppose we run some experiments on entities which are target nodes in the network. Let \(n_p\), \(n_m\) and \(n_z\) denote the set of values which are increased, decreased and remain unchanged in the experimental values respectively. For the source node \(s\), we can tabulate the predictions from the network vs. the experimental values:

\begin{table}[!htbp]
\centering
\begin{tabular}{l l l l l}
\hline
 & Observed \(+\)  & Observed \(-\) & Observed \(0\) &  Total \\ \hline
Predicted \(+\) & \(n_{pp}\) & \(n_{pm}\)  & \(n_{pz}\)  & \(q_p\) \\
Predicted \(-\) & \(n_{mp}\)  & \(n_{mm}\)  & \(n_{mz}\)  & \(q_m\) \\
Predicted \(r\) & \(n_{rp}\)  & \(n_{rm}\)  & \(n_{rz}\)  & \(q_r\) \\
Predicted \(0\) & \(n_{zp}\) & \(n_{zm}\) & \(n_{zz}\) & \(q_z\) \\
Total & \(n_p\)  & \(n_{m}\)   & \(n_z\) & \(T\) \\
\hline
\end{tabular}
\caption{Tabulation of predictions from network edges vs. observations from experimental results.}
\end{table}

\(n_{pp}\) denotes the number of target nodes which \(s\) is predicted to increase by the network and were indeed increased in  experimental values; \(n_{pm}\) the number of target nodes which \(s\) is predicted to increase and were decreased in experimental values; \(n_{pz}\) is the number of target nodes which \(s\) is predicted to increase and were unchanged in experimental values. Similar interpretation follows for all other entries of the table. 
The probability of a table follows the Quaternary Product distribution which is given by:
\begin{align}
P(\text{Table}) &= \frac{ \binom{q_p}{n_{pp},n_{pm},n_{pz}} \binom{q_m}{n_{mp},n_{mm},n_{mz}}\binom{q_z}{n_{zp},n_{zm},n_{zz}}\binom{q_r}{n_{rp},n_{rm},n_{rz}}}{\binom{T}{n_p,n_m,n_z}}.
\end{align}
Note, since the predictions by the network and the experimental values are fixed, then the table has 6 degrees of freedom \(n_{pp}\), \(n_{mm}\), \(n_{rp}\), \(n_{rm}\), \(n_{mp}\) and \(n_{pm}\). The score \(S\) to measure the goodness of fit is given by:
\begin{align}
S(\text{Table}) &= n_{pp} + n_{mm} + n_{rp} + n_{rm} - (n_{mp} + n_{pm}) 
\end{align}
which is the sum of the good predictions (i.e \(n_{pp}\), \(n_{mm}\), \(n_{rp}\) and \(n_{rm}\)) minus the bad predictions (i.e \(n_{mp}\) and \(n_{pm}\)). To compute the probability of a score, we sum the probabilites of all tables with score \(S\) as follows:
\begin{align}
P(S) &= \sum_{P(\text{Table}) = S } P(\text{Table}).
\end{align}

## Functionality
\(\indent\) \(\indent\) \(\indent\) \(\indent\)   **QuaternaryProd** provides different functions for computing the probability of a score, probability mass function, P-value of a score and the domain of the Quaternary Product Scoring Statistic. The probability mass function can be computed if given the margins of the table. 
  
```{R, message=FALSE }
library(QuaternaryProd)

# Compute the probability mass function
pmf <- QP_Pmf(q_p = 20, q_m = 20, q_z = 20, q_r = 0, n_p = 20, n_m = 20, n_z = 20)

# Plot the mass function
plot(names(pmf), pmf, col="blue", xlab = "scores", ylab = "probabilities")
lines(names(pmf), pmf, col = "blue")
```  
  
The package contains optimized functions for computing the P-value of a score. To compute the P-value of score we can use the following:

```{R, message=FALSE }
# Get the P-value of score 5
pval <- QP_Pvalue(score = 5, q_p = 20, q_m = 20, q_z = 20, q_r = 0, 
                                                     n_p = 20, n_m = 20, n_z = 20)
pval

# Compue the P-value only if it is statistically significant otherwise
# return -1
pval <- QP_SigPvalue(score = 5, q_p = 20, q_m = 20, q_z = 20, q_r = 0, 
                                                     n_p = 20, n_m = 20, n_z = 20)
pval
```  

If the user is only interested in obtaining statistically significant P-values, then *QP_SigPvalue* is optimized for this purpose. In either case, the user is advised to compute the P-value of a score using the previous two functions which will be faster than computing the entire probability mass function and then computing the P-value.
Finally, it is possible to also compute the probabilities of scores individually using *QP_Probability* and the support of the distribution using *QP_Support*. 
Since this package is written to the benefit of bioinformaticians, we will provide an example on how to apply this statistic to a publicly available network.
One bioinformatic application is to test how well protein-protein causal networks can predict the results in gene expression data. In the last section of this Vignette, we present an example of computing this statistic over the Stringdb network and given gene expression data.   
  
## Edges from STRINGdb

\(\indent\) \(\indent\) \(\indent\) \(\indent\)  You can use the **STRINGdb** package to interact with the [String database](www.string-db.org). The current release of **STRINGdb** only allows querying of neighbors in the network without the direction of action. For the purposes of the Quaternary Product Scoring Statistic, it is necessary to have the direction of action. It is possible to obtain the signed network by downloading it directly from [String-db network source](http://string-db.org/newstring_cgi/show_download_page.pl), selecting the species of interest, and downloading the protein actions network. Here, we present an example for working with the freely available Homo Sapien protein-protein interaction network from STRINGdb. The network is in tab seperated format which is straightforward to work with. It is possible to get a larger version of network with more relations from STRINGdb by signing a license agreement with the authors of STRINGdb.

### Parse File

\(\indent\) \(\indent\) \(\indent\) \(\indent\)    In this section, we show how to parse the Homo Sapien protein actions network and prepare it to be used with our package. First, we need to upload the network which is attached to **QuaternaryProd** for convenience.

```{r, message=FALSE, eval=FALSE}
library(readr)
library(org.Hs.eg.db)
library(dplyr)
library(stringr)

# Get the full file name containing the STRINGdb relations
ff <- system.file("extdata", "9606.protein.actions.v10.txt.gz", package="QuaternaryProd") 
all_rels <- read_tsv(gzfile(ff), col_names = T)

```

Next, we filter out the important columns and important relations. We remove all rows which do not have a relation *activation*, *inhibition* and *expression*. Moreover, we also consider reverse causality for any relation which has a *direction* value equal to 0.

```{r, message=FALSE, eval=FALSE}
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

# Rename activation as increase, inhibition as decrease, expression and binding
# as regulate
Rels$mode <- sub("activation", "increases", Rels$mode)
Rels$mode <- sub("inhibition", "decreases", Rels$mode)
Rels$mode <- sub("expression", "regulates", Rels$mode)
```

Third, we extract the protein entities from the network, and we map them to their respective genes. Note, the entities could have been possibly a drug or compound, but we are working with this protein interactions network for the purpose of providing a nontrivial example.

```{r, message=FALSE, eval=FALSE}
# Get all unique target protein ensemble ids in the causal network
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

# Add mRNAs to entities
uid <- paste("mRNA_", uid, sep = "")
mRNAs <- data_frame(uid, id, symbol, type="mRNA")
Ents <- bind_rows(Ents, mRNAs)
```

Finally, we filter unique relations in the network, and remove source proteins which do not have more than 10 children in the network.

```{r, message=FALSE, eval=FALSE}
# Get all unique relations
Rels <- Rels[Rels$srcuid %in% Ents$uid & Rels$trguid %in% Ents$uid,]
Rels <- unique(Rels)
Rels$trguid <- paste("mRNA_", Rels$trguid, sep="")

# Leave source proteins which contain at least 10 edges
sufficientRels <- group_by(Rels, srcuid) %>% summarise(count=n()) 
sufficientRels <- sufficientRels %>% filter(count > 10)
Rels <- Rels %>% filter(srcuid %in% sufficientRels$srcuid)
```

### Compute Pvalues Over the Network
\indent  Given new gene expression data, we can compute the scores and P-values for all source nodes in the network. *BioQCREtoNet* is a specialized function for this purpose.

```{r, message=FALSE, eval=FALSE}
# Gene expression data
evidence1 <- system.file("extdata", "e2f3_sig.txt", package = "QuaternaryProd") 
evidence1 <- read.table(evidence1, sep = "\t", header = T, stringsAsFactors = F)
evidence2 <- system.file("extdata", "myc_sig.txt", package = "QuaternaryProd") 
evidence2 <- read.table(evidence2, sep = "\t", header = T, stringsAsFactors = F)
evidence3 <- system.file("extdata", "ras_sig.txt", package = "QuaternaryProd") 
evidence3 <- read.table(evidence3, sep = "\t", header = T, stringsAsFactors = F)

# Remove duplicated entrez ids in evidence and rename column names appropriately
evidence1 <- evidence1[!duplicated(evidence1$entrez),]
names(evidence1) <- c("entrez", "pvalue", "fc")

# Run Quaternary CRE for entire Knowledge base on new evidence
# which computes the statistic for each of the source proteins
CRE_results <- BioQCREtoNet(Rels, evidence1, Ents)
```
*BioQCREtoNet* returns a data frame containing all the source nodes of the causal network, all of which had their respective score P-value computed. The source nodes are ordered in increasing order (Note: details on the columns of the data frame returned can be found in the help page for *BioQCREtoNet*).

## Citation

Please cite:

C. T. Fakhry, P. Choudhary, A. Gutteridge, B. Sidders, P. Chen, D. Ziemek, K. Zarringhalam. Identifying Transcriptional Regulators in Signed and Unsigned Causal Networks, 2015, In submission.

al. FAe (2013). STRING v9.1: protein-protein interaction networks, with increased coverage and integration. Nucleic Acids Research (Database issue), 41.

