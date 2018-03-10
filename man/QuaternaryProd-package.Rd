\name{QuaternaryProd-package}
\alias{QuaternaryProd-package}
%\alias{QuaternaryProd}
\docType{package}
\title{
Computes the Quaternary Dot Product Scoring Statistic for Signed and Unsigned Causal Graphs
}
\description{
	     QuaternaryProd is an R package that performs causal reasoning on biological 
             networks, including publicly available networks such as STRINGdb. QuaternaryProd 
             is an open-source alternative to commercial products such as Inginuity Pathway 
             Analysis. For a given a set of differentially expressed genes, QuaternaryProd 
             computes the significance of upstream regulators in the network by performing causal
             reasoning using the Quaternary Dot Product Scoring Statistic (Quaternary Statistic),
             Ternary Dot product Scoring Statistic (Ternary Statistic) and Fisher's exact test
	     (Enrichment test). The Quaternary Statistic handles signed, unsigned and ambiguous
             edges in the network. Ambiguity arises when the direction of causality is unknown,
             or when the source node (e.g., a protein) has edges with conflicting signs for the 
	     same target gene. On the other hand, the Ternary Statistic provides causal reasoning
             using the signed and unambiguous edges only. The Vignette provides more details on the
             Quaternary Statistic and illustrates an example of how to perform causal reasoning
             using STRINGdb. 
}
\details{
\tabular{ll}{
Package: \tab QuaternaryProd\cr
Type: \tab Package\cr
Version: \tab 1.7.4\cr
Date: \tab 2015-10-22\cr
License: \tab GPL (>= 2)\cr
}

}
\author{
Carl Tony Fakhry, Ping Chen and Kourosh Zarringhalam 

Maintainer: Carl Tony Fakhry <cfakhry@cs.umb.edu>
}
\references{
Carl Tony Fakhry, Parul Choudhary, Alex Gutteridge, Ben Sidders, Ping Chen, Daniel Ziemek, and Kourosh Zarringhalam. Interpreting transcriptional changes using causal graphs: new methods and their practical utility on public networks. BMC Bioinformatics, 17:318, 2016. ISSN 1471-2105. doi: 10.1186/s12859-016-1181-8.

Franceschini, A (2013). STRING v9.1: protein-protein interaction networks, with increased coverage and integration. In:'Nucleic Acids Res. 2013 Jan;41(Database issue):D808-15. doi: 10.1093/nar/gks1094. Epub 2012 Nov 29'.
}
