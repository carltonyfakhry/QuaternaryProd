\name{QuaternaryProd-package}
\alias{QuaternaryProd-package}
%\alias{QuaternaryProd}
\docType{package}
\title{
Computes the Quaternary Product Scoring Statistic for Signed and Unsigned Causal Graphs
}
\description{
One of the most challenging problems in computational
            biology is inference of active regulatory cascades under specific
            molecular and environmental perturbations. As our understanding
            of regulatory mechanisms grows, larger networks of causal and
            non-causal biological interactions become available. However, only
            a fraction of regulatory interactions are tagged as either increasing
            or decreasing a target transcript. This leads to networks with a mix
            of signed and unsigned interactions. In this package, we provide a
            statistical method (Quaternary Product Scoring Statistic) that will infer likely upstream regulators given (1) a
            set of up- and down-regulated transcripts from a specific biological
            experiment and (2) a mixed network of regulatory interactions
            potentially relevant in the current biological context.
}
\details{
\tabular{ll}{
Package: \tab QuaternaryProd\cr
Type: \tab Package\cr
Version: \tab 0.99.2\cr
Date: \tab 2015-10-22\cr
License: \tab GPL (>= 2)\cr
}
%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~
}
\author{
Carl Tony Fakhry, Ping Chen and Kourosh Zarringhalam 

Maintainer: Carl Tony Fakhry <cfakhry@cs.umb.edu>
}
\references{
C. T. Fakhry, P. Choudhary, A. Gutteridge, B. Sidders, P. Chen, D. Ziemek, K. Zarringhalam. Identifying Transcriptional Regulators in Signed and Unsigned Causal Networks, 2015, Submitted.

Franceschini, A (2013). STRING v9.1: protein-protein interaction networks, with increased coverage and integration. In:'Nucleic Acids Res. 2013 Jan;41(Database issue):D808-15. doi: 10.1093/nar/gks1094. Epub 2012 Nov 29'.
}

%~~ Literature or other references for background information ~~
%}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%% ~~ simple examples of the most important functions ~~
%}
