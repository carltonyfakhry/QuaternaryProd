# Get the value of regulation of a trguid from the gene expression
# data (i.e in evidence)
getGeneVals <- function(trguids, gene_expression_data){
  trguids <- as.vector(trguids) 
  val <- rep(0, length(trguids))
  ind <- match(trguids, gene_expression_data[,"uid"])
  ind2 <- na.omit(ind)
  if(length(ind2) > 0){
    val[which(!is.na(ind))] <- gene_expression_data[,"val"][ind2]
  }
  return(val)
}



#' This function runs a causal relation engine by computing the Quaternary Dot
#' Product Scoring Statistic, Ternary Dot Product Scoring Statistic or the Enrichment test over the Homo
#' Sapien STRINGdb causal network.
#' 
#' @description This function runs a causal relation engine by computing the Quaternary Dot
#'              Product Scoring Statistic, Ternary Dot Product Scoring Statistic or the Enrichment test over the Homo
#'              Sapien STRINGdb causal network.
#' 
#' @usage RunCRE_HSAStringDB(gene_expression_data, method = "Quaternary", 
#'                     fc.thresh = log2(1.3), pval.thresh = 0.05) 
#' 
#' @param gene_expression_data A data frame for gene expression data. The \code{gene_expression_data} data frame must have three columns \code{entrez}, 
#'        \code{fc} and \code{pvalue}. \code{entrez} denotes the entrez id of a given gene, \code{fc} denotes
#'        the fold change of a gene, and \code{pvalue} denotes the p-value. The \code{entrez} column must be of type
#'        integer or character, and the \code{fc} and \code{pvalue} columns must be numeric values.
#' @param method Choose one of \code{Quaternary}, \code{Ternary} or \code{Enrichment}. Default is \code{Quaternary}.
#' @param fc.thresh Threshold for fold change in \code{gene_expression_data} data frame. Any row in gene_expression_data with abosolute value of \code{fc}
#'        smaller than \code{fc.thresh} will be ignored. Default value is \code{fc.thresh = log2(1.3)}. 
#' @param pval.thresh Threshold for p-values in \code{gene_expression_data} data frame. All rows in \code{gene_expression_data} with p-values 
#'        greater than \code{pval.thresh} will be ingnored. Default value is \code{pval.thresh = 0.05}. 
#'        
#' @return This function returns a data frame containing parameters concerning the method used. The p-values of each
#'         of the source nodes is also computed, and the data frame
#'         is in increasing order of p-values of the goodness of fit score for the given source nodes. The column
#'         names of the data frame are:
#'         
#' \itemize{        
#' \item  \code{uid} The source node in the STRINGdb network.
#' \item \code{name} Symbol of the source node. 
#' \item \code{regulation} Direction of change of source node.
#' \item \code{correct.pred} Number of correct predictions in \code{gene_expression_data} when compared to predictions made
#'                     by the network.
#' \item \code{incorrect.pred} Number of incorrect predictions in \code{gene_expression_data} when compared to predictions made
#'                     by the network.
#' \item \code{score} The number of correct predictions minus the number of incorrect predictions. 
#' \item \code{total.reachable} Total Number of children of the given source node.
#' \item \code{significant.reachable} Number of children of the given source node that are also present 
#'                                    in \code{gene_expression_data}.
#' \item \code{total.ambiguous} Total number of children of the given source node which are regulated by the given source node without
#'                              knowing the direction of regulation.
#' \item \code{significant.ambiguous} Total number of children of the given source node which are regulated by the given source node without
#'                              knowing the direction of regulation and are also present in \code{gene_expression_data}.  
#' \item \code{unknown} Number of target nodes in the STRINGdb causal network which do not interact with the given source node.
#' \item \code{pvalue} P-value of the score computed according to the selected method.
#' }
#' 
#' @author Carl Tony Fakhry, Ping Chen and Kourosh Zarringhalam
#' 
#' @references Carl Tony Fakhry, Parul Choudhary, Alex Gutteridge, Ben Sidders, Ping Chen, Daniel Ziemek, and
#'             Kourosh Zarringhalam. Interpreting transcriptional changes using causal graphs: new methods and
#'             their practical utility on public networks. BMC Bioinformatics, 17:318, 2016. ISSN 1471-2105.
#'             doi: 10.1186/s12859-016-1181-8.
#'            
#'             Franceschini, A (2013). STRING v9.1: protein-protein interaction networks, with increased coverage 
#'             and integration. In:'Nucleic Acids Res. 2013 Jan;41(Database issue):D808-15. doi: 10.1093/nar/gks1094. 
#'             Epub 2012 Nov 29'.
#'             
#' @examples 
#'
#' # Get gene expression data
#' gene_expression_data1 <- system.file("extdata", "e2f3_sig.txt", package = "QuaternaryProd")
#' gene_expression_data1 <- read.table(gene_expression_data1, sep = "\t", header = TRUE
#'                                                   , stringsAsFactors = FALSE)
#'
#' # Remove duplicated entrez ids in gene_expression_data1 and rename column names appropriately
#' gene_expression_data1 <- gene_expression_data1[!duplicated(gene_expression_data1$entrez),]
#' names(gene_expression_data1) <- c("entrez", "pvalue", "fc")
#'           
#' # Compute the statistic for each source node in the String network
#' results <- RunCRE_HSAStringDB(gene_expression_data1, method = "Quaternary",
#'                               fc.thresh = log2(1.3), pval.thresh = 0.05)
#'
#' @export

RunCRE_HSAStringDB <- function(gene_expression_data, method = "Quaternary", fc.thresh = log2(1.3), pval.thresh = 0.05){
  
  f1 <- system.file("extdata", "newStringRels.dat", package="QuaternaryProd")
  relations <- read.table(f1, header = T, stringsAsFactors = F)
  f2 <- system.file("extdata", "newStringEnts.dat", package="QuaternaryProd")
  entities <- read.table(f2, header = T, stringsAsFactors = F)
  # f1 <- system.file("extdata", "StringRels.dat", package="QuaternaryProd") 
  # relations <- read.table(f1, header = T, stringsAsFactors = F)
  # f2 <- system.file("extdata", "StringEnts.dat", package="QuaternaryProd") 
  # entities <- read.table(f2, header = T, stringsAsFactors = F)
  
  # Check method
  if(!(method %in% c("Quaternary", "Ternary", "Enrichment"))){
    stop("Method must be one of Quaternary, Ternary or Enrichment!")
  }
  
  # Check fc.thresh
  if(length(fc.thresh) != 1 || !is.numeric(fc.thresh)){
    stop("fc.thresh should be a numeric value!")
  }
  
  # Check pval.thresh
  if(length(pval.thresh) != 1 || !is.numeric(pval.thresh) | pval.thresh > 1 | pval.thresh < 0){
    stop("pval.thresh should be a numeric value >= 0 and <= 1!")
  }
  
  if (!is.data.frame(gene_expression_data)){ 
    stop("Please enter a data frame for gene_expression_data!")
  }
  gene_expression_data <- data.frame(gene_expression_data)
  
  
  if (nrow(gene_expression_data) == 0){
    stop("Please enter a non-empty data frame for gene_expression_data!")
  }
  
  # Do more validation for gene_expression_data data frame
  if(ncol(gene_expression_data) == 3){
 
    # Check for NAs in gene_expression_data
    if (length(which(is.na(gene_expression_data))) > 0){
      stop("gene_expression_data data frame contains rows with NAs. Please only provide complete cases!")
    }
        
    # Check if there are duplicated rows in gene_expression_data
    if(nrow(gene_expression_data[duplicated(gene_expression_data),]) > 0){
      stop("Duplicate rows in gene_expression_data found. Please remove duplicate rows!")
    }
    
    # Check for correct column names
    pval.ind <- which(colnames(gene_expression_data) %in% c("pvalue"))
    fc.ind <- which(colnames(gene_expression_data) %in% c("fc") )
    id.ind <- which(colnames(gene_expression_data) %in% c("entrez"))
    
    # Stop if column names are not correct
    if(length(id.ind) == 0 | length(fc.ind) == 0 | length(pval.ind) == 0){
      stop("Please label gene_expression_data column names as entrez, fc, pvalue!")
    }
    
    # Rearrange data
    gene_expression_data <- gene_expression_data[,c(id.ind, fc.ind, pval.ind)]
    names(gene_expression_data) <- c("entrez", "fc", "pvalue")
    
    # Check data types of entities column
    n.e0 <- nrow(gene_expression_data)
    if(is.character(gene_expression_data[,"entrez"])){
      suppressWarnings(gene_expression_data[,"entrez"] <- as.double(gene_expression_data[,"entrez"]))
      gene_expression_data <- gene_expression_data[complete.cases(gene_expression_data),]
      n.e0 <- n.e0 - nrow(gene_expression_data)
      if(any(gene_expression_data[,"entrez"] != floor(gene_expression_data[,"entrez"]))){
        stop("None of the entrez ids in gene_expression_data can be a double or a character representation of a double!")
      }
      gene_expression_data[,"entrez"] <- as.integer(gene_expression_data[,"entrez"])
    }else if(is.integer(gene_expression_data[,"entrez"])){
    }else{
      stop("In gene_expression_data, column entrez must be of type character or integer!")
    }
    
    if(!is.numeric(gene_expression_data[,"fc"])){
      stop("In gene_expression_data, column fc must be numeric!")
    }
    
    if(!is.numeric(gene_expression_data[,"pvalue"])){
      stop("In gene_expression_data, column pvalue must be numeric!")
    }

    # Make sure p-values are positive and less than 1
    if(!all(gene_expression_data[,"pvalue"] >= 0)){
      stop("All P-values in gene_expression_data data frame must be greater or equal to zero!")
    }
    
    if(!all(gene_expression_data[,"pvalue"] <= 1)){
      stop("All P-values in gene_expression_data data frame must be smaller than one!")
    }
    
    # Remove duplicate entrez in gene_expression_data
    if(any(duplicated(gene_expression_data[,"entrez"]))){
      stop("Remove duplicated entrez ids in gene_expression_data!")
    }
    
    # Filter by p-value and fold-change
    gene_expression_data <- gene_expression_data[abs(gene_expression_data[,"fc"]) >= fc.thresh & gene_expression_data[,"pvalue"] <= pval.thresh,]
    gene_expression_data <- gene_expression_data[,c("entrez","fc")]
    
    # Check if any gene_expression_data left
    if(nrow(gene_expression_data) == 0){
      stop("No rows left in gene_expression_data after rows are filtered according to fc.thresh and pval.thresh!")
    }
    
    # Where fc is greater than 0 replace with increase i.e 1 and where fc is negative
    # replace with -1 i.e decreases
    gene_expression_data[,2] <- ifelse(gene_expression_data[,2] > 0, 1, -1)
    names(gene_expression_data) <- c("entrez", "val")
    
  } else{
    
    stop("Please provide a valid gene_expression_data data frame containing columns: entrez, fc, pvalue")
  
  }
  
  # Make sure all gene_expression_data entrez are in entities ids
  n.e1 <- nrow(gene_expression_data)
  gene_expression_data <- gene_expression_data[gene_expression_data[,"entrez"] %in% entities[,"id"],]
  n.e2 <- nrow(gene_expression_data)
  if (n.e2 == 0){
    stop("All entrez ids in gene_expression_data are not present in entities data frame!")
  }
  cat(paste((n.e1-n.e2+n.e0), "rows from gene_expression_data removed due", "\n", "to entrez ids being unrepsented in StringDB entities!", sep = " "))

  # Make sure no entrez id in entities maps to two entities
  gene_expression_data.tmp <- merge(gene_expression_data, entities, by.x = "entrez" , by.y = "id")
  gene_expression_data <- data.frame(uid = gene_expression_data.tmp[,"uid"], val = gene_expression_data.tmp[,"val"], stringsAsFactors = F)
    
  # u.hyps <- list.load(system.file("extdata", "u.hyps.yaml", package="QuaternaryProd"))
  # child.uid <- list.load(system.file("extdata", "child.uid.yaml", package="QuaternaryProd"))
  # child.sgn <- list.load(system.file("extdata", "child.sgn.yaml", package="QuaternaryProd"))
  u.hyps <- list.load(system.file("extdata", "new.u.hyps.yaml", package="QuaternaryProd"))
  child.uid <- list.load(system.file("extdata", "new.child.uid.yaml", package="QuaternaryProd"))
  child.sgn <- list.load(system.file("extdata", "new.child.sgn.yaml", package="QuaternaryProd"))
  
  # Get the value of regulation of the children from the gene expression data
  child.val <- lapply(child.uid, function(x, gene_expression_data) 
                      getGeneVals(x, gene_expression_data), gene_expression_data = gene_expression_data)

  # Get gene values from gene_expression_data for all non children
  unique.children <- unique(relations$trguid)
  non.child.val <- lapply(child.uid, function(x, gene_expression_data, unique.children) 
                          getGeneVals(unique.children[which(!(unique.children %in% x))], gene_expression_data), 
                          gene_expression_data = gene_expression_data, unique.children = unique.children)
  
  results <- data.frame(matrix(0, nrow  = 2 * length(u.hyps), ncol = 12), stringsAsFactors = F)
  colnames(results) <- c('uid', 'name', 'regulation', 'correct.pred', 'incorrect.pred', 'score',
                        'total.reachable', 'significant.reachable', 'total.ambiguous', 'significant.ambiguous',
                        'unknown', 'pvalue')
  
  
  for(p.s in 1:length(u.hyps)){
    results[(2*(p.s-1)+1), 1] <- entities$ensembleid[which(entities$uid == u.hyps[p.s])]
    results[(2*p.s), 1]       <- entities$ensembleid[which(entities$uid == u.hyps[p.s])]
    results[(2*(p.s-1)+1), 2] <- entities$symbol[which(entities$uid == u.hyps[p.s])]
    results[(2*p.s), 2]       <- entities$symbol[which(entities$uid == u.hyps[p.s])]
    # results[(2*(p.s-1)+1), 1] <- u.hyps[p.s]
    # results[(2*p.s), 1]       <- u.hyps[p.s]
    # results[(2*(p.s-1)+1), 2] <- entities$name[which(entities$uid == u.hyps[p.s])]
    # results[(2*p.s), 2]       <- entities$name[which(entities$uid == u.hyps[p.s])]
    results[(2*(p.s-1)+1), 3] <- 'up'
    results[(2*p.s), 3]       <- 'down'
    
    npp = length(which(child.sgn[[p.s]] == 1 &  child.val[[p.s]] == 1))
    npm = length(which(child.sgn[[p.s]] == 1 &  child.val[[p.s]] == -1))
    npz = length(which(child.sgn[[p.s]] == 1 &  child.val[[p.s]] == 0))
    
    nmp = length(which(child.sgn[[p.s]] == -1 &  child.val[[p.s]] == 1))
    nmm = length(which(child.sgn[[p.s]] == -1 &  child.val[[p.s]] == -1))
    nmz = length(which(child.sgn[[p.s]] == -1 &  child.val[[p.s]] == 0))
    
    nrp = length(which(child.sgn[[p.s]] == 0 &  child.val[[p.s]] == 1))
    nrm = length(which(child.sgn[[p.s]] == 0 &  child.val[[p.s]] == -1))
    nrz = length(which(child.sgn[[p.s]] == 0 &  child.val[[p.s]] == 0))
    
    nzp = length(which(non.child.val[[p.s]] == 1))
    nzm = length(which(non.child.val[[p.s]] == -1))
    nzz = length(which(non.child.val[[p.s]] == 0))
    
    if (method == "Quaternary"){
      pval <- runCRE(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, method = 'Quaternary')
    } else if (method == "Ternary"){
      pval <- runCRE(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, method = 'Ternary')
    } else if (method == "Enrichment"){
      pval <- runCRE(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, method = 'Enrichment')
    } else {
      stop("Select one of methods: Quaternary, Ternary or Enrichment")
    }
    
    qPlus  <- npp + npm + npz
    qMinus <- nmp + nmm + nmz
    qR     <- nrp + nrm + nrz
    qZero  <- nzp + nzm + nzz
    
    results[(2*(p.s-1)+1), 4]  <- npp + nmm
    results[(2*(p.s-1)+1), 5]  <- npm + nmp
    results[(2*(p.s-1)+1), 6]  <- npp + nmm - (npm + nmp)
    results[(2*(p.s-1)+1), 7]  <- qPlus + qMinus + qR
    results[(2*(p.s-1)+1), 8]  <- npp + npm + nmp + nmm + nrp + nrm
    results[(2*(p.s-1)+1), 9]  <- qR
    results[(2*(p.s-1)+1), 10] <- nrp + nrm
    results[(2*(p.s-1)+1), 11] <- qZero
    results[(2*(p.s-1)+1), 12] <- pval$pval.up
    
    results[(2*p.s), 4]  <- nmp + npm
    results[(2*p.s), 5]  <- npp + nmm
    results[(2*p.s), 6]  <- nmp + npm - (npp + nmm)
    results[(2*p.s), 7]  <- qPlus + qMinus + qR
    results[(2*p.s), 8]  <- npp + npm + nmp + nmm + nrp + nrm
    results[(2*p.s), 9]  <- qR
    results[(2*p.s), 10] <- nrp + nrm
    results[(2*p.s), 11] <- qZero
    results[(2*p.s), 12] <- pval$pval.down
  }
  
  
  results <- results[order(results$pvalue), ]
  rownames(results) <- 1:nrow(results)
  
  return(results)
  
}  



# This function runs the CRE based on the version specified
runCRE = function(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, method){
  
  if(method == 'Quaternary'){
    
    qR     <- nrp + nrm + nrz
    qZero  <- nzp + nzm + nzz
    nPlus  <- npp + nmp + nrp + nzp
    nMinus <- npm + nmm + nrm + nzm
    nZero  <- npz + nmz + nrz + nzz
    
    ## Assume up-regulated
    qPlus  <- npp + npm + npz
    qMinus <- nmp + nmm + nmz
    score  <- npp + nmm + nrp + nrm - (npm + nmp)
    pval.up <- QP_Pvalue(score = score, q_p = qPlus, q_m = qMinus, q_z = qZero,
                          q_r = qR, n_p = nPlus, n_m = nMinus, n_z = nZero)
    
    ## Assume down-regulated
    qPlus  <- nmp + nmm + nmz
    qMinus <- npp + npm + npz
    score  <- nmp + npm + nrp + nrm - (npp + nmm)
    pval.down <- QP_Pvalue(score = score, q_p = qPlus, q_m = qMinus, q_z = qZero,
                            q_r = qR, n_p = nPlus, n_m = nMinus, n_z = nZero)
    
  } else if(method == 'Ternary'){
    qR     <- 0
    qZero  <- nzp + nzm + nzz
    nPlus  <- npp + nmp + nzp
    nMinus <- npm + nmm + nzm
    nZero  <- npz + nmz + nzz
    
    ## Assume up-regulated
    qPlus  <- npp + npm + npz
    qMinus <- nmp + nmm + nmz
    score  <- npp + nmm - (npm + nmp)
    pval.up <- QP_Pvalue(score = score, q_p = qPlus, q_m = qMinus, q_z = qZero,
                          q_r = qR, n_p = nPlus, n_m = nMinus, n_z = nZero)
    
    ## Assume down-regulated
    qPlus  <- nmp + nmm + nmz
    qMinus <- npp + npm + npz
    score  <- nmp + npm - (npp + nmm)
    pval.down <- QP_Pvalue(score = score, q_p = qPlus, q_m = qMinus, q_z = qZero,
                            q_r = qR, n_p = nPlus, n_m = nMinus, n_z = nZero)
    
    
  } else if(method == 'Enrichment'){
    
    nrp    <- npp + nmp + nrp
    nrm    <- npm + nmm + nrm
    nrz    <- npz + nmz + nrz
    
    qPlus  <- 0
    qMinus <- 0
    qR     <- nrp + nrm + nrz
    qZero  <- nzp + nzm + nzz
    
    nPlus  <- nrp + nzp
    nMinus <- nrm + nzm
    nZero  <- nrz + nzz
    
    score  <- nrp + nrm
    
    pval.up <- QP_Pvalue(score = score, q_p = qPlus, q_m = qMinus, q_z = qZero,
                          q_r = qR, n_p = nPlus, n_m = nMinus, n_z = nZero)
    pval.down <- pval.up
    
  }
  
  pval <- list(pval.up = pval.up, pval.down = pval.down)
  
}

