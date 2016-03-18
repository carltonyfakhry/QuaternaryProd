#' Computes the support for the scores.
#' 
#' @description This function computes the support of the Quaternary Dot Product Scoring distribution
#' for signed causal graphs. This includes all scores which have probabilities strictly greater than 0.
#' 
#' @usage QP_Support(q_p, q_m, q_z, q_r, n_p, n_m, n_z)
#' 
#' @param q_p Expected number of positive predictions.
#' @param q_m Expected number of negative predictions.
#' @param q_z Expected number of nil predictions.
#' @param q_r Expected number of regulated predictions. 
#'
#' @param n_p Number of positive predictions from experiments.
#' @param n_m Number of negative predictions from experiments.
#' @param n_z Number of nil predictions from experiments. 
#'        
#' @return Integer vector of support.
#' 
#' @author Carl Tony Fakhry, Ping Chen and Kourosh Zarringhalam 
#' 
#' @references Fakhry, C. T. et al. (2016). Interpreting transcriptional changes using causal 
#'             graphs: new methods and their practical utility on public networks. submitted.
#'            
#'             Franceschini, A (2013). STRING v9.1: protein-protein interaction networks, with increased coverage 
#'             and integration. In:'Nucleic Acids Res. 2013 Jan;41(Database issue):D808-15. doi: 10.1093/nar/gks1094. 
#'             Epub 2012 Nov 29'.
#'             
#' @examples # Compute the support of the Quaternary Dot Product Scoring distribution with the given margins.
#' QP_Support(50,50,50,0,50,50,50)
#' 
#' @export

QP_Support <- function(q_p, q_m, q_z, q_r, n_p, n_m, n_z){
  
  # All row constraints must be numbers.
  if (!is.numeric(q_p) || !is.numeric(q_m) || !is.numeric(q_r) || !is.numeric(q_z)){
    stop("All constraints must be integers!\n")
  }
  
  # All column constraints must be numbers.
  if (!is.numeric(n_p) || !is.numeric(n_m) || !is.numeric(n_z)){
    stop("All constraints must be integers!\n")
  }
  
  # All constraints must be positive.
  if (q_p < 0 || q_m < 0 || q_r < 0 || q_z < 0 || n_p < 0 || n_m < 0 || n_z < 0){
    stop("All constraints must be positive!\n")
  }
  
  # Sum of row constraints must equal sum of column constraints.
  if ((q_p + q_m + q_r + q_z) != (n_p + n_m + n_z)){
    stop("Sum of row constraints must equal sum of column constraints!\n")
  }
  
  # All row constraints must be integers.
  if (floor(q_p) != q_p || floor(q_m) != q_m || floor(q_r) != q_r || floor(q_z) != q_z){
    stop("All constraints must be integers!\n")
  }
  
  # All column constraints must be integers.
  if (floor(n_p) != n_p || floor(n_m) != n_m || floor(n_z) != n_z){
    stop("All constraints must be integers!\n")
  }
  
  domain <- domain_4by3(q_p, q_m, q_r, q_z, n_p, n_m, n_z)

  return(domain)
}



#' Computes the probability mass function of the scores.
#' 
#' @description This function computes the probability mass function for the Quaternary 
#' Product Scoring Statistic for signed causal graphs. This includes scores with probabilities strictly 
#' greater than zero. 
#' 
#' @usage QP_Pmf(q_p, q_m, q_z, q_r, n_p, n_m, n_z, epsilon = 1e-16)
#' 
#' @param q_p Expected number of positive predictions.
#' @param q_m Expected number of negative predictions.
#' @param q_z Expected number of nil predictions.
#' @param q_r Expected number of regulated predictions. 
#'
#' @param n_p Number of positive predictions from experiments.
#' @param n_m Number of negative predictions from experiments.
#' @param n_z Number of nil predictions from experiments. 
#' @param epsilon parameter for thresholding probabilities of matrices. Default value is 1e-16.
#' 
#' @return Vector of probabilities for scores in the support.
#' 
#' @author Carl Tony Fakhry, Ping Chen and Kourosh Zarringhalam
#' 
#' @references Fakhry, C. T. et al. (2016). Interpreting transcriptional changes using causal 
#'             graphs: new methods and their practical utility on public networks. submitted.
#'            
#'             Franceschini, A (2013). STRING v9.1: protein-protein interaction networks, with increased coverage 
#'             and integration. In:'Nucleic Acids Res. 2013 Jan;41(Database issue):D808-15. doi: 10.1093/nar/gks1094. 
#'             Epub 2012 Nov 29'.
#'             
#' @details This function computes the probability for each score in the support of the distribution.
#'          The returned value is a vector of probabilities where the 
#'          returned vector has names set equal to the corresponding scores.
#' 
#'          Setting epsilon to zero will compute the probability mass function without ignoring any matrices 
#'          with probabilities smaller than epsilon*D_max (D_max is the numerator associated with the matrix of highest
#'          probability for the given constraints). The default value of 1e-16 is experimentally validated to be 
#'          a very reasonable threshold. Setting the threshold to higher values which are smaller than 1 will lead to understimating
#'          the probabilities of each score since more tables will be ignored. 
#' 
#' @examples # Compute the probability mass function of the Quaternary 
#'           # Product Scoring Statistic for the given table margins.
#' pmf <- QP_Pmf(50,50,50,0,50,50,50)
#' 
#' @seealso \code{\link{QP_Pvalue}}, \code{\link{QP_Support}}
#' 
#' @export

QP_Pmf <- function(q_p, q_m, q_z, q_r, n_p, n_m, n_z, epsilon = 1e-16){
  
  # All row constraints must be numbers.
  if (!is.numeric(q_p) || !is.numeric(q_m) || !is.numeric(q_r) || !is.numeric(q_z)){
    stop("All constraints must be integers!\n")
  }
  
  # All column constraints must be numbers.
  if (!is.numeric(n_p) || !is.numeric(n_m) || !is.numeric(n_z)){
    stop("All constraints must be integers!\n")
  }
  
  # All constraints must be positive.
  if (q_p < 0 || q_m < 0 || q_r < 0 || q_z < 0 || n_p < 0 || n_m < 0 || n_z < 0){
    stop("All constraints must be positive!\n")
  }
  
  # Epsilon must be positive between 0 and 1.
  if (epsilon < 0 || epsilon > 1){
    stop("Epsilon must be positive between 0 and 1!\n")
  }
  
  # Sum of row constraints must equal sum of column constraints.
  if ((q_p + q_m + q_r + q_z) != (n_p + n_m + n_z)){
    stop("Sum of row constraints must equal sum of column constraints!\n")
  }
  
  # All row constraints must be integers.
  if (floor(q_p) != q_p || floor(q_m) != q_m || floor(q_r) != q_r || floor(q_z) != q_z){
    stop("All constraints must be integers!\n")
  }
  
  # All column constraints must be integers.
  if (floor(n_p) != n_p || floor(n_m) != n_m || floor(n_z) != n_z){
    stop("All constraints must be integers!\n")
  }
  
  distribution <- pmf_4by3(q_p, q_m, q_r, q_z, n_p, n_m, n_z, epsilon)

  return(distribution)
}



#' Computes the p-value of a score.
#' 
#' @description This function computes the right sided p-value for the Quaternary
#'              Product Scoring Statistic.
#' @usage QP_Pvalue(score, q_p, q_m, q_z, q_r, n_p, n_m, n_z, epsilon = 1e-16)  
#' 
#' @param score The score for which the p-value will be computed.
#' @param q_p Expected number of positive predictions.
#' @param q_m Expected number of negative predictions.
#' @param q_z Expected number of nil predictions.
#' @param q_r Expected number of regulated predictions. 
#' @param n_p Number of positive predictions from experiments.
#' @param n_m Number of negative predictions from experiments.
#' @param n_z Number of nil predictions from experiments. 
#' @param epsilon Threshold for probabilities of matrices. Default value is 1e-16.
#' 
#' @return This function returns a numerical value, where the numerical value is the p-value of the score.
#' 
#' @author Carl Tony Fakhry, Ping Chen and Kourosh Zarringhalam
#' 
#' @references Fakhry, C. T. et al. (2016). Interpreting transcriptional changes using causal 
#'             graphs: new methods and their practical utility on public networks. submitted.
#'            
#'             Franceschini, A (2013). STRING v9.1: protein-protein interaction networks, with increased coverage 
#'             and integration. In:'Nucleic Acids Res. 2013 Jan;41(Database issue):D808-15. doi: 10.1093/nar/gks1094. 
#'             Epub 2012 Nov 29'.
#'             
#' @details Setting epsilon to zero will compute the probability mass function without ignoring any matrices 
#'          with probabilities smaller than epsilon*D_max (D_max is the numerator associated with the matrix of highest
#'          probability for the given constraints). The default value of 1e-16 is experimentally validated to be 
#'          a very reasonable threshold. Setting the threshold to higher values which are smaller than 1 will lead to understimating
#'          the probabilities of each score since more tables will be ignored. 
#' 
#' @examples # Computing The p-value of score 50 
#'           # for the given table margins. 
#' pval <- QP_Pvalue(50,50,50,50,0,50,50,50)
#' 
#' @seealso \code{\link{QP_SigPvalue}}
#' 
#' @export

QP_Pvalue <- function(score, q_p, q_m, q_z, q_r, n_p, n_m, n_z, epsilon = 1e-16){
  
  # All row constraints must be numbers.
  if (!is.numeric(q_p) || !is.numeric(q_m) || !is.numeric(q_r) || !is.numeric(q_z)){
    stop("All constraints must be integers!\n")
  }
  
  # All column constraints must be numbers.
  if (!is.numeric(n_p) || !is.numeric(n_m) || !is.numeric(n_z)){
    stop("All constraints must be integers!\n")
  }
  
  # All constraints must be positive.
  if (q_p < 0 || q_m < 0 || q_r < 0 || q_z < 0 || n_p < 0 || n_m < 0 || n_z < 0){
    stop("All constraints must be positive!\n")
  }
  
  # Epsilon must be positive between 0 and 1.
  if (epsilon < 0 || epsilon > 1){
    stop("Epsilon must be positive between 0 and 1!\n")
  }
  
  # Sum of row constraints must equal sum of column constraints.
  if ((q_p + q_m + q_r + q_z) != (n_p + n_m + n_z)){
    stop("Sum of row constraints must equal sum of column constraints!\n")
  }
  
  # All row constraints must be integers.
  if (floor(q_p) != q_p || floor(q_m) != q_m || floor(q_r) != q_r || floor(q_z) != q_z){
    stop("All constraints must be integers!\n")
  }
  
  # All column constraints must be integers.
  if (floor(n_p) != n_p || floor(n_m) != n_m || floor(n_z) != n_z){
    stop("All constraints must be integers!\n")
  }
  
  pval <- Pvalue_4by3(q_p, q_m, q_r, q_z, n_p, n_m, n_z, score, epsilon)

  return(pval)
}



#' Computes the p-value for a statistically significant score.
#' 
#' @description This function computes the right sided p-value for the Quaternary
#'              Product Scoring Statistic for statistically significant scores.
#' @usage QP_SigPvalue(score, q_p, q_m, q_z, q_r, n_p, n_m, n_z, epsilon = 1e-16, sig_level = 0.05)  
#' 
#' @param score The score for which the p-value will be computed.
#' @param q_p Expected number of positive predictions.
#' @param q_m Expected number of negative predictions.
#' @param q_z Expected number of nil predictions.
#' @param q_r Expected number of regulated predictions. 
#' @param n_p Number of positive predictions from experiments.
#' @param n_m Number of negative predictions from experiments.
#' @param n_z Number of nil predictions from experiments. 
#' @param epsilon Threshold for probabilities of matrices. Default value is 1e-16.
#' @param sig_level Significance level of test hypothesis. Default value is 0.05.
#' 
#' @return This function returns a numerical value, where the numerical value is the p-value of a score if the score
#'         is statistically significant otherwise it returns -1.
#'         
#' 
#' @author Carl Tony Fakhry, Ping Chen and Kourosh Zarringhalam
#' 
#' @references Fakhry, C. T. et al. (2016). Interpreting transcriptional changes using causal 
#'             graphs: new methods and their practical utility on public networks. submitted.
#'            
#'             Franceschini, A (2013). STRING v9.1: protein-protein interaction networks, with increased coverage 
#'             and integration. In:'Nucleic Acids Res. 2013 Jan;41(Database issue):D808-15. doi: 10.1093/nar/gks1094. 
#'             Epub 2012 Nov 29'.
#'              
#' @details Setting epsilon to zero will compute the probability mass function without ignoring any matrices 
#'          with probabilities smaller than epsilon*D_max (D_max is the numerator associated with the matrix of highest
#'          probability for the given constraints). The default value of 1e-16 is experimentally validated to be 
#'          a very reasonable threshold. Setting the threshold to higher values which are smaller than 1 will lead to understimating
#'          the probabilities of each score since more tables will be ignored. If the score is not statistically significant, 
#'          then a value of -1 will be returned.
#' 
#' @examples # Computing The p-value of score 50 
#'           # for the given table margins. 
#' pval <- QP_SigPvalue(50,50,50,50,0,50,50,50)
#' 
#' @seealso \code{\link{QP_Pvalue}}
#' 
#' @export

QP_SigPvalue <- function(score, q_p, q_m, q_z, q_r, n_p, n_m, n_z, epsilon = 1e-16, sig_level = 0.05){
  
  # All row constraints must be numbers.
  if (!is.numeric(q_p) || !is.numeric(q_m) || !is.numeric(q_r) || !is.numeric(q_z)){
    stop("All constraints must be integers!\n")
  }
  
  # All column constraints must be numbers.
  if (!is.numeric(n_p) || !is.numeric(n_m) || !is.numeric(n_z)){
    stop("All constraints must be integers!\n")
  }
  
  # All constraints must be positive.
  if (q_p < 0 || q_m < 0 || q_r < 0 || q_z < 0 || n_p < 0 || n_m < 0 || n_z < 0){
    stop("All constraints must be positive!\n")
  }
  
  # Epsilon must be positive between 0 and 1.
  if (epsilon < 0 || epsilon > 1){
    stop("Epsilon must be positive between 0 and 1!\n")
  }
  
  # Sum of row constraints must equal sum of column constraints.
  if ((q_p + q_m + q_r + q_z) != (n_p + n_m + n_z)){
    stop("Sum of row constraints must equal sum of column constraints!\n")
  }
  
  # All row constraints must be integers.
  if (floor(q_p) != q_p || floor(q_m) != q_m || floor(q_r) != q_r || floor(q_z) != q_z){
    stop("All constraints must be integers!\n")
  }
  
  # All column constraints must be integers.
  if (floor(n_p) != n_p || floor(n_m) != n_m || floor(n_z) != n_z){
    stop("All constraints must be integers!\n")
  }
  
  pval <- SigPvalue_4by3(q_p, q_m, q_r, q_z, n_p, n_m, n_z, score, epsilon, sig_level)
  
  return(pval)
}



#' Computes the probability of a score.
#' 
#' @description This function computes the probability of a score in the Ternary or Quaternary
#'              Product scoring distribution.
#' @usage QP_Probability(score, q_p, q_m, q_z, q_r, n_p, n_m, n_z, epsilon = 1e-16)  
#' 
#' @param score The score for which the probability will be computed.
#' @param q_p Expected number of positive predictions.
#' @param q_m Expected number of negative predictions.
#' @param q_z Expected number of nil predictions.
#' @param q_r Expected number of regulated predictions. 
#' @param n_p Number of positive predictions from experiments.
#' @param n_m Number of negative predictions from experiments.
#' @param n_z Number of nil predictions from experiments. 
#' @param epsilon Threshold for probabilities of matrices. Default value is 1e-16.
#' 
#' @return This function returns a numerical value, where the numerical value is the probability of the score.
#' 
#' @author Carl Tony Fakhry, Ping Chen and Kourosh Zarringhalam
#' 
#' @references Fakhry, C. T. et al. (2016). Interpreting transcriptional changes using causal 
#'             graphs: new methods and their practical utility on public networks. submitted.
#'            
#'             Franceschini, A (2013). STRING v9.1: protein-protein interaction networks, with increased coverage 
#'             and integration. In:'Nucleic Acids Res. 2013 Jan;41(Database issue):D808-15. doi: 10.1093/nar/gks1094. 
#'             Epub 2012 Nov 29'.
#'              
#' @details Setting epsilon to zero will compute the probability mass function without ignoring any matrices 
#'          with probabilities smaller than epsilon*D_max (D_max is the numerator associated with the matrix of highest
#'          probability for the given constraints). The default value of 1e-16 is experimentally validated to be 
#'          a very reasonable threshold. Setting the threshold to higher values which are smaller than 1 will lead to understimating
#'          the probabilities of each score since more tables will be ignored. 
#'          
#'          For computing p-values, the user is advised to use the p-value function which is optimized
#'          for such purposes.
#' 
#' @examples # Computing The probability of score 50 
#'           # for the given table margins. 
#' prob <- QP_Probability(0,50,50,50,0,50,50,50)
#' 
#' @seealso \code{\link{QP_Pmf}}, \code{\link{QP_Pvalue}}, \code{\link{QP_SigPvalue}}
#' 
#' @export

QP_Probability <- function(score, q_p, q_m, q_z, q_r, n_p, n_m, n_z, epsilon = 1e-16){
  
  # All row constraints must be numbers.
  if (!is.numeric(q_p) || !is.numeric(q_m) || !is.numeric(q_r) || !is.numeric(q_z)){
    stop("All constraints must be integers!\n")
  }
  
  # All column constraints must be numbers.
  if (!is.numeric(n_p) || !is.numeric(n_m) || !is.numeric(n_z)){
    stop("All constraints must be integers!\n")
  }
  
  # All constraints must be positive.
  if (q_p < 0 || q_m < 0 || q_r < 0 || q_z < 0 || n_p < 0 || n_m < 0 || n_z < 0){
    stop("All constraints must be positive!\n")
  }
  
  # Epsilon must be positive between 0 and 1.
  if (epsilon < 0 || epsilon > 1){
    stop("Epsilon must be positive between 0 and 1!\n")
  }
  
  # Sum of row constraints must equal sum of column constraints.
  if ((q_p + q_m + q_r + q_z) != (n_p + n_m + n_z)){
    stop("Sum of row constraints must equal sum of column constraints!\n")
  }
  
  # All row constraints must be integers.
  if (floor(q_p) != q_p || floor(q_m) != q_m || floor(q_r) != q_r || floor(q_z) != q_z){
    stop("All constraints must be integers!\n")
  }
  
  # All column constraints must be integers.
  if (floor(n_p) != n_p || floor(n_m) != n_m || floor(n_z) != n_z){
    stop("All constraints must be integers!\n")
  }
  
  prob <- user_probability_of_score_4by3(score, q_p, q_m, q_r, q_z, n_p, n_m, n_z, epsilon)

  return(prob)
}

