// This is the header file for the Quaternary Dot Product Scoring Statistic functions.

#ifndef QUATERNARYPROD
#define QUATERNARYPROD
 
#include <Rcpp.h>
using namespace Rcpp;
 
 

/**
 * 
 * Get the probability of a given score.
 * 
 */
double probability_of_score_4by3(IntegerVector arr, double num, double total, double epsilon = 1e-16);
  
  
  
/**
 * 
 * This function computes the right sided Pvalue of a given score of the
 * Quaternary Dot Product Scoring Statistic.
 * 
 */
double Pvalue_4by3(int q_p, int q_m, int q_r, int q_z, int n_p, int n_m, int n_z, int target_score, double epsilon = 1e-16);
 
 

/**
  * 
  * This function computes the right sided Pvalue of a statistically significant score of the
  * Quaternary Dot Product Scoring Statistic otherwise if the score is not statistically significant,
  * a value of -1 is returned.
  * 
  */
double SigPvalue_4by3(int q_p, int q_m, int q_r, int q_z, int n_p, int n_m, int n_z, int target_score, double epsilon = 1e-16, double significance_level = 0.05);
 
 
  
/**
 * 
 * Computes the support of the Quaternary Dot Product Scoring Statistic.
 * 
 */
IntegerVector domain_4by3(int q_p, int q_m, int q_r, int q_z, int n_p, int n_m, int n_z);
 


/**
 * 
 * Computes the pmf of the Quaternary Dot Product Scoring Statistic.
 * 
 */
NumericVector pmf_4by3(int q_p, int q_m, int q_r, int q_z, int n_p, int n_m, int n_z, double epsilon = 1e-16);
 
 

/**
 * 
 * This is the user's method for computing the probability of a
 * score.
 * 
 */
double user_probability_of_score_4by3(int target_score, int q_p, int q_m, int q_r, int q_z, int n_p, int n_m, int n_z, double epsilon);
  
  
  
#endif