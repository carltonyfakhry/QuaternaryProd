// This is the header file for the Quaternary Product Scoring Statistic utils functions.

#ifndef QUATERNARYPRODUTILS
#define QUATERNARYPRODUTILS
 
#include <Rcpp.h>
using namespace Rcpp;



/**
 *
 * This function gets a vector (i.e matrix in the dictionary ordering
 * of our method) with the maximum score.
 *
 */
IntegerVector max_arr_score_4by3(IntegerVector constraints);



/**
 *
 * This function gets the vector (i.e matrix in the dictionary ordering
 * of our method) with the minimum score.
 *
 */
IntegerVector min_arr_score_4by3(IntegerVector constraints);
 
 

 /** 
  * 
  * This function checks if two vectors are equal.
  * 
  */
bool are_equal(IntegerVector vec1, IntegerVector vec2);
 
 
 
 /**
  * 
  * Computes the measure of evenness.
  * 
  * 
  */
 double get_measure_4by3(IntegerVector vec);

   
   
/**
 *
 * Get the numerator of the numerator of the Quaternary Product 
 * Probability distribution.
 * 
 */
double get_numerator_4by3(IntegerVector constraints);



/**
 * 
 * Get the denominator of the Quaternary Product Probability.
 * 
 */
double get_total_4by3(IntegerVector constraints);



/**
 * 
 * Compute the D-value of the table.
 * 
 */
double get_Dvalue_4by3(IntegerVector vec, double num);
  
  

/**
 * 
 * This function computes the probability of a vector (i.e the matrix 
 * in the dictionary odering of this method).
 * 
 */
double get_probability_4by3(IntegerVector vec, double num, double total);



/**
 * 
 * This function populates all the score increasing moves.
 * 
 */
IntegerMatrix get_score_increase_moves_4by3();



/**
 * 
 * This function populates all the score decreasing moves.
 * 
 */
IntegerMatrix get_score_decrease_moves_4by3();



/**
 * 
 * This function populates all the score preserving moves
 * with decreasing dictionary ordering.
 * 
 */
IntegerMatrix get_dict_decrease_moves_4by3();



/**
 * 
 * This function populates all the score preserving moves
 * with increasing dictionary ordering.
 * 
 */
IntegerMatrix get_dict_increase_moves_4by3();



/**
 * 
 * This function computes the score.
 * 
 */
int get_score_4by3(IntegerVector vec);



/**
 * 
 * This function computes the constraints of a vector(i.e a matrix 
 * in the dictionary ordering of our method).
 * 
 */
IntegerVector get_constraints_4by3(IntegerVector vec);



/** 
 * 
 * Set the maximum dictionary ordering of a vector keeping n_pp, 
 * n_mm, n_rp and n_rm fixed in the dictionary ordering of this method.
 * 
 */ 
IntegerVector set_max_s1_4by3(IntegerVector vec);



/** 
 * 
 * Set the minimum dictionary ordering of a vector keeping n_pp, 
 * n_mm, n_rp and n_rm fixed in the dictionary ordering of this method.
 * 
 */ 
IntegerVector set_min_s1_4by3(IntegerVector vec);



/**
 * 
 * This function makes all moves and takes the one with the greatest Dvalue
 * (Note: log(D-value) >= 0).
 * 
 */
IntegerVector find_move_4by3(IntegerMatrix moves, IntegerVector current_state, double num);



/**
 * 
 * This function finds the matrix of the next bigger score with the highest D-value.
 * 
 */
IntegerVector get_next_max_score_4by3(IntegerMatrix moves, IntegerVector vec, double num);



/**
 * 
 * This function finds the matrix of the next smaller score with the highest D-value.
 * 
 */
IntegerVector get_next_max_score2_4by3(IntegerMatrix moves, IntegerVector vec, double num);



/**
 * 
 * Gets a matrix with the highest Dvalue (i.e highest probability) of a given score if possible else
 * it returns an empty vector.
 * 
 */
IntegerVector get_mat_of_score_4by3(IntegerVector constraints, int score_target, double num);
  
  

/**
 * 
 * This function computes a very accurate estimate of the 
 * largest D-value.
 *
 */
NumericVector max_Dvalue_arr_4by3(IntegerVector constraints);



#endif
