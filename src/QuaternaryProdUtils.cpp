#include <Rcpp.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R.h>
#include <math.h>
#include "QuaternaryProdUtils.h"
using namespace Rcpp;



// Authors: Carl Tony Fakhry, Kourosh Zarringhalam, Ping Chen.

 
 
/**
 *
 * This function gets a vector (i.e matrix in the dictionary ordering
 * of our method) with the maximum score.
 *
 */

IntegerVector max_arr_score_4by3(IntegerVector constraints){
  int q_p = constraints[0];  
  int q_m = constraints[1];  
  int q_r = constraints[2];  
  int q_z = constraints[3];  
  int n_p = constraints[4];  
  int n_m = constraints[5]; 
  int n_pp = min(IntegerVector::create(q_p, n_p));
  int n_mm = min(IntegerVector::create(q_m, n_m));
  int n_rp = min(IntegerVector::create(q_r, n_p - n_pp));
  int n_rm = min(IntegerVector::create(q_r - n_rp, n_m - n_mm));
  int n_zp = min(IntegerVector::create(q_z, n_p - n_pp - n_rp));
  int n_zm = min(IntegerVector::create(q_z - n_zp, n_m - n_mm - n_rm));
  int n_zz = q_z - n_zp - n_zm;
  int n_rz = q_r - n_rp - n_rm;
  int n_mp = n_p - n_pp - n_rp - n_zp;
  int n_pm = n_m - n_mm - n_rm - n_zm;
  int n_pz = q_p - n_pp - n_pm;
  int n_mz = q_m - n_mm - n_mp;
  IntegerVector res = IntegerVector::create(n_rp, n_rm, n_pp, n_mm, n_mp, n_pm, n_pz, n_mz, n_rz, n_zp, n_zm, n_zz);
  return res;  
} 



/**
 *
 * This function gets the vector (i.e matrix in the dictionary ordering
 * of our method) with the minimum score.
 *
 */

IntegerVector min_arr_score_4by3(IntegerVector constraints){
  int q_p = constraints[0];  
  int q_m = constraints[1];  
  int q_r = constraints[2];  
  int q_z = constraints[3];  
  int n_p = constraints[4];  
  int n_m = constraints[5]; 
  int n_z = constraints[6];
  int n_pm = min(IntegerVector::create(q_p, n_m));
  int n_mp = min(IntegerVector::create(q_m, n_p));
  int n_zp = min(IntegerVector::create(q_z, n_p - n_mp));
  int n_zm = min(IntegerVector::create(q_z - n_zp, n_m - n_pm));
  int n_pz = min(IntegerVector::create(q_p - n_pm, n_z));
  int n_mz = min(IntegerVector::create(q_m - n_mp, n_z - n_pz));
  int n_rp = min(IntegerVector::create(q_r, n_p - n_mp - n_zp));
  int n_rm = min(IntegerVector::create(q_r - n_rp, n_m - n_pm - n_zm));
  int n_pp = q_p - n_pm - n_pz;
  int n_mm = q_m - n_mp - n_mz;
  int n_rz = q_r - n_rp - n_rm;
  int n_zz = q_z - n_zm - n_zp;
  IntegerVector res = IntegerVector::create(n_rp, n_rm, n_pp, n_mm, n_mp, n_pm, n_pz, n_mz, n_rz, n_zp, n_zm, n_zz);
  if (q_r > 0 && (n_p > q_m || n_m > q_p)){
    IntegerVector modified;
    IntegerVector s1 = IntegerVector::create(-1,1,0,0,0,0,0,0,0,1,-1,0);
    modified = s1 + res;
    if (is_true(all(modified >= 0))){
      double n = res[0] + res[9];
      double K = res[0] + res[1];
      double N = res[0] + res[1] + res[9] + res[10];
      double hypergeometric_mode = floor(((n+1)*(K+1))/(N+2));
      int max_move = res[0] - hypergeometric_mode;
      modified = s1*max_move + res;
      return modified;
    }
    IntegerVector s2 = IntegerVector::create(-1,0,1,0,0,0,-1,0,1,0,0,0);
    modified = s2 + res;
    if (is_true(all(modified >= 0))){
      double n = res[0] + res[2];
      double K = res[0] + res[8];
      double N = res[0] + res[2] + res[6] + res[8];
      double hypergeometric_mode = floor(((n+1)*(K+1))/(N+2));
      int max_move = res[0] - hypergeometric_mode;
      modified = s2*max_move + res;
      return modified;
    }
    IntegerVector s3 = IntegerVector::create(0,-1,0,1,0,0,0,-1,1,0,0,0);
    modified = s3 + res;
    if (is_true(all(modified >= 0))){
      double n = res[1] + res[3];
      double K = res[1] + res[8];
      double N = res[1] + res[3] + res[7] + res[8];
      double hypergeometric_mode = floor(((n+1)*(K+1))/(N+2));
      int max_move = res[1] - hypergeometric_mode;
      modified = s3*max_move + res;
      return modified;
    }
  }
  return res;  
}  



/**
 * 
 * Computes the measure of evenness.
 * 
 * 
 */

double get_measure_4by3(IntegerVector vec){
  double diff = 0;
  for(int i = 0 ; i < vec.size()-1 ; i++){
    for(int j = i+1; j < vec.size() ; j++){
      diff += pow(vec[i] - vec[j],2);
    }
  }  
  return sqrt(diff);
}



/**
 *
 * Get the numerator of the numerator of the Quaternary Product 
 * Probability distribution.
 * 
 */

double get_numerator_4by3(IntegerVector constraints){
  IntegerVector rows = IntegerVector::create(constraints[0], constraints[1], constraints[2], constraints[3]);
  double num = sum(lfactorial(rows));
  return num;
}



/**
 * 
 * Get the denominator of the Quaternary Product Probability.
 * 
 */

double get_total_4by3(IntegerVector constraints){
  IntegerVector cols = IntegerVector::create(constraints[4], constraints[5], constraints[6]);
  double total1 = Rcpp::internal::lfactorial(sum(cols));
  double total2 = sum(lfactorial(cols));
  double total = total1 - total2;
  return total;
}



/**
 * 
 * Compute the D-value of the table.
 * 
 */

double get_Dvalue_4by3(IntegerVector vec, double num){
  double dvalue = num - sum(lfactorial(vec));
  return dvalue;
}



/**
 * 
 * This function computes the probability of a vector (i.e the matrix 
 * in the dictionary odering of this method).
 * 
 */

double get_probability_4by3(IntegerVector vec, double num, double total){
  double tot = sum(lfactorial(vec));
  double lprobability = exp((num - tot) - total);
  return lprobability;
}



/**
 * 
 * This function populates all the score increasing moves.
 * 
 */

IntegerMatrix get_score_increase_moves_4by3(){
  IntegerMatrix moves(23,12);    
  // Moves that increase the score by 1
  moves.row(0) = NumericVector::create(-1, 1, 1, 0, 0, 0, -1, 0, 0, 0, -1, 1);
  moves.row(1) = NumericVector::create(-1, 1, 1, 0, 1, -1, 0, -1, 0, -1, 0, 1);
  moves.row(2) = NumericVector::create(0, 0, 0, 0, -1, 0, 0, 1, 0, 1, 0, -1);
  moves.row(3) = NumericVector::create(0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 1, -1);
  moves.row(4) = NumericVector::create(0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -2, 1);
  moves.row(5) = NumericVector::create(0, 0, 0, 1, 0, 0, 0, -1, 0, 0, -1, 1);
  moves.row(6) = NumericVector::create(0, 0, 0, 1, 1, -1, 1, -2, 0, -1, 0, 1);
  moves.row(7) = NumericVector::create(0, 0, 1, 0, -1, 1, -2, 1, 0, 0, -1, 1);
  moves.row(8) = NumericVector::create(0, 0, 1, 0, 0, 0, -1, 0, 0, -1, 0, 1);
  moves.row(9) = NumericVector::create(0, 0, 1, 0, 1, -1, 0, -1, 0, -2, 1, 1);
  moves.row(10) = NumericVector::create(0, 1, 0, 0, 0, 0, 0, 0, -1, 0, -1, 1);
  moves.row(11) = NumericVector::create(1, -1, 0, 1, -1, 1, -1, 0, 0, 0, -1, 1);
  moves.row(12) = NumericVector::create(1, -1, 0, 1, 0, 0, 0, -1, 0, -1, 0, 1);
  moves.row(13) = NumericVector::create(1, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 1);
  // Moves that increase the score by 2
  moves.row(14) = NumericVector::create(-1, 1, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0);
  moves.row(15) = NumericVector::create(0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, 0);
  moves.row(16) = NumericVector::create(0, 0, 0, 1, 0, -1, 1, -1, 0, 0, 0, 0);
  moves.row(17) = NumericVector::create(0, 0, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0);
  moves.row(18) = NumericVector::create(0, 0, 1, 0, 0, -1, 0, 0, 0, -1, 1, 0);
  moves.row(19) = NumericVector::create(0, 1, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0);
  moves.row(20) = NumericVector::create(1, -1, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0);
  moves.row(21) = NumericVector::create(1, 0, 0, 0, -1, 0, 0, 1, -1, 0, 0, 0);
  // Moves that increase the score by 4
  moves.row(22) = NumericVector::create(0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 0, 0);
  return moves;
}


/**
 * 
 * This function populates all the score decreasing moves.
 * 
 */

IntegerMatrix get_score_decrease_moves_4by3(){
  IntegerMatrix moves(33,12);
  // Moves that decrease the score by -1
  moves.row(0) = NumericVector::create(-1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1);
  moves.row(1) = NumericVector::create(-1, 0, 0, 1, 1, 0, 0, -2, 1, 0, -1, 1);
  moves.row(2) = NumericVector::create(-1, 1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 1);
  moves.row(3) = NumericVector::create(-1, 1, 1, -1, 1, 0, -1, 0, 0, -1, 0, 1);
  moves.row(4) = NumericVector::create(-1, 2, 0, -1, 1, 0, 0, 0, -1, 0, -1, 1);
  moves.row(5) = NumericVector::create(0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 1, -1);
  moves.row(6) = NumericVector::create(0, -1, 1, 0, 0, 1, -2, 0, 1, -1, 0, 1);
  moves.row(7) = NumericVector::create(0, 0, -1, 0, 0, 0, 1, 0, 0, 1, 0, -1);
  moves.row(8) = NumericVector::create(0, 0, -1, 1, 0, 1, 0, -1, 0, 1, -2, 1);
  moves.row(9) = NumericVector::create(0, 0, -1, 1, 1, 0, 1, -2, 0, 0, -1, 1);
  moves.row(10) = NumericVector::create(0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 1, -1);
  moves.row(11) = NumericVector::create(0, 0, 0, 0, 0, 1, -1, 0, 0, 0, -1, 1);
  moves.row(12) = NumericVector::create(0, 0, 0, 0, 1, 0, 0, -1, 0, -1, 0, 1);
  moves.row(13) = NumericVector::create(0, 0, 1, -1, 0, 1, -2, 1, 0, -1, 0, 1);
  moves.row(14) = NumericVector::create(0, 0, 1, -1, 1, 0, -1, 0, 0, -2, 1, 1);
  moves.row(15) = NumericVector::create(0, 1, -1, 0, 0, 1, 0, 0, -1, 1, -2, 1);
  moves.row(16) = NumericVector::create(0, 1, -1, 0, 1, 0, 1, -1, -1, 0, -1, 1);
  moves.row(17) = NumericVector::create(0, 1, 0, -1, 1, 0, 0, 0, -1, -1, 0, 1);
  moves.row(18) = NumericVector::create(1, -1, -1, 1, 0, 1, 0, -1, 0, 0, -1, 1);
  moves.row(19) = NumericVector::create(1, -1, 0, 0, 0, 1, -1, 0, 0, -1, 0, 1);
  moves.row(20) = NumericVector::create(1, 0, -1, 0, 0, 1, 0, 0, -1, 0, -1, 1);
  moves.row(21) = NumericVector::create(1, 0, 0, -1, 0, 1, -1, 1, -1, -1, 0, 1);
  moves.row(22) = NumericVector::create(2, -1, -1, 0, 0, 1, 0, 0, -1, -1, 0, 1);
  // Moves that decrease the score by -2
  moves.row(23) = NumericVector::create(-1, 0, 0, 0, 1, 0, 0, -1, 1, 0, 0, 0);
  moves.row(24) = NumericVector::create(-1, 1, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0);
  moves.row(25) = NumericVector::create(0, -1, 0, 0, 0, 1, -1, 0, 1, 0, 0, 0);
  moves.row(26) = NumericVector::create(0, 0, -1, 0, 0, 1, 0, 0, 0, 1, -1, 0);
  moves.row(27) = NumericVector::create(0, 0, -1, 0, 1, 0, 1, -1, 0, 0, 0, 0);
  moves.row(28) = NumericVector::create(0, 0, 0, -1, 0, 1, -1, 1, 0, 0, 0, 0);
  moves.row(29) = NumericVector::create(0, 0, 0, -1, 1, 0, 0, 0, 0, -1, 1, 0);
  moves.row(30) = NumericVector::create(1, -1, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0);
  moves.row(31) = NumericVector::create(1, 0, -1, -1, 0, 1, 0, 1, -1, 0, 0, 0);
  // Moves that decrease the score by -4
  moves.row(32) = NumericVector::create(0, 0, -1, -1, 1, 1, 0, 0, 0, 0, 0, 0);
  return moves;
}


/**
 * 
 * This function populates all the score preserving moves
 * with decreasing dictionary ordering.
 * 
 */

IntegerMatrix get_dict_decrease_moves_4by3(){
  IntegerMatrix moves(25,12);     
  moves.row(0) = NumericVector::create(0, 0, 0, 0, -1, 1, -1, 1, 0, 1, -1, 0);
  moves.row(1) = NumericVector::create(0, 0, 0, -1, -1, 0, 0, 2, 0, 1, 1, -2);
  moves.row(2) = NumericVector::create(0, 0, 0, -1, 0, -1, 1, 1, 0, 0, 2, -2);
  moves.row(3) = NumericVector::create(0, 0, -1, 0, -1, 0, 1, 1, 0, 2, 0, -2);
  moves.row(4) = NumericVector::create(0, 0, -1, 0, 0, -1, 2, 0, 0, 1, 1, -2);
  moves.row(5) = NumericVector::create(0, 0, -1, 1, -1, 1, 0, 0, 0, 2, -2, 0);
  moves.row(6) = NumericVector::create(0, 0, -1, 1, 0, 0, 1, -1, 0, 1, -1, 0);
  moves.row(7) = NumericVector::create(0, 0, -1, 1, 1, -1, 2, -2, 0, 0, 0, 0);
  moves.row(8) = NumericVector::create(0, -1, 0, 0, -1, 0, 0, 1, 1, 1, 1, -2);
  moves.row(9) = NumericVector::create(0, -1, 0, 0, 0, -1, 1, 0, 1, 0, 2, -2);
  moves.row(10) = NumericVector::create(0, -1, 0, 1, -1, 1, -1, 0, 1, 1, -1, 0);
  moves.row(11) = NumericVector::create(0, -1, 0, 1, 0, 0, 0, -1, 1, 0, 0, 0);
  moves.row(12) = NumericVector::create(0, -1, 1, 0, -1, 1, -2, 1, 1, 0, 0, 0);
  moves.row(13) = NumericVector::create(0, -1, 1, 0, 0, 0, -1, 0, 1, -1, 1, 0);
  moves.row(14) = NumericVector::create(-1, 2, 0, -1, 1, -1, 1, 0, -1, 0, 0, 0);
  moves.row(15) = NumericVector::create(-1, 1, 1, -1, 1, -1, 0, 0, 0, -1, 1, 0);
  moves.row(16) = NumericVector::create(-1, 1, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0);
  moves.row(17) = NumericVector::create(-1, 1, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0);
  moves.row(18) = NumericVector::create(-1, 0, 1, 0, 1, -1, 0, -1, 1, -1, 1, 0);
  moves.row(19) = NumericVector::create(-1, 0, 1, 0, 0, 0, -1, 0, 1, 0, 0, 0);
  moves.row(20) = NumericVector::create(-1, 0, 0, 1, 1, -1, 1, -2, 1, 0, 0, 0);
  moves.row(21) = NumericVector::create(-1, 0, 0, 0, 0, -1, 1, 0, 1, 1, 1, -2);
  moves.row(22) = NumericVector::create(-1, 0, 0, 0, -1, 0, 0, 1, 1, 2, 0, -2);
  moves.row(23) = NumericVector::create(-2, 2, 1, -1, 1, -1, 0, 0, 0, 0, 0, 0);
  moves.row(24) = NumericVector::create(-2, 1, 1, 0, 1, -1, 0, -1, 1, 0, 0, 0);
  return moves;
}


/**
 * 
 * This function populates all the score preserving moves
 * with increasing dictionary ordering.
 * 
 */

IntegerMatrix get_dict_increase_moves_4by3(){
  IntegerMatrix moves(38,12);
  moves.row(0) = NumericVector::create(0, 0, 0, 0, 1, -1, 1, -1, 0, -1, 1, 0);
  moves.row(1) = NumericVector::create(0, 0, 0, 1, -1, 2, -2, 0, 0, 1, -3, 2);
  moves.row(2) = NumericVector::create(0, 0, 0, 1, 0, 1, -1, -1, 0, 0, -2, 2);
  moves.row(3) = NumericVector::create(0, 0, 0, 1, 1, 0, 0, -2, 0, -1, -1, 2);
  moves.row(4) = NumericVector::create(0, 0, 0, 1, 2, -1, 1, -3, 0, -2, 0, 2);
  moves.row(5) = NumericVector::create(0, 0, 1, -3, -1, -1, 0, 4, 0, 0, 4, -4);
  moves.row(6) = NumericVector::create(0, 0, 1, -2, -1, 0, -1, 3, 0, 0, 2, -2);
  moves.row(7) = NumericVector::create(0, 0, 1, -2, 0, -1, 0, 2, 0, -1, 3, -2);
  moves.row(8) = NumericVector::create(0, 0, 1, -1, -1, 1, -2, 2, 0, 0, 0, 0);
  moves.row(9) = NumericVector::create(0, 0, 1, -1, 0, 0, -1, 1, 0, -1, 1, 0);
  moves.row(10) = NumericVector::create(0, 0, 1, -1, 1, -1, 0, 0, 0, -2, 2, 0);
  moves.row(11) = NumericVector::create(0, 0, 1, 0, -1, 2, -3, 1, 0, 0, -2, 2);
  moves.row(12) = NumericVector::create(0, 0, 1, 0, 0, 1, -2, 0, 0, -1, -1, 2);
  moves.row(13) = NumericVector::create(0, 0, 1, 0, 1, 0, -1, -1, 0, -2, 0, 2);
  moves.row(14) = NumericVector::create(0, 0, 1, 0, 2, -1, 0, -2, 0, -3, 1, 2);
  moves.row(15) = NumericVector::create(0, 1, -2, 0, 0, -1, 3, 0, -1, 2, 0, -2);
  moves.row(16) = NumericVector::create(0, 1, -1, 0, 0, 0, 1, 0, -1, 1, -1, 0);
  moves.row(17) = NumericVector::create(0, 1, -1, 0, 1, -1, 2, -1, -1, 0, 0, 0);
  moves.row(18) = NumericVector::create(0, 1, 0, -1, 0, 0, 0, 1, -1, 0, 0, 0);
  moves.row(19) = NumericVector::create(0, 1, 0, 0, 0, 1, -1, 0, -1, 0, -2, 2);
  moves.row(20) = NumericVector::create(0, 1, 0, 0, 1, 0, 0, -1, -1, -1, -1, 2);
  moves.row(21) = NumericVector::create(0, 1, 0, 0, 2, -1, 1, -2, -1, -2, 0, 2);
  moves.row(22) = NumericVector::create(1, -2, 0, 0, -1, 0, 0, 1, 1, 0, 2, -2);
  moves.row(23) = NumericVector::create(1, -1, -1, 0, -1, 0, 1, 1, 0, 1, 1, -2);
  moves.row(24) = NumericVector::create(1, -1, -1, 0, 0, -1, 2, 0, 0, 0, 2, -2);
  moves.row(25) = NumericVector::create(1, -1, -1, 1, -1, 1, 0, 0, 0, 1, -1, 0);
  moves.row(26) = NumericVector::create(1, -1, -1, 1, 0, 0, 1, -1, 0, 0, 0, 0);
  moves.row(27) = NumericVector::create(1, -1, 0, -1, -1, 0, 0, 2, 0, 0, 2, -2);
  moves.row(28) = NumericVector::create(1, -1, 0, 0, -1, 1, -1, 1, 0, 0, 0, 0);
  moves.row(29) = NumericVector::create(1, -1, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0);
  moves.row(30) = NumericVector::create(1, 0, -1, 0, 0, 0, 1, 0, -1, 0, 0, 0);
  moves.row(31) = NumericVector::create(1, 0, 0, -2, -1, 0, 0, 3, -1, 0, 2, -2);
  moves.row(32) = NumericVector::create(1, 0, 0, -1, -1, 1, -1, 2, -1, 0, 0, 0);
  moves.row(33) = NumericVector::create(1, 0, 0, -1, 0, 0, 0, 1, -1, -1, 1, 0);
  moves.row(34) = NumericVector::create(1, 0, 0, 0, -1, 2, -2, 1, -1, 0, -2, 2);
  moves.row(35) = NumericVector::create(1, 0, 0, 0, 0, 1, -1, 0, -1, -1, -1, 2);
  moves.row(36) = NumericVector::create(1, 0, 0, 0, 1, 0, 0, -1, -1, -2, 0, 2);
  moves.row(37) = NumericVector::create(2, -2, -1, 1, -1, 1, 0, 0, 0, 0, 0, 0);
  return moves;
}



/**
 * 
 * This function computes the score.
 * 
 */

int get_score_4by3(IntegerVector vec){
  return vec[0] + vec[1] + vec[2] + vec[3] - (vec[4] + vec[5]);
}


/**
 * 
 * This function computes the constraints of a vector(i.e a matrix 
 * in the dictionary ordering of our method).
 * 
 */

IntegerVector get_constraints_4by3(IntegerVector vec){
  int q_p = vec[2] + vec[5] + vec[6];
  int q_m = vec[3] + vec[4] + vec[7];
  int q_r = vec[0] + vec[1] + vec[8]; 
  int q_z = vec[9] + vec[10] + vec[11];
  int n_p = vec[0] + vec[2] + vec[4] + vec[9];
  int n_m = vec[1] + vec[3] + vec[5] + vec[10];
  int n_z = vec[6] + vec[7] + vec[8] + vec[11];
  IntegerVector res = IntegerVector::create(q_p,q_m,q_r,q_z,n_p,n_m,n_z);
  return res;  
}


/** 
 * 
 * Set the maximum dictionary ordering of a vector keeping n_pp, 
 * n_mm, n_rp and n_rm fixed in the dictionary ordering of this method.
 * 
 */ 

IntegerVector set_max_s1_4by3(IntegerVector vec){
  IntegerVector maximum = clone(vec);
  IntegerVector inv_1 = IntegerVector::create(0, 0, 0, 0, 1, -1, 1, -1, 0, -1, 1, 0);
  int k = min(IntegerVector::create(maximum[5],maximum[7],maximum[9]));
  if(k > 0){
    IntegerVector move = k*inv_1;
    maximum = move + maximum;
  }
  return maximum;
}


/** 
 * 
 * Set the minimum dictionary ordering of a vector keeping n_pp, 
 * n_mm, n_rp and n_rm fixed in the dictionary ordering of this method.
 * 
 */ 

IntegerVector set_min_s1_4by3(IntegerVector vec){
  IntegerVector minimum = clone(vec) ;
  IntegerVector inv_1 = IntegerVector::create(0, 0, 0, 0, -1, 1, -1, 1, 0, 1, -1, 0);
  int k = min(IntegerVector::create(minimum[4],minimum[6],minimum[10]));
  if(k > 0){
    IntegerVector move = k*inv_1;
    minimum = move + minimum;
  }
  return minimum;
}


/**
 * 
 * This function makes all moves and takes the one with the greatest Dvalue
 * (Note: log(D-value) >= 0).
 * 
 */

IntegerVector find_move_4by3(IntegerMatrix moves, IntegerVector current_state, double num){
  double max_prob = -1*INFINITY;
  IntegerVector maximum;
  for(int k = 0 ; k < moves.nrow() ; k++){
    IntegerVector next_move = current_state + moves.row(k);
    if (is_false(any(next_move < 0))){
      if (abs(next_move[0] - current_state[0]) <= 1){
        double temp_prob = get_Dvalue_4by3(next_move, num);
        if (temp_prob > max_prob){
          maximum = next_move;
          max_prob = temp_prob;
        }
      } 
      // Handle n_rp moves which change n_rp by more than 1.
      else if (maximum.size() == 0 && abs(next_move[0] - current_state[0]) == 2){
        double temp_prob = get_Dvalue_4by3(next_move, num);
        if (temp_prob > max_prob){
          maximum = next_move;
          max_prob = temp_prob;
        }
      }
    } 
  }
  return maximum;
}



/**
 * 
 * This function finds the matrix of the next bigger score with the highest D-value.
 * 
 */

IntegerVector get_next_max_score_4by3(IntegerMatrix moves, IntegerVector vec, double num){
  IntegerVector maximum;
  IntegerMatrix moves2 = moves(Range(0,13), Range(0,11));
  maximum = find_move_4by3(moves2, vec, num);
  if (maximum.size() == 0){
    moves2 = moves(Range(14,21), Range(0,11));
    maximum = find_move_4by3(moves2, vec, num);
  }
  if (maximum.size() == 0){
    moves2 = moves(Range(22,22), Range(0,11));
    maximum = find_move_4by3(moves2, vec, num);
  }
  return maximum;
}



/**
 * 
 * This function finds the matrix of the next smaller score with the highest D-value.
 * 
 */
IntegerVector get_next_max_score2_4by3(IntegerMatrix moves, IntegerVector vec, double num){
  IntegerVector maximum;
  IntegerMatrix moves2 = moves(Range(0,22), Range(0,11));
  maximum = find_move_4by3(moves2, vec, num);
  if (maximum.size() == 0){
    moves2 = moves(Range(23,31), Range(0,11));
    maximum = find_move_4by3(moves2, vec, num);
  }
  if (maximum.size() == 0){
    moves2 = moves(Range(32,32), Range(0,11));
    maximum = find_move_4by3(moves2, vec, num);
  }
  return maximum;
}



/**
 * 
 * Gets a matrix with the highest Dvalue (i.e highest probability) of a given score or the first
 * score greater than the score_target.
 * 
 */

IntegerVector get_mat_of_score_4by3(IntegerVector constraints, int score_target, double num){
  IntegerMatrix moves = get_score_increase_moves_4by3();
  IntegerVector start = min_arr_score_4by3(constraints);
  int score = get_score_4by3(start);
  while(true){
    if (score >= score_target){
      break;
    }
    
    // Start is never supposed to be empty.
    if (start.size() == 0){
      Rcpp::stop("Something went wrong! Report error to package maintainer.\n");
    }
    
    start = get_next_max_score_4by3(moves, start, num);
    score = get_score_4by3(start);
  }
  return start;
}


/**
 * 
 * This function computes a very accurate estimate of the 
 * largest D-value.
 *
 */

NumericVector max_Dvalue_arr_4by3(IntegerVector constraints){
  long q_p = constraints[0];
  long q_m = constraints[1];
  long q_r = constraints[2];
  long q_z = constraints[3];
  long n_p = constraints[4];
  long n_m = constraints[5];
  long n_z = constraints[6];
  double tot = q_p + q_m + q_z + q_r;
  double n_pp = (q_p*n_p)/tot;
  double n_mm = (q_m*n_m)/tot;
  double n_rp = (q_r*n_p)/tot;
  double n_rm = (q_r*n_m)/tot;
  double n_pm = (q_p*n_m)/tot;
  double n_mp = (q_m*n_p)/tot;
  double n_zm = (q_z*n_m)/tot;
  double n_zp = (q_z*n_p)/tot;
  double n_zz = (q_z*n_z)/tot;
  double n_rz = (q_r*n_z)/tot;
  double n_pz = (q_p*n_z)/tot;
  double n_mz = (q_m*n_z)/tot;
  NumericVector max_mat = NumericVector::create(n_rp, n_rm, n_pp, n_mm, n_mp, n_pm, n_pz, n_mz, n_rz, n_zp, n_zm, n_zz);
  return max_mat;
}
