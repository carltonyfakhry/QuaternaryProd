#include <Rcpp.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R.h>
#include <math.h>
#include "QuaternaryProd.h"
#include "QuaternaryProdUtils.h"
using namespace Rcpp;



// Authors: Carl Tony Fakhry, Kourosh Zarringhalam, Ping Chen.

 

/**
 * 
 * Get the probability of a given score.
 * 
 */
 
double probability_of_score_4by3(IntegerVector arr, double num, double total, double epsilon){
  IntegerVector start(clone(arr));
  IntegerMatrix dict_dec = get_dict_decrease_moves_4by3();
  IntegerMatrix dict_inc = get_dict_increase_moves_4by3();
  IntegerMatrix first_col = dict_dec(Range(14,24), Range(0,11));
  IntegerMatrix inv_first_col = dict_inc(Range(22,37), Range(0,11)) ;
  IntegerMatrix second_col = dict_dec(Range(8,13), Range(0,11));
  IntegerMatrix inv_second_col = dict_inc(Range(15,21), Range(0,11)) ;
  IntegerMatrix third_col = dict_dec(Range(3,7), Range(0,11));
  IntegerMatrix inv_third_col = dict_inc(Range(5,14), Range(0,11)) ;
  IntegerMatrix fourth_col = dict_dec(Range(1,2), Range(0,11));
  IntegerMatrix inv_fourth_col = dict_inc(Range(1,4), Range(0,11));
  IntegerVector s_1 = dict_dec.row(0);
  IntegerVector is_1 = dict_inc.row(0);
  IntegerVector constraints = get_constraints_4by3(start) ;
  IntegerVector move_first_up, move_first_down; 
  IntegerVector move_second_up, move_second_down;
  IntegerVector move_third_up, move_third_down;
  IntegerVector move_fourth_up, move_fourth_down;
  int direction4 = 0;
  int direction3 = 0;
  int direction2 = 0; 
  int direction1 = 0;
  bool first_up = false ;
  bool second_up = false ;
  bool third_up = false;
  bool fourth_up = false;
  bool first_down = false ;
  bool second_down = false ;
  bool third_down = false;
  bool fourth_down = false;
  double probability = 0.0; 
  double dval1 = 0.0;
  double dval2 = 0.0;
  double dval3 = 0.0;
  double dval4 = 0.0;
  IntegerVector start2, start3, start4, min_s1, max_s1;
  while (true){
    
    Rcpp::checkUserInterrupt();
    
    dval1 = get_Dvalue_4by3(start, num);
    start2 = clone(start);
    
    // All given attempts to find a move, will find a move to a 
    // matrix with highest probability in its given class while
    // either increasing or decreasing dictionary ordering of
    // n_pp and n_mm.
    
    // Find a move that increases the dictionary ordering of n_rp.
    if (!first_up && direction1 >= 0){
      move_first_up = find_move_4by3(inv_first_col, start, num);
      if (move_first_up.size() > 0){
        first_up = true;
      }
    }
    
    // Find a move that decreases the dictionary ordering of n_rp.
    if (!first_down && direction1 <= 0){
      move_first_down = find_move_4by3(first_col, start, num);
      if (move_first_down.size() > 0){
        first_down = true;
      }
    }
    
    while (true){
    
      dval2 = get_Dvalue_4by3(start2, num);
      start3 = clone(start2);
      
      // Find a move that increases the dictionary ordering of n_rm.
      if (!second_up && direction2 >= 0){
        move_second_up = find_move_4by3(inv_second_col, start2, num);
        if (move_second_up.size() > 0){
          second_up = true;
        }
      }
      
      // Find a move that decreases the dictionary ordering of n_rm.
      if (!second_down && direction2 <= 0){
        move_second_down = find_move_4by3(second_col, start2, num);
        if (move_second_down.size() > 0){
          second_down = true;
        }
      }
      
      while (true){
        
        dval3 = get_Dvalue_4by3(start3, num);
        start4 = clone(start3);
        
        // Find a move that increases the dictionary ordering of n_pp.
        if (!third_up && direction3 >= 0){
          move_third_up = find_move_4by3(inv_third_col, start3, num);
          if (move_third_up.size() > 0){
            third_up = true;
          }
        }
        
        // Find a move that decreases the dictionary ordering of n_pp.
        if (!third_down && direction3 <= 0){
          move_third_down = find_move_4by3(third_col, start3, num);
          if (move_third_down.size() > 0){
            third_down = true;
          }
        }
        
        while (true){
          
          // Find a move that increases the dictionary ordering of n_mm.
          if (!fourth_up && direction4 >= 0){
            move_fourth_up = find_move_4by3(inv_fourth_col, start4, num);
            if (move_fourth_up.size() > 0){
              fourth_up = true;
            }
          }
          
          // Find a move that decreases the dictionary ordering of n_mm.
          if (!fourth_down && direction4 <= 0){
            move_fourth_down = find_move_4by3(fourth_col, start4, num);
            if (move_fourth_down.size() > 0){
              fourth_down = true;
            }
          }
                
          dval4 = get_Dvalue_4by3(start4, num);
          probability += exp(dval4 - total);
                
          // Get the matrices with min and max dictionary ordering
          // for a given fixed n_rp, n_rm, n_pp and n_mm.
          min_s1 = set_min_s1_4by3(start4);
          max_s1 = set_max_s1_4by3(start4);
                                            
          // Upper and lower bounds on n_mp.
          int u_nmp = max_s1[4];
          int l_nmp = min_s1[4];
          
          // This is the matrix with max probability (i.e max D-value)
          // for a set of tables with fixed n_pp and n_mm.
          IntegerVector max_npp_nmm = clone(start4);
            
          int n_mp = max_npp_nmm[4] + 1;  
          double prev_dval = dval4;
          double i = 1.0;
          IntegerVector prev_vec = clone(start4);

          // Add up the probabilities of all tables with higher dictionary
          // ordering than the matrix with highest probability for a set of
          // tables with fixed n_rp, n_rm, n_pp and n_mm.
          while (n_mp <= u_nmp){
            prev_dval = prev_dval + log((((prev_vec[5]-i+1)*(prev_vec[7]-i+1)*(prev_vec[9]-i+1))/((prev_vec[4]+i)*(prev_vec[6]+i)*(prev_vec[10]+i))));
            probability += exp(prev_dval - total);
            n_mp += 1;
            i += 1;
            if (prev_dval < epsilon){
              break;
            }
          }
          
          n_mp = max_npp_nmm[4] - 1;  
          prev_dval = dval4;
          i = 1.0;
          
          // Add up the probabilities of all tables with lower dictionary
          // ordering than the matrix with highest probability for a set of
          // tables with fixed n_pp, n_mm, n_rp and n_rm.
          while (n_mp >= l_nmp){
            prev_dval = prev_dval + log((((prev_vec[4]-i+1)*(prev_vec[6]-i+1)*(prev_vec[10]-i+1))/((prev_vec[5]+i)*(prev_vec[7]+i)*(prev_vec[9]+i))));
            probability += exp(prev_dval - total);
            n_mp -= 1;
            i += 1;
            if (prev_dval < epsilon){
              break;
            }
          }
                
          // If fourth inner while loop starting at matrix with
          // D-value lower than epsilon.
          if (dval4 < epsilon && direction4 == 0){
            fourth_down = false;
            fourth_up = false;
            break;
          }
          
          // If increasing n_mm, but new matrix of highest probability has 
          // D-value smaller than epsilon, and there exists a move to 
          // decrease orginal n_mm (i.e n_rm where the fourth inner for loop 
          // started).
          else if (dval4 < epsilon && direction4 > 0 && fourth_down){
            direction4 = -1;
            start4 = move_fourth_down;
            fourth_down = false;
            fourth_up = false;
          } 
          
          // If increasing n_mm, but new matrix has D-value less than 
          // epsilon.
          else if (dval4 < epsilon && direction4 > 0) {
            direction4 = 0;
            fourth_down = false;
            fourth_up = false;
            break;
          }
          
          // If decreasing n_mm, but new matrix has D-value less than 
          // epsilon.
          else if (dval4 < epsilon && direction4 < 0){
            direction4 = 0;
            fourth_down = false;
            fourth_up = false;
            break;
          }
          
          // If increasing n_mm.
          else if (fourth_up){
            start4 = move_fourth_up;
            direction4 = 1;
            fourth_up = false;
          }
          
          // If decreasing n_mm.
          else if (fourth_down){
            start4 = move_fourth_down;
            direction4 = -1;
            fourth_down = false;
          }
                                       
          // If neither decreasing or increasing n_mm.
          else{
            direction4 = 0;
            fourth_down = false;
            fourth_up = false;
            break;
          }     
        }
        
        // If third inner while loop starting at matrix with
        // D-value lower than epsilon.
        if (dval3 < epsilon && direction3 == 0){
          third_down = false;
          third_up = false;
          break;
        }
        
        // If increasing n_pp, but new matrix of highest probability has 
        // D-value smaller than epsilon, and there exists a move to 
        // decrease orginal n_pp (i.e n_rp where the third for loop 
        // started).
        else if (dval3 < epsilon && direction3 > 0 && third_down){
          direction3 = -1;
          start3 = move_third_down;
          third_down = false;
          third_up = false;
        } 
        
        // If increasing n_pp, but new matrix has D-value less than 
        // epsilon.
        else if (dval3 < epsilon && direction3 > 0) {
          direction3 = 0;
          third_down = false;
          third_up = false;
          break;
        }
        
        // If decreasing n_pp, but new matrix has D-value less than 
        // epsilon.
        else if (dval3 < epsilon && direction3 < 0){
          direction3 = 0;
          third_down = false;
          third_up = false;
          break;
        }
        
        // If increasing n_pp.
        else if (third_up){
          start3 = move_third_up;
          direction3 = 1;
          third_up = false;
        }
        
        // If decreasing n_pp.
        else if (third_down){
          start3 = move_third_down;
          direction3 = -1;
          third_down = false;
        }
        
        // If no direction left to increase or decrease n_pp.
        else{
          direction3 = 0;
          third_down = false;
          third_up = false;
          break;
        }     
      }

      // If second inner while loop starting at matrix with
      // D-value lower than epsilon.
      if (dval2 < epsilon && direction2 == 0){
        second_down = false;
        second_up = false;
        break;
      }
      
      // If increasing n_rm, but new matrix of highest probability has 
      // D-value smaller than epsilon, and there exists a move to 
      // decrease orginal n_rm (i.e n_mm where the second for loop 
      // started).
      else if (dval2 < epsilon && direction2 > 0 && second_down){
        direction2 = -1;
        start2 = move_second_down;
        second_down = false;
        second_up = false;
      } 
      
      // If increasing n_rm, but new matrix has D-value less than 
      // epsilon.
      else if (dval2 < epsilon && direction2 > 0) {
        direction2 = 0;
        second_down = false;
        second_up = false;
        break;
      }
      
      // If decreasing n_rm, but new matrix has D-value less than 
      // epsilon.
      else if (dval2 < epsilon && direction2 < 0){
        direction2 = 0;
        second_down = false;
        second_up = false;
        break;
      }
      
      // If increasing n_rm.
      else if (second_up){
        start2 = move_second_up;
        direction2 = 1;
        second_up = false;
      }
      
      // If decreasing n_rm.
      else if (second_down){
        start2 = move_second_down;
        direction2 = -1;
        second_down = false;
      }
      
      // If no direction left to increase or decrease n_rm.
      else{
        direction2 = 0;
        second_down = false;
        second_up = false;
        break;
      }     
    }    
    
    // If first while loop starting at matrix with
    // D-value lower than epsilon.
    if (dval1 < epsilon && direction1 == 0){
      first_down = false;
      first_up = false;
      break;
    }
    
    // If increasing n_rp, but new matrix of highest probability has 
    // D-value smaller than epsilon, and there exists a move to 
    // decrease orginal n_rp (i.e n_pp where the firstfor loop 
    // started).
    else if (dval1 < epsilon && direction1 > 0 && first_down){
      direction1 = -1;
      start = move_first_down;
      first_down = false;
      first_up = false;
    } 
    
    // If increasing n_rp, but new matrix has D-value less than 
    // epsilon.
    else if (dval1 < epsilon && direction1 > 0) {
      direction1 = 0;
      first_down = false;
      first_up = false;
      break;
    }
    
    // If decreasing n_rp, but new matrix has D-value less than 
    // epsilon.
    else if (dval1 < epsilon && direction1 < 0){
      direction1 = 0;
      first_down = false;
      first_up = false;
      break;
    }
    
    // If increasing n_rp.
    else if (first_up){
      start = move_first_up;
      direction1 = 1;
      first_up = false;
    }
    
    // If decreasing n_rp.
    else if (first_down){
      start = move_first_down;
      direction1 = -1;
      first_down = false;
    }
    
    // If no direction left to increase or decrease n_rp.
    else{
      first_up = false;
      first_down = false;
      break;
    }   
  }
  
  return probability;
}



/**
 * 
 * This function computes the right sided Pvalue of a given score of the
 * Quaternary Product Scoring Statistic.
 * 
 */

// [[Rcpp::export]]
double Pvalue_4by3(int q_p, int q_m, int q_r, int q_z, int n_p, int n_m, int n_z, int target_score, double epsilon){
  
  IntegerVector constraints = IntegerVector::create(q_p, q_m, q_r, q_z, n_p, n_m, n_z);
  
  // Get matrices with min and max scores.
  IntegerVector max_arr = max_arr_score_4by3(constraints);
  IntegerVector min_arr = min_arr_score_4by3(constraints);
  
  // Get min and max scores.
  int max_score =  get_score_4by3(max_arr);
  int min_score = get_score_4by3(min_arr);
  
  // Get denominator (i.e total) for probability.
  double total = get_total_4by3(constraints);
  
  // Get numerator for D-value.
  double num = get_numerator_4by3(constraints);
  
  // Check if score is outside boundary of support.
  if (target_score > max_score || target_score < min_score){
    Rcpp::stop("Score outside boundaries of support!\n");
  }
  
  IntegerVector start;
  IntegerMatrix moves;
  
  // Get the matrix with highest D-value.
  NumericVector max_dvalue_arr = max_Dvalue_arr_4by3(constraints);
  
  // Center is an pproximate of the score with highest probability.
  int center = max_dvalue_arr[0] + max_dvalue_arr[1] + max_dvalue_arr[2] + max_dvalue_arr[3] - max_dvalue_arr[4] - max_dvalue_arr[5];
  
  if (center >= min_score && center <= max_score){
    center = floor(center);
  }
  else if (center > max_score){
    center = max_score;
  }
  else{
    center = min_score;
  }
  
  int score;
  int stop_score;
  double probability = 0.0 ;
  double temp_prob;
  bool complement = false;
  
  // Get log of Dvalue.
  double max_Dval = num - sum(lfactorial(max_dvalue_arr));

  // Update log scaled epsilon, epsilon = D-value*epsilon.
  epsilon = max_Dval + log(epsilon);
  
  // In case there is only onle element in the domain,
  // then that element has probability 1.
  if (min_score == max_score){
    return 1;
  }
  
  // Start from min score if target score is less than center.
  if (target_score < center){
    complement = true;
    score = min_score;
    stop_score = target_score;
    start = min_arr;
    moves = get_score_increase_moves_4by3();
    while(score < stop_score){  
      temp_prob = probability_of_score_4by3(start, num, total, epsilon);
      probability += temp_prob;  
      start = get_next_max_score_4by3(moves, start, num);
      score = get_score_4by3(start);
    }
  }
  // Else start from target score if score is greater or equal to center.
  else {  
    stop_score = max_score;
    moves = get_score_increase_moves_4by3();
    // Get table of target_score or first score right after it
    start = get_mat_of_score_4by3(constraints, target_score, num);
    while(true){  
      temp_prob = probability_of_score_4by3(start, num, total, epsilon);
      probability += temp_prob; 
      score = get_score_4by3(start);
      if (score == stop_score){
        break;
      }
      start = get_next_max_score_4by3(moves, start, num);
    }
  }
  
  if (complement){
    
    // Handle small floating point arithmetic error.
    // If probability slightly smaller than 0, then
    // compute the right tail instead of just taking 
    // the compliment.
    if (1 - probability < 0){
      probability = 0.0;
      stop_score = max_score;
      moves = get_score_increase_moves_4by3();
      // Get table of target_score or first score right after it
      start = get_mat_of_score_4by3(constraints, target_score, num);
      while(true){  
        temp_prob = probability_of_score_4by3(start, num, total, epsilon);
        probability += temp_prob; 
        score = get_score_4by3(start);
        if (score == stop_score){
          break;
        }
        start = get_next_max_score_4by3(moves, start, num);
      }
      return probability;
    }
    return 1 - probability;
  }
  
  // Handle small floating point arithmetic error.
  if (probability > 1){
    return 1;
  }
  return probability;
}



/**
 * 
 * This function computes the right sided Pvalue of a statistically significant score of the
 * Quaternary Product Scoring Statistic otherwise if the score is not statistically significant,
 * a value of -1 is returned.
 * 
 */

// [[Rcpp::export]]
double SigPvalue_4by3(int q_p, int q_m, int q_r, int q_z, int n_p, int n_m, int n_z, int target_score, double epsilon, double significance_level){
 
  IntegerVector constraints = IntegerVector::create(q_p, q_m, q_r, q_z, n_p, n_m, n_z);
  
  // Get matrices with min and max scores.
  IntegerVector max_arr = max_arr_score_4by3(constraints);
  IntegerVector min_arr = min_arr_score_4by3(constraints);
  
  // Get min and max scores.
  int max_score =  get_score_4by3(max_arr);
  int min_score = get_score_4by3(min_arr);
  
  // Get denominator (i.e total) for probability.
  double total = get_total_4by3(constraints);
  
  // Get numerator for D-value.
  double num = get_numerator_4by3(constraints);
  
  // Check if score is outside boundary of support.
  if (target_score > max_score || target_score < min_score){
    Rcpp::stop("Score outside boundaries of support!\n");
  }
  
  IntegerVector start;
  IntegerMatrix moves;
  
  // Get the matrix with highest D-value.
  NumericVector max_dvalue_arr = max_Dvalue_arr_4by3(constraints);
  
  // Center is an pproximate of the score with highest probability.
  int center = max_dvalue_arr[0] + max_dvalue_arr[1] + max_dvalue_arr[2] + max_dvalue_arr[3] - max_dvalue_arr[4] - max_dvalue_arr[5];
  
  if (center >= min_score && center <= max_score){
    center = floor(center);
  }
  else if (center > max_score){
    center = max_score;
  }
  else{
    center = min_score;
  }
  
  int score;
  int stop_score;
  double probability = 0.0 ;
  double temp_prob;
  bool complement = false;
  
  // Get log of Dvalue.
  double max_Dval = num - sum(lfactorial(max_dvalue_arr));
  
  // Update log scaled epsilon, epsilon = D-value*epsilon.
  epsilon = max_Dval + log(epsilon);
  
  // In case there is only onle element in the domain,
  // then that element has probability 1.
  if (min_score == max_score){
    return -1;
  }

  // Start from min score if target score is less than center.
  if (target_score < center){
    complement = true;
    score = min_score;
    stop_score = target_score;
    start = min_arr;
    moves = get_score_increase_moves_4by3();
    while(score < stop_score){  
      temp_prob = probability_of_score_4by3(start, num, total, epsilon);
      probability += temp_prob;  
      start = get_next_max_score_4by3(moves, start, num);
      score = get_score_4by3(start);
    }
  }
  // Else start from target score if score is greater or equal to center.
  else {
    moves = get_score_decrease_moves_4by3();
    start = max_arr;
    score = max_score;
    while(true){  
      temp_prob = probability_of_score_4by3(start, num, total, epsilon);
      probability += temp_prob;  
      if (probability > significance_level){
        return -1;
      }
      start = get_next_max_score2_4by3(moves, start, num);
      score = get_score_4by3(start);
      if (score < target_score){
        break;
      }
      if (score == min_score){
        temp_prob = probability_of_score_4by3(start, num, total, epsilon);
        probability += temp_prob;  
        if (probability > significance_level){
          return -1;
        }
      }
    }
  } 
  
  if (complement){
    
    // Handle small floating point arithmetic error.
    // If probability slightly smaller than 0, then
    // compute the right tail instead of just taking 
    // the compliment.
    if (1 - probability < 0){
      probability = 0.0;
      stop_score = max_score;
      moves = get_score_increase_moves_4by3();
      start = get_mat_of_score_4by3(constraints, target_score, num);
      while(true){  
        temp_prob = probability_of_score_4by3(start, num, total, epsilon);
        probability += temp_prob; 
        score = get_score_4by3(start);
        if (score == stop_score){
          break;
        }
        start = get_next_max_score_4by3(moves, start, num);
      }
      if (probability > significance_level){
        return -1;
      } 
      else {
        return probability;
      }
    }
    
    if (1 - probability > significance_level){
      return -1;
    } 
    else {
      return 1 - probability;
    }
  }

  return probability;  
  
}



/**
 * 
 * Computes the support of the Quaternary Product Scoring Statistic.
 * 
 */

// [[Rcpp::export]]
IntegerVector domain_4by3(int q_p, int q_m, int q_r, int q_z, int n_p, int n_m, int n_z){

  IntegerVector constraints = IntegerVector::create(q_p, q_m, q_r, q_z, n_p, n_m, n_z);
  
  // Create an integer vector for the domain.
  IntegerVector domain = IntegerVector::create();
  
  // Get the score increase moves for the 4by3 tables.
  IntegerMatrix moves = get_score_increase_moves_4by3();
  
  // Get the table with maximum score.
  IntegerVector max_arr = max_arr_score_4by3(constraints);
  
  // Get the table with minimum score.
  IntegerVector start = min_arr_score_4by3(constraints);
  
  // Get the max and min scores.
  int max_score =  get_score_4by3(max_arr);
  int score = get_score_4by3(start);
  
  // Get numerator of the D-value.
  double num = get_numerator_4by3(constraints);
  
  // Enter score into domain.
  domain.push_back(score);
  
  // Special case if there is only one element in the support.
  if (score == max_score){
    return domain;
  }
  while(true){
    
    // Get table of max probability for next score.
    start = get_next_max_score_4by3(moves, start, num);

    // Get score of table with max probability.
    score = get_score_4by3(start);

    // Add score to domain vector.
    domain.push_back(score);
    
    // Start is never supposed to be empty.
    if (start.size() == 0){
      Rcpp::stop("Something went wrong! Report error to package maintainer.\n");
    }
    
    // Break out of loop if max score is reached.
    if (score == max_score){
      break;
    } 
    
    Rcpp::checkUserInterrupt();
    
  }
  
  return domain;
}



/**
 * 
 * Computes the pmf of the Quaternary Product Scoring Statistic.
 * 
 */

// [[Rcpp::export]]
NumericVector pmf_4by3(int q_p, int q_m, int q_r, int q_z, int n_p, int n_m, int n_z, double epsilon){

  IntegerVector constraints = IntegerVector::create(q_p, q_m, q_r, q_z, n_p, n_m, n_z);
  
  // Create an integer vector for the domain.
  IntegerVector domain = IntegerVector::create();
  
  // Create a vector to hold the probabilities.
  NumericVector probabilities = NumericVector::create();
  
  // Get the score increase moves for the 4by3 tables.
  IntegerMatrix moves = get_score_increase_moves_4by3();
  
  // Get the table with maximum score.
  IntegerVector max_arr = max_arr_score_4by3(constraints);
  
  // Get the table with minimum score.
  IntegerVector start = min_arr_score_4by3(constraints);
  
  // Get the max and min scores.
  int max_score =  get_score_4by3(max_arr);
  int score = get_score_4by3(start);
  
  // Get numerator of the D-value.
  double num = get_numerator_4by3(constraints);
  
  // Get the Denominator of the Quaternary Product probability.
  double total = get_total_4by3(constraints);
  
  // Get the matrix with highest D-value.
  NumericVector max_dvalue_arr = max_Dvalue_arr_4by3(constraints);
  
  // Get log of Dvalue.
  double max_Dval = num - sum(lfactorial(max_dvalue_arr));
  
  // Update log scaled epsilon, epsilon = D-value*epsilon.
  epsilon = max_Dval + log(epsilon);
  
  // Enter score into domain.
  domain.push_back(score);
  
  // Compute probability of min score.
  double temp_prob = probability_of_score_4by3(start, num, total, epsilon);
  probabilities.push_back(temp_prob);
  
  // Special case if there is only one element in the support.
  if (score == max_score){
    probabilities.names() = domain;
    return probabilities;
  }
  
  
  while(true){
    
    // Get table of max probability for next score.
    start = get_next_max_score_4by3(moves, start, num);
    
    // Get score of table with max probability.
    score = get_score_4by3(start);
    
    // Add score to domain vector.
    domain.push_back(score);
    
    // Get probability of score.
    temp_prob = probability_of_score_4by3(start, num, total, epsilon);
    
    // Add probability to probabilities vector.
    probabilities.push_back(temp_prob);
    
    // Start is never supposed to be empty.
    if (start.size() == 0){
      Rcpp::stop("Something went wrong! Report error to package maintainer.\n");
    }
    
    // Break out of loop if max score is reached.
    if (score == max_score){
      break;
    }
  }

  probabilities.names() = domain;
  return probabilities;
}


/**
 * 
 * This is the user's method for computing the probability of a
 * score.
 * 
 */

// [[Rcpp::export]]
double user_probability_of_score_4by3(int target_score, int q_p, int q_m, int q_r, int q_z, int n_p, int n_m, int n_z, double epsilon){
  
  IntegerVector constraints = IntegerVector::create(q_p, q_m, q_r, q_z, n_p, n_m, n_z);

  // Get numerator of the D-value.
  double num = get_numerator_4by3(constraints);
  
  // Get the Denominator of the Quaternary Product probability.
  double total = get_total_4by3(constraints);
  
  // Get matrices with min and max scores.
  IntegerVector max_arr = max_arr_score_4by3(constraints);
  IntegerVector min_arr = min_arr_score_4by3(constraints);
  
  // Get min and max scores.
  int max_score =  get_score_4by3(max_arr);
  int min_score = get_score_4by3(min_arr);
  
  // Check if score is outside boundary of support.
  if (target_score > max_score || target_score < min_score){
    Rcpp::stop("Score outside boundaries of support!\n");
  }
  
  // Get the matrix with highest D-value.
  NumericVector max_dvalue_arr = max_Dvalue_arr_4by3(constraints);
  
  // Get log of Dvalue.
  double max_Dval = num - sum(lfactorial(max_dvalue_arr));
  
  // Update log scaled epsilon, epsilon = D-value*epsilon.
  epsilon = max_Dval + log(epsilon);
  
  // Get table of max probability of a given score.
  IntegerVector mat_of_score = get_mat_of_score_4by3(constraints, target_score, num);

  // Make sure we found a matrix of a the target_score
  int score = get_score_4by3(mat_of_score);
  
  if (score != target_score){
    Rcpp::stop("The target score is not present in the support of the distribution therefore it has zero probability!");
  }
  
  // Compute probability of a given score.
  double probability = probability_of_score_4by3(mat_of_score, num, total, epsilon);
  

  return probability;
}
