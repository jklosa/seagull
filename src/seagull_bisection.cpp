#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

inline static double sqrt_double(double x) { return ::sqrt(x); }

using namespace Rcpp;
using namespace arma;

//' Internal bisection algorithm
//' 
//' @description This algorithm finds the smallest positive root of a polynomial
//' of second degree in \eqn{\lambda}. Bisection is an implicit algorithm, i.e.,
//' it calls itself until a certain precision is reached.
//' 
//' @name seagull_bisection
//' 
//' @aliases seagull_bisection
//' 
//' @param ROWS the length of the input vectors.
//' 
//' @param ALPHA mixing parameter of the penalty terms. Satisfies: \eqn{0 <
//' \alpha < 1}. The penalty term looks as follows: \deqn{\alpha *
//' "lasso penalty" + (1-\alpha) * "group lasso penalty".}
//' 
//' @param LEFT_BORDER value of the left border of the current interval that
//' for sure harbors a root.
//' 
//' @param RIGHT_BORDER value of the right border of the current interval that
//' for sure harbors a root.
//' 
//' @param GROUP_WEIGHT a multiplicative scalar which is part of the polynomial.
//' 
//' @param VECTOR_WEIGHTS an input vector of multiplicative scalars which are
//' part of the polynomial. This vector is a subset of the vector of weights for
//' features.
//' 
//' @param VECTOR_IN another input vector which is required to compute the value
//' of the polynomial.
//' 
//' @details The polynomial has the following form:
//' \deqn{\sum_j (|vector_j| - \alpha weight_j \lambda )^2_+ - (1 - \alpha)^2
//' weight^2 \lambda^2.} The polynomial is non-trivial, because summands are
//' part of the sum if and only if the terms are non-negative.
//' 
//' @return If a certain precision (\code{TOLERANCE}) is reached, this algorithm
//' returns the center point of the current interval, in which the root is
//' located. Otherwise, the function calls itself using half of the initial
//' interval, in which the root is surely located.
//' 
// [[Rcpp::export]]
double seagull_bisection(
  int ROWS,
  double ALPHA,
  double LEFT_BORDER,
  double RIGHT_BORDER,
  double GROUP_WEIGHT,
  arma::colvec VECTOR_WEIGHTS,
  arma::colvec &VECTOR_IN) {
  
  int index_i       = 0;
  double TOLERANCE  = 0.0000000000001;
  double MID_POINT  = 0.5 * (LEFT_BORDER + RIGHT_BORDER);
  double FUNC_LEFT  = 0.0;
  double FUNC_MID   = 0.0;
  double FUNC_RIGHT = 0.0;
  double TEMP_LEFT  = 0.0;
  double TEMP_MID   = 0.0;
  double TEMP_RIGHT = 0.0;
  
  
  /*********************************************************
   **     Calculate the value of the function at the      **
   **     left border, the mid point, and the right       **
   **     border of the interval:                         **
   *********************************************************/
  for (index_i = 0; index_i < ROWS; index_i++) {
    if (VECTOR_IN(index_i) < 0.0) {
      TEMP_LEFT  = -1.0 * VECTOR_IN(index_i) - ALPHA * LEFT_BORDER * VECTOR_WEIGHTS(index_i);
      TEMP_MID   = -1.0 * VECTOR_IN(index_i) - ALPHA * MID_POINT * VECTOR_WEIGHTS(index_i);
      TEMP_RIGHT = -1.0 * VECTOR_IN(index_i) - ALPHA * RIGHT_BORDER * VECTOR_WEIGHTS(index_i);
    } else {
      TEMP_LEFT  = VECTOR_IN(index_i) - ALPHA * LEFT_BORDER * VECTOR_WEIGHTS(index_i);
      TEMP_MID   = VECTOR_IN(index_i) - ALPHA * MID_POINT * VECTOR_WEIGHTS(index_i);
      TEMP_RIGHT = VECTOR_IN(index_i) - ALPHA * RIGHT_BORDER * VECTOR_WEIGHTS(index_i);
    }
    
    if (TEMP_LEFT > 0.0) {
      FUNC_LEFT = FUNC_LEFT + TEMP_LEFT * TEMP_LEFT;
    }
    
    if (TEMP_MID > 0.0) {
      FUNC_MID = FUNC_MID + TEMP_MID * TEMP_MID;
    }
    
    if (TEMP_RIGHT > 0.0) {
      FUNC_RIGHT = FUNC_RIGHT + TEMP_RIGHT * TEMP_RIGHT;
    }
  }
  FUNC_LEFT  = FUNC_LEFT - (1.0 - ALPHA) * (1.0 - ALPHA) * LEFT_BORDER * LEFT_BORDER * GROUP_WEIGHT;
  FUNC_MID   = FUNC_MID - (1.0 - ALPHA) * (1.0 - ALPHA) * MID_POINT * MID_POINT * GROUP_WEIGHT;
  FUNC_RIGHT = FUNC_RIGHT - (1.0 - ALPHA) * (1.0 - ALPHA) * RIGHT_BORDER * RIGHT_BORDER * GROUP_WEIGHT;
  
  
  /*********************************************************
   **     Check for change of sign within sub-intervals   **
   **     and redo bisection:                             **
   *********************************************************/
  if (FUNC_LEFT * FUNC_MID < 0.0) {
    if (std::abs (LEFT_BORDER - MID_POINT) > TOLERANCE) {
      return seagull_bisection(ROWS, ALPHA, LEFT_BORDER, MID_POINT, GROUP_WEIGHT, VECTOR_WEIGHTS, VECTOR_IN);
    } else {
      return MID_POINT;
    }
  } else if (FUNC_MID * FUNC_RIGHT < 0.0) {
    if (std::abs (MID_POINT - RIGHT_BORDER) > TOLERANCE) {
      return seagull_bisection(ROWS, ALPHA, MID_POINT, RIGHT_BORDER, GROUP_WEIGHT, VECTOR_WEIGHTS, VECTOR_IN);
    } else {
      return MID_POINT;
    }
  } else {
    return MID_POINT;
  }
}
