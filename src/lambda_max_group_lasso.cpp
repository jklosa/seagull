#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

inline static double sqrt_double(double x) { return ::sqrt(x); }

using namespace Rcpp;
using namespace arma;

//' Maximal \eqn{\lambda}
//' 
//' @name lambda_max
//' 
//' @aliases lambda_max_group_lasso
//' 
//' @param VECTOR_Y numeric vector of observations.
//' 
//' @param VECTOR_GROUPS integer vector specifying which effect (fixed and
//' random) belongs to which group.
//' 
//' @param VECTOR_WEIGHTS_FEATURES numeric vector of weights for the vectors of
//' fixed and random effects \eqn{[b^T, u^T]^T}. The entries may be permuted
//' corresponding to their group assignments.
//' 
//' @param VECTOR_BETA numeric vector of features. At the end of this function,
//' the random effects are initialized with zero, but the fixed effects are
//' initialized via a least squares procedure.
//' 
//' @param MATRIX_X numeric design matrix relating y to fixed and random
//' effects \eqn{[X Z]}.
//' 
//' @export
// [[Rcpp::export]]
double lambda_max_group_lasso(
  arma::colvec &VECTOR_Y,
  IntegerVector VECTOR_GROUPS,
  arma::colvec &VECTOR_WEIGHTS_FEATURES,
  arma::colvec &VECTOR_BETA,
  arma::mat &MATRIX_X) {
  
  int n                    = MATRIX_X.n_rows;
  int p                    = MATRIX_X.n_cols;
  int index_i              = 0;
  int index_j              = 0;
  int COUNTER              = 0;
  int COUNTER_GROUP_SIZE   = 0;
  int NUMBER_ZEROS_WEIGHTS = 0;
  int NUMBER_GROUPS        = max(VECTOR_GROUPS);
  double LAMBDA_MAX        = 0.0;
  double TEMP              = 0.0;
  
  //Determine the number of weights equal to zero:
  for (index_j = 0; index_j < p; index_j++) {
    if (VECTOR_WEIGHTS_FEATURES(index_j) == 0.0) {
      NUMBER_ZEROS_WEIGHTS = NUMBER_ZEROS_WEIGHTS + 1;
    }
  }
  
  IntegerVector VECTOR_INDEX_START (NUMBER_GROUPS);
  IntegerVector VECTOR_INDEX_END (NUMBER_GROUPS);
  NumericVector VECTOR_X_TRANSP_RESIDUAL_ACTIVEc (p);
  NumericVector VECTOR_WEIGHTS_GROUPSc (NUMBER_GROUPS);
  NumericVector VECTOR_L2_NORM_GROUPSc (NUMBER_GROUPS);
  
  colvec VECTOR_X_TRANSP_RESIDUAL_ACTIVE(VECTOR_X_TRANSP_RESIDUAL_ACTIVEc.begin(), p, false);
  colvec VECTOR_WEIGHTS_GROUPS(VECTOR_WEIGHTS_GROUPSc.begin(), NUMBER_GROUPS, false);
  colvec VECTOR_L2_NORM_GROUPS(VECTOR_L2_NORM_GROUPSc.begin(), NUMBER_GROUPS, false);
  
  
  /*********************************************************
   **     Create a vector of group sizes, a vector of     **
   **     start indices, and a vector of end indices      **
   **     from the vector of groups. And also create a    **
   **     vector of group weights as group means from     **
   **     the vector of feature weights:                  **
   *********************************************************/
  COUNTER = 1;
  for (index_i = 0; index_i < NUMBER_GROUPS; index_i++) {
    COUNTER_GROUP_SIZE = 0;
    for (index_j = 0; index_j < p; index_j++) {
      if (VECTOR_GROUPS(index_j) == (index_i + 1)) {
        COUNTER_GROUP_SIZE = COUNTER_GROUP_SIZE + 1;
        VECTOR_INDEX_START(index_i)    = index_j - COUNTER_GROUP_SIZE + 1;
        VECTOR_INDEX_END(index_i)      = index_j;
        VECTOR_WEIGHTS_GROUPS(index_i) = VECTOR_WEIGHTS_GROUPS(index_i) + VECTOR_WEIGHTS_FEATURES(index_j);
      }
    }
    VECTOR_WEIGHTS_GROUPS(index_i) = sqrt_double(VECTOR_WEIGHTS_GROUPS(index_i));
  }
  
  
  /*********************************************************
   **     Treatment, if unpenalized features are          **
   **     involved:                                       **
   *********************************************************/
  if (NUMBER_ZEROS_WEIGHTS > 0) {
    NumericVector VECTOR_RESIDUAL_ACTIVEc (n);
    NumericVector VECTOR_BETA_ACTIVEc (NUMBER_ZEROS_WEIGHTS);
    NumericMatrix MATRIX_X_ACTIVEc (n, NUMBER_ZEROS_WEIGHTS);
    
    colvec VECTOR_RESIDUAL_ACTIVE(VECTOR_RESIDUAL_ACTIVEc.begin(), n, false);
    colvec VECTOR_BETA_ACTIVE(VECTOR_BETA_ACTIVEc.begin(), NUMBER_ZEROS_WEIGHTS, false);
    mat MATRIX_X_ACTIVE(MATRIX_X_ACTIVEc.begin(), n, NUMBER_ZEROS_WEIGHTS, false);
    
    
    /*******************************************************
     **     Calculations with "active" set:               **
     *******************************************************/
    //Determine the "active" set and create X_A = X_active:
    COUNTER = 0;
    for (index_j = 0; index_j < p; index_j++) {
      if (VECTOR_WEIGHTS_FEATURES(index_j) == 0.0) {
        for (index_i = 0; index_i < n; index_i++) {
          MATRIX_X_ACTIVE(index_i, COUNTER) = MATRIX_X(index_i, index_j);
        }
        COUNTER = COUNTER + 1;
      }
    }
    
    //Solve for beta_A in y = X_A * beta_A:
    VECTOR_BETA_ACTIVE = solve(MATRIX_X_ACTIVE, VECTOR_Y);
    
    //Create beta with beta_A:
    COUNTER = 0;
    for (index_j = 0; index_j < p; index_j++) {
      if (VECTOR_WEIGHTS_FEATURES(index_j) == 0.0) {
        VECTOR_BETA(index_j) = VECTOR_BETA_ACTIVE(COUNTER);
        COUNTER = COUNTER + 1;
      }
    }
    
    //Calculate res_A = y - X_A*beta_A:
    VECTOR_RESIDUAL_ACTIVE = VECTOR_Y - (MATRIX_X_ACTIVE * VECTOR_BETA_ACTIVE);
    
    //Calculate t(X)*res_A:
    VECTOR_X_TRANSP_RESIDUAL_ACTIVE = MATRIX_X.t() * VECTOR_RESIDUAL_ACTIVE;
    
    
  /*********************************************************
   **     Treatment, if only penalized features are       **
   **     involved:                                       **
   *********************************************************/
  } else {
    //Calculate t(X)*y:
    VECTOR_X_TRANSP_RESIDUAL_ACTIVE = MATRIX_X.t() * VECTOR_Y;
  }
  
  //Scale t(X)*res_A in groups with n*weight*sqrt(group_size), if weight>0:
  // Calculate l2-norm in groups:
  for (index_i = 0; index_i < NUMBER_GROUPS; index_i++) {
    if (VECTOR_WEIGHTS_GROUPS(index_i) == 0.0) {
      VECTOR_L2_NORM_GROUPS(index_i) = 0.0;
    } else {
      TEMP = static_cast<double>(n) * VECTOR_WEIGHTS_GROUPS(index_i);
      for (index_j = VECTOR_INDEX_START(index_i); index_j <= VECTOR_INDEX_END(index_i); index_j++) {
        VECTOR_L2_NORM_GROUPS(index_i) = VECTOR_L2_NORM_GROUPS(index_i) + VECTOR_X_TRANSP_RESIDUAL_ACTIVE(index_j) * 
          VECTOR_X_TRANSP_RESIDUAL_ACTIVE(index_j);
      }
      VECTOR_L2_NORM_GROUPS(index_i) = sqrt_double(VECTOR_L2_NORM_GROUPS(index_i)) / TEMP;
    }
  }
  
  //Determine lambda_max and perform numeric correction:
  for (index_i = 0; index_i < NUMBER_GROUPS; index_i++) {
    if (VECTOR_L2_NORM_GROUPS(index_i) < 0.0) {
      VECTOR_L2_NORM_GROUPS(index_i) = -1.0 * VECTOR_L2_NORM_GROUPS(index_i);
    }
  }
  LAMBDA_MAX = max(VECTOR_L2_NORM_GROUPS);
  return LAMBDA_MAX*1.00001;
}
