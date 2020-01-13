#ifndef BISECTION_SEAGULL_H
#define BISECTION_SEAGULL_H

#include "RcppArmadillo.h"
double seagull_bisection(
  int ROWS,
  double ALPHA,
  double LEFT_BORDER,
  double RIGHT_BORDER,
  double GROUP_WEIGHT,
  arma::colvec VECTOR_WEIGHTS,
  arma::colvec &VECTOR_IN);

#endif