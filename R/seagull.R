#' @title Mixed model fitting with lasso, group lasso, or sparse-group lasso
#' regularization
#' 
#' @description Fit a mixed model with lasso, group lasso, or sparse-group lasso
#' via proximal gradient descent. As this is an iterative algorithm, the step
#' size for each iteration is determined via backtracking line search. A grid
#' search for the regularization parameter \eqn{\lambda} is performed using warm
#' starts. Depending on the input parameter \code{alpha} this function
#' subsequently calls one of the three implemented
#' \code{\link[seagull]{lasso_variants}}. None of the input variables will be
#' centered or standardized throughout any of the calculations of this package.
#' 
#' @docType package
#' 
#' @author Jan Klosa
#' 
#' @import Rcpp
#' 
#' @importFrom Rcpp evalCpp
#' 
#' @useDynLib seagull
#' 
#' @name seagull
#' 
#' @keywords models regression
#' 
#' @usage seagull(y, X, Z, weights_u, groups, alpha, rel_acc, max_lambda, xi,
#'         loops_lambda, max_iter, gamma_bls, trace_progress)
#' 
#' @param y numeric vector of observations.
#' 
#' @param X (optional) numeric design matrix relating y to fixed effects b.
#' 
#' @param Z numeric design matrix relating y to random effects u.
#' 
#' @param weights_u (optional) numeric vector of weights for the vector of
#' random effects u. If no weights are passed by the user, all weights for u are
#' set to \code{1}.
#' 
#' @param groups (not required for the lasso) integer vector specifying which
#' effect belongs to which group. If the lasso is called, this vector is without
#' use and may remain empty. For group and sparse-group lasso, at least each
#' random effect variable needs to be assigned to one group. If in this case
#' also fixed effect variables are present in the model, but without assigned
#' group values, then the same value will be assigned to all fixed effect
#' variables automatically. Float or double input values will be truncated.
#' 
#' @param alpha mixing parameter for the sparse-group lasso. Has to satisfy:
#' \eqn{0 \le \alpha \le 1}. The penalty term looks as follows: \deqn{\alpha *
#' '"lasso penalty" + (1-\alpha) * "group lasso penalty".} If \code{alpha=1},
#' the lasso is called. If \code{alpha=0}, the group lasso is called. Default
#' value is \code{0.9}.
#' 
#' @param rel_acc (optional) relative accuracy of the solution to stop the
#' algorithm for the current value of \eqn{\lambda}. The algorithm stops after
#' iteration m, if: \deqn{||sol^{(m)} - sol^{(m-1)}||_\infty < rel_acc *
#' ||sol1{(m-1)}||_2.} Default value is \code{0.0001}.
#' 
#' @param max_lambda (optional) maximum value for the penalty parameter. Default
#' option is an integrated algorithm to calculate this value. This is the start
#' value for the grid search of the penalty parameter \eqn{\lambda}. (More
#' details about the integrated algorithms are available here:
#' \code{\link[seagull]{lambda_max}}.)
#' 
#' @param xi (optional) multiplicative parameter to determine the minimum value
#' of \eqn{\lambda} for the grid search, i.e. \eqn{\lambda_{min} = \xi *
#' \lambda_{max}}. Has to satisfy: \eqn{0 < \xi \le 1}. If \code{xi=1}, only a
#' single solution for \eqn{\lambda = \lambda_{max}} is calculated. Default
#' value is \code{0.01}.
#' 
#' @param loops_lambda (optional) number of lambdas for the grid search between
#' \eqn{\lambda_{max}} and \eqn{\xi * \lambda_{max}}. Loops are performed on a 
#' logarithmic grid. Float or double input values will be truncated. Default
#' value is \code{50}.
#' 
#' @param max_iter (optional) maximum number of iterations for each value of the
#' penalty parameter \eqn{\lambda}. Determines the end of the calculation if the
#' algorithm didn't converge according to \code{rel_acc} before. Float or double
#' input values will be truncated. Default value is \code{1000}.
#' 
#' @param gamma_bls (optional) multiplicative parameter to decrease the step
#' size during backtracking line search. Has to satisfy: \eqn{0 < \gamma_{bls}
#' < 1}. Default value is \code{0.8}.
#' 
#' @param trace_progress (optional) if \code{TRUE}, a message will occur on the
#' screen after each finished loop of the \eqn{\lambda} grid. This is
#' particularly useful for larger data sets.
#' 
#' @details The underlying mixed model has the form: \deqn{y = X b + Z u +
#' residual,} where \eqn{b} is the vector of fixed effects and \eqn{u} is the
#' vector of random effects. The penalty of the sparse-group lasso (without
#' additional weights for features) is then: \deqn{\alpha \lambda ||u||_1 + (1 -
#' \alpha) \lambda \sum_l \omega^G_l ||u^{(l)}||_2.} The variable
#' \eqn{\omega^G_l} is a weight for group \eqn{l}. If \eqn{\alpha = 1}, this
#' leads to the lasso. If \eqn{\alpha = 0}, this leads to the group lasso.
#' 
#' The above penalty can be enhanced by introducing additional weights
#' \eqn{\omega^F} for features in the lasso penalty: \deqn{\alpha \lambda
#' \sum_j | \omega^F_j u_j | + (1 - \alpha) \lambda \sum_l \omega^G_l
#' ||u^{(l)}||_2.}
#' 
#' The group weights \eqn{\omega^G_l} are implemented as follows:
#' \eqn{\omega^G_l = \sum_{j, j \in G} \sqrt{\omega^F_{j}}}. So, in the case
#' that the weights for features in a particular group are all equal to one,
#' this term collapses to the square root of the group size.
#' 
#' In addition, the above penalty can be formally rewritten by including the
#' fixed effects and assigning weights equal to zero to all these features. This
#' is how the algorithms are implemented. Consequently, if a weight for any
#' random effect is set to zero by the user, the function drops an error and the
#' calculation will not start.
#' 
#' @return A list of estimates and parameters relevant for the computation:
#' \describe{
#'   \item{fixed_effects}{estimates for the fixed effects b, if present in the
#'   model. Each row corresponds to a particular value of \eqn{\lambda}.}
#'   \item{random_effects}{predictions for the random effects u. Each row
#'   corresponds to a particular value of \eqn{\lambda}.}
#'   \item{lambda}{all values for \eqn{\lambda} which were used during the grid
#'   search.}
#'   \item{iterations}{a sequence of actual iterations for each value of
#'   \eqn{\lambda}. If an occurring number is equal to \code{max_iter}, then the
#'   algorithm most likely did not converge to \code{rel_acc} during the
#'   corresponding run of the grid search.}
#' }
#' The following parameters are also returned. But primarily for the purpose of
#' comparison and repetition: \code{alpha} (only for the sparse-group lasso),
#' \code{max_iter}, \code{gamma_bls}, \code{xi}, and \code{loops_lambda}.
#' 
#' @seealso \code{\link{plot}}
#' 
#' @examples
#' set.seed(62)
#' n <- 50  ## observations
#' p <- 8   ## variables
#' 
#' ## Create a deisng matrix X for fixed effects to model the intercept:
#' X <- matrix(1, nrow = n, ncol = 1)
#' 
#' ## Create a design matrix Z for random effects:
#' Z <- matrix(rnorm(p * n, mean = 0, sd = 1), nrow = n)
#' 
#' ## Intercept b, random effect vector u, and response y:
#' b <- 0.2
#' u <- c(0, 1.5, 0, 0.5, 0, 0, -2, 1)
#' y <- X%*%b + Z%*%u + rnorm(n, mean = 0, sd = 1)
#' 
#' ## Create a vector of three groups corresponding to vector u:
#' group_indices <- c(1L, 2L, 2L, 2L, 1L, 1L, 3L, 1L)
#' 
#' ## Calculate the solution:
#' fit_l   <- seagull(y = y, X = X, Z = Z, alpha = 1.0)
#' fit_gl  <- seagull(y = y, X = X, Z = Z, alpha = 0.0, groups = group_indices)
#' fit_sgl <- seagull(y = y, X = X, Z = Z, groups = group_indices)
#' 
#' ## Combine the estimates for fixed and random effects:
#' estimates_l   <- cbind(fit_l$fixed_effects, fit_l$random_effects)
#' estimates_gl  <- cbind(fit_gl$fixed_effects, fit_gl$random_effects)
#' estimates_sgl <- cbind(fit_sgl$fixed_effects, fit_sgl$random_effects)
#' true_effects  <- c(b, u)
#' 
#' ## Calculate mean squared errors along the solution paths:
#' MSE_l   <- rep(as.numeric(NA), fit_l$loops_lambda)
#' MSE_gl  <- rep(as.numeric(NA), fit_l$loops_lambda)
#' MSE_sgl <- rep(as.numeric(NA), fit_l$loops_lambda)
#' 
#' for (i in 1:fit_l$loops_lambda) {
#'   MSE_l[i] <- t(estimates_l[i,] - true_effects)%*%(estimates_l[i,] - true_effects)/(p+1)
#'   MSE_gl[i] <- t(estimates_gl[i,] - true_effects)%*%(estimates_gl[i,] - true_effects)/(p+1)
#'   MSE_sgl[i] <- t(estimates_sgl[i,] - true_effects)%*%(estimates_sgl[i,] - true_effects)/(p+1)
#' }
#' 
#' ## Plot a fraction of the results of the MSEs:
#' plot(x = seq(1, fit_l$loops_lambda, 1)[25:50], MSE_l[25:50], type = "l", lwd = 2)
#' points(x = seq(1, fit_l$loops_lambda, 1)[25:50], MSE_gl[25:50], type = "l", lwd = 2, col = "blue")
#' points(x = seq(1, fit_l$loops_lambda, 1)[25:50], MSE_sgl[25:50], type = "l", lwd = 2, col = "red")
#' 
#' \dontrun{
#'   
#'   ## A larger example with simulated genetic data:
#'   data(seagull_data)
#'   
#'   fit_l1 <- seagull(y = phenotypes[,1], Z = genotypes, alpha = 1.0)
#'   fit_l2 <- seagull(y = phenotypes[,2], Z = genotypes, alpha = 1.0)
#'   fit_l3 <- seagull(y = phenotypes[,3], Z = genotypes, alpha = 1.0,
#'             trace_progress = T)
#'   
#'   fit_gl1 <- seagull(y = phenotypes[,1], Z = genotypes, alpha = 0.0,
#'              groups = groups)
#'   fit_gl2 <- seagull(y = phenotypes[,2], Z = genotypes, alpha = 0.0,
#'              groups = groups)
#'   fit_gl3 <- seagull(y = phenotypes[,3], Z = genotypes, alpha = 0.0,
#'              groups = groups, trace_progress = T)
#'   
#'   fit_sgl1 <- seagull(y = phenotypes[,1], Z = genotypes, groups = groups)
#'   fit_sgl2 <- seagull(y = phenotypes[,2], Z = genotypes, groups = groups)
#'   fit_sgl3 <- seagull(y = phenotypes[,3], Z = genotypes, groups = groups,
#'               trace_progress = T)
#' }
#' 
#' @export
seagull <- function(
  y,
  X,
  Z,
  weights_u,
  groups,
  alpha,
  rel_acc,
  max_lambda,
  xi,
  loops_lambda,
  max_iter,
  gamma_bls,
  trace_progress) {
  
  
  # Check vector y and give feedback to the user, if necessary.
  # This vector shall not be empty and it shall not be a matrix.
  # If the data includes NA, NaN, or Inf, then stop:
  if (missing(y)) {
    stop("Vector y is missing. Please restart when solved.")
  }
  if (is.null(y)) {
    stop("Vector y is empty. Please restart when solved.")
  }
  if (!is.numeric(y)) {
    stop("Non-numeric values in vector y detected. Please restart when solved.")
  }
  if (!is.null(dim(y)) && dim(y)[1] != 1 && dim(y)[2] != 1) {
    stop("Input y shall not be a matrix. Please restart when solved.")
  }
  if (any(is.na(y)) || any(is.infinite(y))) {
    stop("NA, NaN, or Inf detected in vector y. Please restart when solved.")
  }
  
  
  # Check matrix X and give feedback to the user, if necessary.
  # This matrix is allowed to be empty.
  # If data is a vector, transform it.
  # If X has less rows than columns and no proper max_lambda is provided, then stop.
  # If the data includes NA, NaN, or Inf, then stop:
  if (missing(X)) {
    X <- NULL
    p1 <- 0L
  }
  if (is.null(X)) {
    p1 <- 0L
  }
  if (!is.numeric(X) && !is.null(X)) {
    stop("Non-numeric values in matrix X detected. Please restart when solved.")
  }
  if (!is.null(X)) {
    if (is.null(dim(X))) {
      if (length(y) == 1) {
        X <- matrix(X, nrow=1, ncol=length(X))
      } else {
        X <- as.matrix(X)
      }
    }
    p1 <- dim(X)[2]
  }
  if (!is.null(X) && (dim(X)[1] < dim(X)[2])) {
    # X may have less rows than columns, if max_lambda is provided. So, check this:
    if (missing(max_lambda)) {
      stop("X has less rows than columns and max_lambda is empty. Please provide a single positive value for max_lambda before restarting.")
    }
    if (is.null(max_lambda)) {
      stop("X has less rows than columns and max_lambda is NULL. Please provide a single positive value for max_lambda before restarting.")
    } else if (!is.numeric(max_lambda)) {
      stop("X has less rows than columns and max_lambda is non-numeric. Please provide a single positive value for max_lambda before restarting.")
    } else if ((length(max_lambda) > 1)) {
      stop("X has less rows than columns and the length of max_lambda is greater than 1. Please provide a single positive value for max_lambda before restarting.")
    } else if (is.na(max_lambda) || is.infinite(max_lambda)) {
      stop("X has less rows than columns and max_lambda is either NA, NaN, or Inf. Please provide a single positive value for max_lambda before restarting.")
    } else if (max_lambda <= 0) {
      stop("X has less rows than columns and max_lambda is non-positive. Please provide a single positive value for max_lambda before restarting.")
    }
  }
  if (any(is.na(X)) || any(is.infinite(X))) {
    stop("NA, NaN, or Inf detected in matrix X. Please restart when solved.")
  }
  
  
  # Check matrix Z and give feedback to the user, if necessary.
  # This matrix shall not be empty.
  # If data is a vector, transform it.
  # If the data includes NA, NaN, or Inf, then stop:
  if (missing(Z)) {
    stop("Matrix Z is missing. Please restart when solved.")
  }
  if (is.null(Z)) {
    stop("Matrix Z is empty. Please restart when solved.")
  }
  if (!is.numeric(Z)) {
    stop("Non-numeric values in matrix Z detected. Please restart when solved.")
  }
  if (is.null(dim(Z))) {
    if (length(y) == 1) {
      Z <- matrix(Z, nrow=1, ncol=length(Z))
    } else {
      Z <- as.matrix(Z)
    }
  }
  p2 <- dim(Z)[2]
  if (any(is.na(Z)) || any(is.infinite(Z))) {
    stop("NA, NaN, or Inf detected in matrix Z. Please restart when solved.")
  }
  
  
  # Check vector weights_u and give feedback to the user, if necessary.
  # This vector is allowed to be empty, but it shall not be a matrix.
  # If missing or empty, fill with 1's.
  # If the data includes NA, NaN, or Inf, then stop.
  # If the data includes values <= 0, then stop:
  if (missing(weights_u)) {
    weights_u <- rep(1.0, dim(Z)[2])
  }
  if (is.null(weights_u)) {
    weights_u <- rep(1.0, dim(Z)[2])
  }
  if (!is.numeric(weights_u)) {
    stop("Non-numeric values in vector weights_u detected. Please restart when solved.")
  }
  if (!is.null(dim(weights_u)) && dim(weights_u)[1] != 1 && dim(weights_u)[2] != 1) {
    stop("Input weights_u shall not be a matrix. Please restart when solved.")
  }
  if (any(is.na(weights_u)) || any(is.infinite(weights_u))) {
    stop("NA, NaN, or Inf detected in vector weights_u. Please restart when solved.")
  }
  if (any(weights_u <= 0.0)) {
    stop("Weights <= 0 detected. Please restart when solved.")
  }
  
  
  # Check alpha and give feedback to the user, if necessary:
  if (missing(alpha)) {
    alpha <- 0.9
  }
  if (is.null(alpha)) {
    warning("The parameter alpha is equal to NULL. Reset to default value (=0.9).")
    alpha <- 0.9
  } else if (!is.numeric(alpha)) {
    warning("The parameter alpha is non-numeric. Reset to default value (=0.9).")
    alpha <- 0.9
  } else if (length(alpha) > 1) {
    warning("The length of the parameter alpha is greater than 1. Reset to default value (=0.9).")
    alpha <- 0.9
  } else if (is.na(alpha) || is.infinite(alpha)) {
    warning("The parameter alpha is either NA, NaN, or Inf. Reset to default value (=0.9).")
    alpha <- 0.9
  } else if (alpha < 0.0 || alpha > 1.0) {
    warning("The parameter alpha is out of range. Reset to default value (=0.9).")
    alpha <- 0.9
  }
  alpha <- as.double(alpha)
  
  
  # Check vector groups and give feedback to the user, if necessary.
  # This vector shall not be empty and it shall not be a matrix.
  # If the data includes NA, NaN, or Inf, then stop:
  if (alpha < 1.0) {
    if (missing(groups)) {
      stop("Vector groups is missing. Please restart when solved.")
    }
    if (is.null(groups)) {
      stop("Vector groups is empty. Please restart when solved.")
    }
    if (!is.numeric(groups)) {
      stop("Non-numeric values in vector groups detected. Please restart when solved.")
    }
    if (!is.null(dim(groups)) && dim(groups)[1] != 1 && dim(groups)[2] != 1) {
      stop("Input groups shall not be a matrix. Please restart when solved.")
    }
    if (any(is.na(groups)) || any(is.infinite(groups))) {
      stop("NA, NaN, or Inf detected in vector groups. Please restart when solved.")
    }
    groups <- as.integer(groups)
  }
  
  
  # Check for mismatching dimensions:
  if ((length(y) != dim(X)[1]) && !is.null(X)) {
    stop("Mismatching dimensions of vector y and matrix X. Please restart when solved.")
  }
  if (length(y) != dim(Z)[1]) {
    stop("Mismatching dimensions of vector y and matrix Z. Please restart when solved.")
  }
  if (dim(Z)[2] != length(weights_u)) {
    stop("Mismatching dimensions of matrix Z and vector weights_u. Please restart when solved.")
  }
  if (alpha < 1.0 && (length(groups) != (p1 + p2)) && (length(groups) != p2)) {
    stop("Mismatching dimension of vector groups. Please assign one group value to each random effect, and either to all or to none of the fixed effects.")
  }
  
  
  # Create correct input:
  if (p1 == 0) {
    weights_u_tilde <- weights_u
    p <- p2
  } else {
    p <- p1 + p2
    weights_u_tilde <- rep(0.0, p)
    for (i in 1:p2) {
      weights_u_tilde[p1 + i] <- weights_u[i]
    }
  }
  X_tilde <- cbind(X, Z)
  b_tilde <- rep(0, p)
  
  
  if (alpha < 1.0) {
    # Assign all fixed effects to one group, if fixed effects were not assigned to any group by the user:
    if ((p1 > 0) && (length(groups) == p2)) {
      groups_temp <- rep(as.integer(min(groups)) - 1L, p1)
      groups      <- c(groups_temp, groups)
    }
    
    
    # Ensure positivity of group assignments:
    temp <- as.integer(min(groups))
    if (temp <= 0) {
      groups <- groups - temp + 1L
    }
    
    
    # Sort the vector "groups", the matrix "X_tilde", and the vector "weights_u_tilde" by ascending order of group number:
    index_permutation <- order(groups)
    if (p > 1) {
      X_tilde         <- X_tilde[, index_permutation]
      groups          <- groups[index_permutation]
      weights_u_tilde <- weights_u_tilde[index_permutation]
    }
    
    
    # Renumber the vector of groups:
    temp_diff_groups  <- sort(unique(groups))
    for (i in 1:length(temp_diff_groups)) {
      groups[groups==temp_diff_groups[i]] = i
    }
  }
  
  
  # Check rel_acc and give feedback to the user, if necessary:
  if (missing(rel_acc)) {
    rel_acc <- 0.0001
  }
  if (is.null(rel_acc)) {
    warning("The parameter rel_acc is equal to NULL. Reset to default value (=0.0001).")
    rel_acc <- 0.0001
  } else if (!is.numeric(rel_acc)) {
    warning("The parameter rel_acc is non-numeric. Reset to default value (=0.0001).")
    rel_acc <- 0.0001
  } else if (length(rel_acc) > 1) {
    warning("The length of the parameter rel_acc is greater than 1. Reset to default value (=0.0001).")
    rel_acc <- 0.0001
  } else if (is.na(rel_acc) || is.infinite(rel_acc)) {
    warning("The parameter rel_acc is either NA, NaN, or Inf. Reset to default value (=0.0001).")
    rel_acc <- 0.0001
  } else if (rel_acc <= 0.0) {
    warning("The parameter rel_acc is non-positive. Reset to default value (=0.0001).")
    rel_acc <- 0.0001
  }
  rel_acc <- as.double(rel_acc)
  
  
  # Check max_iter and give feedback to the user, if necessary:
  if (missing(max_iter)) {
    max_iter <- 1000L
  }
  if (is.null(max_iter)) {
    warning("The parameter max_iter is equal to NULL. Reset to default value (=1000).")
    max_iter <- 1000L
  } else if (!is.numeric(max_iter)) {
    warning("The parameter max_iter is non-numeric. Reset to default value (=1000).")
    max_iter <- 1000L
  } else if (length(max_iter) > 1) {
    warning("The length of the parameter max_iter is greater than 1. Reset to default value (=1000).")
    max_iter <- 1000L
  } else if (is.na(max_iter) || is.infinite(max_iter)) {
    warning("The parameter max_iter is either NA, NaN, or Inf. Reset to default value (=1000).")
    max_iter <- 1000L
  } else if (max_iter <= 0) {
    warning("The parameter max_iter is non-positive. Reset to default value (=1000).")
    max_iter <- 1000L
  }
  max_iter <- as.integer(max_iter)
  
  
  # Check gamma_bls and give feedback to the user, if necessary:
  if (missing(gamma_bls)) {
    gamma_bls <- 0.8
  }
  if (is.null(gamma_bls)) {
    warning("The parameter gamma_bls is equal to NULL. Reset to default value (=0.8).")
    gamma_bls <- 0.8
  } else if (!is.numeric(gamma_bls)) {
    warning("The parameter gamma_bls is non-numeric. Reset to default value (=0.8).")
    gamma_bls <- 0.8
  } else if (length(gamma_bls) > 1) {
    warning("The length of the parameter gamma_bls is greater than 1. Reset to default value (=0.8).")
    gamma_bls <- 0.8
  } else if (is.na(gamma_bls) || is.infinite(gamma_bls)) {
    warning("The parameter gamma_bls is either NA, NaN, or Inf. Reset to default value (=0.8).")
    gamma_bls <- 0.8
  } else if ((gamma_bls <= 0.0) || (gamma_bls >= 1.0)) {
    warning("The parameter gamma_bls is out of range. Reset to default value (=0.8).")
    gamma_bls <- 0.8
  }
  gamma_bls <- as.double(gamma_bls)
  
  
  # Check max_lambda and give feedback to the user, if necessary:
  if (missing(max_lambda)) {
    if (alpha == 1.0) {
      max_lambda <- lambda_max_lasso(y, weights_u_tilde, b_tilde, X_tilde)
    } else if (alpha == 0.0) {
      max_lambda <- lambda_max_group_lasso(y, groups, weights_u_tilde, b_tilde, X_tilde)
    } else {
      max_lambda <- lambda_max_sparse_group_lasso(alpha, y, groups, weights_u_tilde, b_tilde, X_tilde)
    }
  }
  if (is.null(max_lambda)) {
    warning("The parameter max_lambda is equal to NULL. Use default algorithm instead.")
    if (alpha == 1.0) {
      max_lambda <- lambda_max_lasso(y, weights_u_tilde, b_tilde, X_tilde)
    } else if (alpha == 0.0) {
      max_lambda <- lambda_max_group_lasso(y, groups, weights_u_tilde, b_tilde, X_tilde)
    } else {
      max_lambda <- lambda_max_sparse_group_lasso(alpha, y, groups, weights_u_tilde, b_tilde, X_tilde)
    }
  } else if (!is.numeric(max_lambda)) {
    warning("The parameter max_lambda is non-numeric. Use default algorithm instead.")
    if (alpha == 1.0) {
      max_lambda <- lambda_max_lasso(y, weights_u_tilde, b_tilde, X_tilde)
    } else if (alpha == 0.0) {
      max_lambda <- lambda_max_group_lasso(y, groups, weights_u_tilde, b_tilde, X_tilde)
    } else {
      max_lambda <- lambda_max_sparse_group_lasso(alpha, y, groups, weights_u_tilde, b_tilde, X_tilde)
    }
  } else if (length(max_lambda) > 1) {
    warning("The length of the parameter max_lambda is greater than 1. Use default algorithm instead.")
    if (alpha == 1.0) {
      max_lambda <- lambda_max_lasso(y, weights_u_tilde, b_tilde, X_tilde)
    } else if (alpha == 0.0) {
      max_lambda <- lambda_max_group_lasso(y, groups, weights_u_tilde, b_tilde, X_tilde)
    } else {
      max_lambda <- lambda_max_sparse_group_lasso(alpha, y, groups, weights_u_tilde, b_tilde, X_tilde)
    }
  } else if (is.na(max_lambda) || is.infinite(max_lambda)) {
    warning("The parameter max_lambda is either NA, NaN, or Inf. Use default algorithm instead.")
    if (alpha == 1.0) {
      max_lambda <- lambda_max_lasso(y, weights_u_tilde, b_tilde, X_tilde)
    } else if (alpha == 0.0) {
      max_lambda <- lambda_max_group_lasso(y, groups, weights_u_tilde, b_tilde, X_tilde)
    } else {
      max_lambda <- lambda_max_sparse_group_lasso(alpha, y, groups, weights_u_tilde, b_tilde, X_tilde)
    }
  } else if (max_lambda <= 0) {
    warning("The parameter max_lambda is non-positive. Use default algorithm instead.")
    if (alpha == 1.0) {
      max_lambda <- lambda_max_lasso(y, weights_u_tilde, b_tilde, X_tilde)
    } else if (alpha == 0.0) {
      max_lambda <- lambda_max_group_lasso(y, groups, weights_u_tilde, b_tilde, X_tilde)
    } else {
      max_lambda <- lambda_max_sparse_group_lasso(alpha, y, groups, weights_u_tilde, b_tilde, X_tilde)
    }
  }
  max_lambda <- as.double(max_lambda)
  
  
  # Check xi and give feedback to the user, if necessary:
  if (missing(xi)) {
    xi <- 0.01
  }
  if (is.null(xi)) {
    warning("The parameter xi is equal to NULL. Reset to default value (=0.01).")
    xi <- 0.01
  } else if (!is.numeric(xi)) {
    warning("The parameter xi is non-numeric. Reset to default value (=0.01).")
    xi <- 0.01
  } else if (length(xi) > 1) {
    warning("The length of the parameter xi is greater than 1. Reset to default value (=0.01).")
    xi <- 0.01
  } else if (is.na(xi) || is.infinite(xi)) {
    warning("The parameter xi is either NA, NaN, or Inf. Reset to default value (=0.01).")
    xi <- 0.01
  } else if ((xi <= 0.0) || (xi > 1.0)) {
    warning("The parameter xi is out of range. Reset to default value (=0.01).")
    xi <- 0.01
  }
  xi <- as.double(xi)
  
  
  # Check loops_lambda and give feedback to the user, if necessary:
  if (missing(loops_lambda)) {
    loops_lambda <- 50L
  }
  if (is.null(loops_lambda)) {
    warning("The parameter loops_lambda is equal to NULL. Reset to default value (=50).")
    loops_lambda <- 50L
  } else if (!is.numeric(loops_lambda)) {
    warning("The parameter loops_lambda is non-numeric. Reset to default value (=50).")
    loops_lambda <- 50L
  } else if (length(loops_lambda) > 1) {
    warning("The length of the parameter loops_lambda is greater than 1. Reset to default value (=50).")
    loops_lambda <- 50L
  } else if (is.na(loops_lambda) || is.infinite(loops_lambda)) {
    warning("The parameter loops_lambda is either NA, NaN, or Inf. Reset to default value (=50).")
    loops_lambda <- 50L
  } else if (loops_lambda <= 0) {
    warning("The parameter loops_lambda is non-positive. Reset to default value (=50).")
    loops_lambda <- 50L
  } else if (xi == 1.0) {
    warning("Since the parameter xi = 1, the parameter loops_lambda will be set to 1.")
    loops_lambda <- 1L
  }
  loops_lambda <- as.integer(loops_lambda)
  
  
  # Check trace_progress and give feedback to the user, if necessary:
  if (missing(trace_progress)) {
    trace_progress <- FALSE
  }
  if (is.null(trace_progress)) {
    warning("The parameter trace_progress is equal to NULL. Reset to default value (=FALSE).")
    trace_progress <- FALSE
  } else if (is.numeric(trace_progress) || is.character(trace_progress)) {
    warning("The parameter trace_progress is not Boolean. Reset to default value (=FALSE).")
    trace_progress <- FALSE
  } else if (length(trace_progress) > 1) {
    warning("The length of the parameter trace_progress is greater than 1. Reset to default value (=FALSE).")
    trace_progress <- FALSE
  } else if (is.na(trace_progress)) {
    warning("The parameter trace_progress is NA. Reset to default value (=FALSE).")
    trace_progress <- FALSE
  }
  
  
  # Calculate the solution. Choose the algorithm according to alpha:
  if (alpha == 1.0) {
    res <- seagull_lasso(y, X_tilde, weights_u_tilde, b_tilde, rel_acc,
                         max_iter, gamma_bls, max_lambda, xi, loops_lambda, p1,
                         trace_progress)
  } else if (alpha == 0.0) {
    res <- seagull_group_lasso(y, X_tilde, weights_u_tilde, groups, b_tilde,
                               index_permutation, rel_acc, max_iter, gamma_bls,
                               max_lambda, xi, loops_lambda, p1, trace_progress)
  } else {
    res <- seagull_sparse_group_lasso(y, X_tilde, weights_u_tilde, groups,
                                      b_tilde, index_permutation, alpha,
                                      rel_acc, max_iter, gamma_bls, max_lambda,
                                      xi, loops_lambda, p1, trace_progress)
  }
}