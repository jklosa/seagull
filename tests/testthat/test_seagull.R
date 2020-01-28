context("seagull: wrapper")

test_that("Test for input y works.", {
  expect_error(seagull(), "Vector y is missing. Please restart when solved.")
  expect_error(seagull(y = NULL), "Vector y is empty. Please restart when solved.")
  expect_error(seagull(y = TRUE), "Non-numeric values in vector y detected. Please restart when solved.")
  expect_error(seagull(y = matrix(1, 2, 2)), "Input y shall not be a matrix. Please restart when solved.")
  expect_error(seagull(y = as.numeric(NA)), "NA, NaN, or Inf detected in vector y. Please restart when solved.")
  expect_error(seagull(y = as.numeric(Inf)), "NA, NaN, or Inf detected in vector y. Please restart when solved.")
})

test_that("Test for input X works.", {
  expect_error(seagull(y = 1, X = TRUE), "Non-numeric values in matrix X detected. Please restart when solved.")
  expect_error(seagull(y = 1, X = matrix(1, 1, 2)), "X has less rows than columns and max_lambda is empty. Please provide a single positive value for max_lambda before restarting.")
  expect_error(seagull(y = 1, X = matrix(1, 1, 2), max_lambda = c()), "X has less rows than columns and max_lambda is NULL. Please provide a single positive value for max_lambda before restarting.")
  expect_error(seagull(y = 1, X = matrix(1, 1, 2), max_lambda = TRUE), "X has less rows than columns and max_lambda is non-numeric. Please provide a single positive value for max_lambda before restarting.")
  expect_error(seagull(y = 1, X = matrix(1, 1, 2), max_lambda = c(1, 1)), "X has less rows than columns and the length of max_lambda is greater than 1. Please provide a single positive value for max_lambda before restarting.")
  expect_error(seagull(y = 1, X = matrix(1, 1, 2), max_lambda = as.numeric(NA)), "X has less rows than columns and max_lambda is either NA, NaN, or Inf. Please provide a single positive value for max_lambda before restarting.")
  expect_error(seagull(y = 1, X = matrix(1, 1, 2), max_lambda = as.numeric(Inf)), "X has less rows than columns and max_lambda is either NA, NaN, or Inf. Please provide a single positive value for max_lambda before restarting.")
  expect_error(seagull(y = 1, X = matrix(1, 1, 2), max_lambda = 0), "X has less rows than columns and max_lambda is non-positive. Please provide a single positive value for max_lambda before restarting.")
  expect_error(seagull(y = 1, X = as.numeric(NA)), "NA, NaN, or Inf detected in matrix X. Please restart when solved.")
  expect_error(seagull(y = 1, X = as.numeric(Inf)), "NA, NaN, or Inf detected in matrix X. Please restart when solved.")
})

test_that("Test for input Z works.", {
  expect_error(seagull(y = 1, X = 1), "Matrix Z is missing. Please restart when solved.")
  expect_error(seagull(y = 1, X = 1, Z = NULL), "Matrix Z is empty. Please restart when solved.")
  expect_error(seagull(y = 1, X = 1, Z = TRUE), "Non-numeric values in matrix Z detected. Please restart when solved.")
  expect_error(seagull(y = 1, X = 1, Z = as.numeric(NA)), "NA, NaN, or Inf detected in matrix Z. Please restart when solved.")
  expect_error(seagull(y = 1, X = 1, Z = as.numeric(Inf)), "NA, NaN, or Inf detected in matrix Z. Please restart when solved.")
})

test_that("Test for input weights_u works.", {
  expect_error(seagull(y = 1, X = 1, Z = 1, weights_u = TRUE), "Non-numeric values in vector weights_u detected. Please restart when solved.")
  expect_error(seagull(y = 1, X = 1, Z = 1, weights_u = matrix(2, 2, 2)), "Input weights_u shall not be a matrix. Please restart when solved.")
  expect_error(seagull(y = 1, X = 1, Z = 1, weights_u = as.numeric(NA)), "NA, NaN, or Inf detected in vector weights_u. Please restart when solved.")
  expect_error(seagull(y = 1, X = 1, Z = 1, weights_u = as.numeric(Inf)), "NA, NaN, or Inf detected in vector weights_u. Please restart when solved.")
  expect_error(seagull(y = 1, X = 1, Z = 1, weights_u = 0), "Weights <= 0 detected. Please restart when solved.")
})

test_that("Test for input group works.", {
  expect_error(seagull(y = 1, X = 1, Z = 1), "Vector groups is missing. Please restart when solved.")
  expect_error(seagull(y = 1, X = 1, Z = 1, groups = NULL), "Vector groups is empty. Please restart when solved.")
  expect_error(seagull(y = 1, X = 1, Z = 1, groups = TRUE), "Non-numeric values in vector groups detected. Please restart when solved.")
  expect_error(seagull(y = 1, X = 1, Z = 1, groups = matrix(1, 2, 2)), "Input groups shall not be a matrix. Please restart when solved.")
  expect_error(seagull(y = 1, X = 1, Z = 1, groups = as.numeric(NA)), "NA, NaN, or Inf detected in vector groups. Please restart when solved.")
  expect_error(seagull(y = 1, X = 1, Z = 1, groups = as.numeric(Inf)), "NA, NaN, or Inf detected in vector groups. Please restart when solved.")
})

test_that("Test for dimensionality works.", {
  expect_error(seagull(y = c(1, 1), X = 1, Z = matrix(1, 2, 3), groups = 1), "Mismatching dimensions of vector y and matrix X. Please restart when solved.")
  expect_error(seagull(y = c(1, 1), X = matrix(1, 2, 2), Z = 1, groups = 1), "Mismatching dimensions of vector y and matrix Z. Please restart when solved.")
  expect_error(seagull(y = c(1, 1), X = matrix(1, 2, 2), Z = matrix(1, 2, 3), weights_u = 1, groups = 1), "Mismatching dimensions of matrix Z and vector weights_u. Please restart when solved.")
  expect_error(seagull(y = c(1, 1), X = matrix(1, 2, 2), Z = matrix(1, 2, 3), groups = 1), "Mismatching dimension of vector groups. Please assign one group value to each random effect, and either to all or to none of the fixed effects.")
})

# test_that("Set alpha to default.", {
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, alpha = NULL))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, alpha = "a"))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, alpha = c(1, 1)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, alpha = as.numeric(NA)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, alpha = as.numeric(Inf)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, alpha = -0.1))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, alpha = 1.1))
# })
# 
# test_that("Set relative accuracy to default.", {
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, rel_acc = NULL))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, rel_acc = "a"))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, rel_acc = c(1, 1)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, rel_acc = as.numeric(NA)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, rel_acc = as.numeric(Inf)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, rel_acc = 0.0))
# })
# 
# test_that("Set number of maximal iterations to default.", {
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, max_iter = NULL))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, max_iter = "a"))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, max_iter = c(1, 1)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, max_iter = as.numeric(NA)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, max_iter = as.numeric(Inf)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, max_iter = 0))
# })
# 
# test_that("Set gamma for backtracking line search to default.", {
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, gamma_bls = NULL))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, gamma_bls = "a"))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, gamma_bls = c(1, 1)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, gamma_bls = as.numeric(NA)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, gamma_bls = as.numeric(Inf)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, gamma_bls = 0))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, gamma_bls = 1))
# })
# 
# test_that("Calculate maximal lambda with default algorithm.", {
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, max_lambda = NULL))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, max_lambda = "a"))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, max_lambda = c(1, 1)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, max_lambda = as.numeric(NA)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, max_lambda = as.numeric(Inf)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, max_lambda = 0))
# })
# 
# test_that("Set xi to default.", {
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, xi = NULL))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, xi = "a"))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, xi = c(1, 1)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, xi = as.numeric(NA)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, xi = as.numeric(Inf)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, xi = 0))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, xi = 1.1))
# })
# 
# test_that("Set number of loops for lambda to default.", {
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, loops_lambda = NULL))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, loops_lambda = "a"))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, loops_lambda = c(1, 1)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, loops_lambda = as.numeric(NA)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, loops_lambda = as.numeric(Inf)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, loops_lambda = 0))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, xi = 1.0))
# })
# 
# test_that("Set trace_progress to default.", {
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, trace_progress = NULL))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, trace_progress = "a"))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, trace_progress = 1))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, trace_progress = c(TRUE, TRUE)))
#   expect_warning(seagull(y = 1, Z = 1, groups = 1, trace_progress = as.logical(NA)))
# })
