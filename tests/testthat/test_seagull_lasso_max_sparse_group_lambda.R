context("seagull: sparse-group lasso max_lambda")

vec1   <- c(3, 3, 3)
matX   <- matrix(c(1, 0, 1, 0, 2, 2), 3, 2)
matZ1  <- matrix(c(3, 0, 0, 0, 6, 0, 0, 0, 9), 3, 3)

alpha1 <- 0.1
alpha2 <- 0.5
alpha3 <- 0.9

g1_1   <- c(1L, 2L, 3L, 4L, 5L)
g1_2   <- c(1L, 1L, 2L, 3L, 4L)
g2_1   <- c(1L, 2L, 3L)

w_u1   <- c(0, 0, 1, 1, 1)
w_u2   <- c(1, 1, 1)

beta1  <- rep(0, 5)
beta2  <- rep(0, 3)

delta  <- 1.00001

res1   <- 3 * delta
res2   <- 9 * delta

test_that("The parameter max_lambda for lasso is correctly calculated.", {
  expect_equal(res1, lambda_max_sparse_group_lasso(alpha1, vec1, g1_1, w_u1, beta1, cbind(matX, matZ1)))
  expect_equal(res1, lambda_max_sparse_group_lasso(alpha2, vec1, g1_1, w_u1, beta1, cbind(matX, matZ1)))
  expect_equal(res1, lambda_max_sparse_group_lasso(alpha3, vec1, g1_1, w_u1, beta1, cbind(matX, matZ1)))
  expect_equal(res1, lambda_max_sparse_group_lasso(alpha1, vec1, g1_2, w_u1, beta1, cbind(matX, matZ1)))
  expect_equal(res1, lambda_max_sparse_group_lasso(alpha2, vec1, g1_2, w_u1, beta1, cbind(matX, matZ1)))
  expect_equal(res1, lambda_max_sparse_group_lasso(alpha3, vec1, g1_2, w_u1, beta1, cbind(matX, matZ1)))
  expect_equal(res2, lambda_max_sparse_group_lasso(alpha1, vec1, g2_1, w_u2, beta2, matZ1))
  expect_equal(res2, lambda_max_sparse_group_lasso(alpha2, vec1, g2_1, w_u2, beta2, matZ1))
  expect_equal(res2, lambda_max_sparse_group_lasso(alpha3, vec1, g2_1, w_u2, beta2, matZ1))
})