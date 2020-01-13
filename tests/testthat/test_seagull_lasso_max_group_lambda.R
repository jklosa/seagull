context("seagull: group lasso max_lambda")

vec1   <- c(3, 3, 3)
matX   <- matrix(c(1, 0, 1, 0, 2, 2), 3, 2)
matZ1  <- matrix(c(3, 0, 0, 0, 6, 0, 0, 0, 9), 3, 3)
matZ2  <- matrix(c(3, 0, 0, 0, 6, 0), 3, 2)
matZ3  <- matrix(c(3, 0, 0), 3, 1)

g1_1   <- c(1L, 2L, 3L, 4L, 5L)
g1_2   <- c(1L, 1L, 2L, 3L, 4L)
g1_3   <- c(1L, 1L, 2L, 3L, 3L)
g1_4   <- c(1L, 1L, 2L, 2L, 2L)
g2_1   <- c(1L, 1L, 2L, 3L)
g2_2   <- c(1L, 1L, 2L, 2L)
g3_1   <- c(1L, 1L, 2L)
g4_1   <- c(1L, 1L, 1L)

w_u1   <- c(0, 0, 1, 1, 1)
w_u2   <- c(0, 0, 1, 1)
w_u3   <- c(0, 0, 1)
w_u4   <- c(1, 1, 1)

beta1  <- rep(0, 5)
beta2  <- rep(0, 4)
beta3  <- rep(0, 3)

delta  <- 1.00001

res1_1 <- 3 * delta
res1_2 <- 3 * delta
res1_3 <- sqrt(13/2) * delta
res1_4 <- sqrt(14/3) * delta
res2_1 <- 2 * delta
res2_2 <- sqrt(5/2) * delta
res3_1 <- 1 * delta
res4_1 <- sqrt(42) * delta

test_that("The parameter max_lambda for lasso is correctly calculated.", {
  expect_equal(res1_1, lambda_max_group_lasso(vec1, g1_1, w_u1, beta1, cbind(matX, matZ1)))
  expect_equal(res1_2, lambda_max_group_lasso(vec1, g1_2, w_u1, beta1, cbind(matX, matZ1)))
  expect_equal(res1_3, lambda_max_group_lasso(vec1, g1_3, w_u1, beta1, cbind(matX, matZ1)))
  expect_equal(res1_4, lambda_max_group_lasso(vec1, g1_4, w_u1, beta1, cbind(matX, matZ1)))
  expect_equal(res2_1, lambda_max_group_lasso(vec1, g2_1, w_u2, beta2, cbind(matX, matZ2)))
  expect_equal(res2_2, lambda_max_group_lasso(vec1, g2_2, w_u2, beta2, cbind(matX, matZ2)))
  expect_equal(res3_1, lambda_max_group_lasso(vec1, g3_1, w_u3, beta3, cbind(matX, matZ3)))
  expect_equal(res4_1, lambda_max_group_lasso(vec1, g4_1, w_u4, beta3, matZ1))
})