context("seagull: lasso max_lambda")

vec1  <- c(3, 3, 3)
matX  <- matrix(c(1, 0, 1, 0, 2, 2), 3, 2)
matZ1 <- matrix(c(3, 0, 0, 0, 6, 0, 0, 0, 9), 3, 3)
matZ2 <- matrix(c(3, 0, 0, 0, 6, 0), 3, 2)
matZ3 <- matrix(c(3, 0, 0), 3, 1)

w_u1  <- c(0, 0, 1, 1, 1)
w_u2  <- c(0, 0, 1, 1)
w_u3  <- c(0, 0, 1)
w_u4  <- c(1, 1, 1)

beta1 <- rep(0, 5)
beta2 <- rep(0, 4)
beta3 <- rep(0, 3)

delta <- 1.00001

res1  <- 3 * delta
res2  <- 2 * delta
res3  <- 1 * delta
res4  <- 9 * delta

test_that("The parameter max_lambda for lasso is correctly calculated.", {
  expect_equal(res1, lambda_max_lasso(vec1, w_u1, beta1, cbind(matX, matZ1)))
  expect_equal(res2, lambda_max_lasso(vec1, w_u2, beta2, cbind(matX, matZ2)))
  expect_equal(res3, lambda_max_lasso(vec1, w_u3, beta3, cbind(matX, matZ3)))
  expect_equal(res4, lambda_max_lasso(vec1, w_u4, beta3, matZ1))
})