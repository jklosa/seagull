## ----gobal_options, include=FALSE---------------------------------------------
knitr::opts_chunk$set(fig.width=10, fig.height=4.5, 
                      echo=FALSE)

## ---- eval = TRUE, echo = FALSE-----------------------------------------------
library("seagull")
data("seagull_data")

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("jklosa/seagull")
#  library("seagull")
#  data("seagull_data")

## ---- eval = TRUE, echo = TRUE------------------------------------------------
fit_sgl1 <- seagull(y = phenotypes[, 1], Z = genotypes, groups = groups)

## ---- eval = TRUE, echo = TRUE------------------------------------------------
attributes(fit_sgl1)

## ---- eval = TRUE, echo = TRUE------------------------------------------------
last_solution <- fit_sgl1$loops_lambda
plot(x = seq(1, dim(genotypes)[2], 1),
     y = fit_sgl1$random_effects[last_solution,],
     xlab = "position", ylab = "effect estimate",
     col = "gray80", pch = 16)

## ---- eval = TRUE, echo = FALSE-----------------------------------------------
print(c("The number of ZEROS in the last solution is: ", length(which(fit_sgl1$random_effects[last_solution,]==0))))

## ---- eval = TRUE, echo = 3---------------------------------------------------
last_solution <- fit_sgl1$loops_lambda
plot(x = seq(1, dim(genotypes)[2], 1),
     y = fit_sgl1$random_effects[last_solution,],
     xlab = "position", ylab = "effect estimate",
     col = "gray80", pch = 16)
points(x = seq(1, dim(genotypes)[2], 1),
       y = fit_sgl1$random_effects[20,],
       pch = 16)

## ---- eval = TRUE, echo = FALSE-----------------------------------------------
print(c("The number of ZEROS in the solution in line 20 is: ", length(which(fit_sgl1$random_effects[20,]==0))))

