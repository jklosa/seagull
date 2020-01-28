
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SeaGuLl

This package provides regularization paths for the lasso, group lasso,
and sparse-group lasso. The underlying mathematical model is a mixed
model, i.e., a model with fixed and random effects. (Whereas it is
actually optional to include any fixed effect.)

The sparse-group lasso contains two penalty terms, which are combined
via a mixing parameter `0 <= alpha <= 1`. Thus, if the parameter is set
to either `1` or `0`, the resulting regularization operator is the lasso
or the the group lasso, respectively.

Key features:

  - The lasso, group lasso, and sparse-group lasso are implemented via
    *proximal gradient descent*

  - By default, a grid search for the penalty parameter `lambda` is
    performed. *Warm starts* are implemented to effectively accelerate
    this procedure.

  - The step size between consecutive iterations is automatically
    determined via *backtracking line search*.

# Installation

To get the current release version from CRAN, please type:

``` r
install.packages("seagull")
```

To get the current development version from github, please type:

``` r
# install.packages("devtools")
devtools::install_github("jklosa/seagull")
```

# Components

A data set is included and can be loaded:

``` r
data("seagull_data")
```

Furthermore, the following functions are available to the user:

  - `seagull`

  - `lambda_max__lasso`

  - `lambda_max_group_lasso`

  - `lambda_max_sparse_group_lasso`

# Example

Please load the data as shown in the section above and get started:

``` r
## Call the lasso:
fit_l <- seagull(y = phenotypes[, 1], Z = genotypes, alpha = 1)

## Call the group lasso:
fit_gl <- seagull(y = phenotypes[, 1], Z = genotypes, groups = groups, alpha = 0)

## Call the sparse-group lasso:
fit_sgl <- seagull(y = phenotypes[, 1], Z = genotypes, groups = groups)
```
