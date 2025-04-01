Download the latest version:  
[![Download](https://img.shields.io/badge/Download-ZIP-blue.svg)](https://github.com/Wei-M-Wei/Survivalml/blob/master/Survivalml_0.1.0.tar.gz)
## Survivalml

This is the first version of the R package 'Survivalml'. Please note that it does not include the raw dataset used in the empirical section. The dataset is available upon request.

To install and use the package, we recommend downloading 'Survivalml_0.1.0.tar.gz' and installing the package locally by
```{r }
install.packages('your path/Survivalml_0.1.0.tar.gz')
```
Another possible way to install
```{r }
install.packages("devtools")  # If not already installed
library(devtools)
install_github("Wei-M-Wei/Survivalml")
```
Once installed, load the package with
```{r }
library(Survivalml)
```
A CRAN release is coming soon.

Several other packages are needed
```{r }
packages <- unique(c("midasml", "dplyr", "RSpectra", "pROC", "openxlsx", "ggplot2", 
              "xtable", "caret", "survival", "parallel", "foreach", "doParallel", 
              "PRROC", "MLmetrics", "survivalROC", "pracma", "dotCall64", 
              "rlang", "readxl", "lubridate", "timeROC"))

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

sapply(packages, install_if_missing)

lapply(packages, library, character.only = TRUE)
```

## Description

This package is based on code from the following sources:
- **sparsegl** (https://cran.r-project.org/web/packages/sparsegl/index.html): Provides a fast implementation of the sparse-group LASSO estimator [^1].
- **midasml** (https://cran.r-project.org/web/packages/midasml/index.html): Applies MIDAS into logistic regression [^2].

## Features
- **Main functionality**: The primary function of this package estimates outcome weighted logistic regression with sparse group LASSO penalty, which can also be reduced to a standard logistic regression model and allows for censored data. For more details, see the paper [^3].
- **Validation example**: An example is included to compare this package with:
  - **glmnet** for LASSO.
  - **optim** for incorporating weights into the logistic model.
    
The example can be used to verify the correctness of the package, see more details in 'check the correctness of this package.R'. Another example of recovering the true parameters in the weighted logistic regression is 'correctness check for weighted logistic regression without penalty and MIDAS.R'.

## Additional resources
- **Replication code**: The repository includes replication code for all simulations and empirical applications.
- **Suggestions welcome**: Further improvements are planned, and we encourage feedback and suggestions to enhance the package.

## An example:
```{r }
library(Survivalml)
n <- 2000
p <- 10
x1 <- rnorm(n)
x2 <- 0.2*x1 + rnorm(n)
x3 <- matrix(rnorm(n*p), nrow = n, ncol = p)

# Define parameters for generating y
intercept <- 5
coef_x1 <- 2
coef_x2 <- 0
coef_x3 <- rep(0,p)

# Generate y based on logistic function
logistic_function <- intercept + coef_x1 * x1 + coef_x2 * x2 + x3 %*% coef_x3
probabilities <- plogis(logistic_function)
y <- rbinom(n, 1, probabilities)
X <- data.frame(x1, x2, x3)
X <- data.matrix(X)
index <- seq(p+2)

# weight is specified randomly
weight <- c(rep(0.2,100), rep(1,100), rep(1,n - 200))

# LASSO when alpha = 1, Group LASSO when alpha = 0
# intercept_zero could be set arbitrarily sicne it's a starting point of the block coordinate descent algorithm, we use 'intercept_zero = 0' through all the simulations and empirical applications
fit = survival_sparsegl(X, y, group = index, nlambda = 100, asparse = 1, weight = weight, intercept_zero = 0, standardize = TRUE)
fit$beta

# Group LASSO
fit = survival_sparsegl(X, y, group = index, nlambda = 100, asparse = 0, weight = weight, intercept_zero = 0, standardize = TRUE)
fit$beta

# Sparse group LASSO
fit = survival_sparsegl(X, y, group = index, nlambda = 100, asparse = 0.5, weight = weight, intercept_zero = 0, standardize = TRUE)
fit$beta

# Cross-validation without censored data, where ' pred.loss = 'censor' ' allows maximizing the weighted log-likelihood. If the data is not censored, ' pred.loss = 'censor' ' maximizes the log-likelihood.
fit_cv = cv.survival_sparsegl(X, y, group = index, nlambda = 100, asparse = 0.5, weight = weight, intercept_zero = 0, standardize = TRUE, asparse = 0.5, nfolds = 5, pred.loss = 'censor', intercept_zero = 0, standardize = TRUE)

# In the case of censored data, please check the 'survival-estimate.r' file, where we describe functions in detail.
    
    
```

## Reference

[^1]: Liang, X., Cohen, A., Sólon Heinsfeld, A., Pestilli, F., & McDonald, D. J. (2024). sparsegl: An R Package for Estimating Sparse Group Lasso. Journal of Statistical Software, 110(6), 1–23. https://doi.org/10.18637/jss.v110.i06

[^2]: Babii, A., Ghysels, E., & Striaukas, J. (2022). Machine learning time series regressions with an application to nowcasting. Journal of Business & Economic Statistics, 40(3), 1094-1106.

[^3]: Miao, W., Beyhum, J., Striaukas, J., & Van Keilegom, I. High-dimensional censored MIDAS logistic regression for corporate survival forecasting. [arXiv:2502.09740
Search](https://arxiv.org/abs/2502.09740).
