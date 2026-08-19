## Survivalml [![Download](https://img.shields.io/badge/Download-ZIP-blue.svg)](https://github.com/Wei-M-Wei/Survivalml/raw/master/Survivalml_0.1.0.tar.gz)

Here is R package 'Survivalml'. Please note that the package itself does not include the raw dataset used in the empirical section. The dataset is available in 'reproduce package'.

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

# help() to check the details of the main function
help(survival_sparsegl)
```

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

This package is developed based on code from the following sources:
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

# In the case of censored data, please check the 'survival_estimate.r' file, where we describe functions in detail.
    
    
```

## Replication instructions

The replication scripts and processed datasets are in the `reproduce package/`
directory. Start R in that directory (or set it as the working directory)
because the scripts use relative paths and several of them source
`import functions for empirical application.R`.

Install `Survivalml` and the packages listed above before running the scripts.
The empirical datasets are available upon request. The simulation and empirical applications can be computationally
intensive because they use bootstrap repetitions, cross-validation, and
parallel processing. Adjust the number of workers to suit your computer.

Use the exact filenames below.

### Main figures and tables

- **Figure 2 and Table 2 (variance ratios and correlation heatmap):**
  `Figure 2 and Table 2 (variance ratio and correlation heatmap).R`
- **Table 4 (prediction simulations):**
  `simulation scenario 1 revised.R`,
  `Simulation scenario 2 revised.R`,
  `Simulation scenario 3 revised.R`,
  `Simulation scenario 4 revised.R`, and
  `simulation scenario 5 revised.R`
- **Table 5 (LASSO-U signal-to-noise-ratio simulations):**
  `SNR LASSOU scenario 1.R`, `SNR LASSOU scenario 2.R`, and
  `SNR LASSOU scenario 5.R`
- **Table 6:** generated by the five prediction-simulation scripts listed for
  Table 4
- **Table 7 (inference simulations):** run both the size and power scripts for
  each scenario:
  `inference test scenario 1.R`,
  `inference test scenario 1 power.R`,
  `inference test scenario 2.R`,
  `inference test scenario 2 power.R`,
  `inference test scenario 3.R`,
  `inference test scenario 3 power.R`,
  `inference test scenario 4.R`,
  `inference test scenario 4 power.R`,
  `inference test scenario 5.R`, and
  `inference test scenario 5 power.R`
- **Tables 9 and 10 (main empirical application):**
  `application 1 for github (s = 6).R`,
  `application 1 for github (s = 10).R`,
  `Benchmark for github (s = 6).R`,
  `Benchmark for github (s = 10).R`,
  `Macro_augmented for github (s = 6).R`,
  `Macro_augmented for github (s = 10).R`,
  `Oversample for github (s = 6).R`, and
  `Oversample for github (s = 10).R`
- **Table 11 (pairwise significance tests):** first generate the empirical
  bootstrap output files, then run `siginificance test.R` (the filename uses
  this spelling)

### Supplementary results

- **Tables S2-S6 (arXiv Tables 13-17):** generated by the five
  prediction-simulation scripts listed for Table 4
- **Figures S1-S3 (arXiv Figures 4-6):** first generate the coefficient files
  with the empirical-application scripts, then run
  `Fig 3 and Fig 6_coefficient analysis for s=6 years.R` and
  `Fig 4 and Fig 6_cofficient analysis for s=10 years.R`
- **Figures S4-S5 (arXiv Figures 7-8):**
  `inference in empirical application s = 6.R` and
  `inference in empirical application s = 10.R`
- **Tables S7-S8 (arXiv Tables 18-19):**
  `application 1 for github (s = 6) another group structure 1 (grouped by covariate type).R`,
  `application 1 for github (s = 6) another group structure 2 (grouped by correlation matrix).R`,
  `application 1 for github (s = 10) another group structure 1 ((grouped by covariate type)).R`, and
  `application 1 for github (s = 10) another group structure 2 ((grouped by correlation matrix)).R`
  <!--
- **Tables S9-S10 (arXiv Tables 20-21):**
  `application 2 for github (s = 6).R`,
  `application 2 for github (s = 10).R`,
  `application 2 Benchmark for github (s = 6).R`,
  `application 2 Benchmark for github (s = 10).R`,
  `application 2 Macro for github (s = 6).R`, and
  `application 2 Macro for github (s = 10).R`
-->

## Reference

[^1]: Liang, X., Cohen, A., Sólon Heinsfeld, A., Pestilli, F., & McDonald, D. J. (2024). sparsegl: An R Package for Estimating Sparse Group Lasso. Journal of Statistical Software, 110(6), 1–23. https://doi.org/10.18637/jss.v110.i06

[^2]: Babii, A., Ghysels, E., & Striaukas, J. (2022). Machine learning time series regressions with an application to nowcasting. Journal of Business & Economic Statistics, 40(3), 1094-1106.

[^3]: Miao, W., Beyhum, J., Striaukas, J., & Van Keilegom, I. High-dimensional censored MIDAS logistic regression for corporate survival forecasting. [arXiv:2502.09740
Search](https://arxiv.org/abs/2502.09740).
