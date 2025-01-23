## Survivalml

This is the first version of the R package survivalml. Please note that it does not include the raw dataset used in the empirical section. The dataset is available upon request.

To install and use the package, we recommend downloading 'Survivalml_0.1.0.tar.gz' and installing the package locally by
```{r }
install.packages('your path/Survivalml_0.1.0.tar.gz')
```
Once installed, load the package with
```{r }
library(Survivalml)
```
A CRAN release is coming soon.


## Description

This package is based on code from the following sources:
- **sparsegl** (https://cran.r-project.org/web/packages/sparsegl/index.html): Provides a fast implementation of the sparse-group LASSO estimator [^1].
- **midasml** (https://cran.r-project.org/web/packages/midasml/index.html): Applies MIDAS into logistic regression [^2].

## Features
- **Main functionality**: The primary function of this package estimates outcome weighted logistic regression with sparse group LASSO penalty, which can also be reduced to a standard logistic regression model and allows for censored data. For more details, see [^3].
- **Validation example**: An example is included to compare this package with:
  - **glmnet** for LASSO
  - **optim** for incorporating weights into the logistic model.
    
This example can be used to verify the correctness of the package (see 'check the correctness of this package.R').

## Additional resources
- **Replication code**: The repository includes replication code for all simulations and empirical applications.
- **Suggestions welcome**: Further improvements are planned, and we encourage feedback and suggestions to enhance the package.

## An example:
```{r }
library(Survivalml)
n <- 2000
p=10
x1 <- rnorm(n)
x2 <- 0.2*x1 + rnorm(n)
x3= matrix(rnorm(n*p), nrow = n, ncol = p)

# Define parameters for generating y
intercept <- 5
coef_x1 <- 2
coef_x2 <- 0
coef_x3 <- rep(0,p)
#Generate y based on logistic function
logistic_function <- intercept + coef_x1 * x1 + coef_x2 * x2 + x3 %*% coef_x3
probabilities <- plogis(logistic_function)
y <- rbinom(n, 1, probabilities)
X <- data.frame(x1, x2, x3)
X = data.matrix(X)
index = seq(p+2)
#weight is specified randomly
weight = c(rep(0.2,100), rep(1,100), rep(1,n - 200))

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

```

## Reference

[^1]: Liang, X., Cohen, A., Sólon Heinsfeld, A., Pestilli, F., & McDonald, D. J. (2024). sparsegl: An R Package for Estimating Sparse Group Lasso. Journal of Statistical Software, 110(6), 1–23. https://doi.org/10.18637/jss.v110.i06

[^2]: Babii, A., Ghysels, E., & Striaukas, J. (2022). Machine learning time series regressions with an application to nowcasting. Journal of Business & Economic Statistics, 40(3), 1094-1106.

[^3]: Wei, M., Jad, B., Jonas, S., & Ingrid Van, K. High-dimensional censored MIDAS logistic regression for corporate survival forecasting.
