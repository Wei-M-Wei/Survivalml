# Script: SNR LASSOU scenario 1.R
# Corporate Survival Analysis with Machine Learning Methods
# Purpose: Run the LASSO-U signal-to-noise-ratio simulation for scenario 1.
#============================================================


#============================================================
# Reproduce simulation scenario 1
#============================================================



rm(list = ls())

source('import functions for empirical application.R')



# Load the required packages
# library(midasml)
library(dplyr)
library(RSpectra)
library(pROC)
library(openxlsx)
library(ggplot2)
library(xtable)
library(caret)
library(survival)
library (parallel)
library(foreach)
library(doParallel)
library (parallel)
library(foreach)
library(doParallel)
library(PRROC)
library(MLmetrics)
library(survivalROC)
library(pracma)
library(dotCall64)
library(rlang)
library(readxl)
library(Survivalml)
library(lubridate)
library(timeROC)

# Data-generating process: scenario 1
generateData_AR <-
  function(s = s, numberOfObservations, numhv, numtrue, degree, jmax, parameters = NA, censor_strength){
    TN = numberOfObservations
    p <- numhv # Set the matrix dimensions
    degree <- degree # Polynomial degree
    jmax <- jmax # Number of high-frequency observations per low-frequency observation
    # Quarterly covariates
    Xdd <- array(NA, dim = c(TN, jmax, numhv))
    # Storage for Legendre-weighted quarterly aggregates
    Xdw <- matrix(NA, nrow = TN, ncol = (degree+1)*numhv)
    # Storage for beta-weighted quarterly aggregates
    Xdw_true <- matrix(NA, nrow = TN, ncol = numtrue)
    # Storage for quarterly data
    Xdd_g <- array(NA, dim = c(TN, jmax, numhv))
    # Storage for Legendre-weighted quarterly aggregates
    Xdw_g <- matrix(NA, nrow = TN, ncol = (degree+1)*numhv)
    # Storage for beta-weighted quarterly aggregates
    Xdw_true_g <- matrix(NA, nrow = TN, ncol = numtrue)
    ##########################################################

    Xd <- matrix(0, nrow = TN*jmax, ncol = p)

    # Level of cross-sectional dependence across variables
    phi = 0.1

    # Degree of time-series dependence across lags
    rho = 0.1
    corr_matrix <- diag(1, p)

    for (v in 1:(p)) {
      for (u in 1:(p)) {
        corr_matrix[v, u] <- (phi)^abs(v - u)
      }
    }

    for ( t in 1:(jmax)){
      if (t == 1){
        Xdd[,t,(1:p)] <- MASS::mvrnorm(TN, mu = rep(0, p), corr_matrix) # Draw the initial innovations
      }
      else{
        Xdd[,t,(1:p)] <- rho * Xdd[,t-1,(1:p)] + MASS::mvrnorm(TN, mu = rep(0, p), corr_matrix*(1 - rho^2))
      }
    }

    Xdd = abs(Xdd)

    # Storage for all lags used in LASSO-MIDAS
    X_orginal = matrix(1, nrow = TN)
    for (i in seq(numhv)){
      X_orginal = cbind(X_orginal, Xdd[,,i])
    }

    # Specify the underlying weight functions w1 and w2, representing Beta(1,3) and Beta(2,3), respectively.
    idx <- 1:(degree+1)
    gindex <- NULL

    for (z in seq(p)){
      z_idx <- idx + (z-1)*(degree+1)
      Xdw[,z_idx] <- Xdd[,,z]%*%w
      if (z == 1){
        Xdw_true[,z] <- Xdd[,,z]%*%w1
      }
      if (z == 2){
        Xdw_true[,z] <- Xdd[,,z]%*%w2
      }
      gindex <- c(gindex, rep(z, times = degree + 1))
    }

    #
    dataset = cbind(rep(1, numberOfObservations),  Xdw)

    # Logit function (probability to log-odds)
    logit_inv <- function(p) log( p / (1 - p))  # Logit function (probability to log-odds)
    uniform_rand <- runif(numberOfObservations)

    # X^T%*%l
    div = cbind(rep(1, numberOfObservations) ,  Xdw_true) %*% c(1, rep(1, length(parameters) - 1))

    dataset %>%
      as_tibble() %>%
      mutate( Ts = s + exp( (logit_inv(uniform_rand) - (cbind(rep(1,numberOfObservations), Xdw_true) %*% parameters)) / (div)), # Equation (10) in the paper
              censoringtime = s + rexp(numberOfObservations, censor_strength),
              time = pmin(Ts, censoringtime),
              status = as.numeric(Ts <= censoringtime),
              X_orginal = X_orginal)
  }


# Check the censoring rate
test_censoring_AR = function(s = s, n, censor_strength){
  res = 0
  for (i in seq(100)){
    data =generateData_AR( s = s, numberOfObservations = n, numhv= numtrue, numtrue=numtrue, degree=degree, jmax=jmax, parameters = beta_true, censor_strength = censor_strength) # 3, 500, 0.7
    res = res + length(which(data$status==0))/n
  }
  return(res/100)
}

compute_ic_margins <- function(X, beta, eps = 1e-8) {
  N <- nrow(X)
  p <- ncol(X)
  S <- which(beta != 0)
  Sc <- setdiff(1:p, S)

  if (length(S) == 0) {
    return(list(
      margin_specific = NA,
      margin_worst    = NA,
      t_specific      = NA,
      t_worst         = NA
    ))
  }

  # Linear predictor and probabilities
  eta <- X %*% beta
  p_hat <- 1 / (1 + exp(-eta))
  W <- p_hat * (1 - p_hat)

  # Weighted Hessian
  XW <- sweep(X, 1, W, "*")       # Multiply each row by its corresponding weight
  Sigma <- t(X) %*% XW / N

  # Partition Sigma
  Sigma_SS <- Sigma[S, S, drop = FALSE]
  Sigma_ScS <- Sigma[Sc, S, drop = FALSE]

  # Stable inverse
  inv_Sigma_SS <- tryCatch(
    solve(Sigma_SS),
    error = function(e) {
      solve(Sigma_SS + eps * diag(length(S)))
    }
  )

  A <- Sigma_ScS %*% inv_Sigma_SS
  sgn <- sign(beta[S])
  t_specific <- max(abs(A %*% sgn))
  t_worst <- max(rowSums(abs(A)))

  return(list(
    margin_specific = 1 - t_specific,
    margin_worst    = 1 - t_worst,
    t_specific      = t_specific,
    t_worst         = t_worst
  ))
}

compute_logit_snr <- function(X, beta_true) {
  # Linear predictor
  eta <- as.numeric(X %*% beta_true)

  # Predicted probabilities
  # p_hat <- 1 / (1 + exp(-eta))

  # Signal variance
  signal_var <- var(eta)

  # Noise variance: mean Bernoulli variance
  noise_var <-  pi^2 / 3 # mean(p_hat * (1 - p_hat))

  # SNR
  snr <- signal_var / noise_var
  return(snr)
}

calculate_theoretical_logit_lasso_margin <- function(X, b_true) {

  # 0. Setup and Checks
  n <- nrow(X)
  p_full <- ncol(X) # p_full = p (covariates) + 1 (intercept)

  if (length(b_true) != p_full) {
    stop(paste("Error: b_true length (", length(b_true),
               ") must match the number of columns in X (", p_full, ")."))
  }

  # 1. Define Indices and Partition Betas
  covariate_indices <- 2:p_full
  beta_covariates <- b_true[covariate_indices]

  # 2. Identify the True Active Set (S) and Inactive Set (S_c) among covariates
  active_set_indices <- which(abs(beta_covariates) > 1e-6)
  inactive_set_indices <- which(abs(beta_covariates) <= 1e-6)

  if (length(active_set_indices) == 0 || length(inactive_set_indices) == 0) {
    warning("Active or Inactive set is empty among covariates. Margin is trivially 1.0.")
    return(1.0)
  }

  # 3. Calculate the full Fisher information matrix (I_full)

  # 3.1. Calculate the true linear predictor (eta_0) and probabilities (pi_0)
  eta_0 <- X %*% b_true
  pi_0 <- as.vector(1 / (1 + exp(-eta_0)))

  # 3.2. Calculate the diagonal weight matrix W_0 (n x n)
  W_0 <- diag(pi_0 * (1 - pi_0))

  # 3.3. Calculate I_full ((p+1) x (p+1))
  # This is the line that defines I_full!
  I_full <- t(X) %*% W_0 %*% X / n

  # 4. Partition the Fisher information matrix (remove the intercept)

  # Drop the first row and column (index 1), which correspond to the intercept.
  I_covariates <- I_full[covariate_indices, covariate_indices]

  # 5. Partition I_covariates into I_11 and I_21
  I_11 <- I_covariates[active_set_indices, active_set_indices]
  I_21 <- I_covariates[inactive_set_indices, active_set_indices]

  # 6. Calculate the maximum-correlation bound using the pseudoinverse (ginv)

  # Use ginv(), which handles singular matrices
  I_11_inverse <- ginv(I_11)

  # Get the sign vector of the active covariate coefficients
  sign_beta_1 <- sign(beta_covariates[active_set_indices])

  # Calculate the critical term (I_21 * I_11^-1 * sign(beta_1))
  bound_vector <- I_21 %*% I_11_inverse %*% sign_beta_1

  # Take the L-infinity norm
  max_correlation_bound <- max(abs(bound_vector))

  # 7. Calculate the Margin
  margin <- 1 - max_correlation_bound

  return(margin)
}

# Example usage (ensure that X_orginal and b_true are defined):
# margin_value <- calculate_theoretical_logit_lasso_margin_robust(X_orginal, b_true)
# Initial settings
s = 6

# Quarterly frequency
high_frequency = 4

# Polynomial degree
degree = 2
jmax <- s * high_frequency

# Number of covariates, K
numhv = 50

# Number of relevant covariates
numtrue = 2

# Gegenbauer polynomials W: L = 3 implies degree = 2
w <- gb(degree = degree, alpha = -1/2, a = 0, b = 1, jmax = jmax)/jmax

# Underlying weight function
w1 <- dbeta((1:jmax)/(jmax), shape1 = 1, shape2 = 3)
w2 <- dbeta((1:jmax)/(jmax), shape1 = 2, shape2 = 3)

# True parameters, including the intercept; see Example 5.1 in the paper
beta_true = c(1, 1, -1)

# Group structure
gindex = NULL
for (z in seq(numhv)){
  gindex <- c(gindex, rep(z, times = degree + 1))
}

# Starting point for the gradient descent algorithm
intercept_zero = 0

# Cross-validation
alpha = c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1)

# Censoring-strength parameter gamma, which yields approximately 81% censoring
censor = 1.68
test_censoring_AR( s = s, n = 1000, censor_strength = censor)

# Generate the prediction horizon t
# data = generateData_AR( s = s, n = 1200, numhv = numhv, numtrue = numtrue, degree = degree, jmax = jmax, parameters = beta_true, censor_strength = censor) # 3, 500, 0.7
# ind = which(data$Ts <=  data$censoringtime)
# t = quantile(data$Ts[ind], probs = t_quan)
t_quantile = 0
for (k in 1:100) {
  set.seed(k)
  print(k)
  data = generateData_AR( s = s, n = 1200, numhv = numhv, numtrue = numtrue, degree = degree, jmax = jmax, parameters = beta_true, censor_strength = censor) # 3, 500, 0.7
  ind = which(data$Ts <=  data$censoringtime)
  t_quantile = t_quantile + c(quantile(data$Ts[ind], probs = 0.1), quantile(data$Ts[ind], probs = 0.3), quantile(data$Ts[ind], probs = 0.5))
}


#
t_quantile = t_quantile/100 # c(quantile(data$Ts[ind], probs = 0.1), quantile(data$Ts[ind], probs = 0.3), quantile(data$Ts[ind], probs = 0.5))
censor_strength = c(censor)
t = t_quantile[1]

# Set up parallel processing
num_cores <- detectCores()
cl <- makeCluster(num_cores-7)
registerDoParallel(cl)
packages_to_export <- c('MASS' ,'timeROC', 'mvtnorm', 'Survivalml', "dplyr", "survival", "MLmetrics", "pROC", "PRROC", 'survivalROC', 'dotCall64', 'glmnet', 'rlang', 'pracma')

# Simulation repetitions
it = 100

# Scenario 1: N = 800
for (censor in censor_strength){
  n = 800
  i=0
  pa = matrix(0, nrow = length(t_quantile), ncol = 11)
  for (t in t_quantile){
    i = i + 1
    it = it

    result = foreach(k = 1:it, .errorhandling= 'remove', .packages = packages_to_export, .export = export_env) %dopar% {
      set.seed(k)

      # Data preparation
      data = generateData_AR( s = s, n = n, numhv = numhv, numtrue = numtrue, degree = degree, jmax = jmax, parameters = beta_true, censor_strength = censor) # 3, 500, 0.7
      dataset_p = data[ which( (1*( data$time <= t ))*(1*( data$status == 1 )) == 1), ]
      dataset_n = data[ which( (1*( data$time <= t ))*(1*( data$status == 1 )) == 0), ]
      indices_p <- sample( nrow( dataset_p ), nrow( dataset_p ) * 0.8 )
      indices_n = sample( nrow( dataset_n ), nrow( dataset_n ) * 0.8 )
      train_dataset_p <- dataset_p[ indices_p, ]
      train_dataset_n <- dataset_n[ indices_n, ]
      train_dataset = rbind( train_dataset_p, train_dataset_n )
      test_dataset_p <- dataset_p[ -indices_p, ]
      test_dataset_n <- dataset_n[ -indices_n, ]
      test_dataset = rbind(test_dataset_p, test_dataset_n)
      train_dataset = testRandomLogitDataset( train_dataset, t = t )
      train_dataset = KM_estimate(train_dataset)
      test_dataset = testRandomLogitDataset( test_dataset, t = t )
      test_dataset = KM_estimate(test_dataset)

      #5-fold cross-validation
      nfold = 5
      index_fold = which(1*(train_dataset$time <= t) * 1*(train_dataset$status == 1) == 1)
      if (length(index_fold) < 5){
        foldid = NULL
      }else{
        foldid = c( form_folds(nrow(train_dataset[index_fold,]), nfold), form_folds(nrow(train_dataset[-index_fold,]), nfold) )
      }

      # Dataset for LASSO-UMIDAS
      X_orginal_train = train_dataset$X_orginal[ , 1:(1+jmax*numhv) ]
      X_orginal_test = test_dataset$X_orginal[ , 1:(1+jmax*numhv) ]

      # Dataset for sg-LASSO-MIDAS
      X_train = train_dataset[ , 1:(numhv*(degree+1)+1) ]
      y_train = 1*(train_dataset$time <= t)
      X_test = as.matrix( test_dataset[, 1:(numhv*(degree+1)+1)] )
      y_test = 1*( test_dataset$time <= t )

      X = rbind(X_train, X_test)
      X_orginal = rbind(X_orginal_train, X_orginal_test)

      # Oracle AUC
      b_truec = c(beta_true[1] + log((t - s)), beta_true[-1]  +  log((t - s)))
      b_truec = c(b_truec[1:(1)], w1*b_truec[1+1], w2*b_truec[1+1+1])
      b_true = c(b_truec, rep(0,(numhv-2)*jmax))
      margin_ic = calculate_theoretical_logit_lasso_margin(X_orginal, b_true)
      snr = compute_logit_snr(X_orginal, b_true)
      # calculate_theoretical_logit_lasso_margin(X_orginal, b_true)
      var_MIDAS = apply(X[,2:7],2,var)

      results = list(var = var_MIDAS, margin_ic = margin_ic, snr = snr)
      return(results)
    }

    var_scenario1 = c()
    for (k in 1:it) {
      var_scenario1 = rbind(var_scenario1, result[[k]]$var)
    }
    v1 = apply(var_scenario1, 2, mean)

    IC_scenario1 = c()
    for (k in 1:it) {
      IC_scenario1 = rbind(IC_scenario1, result[[k]]$margin_ic)
    }
    snr_s = c()
    for (k in 1:it) {
      snr_s = rbind(snr_s, result[[k]]$snr)
    }
    pa[i, ] = c('variance ratio', v1, 'IC_margin', mean(IC_scenario1), 'SNR', mean(snr_s))
    print(pa[i, ])
  }
  table_name <- paste0("LASSOUMIDAS_scenario1_SNR", censor, 800, s,  ".csv")
  write.table(pa, table_name, sep = ";", row.names = FALSE, quote = FALSE)
}

