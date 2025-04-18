#============================================================
# Corporate Survival Analysis with Machine Learning Methods
# Import functions for simulation and empirical application
#============================================================


#============================================================
# Correctness check, weighted logistic regression without penaly and MIDAS
# The goal is to check the MSE of estimated parameters
#============================================================



rm(list = ls())

source('import functions for empirical application.R')



#load necessary packages
library(midasml)
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

# data generating process, scenario 1
generateData_AR <-
  function(s = s, numberOfObservations, numhv, numtrue, degree, jmax, parameters = NA, censor_strength){
    TN = numberOfObservations
    p <- numhv # Set the dimension of the matrix
    degree <- degree # degree of the polynomials
    jmax <- jmax # number of hf obs per lf obs
    # quarterly covariates
    Xdd <- array(NA, dim = c(TN, jmax, numhv))
    # storage for Legendre weighted quarterly
    Xdw <- matrix(NA, nrow = TN, ncol = (degree+1)*numhv)
    # storage for beta weighted quarterly
    Xdw_true <- matrix(NA, nrow = TN, ncol = numtrue)
    # storage for quarterly
    Xdd_g <- array(NA, dim = c(TN, jmax, numhv))
    # storage for Legendre weighted quarterly
    Xdw_g <- matrix(NA, nrow = TN, ncol = (degree+1)*numhv)
    # storage for beta weighted quarterly
    Xdw_true_g <- matrix(NA, nrow = TN, ncol = numtrue)
    ##########################################################

    Xd <- matrix(0, nrow = TN*jmax, ncol = p)

    # level of cross-sectional dependence across all variables
    phi = 0.1

    # degree of time series dependence among its lag
    rho = 0.1
    corr_matrix <- diag(1, p)

    for (v in 1:(p)) {
      for (u in 1:(p)) {
        corr_matrix[v, u] <- (phi)^abs(v - u)
      }
    }

    for ( t in 1:(jmax)){
      if (t == 1){
        Xdd[,t,(1:p)] <- MASS::mvrnorm(TN, mu = rep(0, p), corr_matrix*(1 - rho^2)) #scenario 1
      }
      else{
        Xdd[,t,(1:p)] <- rho * Xdd[,t-1,(1:p)] + MASS::mvrnorm(TN, mu = rep(0, p), corr_matrix*(1 - rho^2))
      }
    }

    Xdd = abs(Xdd)


    # storage for all lags, which is ued for LASSO-MIDAS
    X_orginal = matrix(1, nrow = TN)
    for (i in seq(numhv)){
      X_orginal = cbind(X_orginal, Xdd[,,i])
    }

    # specify the underlying weight function w1 and w2, which represent Beta(1,3) and Beta(2,3) respectively.
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

    # Inverse logit function
    logit_inv <- function(p) log( p / (1 - p))  # Inverse logit function
    uniform_rand <- runif(numberOfObservations)

    # X^T%*%l
    div = cbind(rep(1, numberOfObservations) ,  Xdw_true) %*% c(1, rep(1, length(parameters) - 1))

    dataset %>%
      as_tibble() %>%
      mutate( Ts = s + exp( (logit_inv(uniform_rand) - (cbind(rep(1,numberOfObservations), Xdw_true) %*% parameters)) / (div)), # equation (10) in the paper
              censoringtime = s + rexp(numberOfObservations, censor_strength),
              time = pmin(Ts, censoringtime),
              status = as.numeric(Ts <= censoringtime),
              X_orginal = X_orginal)
  }

# grab the prediction metric
paral_independent = function(result,it){
  res_fit_MIDAS = NULL
  res_fit_MIDAS_LASSO = NULL
  res_fit_LASSO = NULL
  AUC_true_N = NULL
  AUC_LASSO_N = NULL
  AUC_MIDAS_N = NULL
  AUC_MIDAS_LASSO_N = NULL
  balance = NULL
  censoring_proportions = NULL
  it = length(result)
  for (i in seq(it)){
    res_fit_MIDAS= rbind(res_fit_MIDAS , result[[i]]$fit_MIDAS)
    res_fit_MIDAS_LASSO= rbind(res_fit_MIDAS_LASSO , result[[i]]$fit_MIDAS_LASSO)
    res_fit_LASSO= rbind(res_fit_LASSO , result[[i]]$fit_LASSO)
    AUC_true_N = rbind(AUC_true_N , as.numeric( result[[i]]$AUC_true_N ))
    AUC_MIDAS_N = rbind(AUC_MIDAS_N , as.numeric( result[[i]]$AUC_MIDAS_N ))
    AUC_MIDAS_LASSO_N = rbind(AUC_MIDAS_LASSO_N , as.numeric( result[[i]]$AUC_MIDAS_LASSO_N ))
    AUC_LASSO_N = rbind(AUC_LASSO_N , as.numeric( result[[i]]$AUC_LASSO_N ))
    censoring_proportions = rbind(censoring_proportions,  as.numeric( result[[i]]$censoring_proportions ))
    balance  = rbind(balance , as.numeric( result[[i]]$balance ))
  }
  res_fit_av_MIDAS = colMeans(res_fit_MIDAS)
  res_fit_av_LASSO = colMeans(res_fit_LASSO)
  res_fit_av_MIDAS_LASSO = colMeans(res_fit_MIDAS_LASSO)
  res_fit_var_LASSO = mean(apply(res_fit_LASSO, 2, var))
  b_true = c(beta_true[1] + log((t - s)), beta_true[-1]  +  log((t - s)))
  ##########################################
  w1_tr = w1 * b_true[1+1]
  w1_es = w%*%as.numeric(res_fit_av_MIDAS[(1+1):(1+1+degree)])
  w1_es_MIDAS_LASSO = w%*%as.numeric(res_fit_av_MIDAS_LASSO[(1+1):(1+1+degree)])
  w1_es_LASSO = as.numeric(res_fit_av_LASSO[(1+1):(1+1+jmax-1)])
  w1_es_var = round(mean(apply(as.matrix(res_fit_MIDAS[, (1  + 1):(1  + 1 + degree)])%*%t(w), 2, var)), 5)
  w1_es_MSE = round(mean(as.numeric((w1_tr - w1_es)^2) + apply(as.matrix(res_fit_MIDAS[, (1  + 1):(1  + 1 + degree)])%*%t(w), 2, var)), 5)
  w1_es_var_LASSO = round(mean(apply(as.matrix(res_fit_LASSO[, (1  + 1):(1  + 1 + jmax - 1)]), 2, var)), 5)
  w1_es_MSE_LASSO = round(mean(as.numeric((w1_tr - w1_es_LASSO)^2) + apply(as.matrix(res_fit_LASSO[, (1  + 1):(1  + 1 + jmax - 1)]), 2, var)), 5)
  w1_es_var_MIDAS_LASSO = round(mean(apply(as.matrix(res_fit_MIDAS_LASSO[, (1  + 1):(1  + 1 + degree)])%*%t(w), 2, var)), 5)
  w1_es_MSE_MIDAS_LASSO = round(mean(as.numeric((w1_tr - w1_es_MIDAS_LASSO)^2) + apply(as.matrix(res_fit_MIDAS_LASSO[, (1  + 1):(1  + 1 + degree)])%*%t(w), 2, var)), 5)
  #########################################
  w2_tr = w2 * b_true[1+1+1]
  w2_es = w%*%as.numeric(res_fit_av_MIDAS[(1+1+degree+1):(1+1+degree+1+degree)])
  w2_es_LASSO = as.numeric(res_fit_av_LASSO[(1+1+jmax):(1+1+jmax*2-1)])
  w2_es_MIDAS_LASSO = w%*%as.numeric(res_fit_av_MIDAS_LASSO[(1+1+degree+1):(1+1+degree+1+degree)])
  w2_es_var = round(mean(apply(as.matrix(res_fit_MIDAS[, (1+1+degree+1):(1+1+degree+1+degree)])%*%t(w), 2, var)), 5)
  w2_es_var_LASSO = round(mean(apply(as.matrix(res_fit_LASSO[, (1+1+jmax):(1+1+jmax*2-1)]), 2, var)), 5)
  w2_es_var_MIDAS_LASSO = round(mean(apply(as.matrix(res_fit_MIDAS_LASSO[, (1+1+degree+1):(1+1+degree+1+degree)])%*%t(w), 2, var)), 5)
  w2_es_MSE = round(mean(as.numeric((w2_tr - w2_es)^2) + apply(as.matrix(res_fit_MIDAS[, (1+1+degree+1):(1+1+degree+1+degree)])%*%t(w), 2, var)), 5)
  w2_es_MSE_LASSO = round(mean(as.numeric((w2_tr - w2_es_LASSO)^2) + apply(as.matrix(res_fit_LASSO[, (1+1+jmax):(1+1+jmax*2-1)]), 2, var)), 5)
  w2_es_MSE_MIDAS_LASSO = round(mean(as.numeric((w2_tr - w2_es_MIDAS_LASSO)^2) + apply(as.matrix(res_fit_MIDAS_LASSO[, (1+1+degree+1):(1+1+degree+1+degree)])%*%t(w), 2, var)), 5)
  ############################################################
  #########################################################################
  res = list(AUC_true_N = round(colMeans(AUC_true_N),3), AUC_true_N_VAR = round(apply(AUC_true_N, 2, var),3),  AUC_MIDAS_N = round(colMeans(AUC_MIDAS_N),3), AUC_MIDAS_N_VAR = round(apply(AUC_MIDAS_N, 2, var),3), AUC_LASSO_N = round(colMeans(AUC_LASSO_N),3), AUC_LASSO_N_VAR = round(apply(AUC_LASSO_N, 2, var),3), AUC_MIDAS_LASSO_N = round(colMeans(AUC_MIDAS_LASSO_N),3), AUC_MIDAS_LASSO_N_VAR = round(apply(AUC_MIDAS_LASSO_N, 2, var),3),
             w1_es_MSE = w1_es_MSE, w1_es_var = w1_es_var,
             w2_es_MSE = w2_es_MSE, w2_es_var = w2_es_var,
             w1_es_MSE_LASSO = w1_es_MSE_LASSO, w1_es_var_LASSO = w1_es_var_LASSO,
             w2_es_MSE_LASSO = w2_es_MSE_LASSO, w2_es_var_LASSO = w2_es_var_LASSO,
             w1_es_MSE_MIDAS_LASSO = w1_es_MSE_MIDAS_LASSO, w1_es_var_MIDAS_LASSO = w1_es_var_MIDAS_LASSO,
             w2_es_MSE_MIDAS_LASSO = w2_es_MSE_MIDAS_LASSO, w2_es_var_MIDAS_LASSO = w2_es_var_MIDAS_LASSO,
             censoring_proportions = colMeans(censoring_proportions), balance = colMeans(balance))
  return(res)
}

test_censoring_AR = function(s = s, n, censor_strength){
  res = 0
  for (i in seq(100)){
    data =generateData_AR( s = s, numberOfObservations = n, numhv= numtrue, numtrue=numtrue, degree=degree, jmax=jmax, parameters = beta_true, censor_strength = censor_strength) # 3, 500, 0.7
    res = res + length(which(data$status==0))/n
  }
  return(res/100)
}

# initial setting
s = 1

#quratker/year frequency
high_frequency = 2

# degree of the polynomials
degree = 2
jmax <- s * high_frequency

# Number of covariates K
numhv = 1

# Number of relevant covariates
numtrue = 1

# without using MIDAS
w = rep(1, jmax)
w1 = rep(1, jmax)
w2 = rep(1, jmax)

#true parameters including the intercept, see Example 5.1 in the paper
beta_true = c(1, 1)

#group structure
gindex = NULL
for (z in seq(numhv)){
  gindex <- c(gindex, rep(z, times = degree + 1))
}

# starting point of the gradient decent algorithm
intercept_zero = 0

#censoring strength gamma, which can generate approximately 81% censoring
censor = 0.15
test_censoring_AR( s = s, n = 1000, censor_strength = censor)

# Generate the prediction horizon t
set.seed(1)
data = generateData_AR( s = s, n = 1200, numhv = numhv, numtrue = numtrue, degree = degree, jmax = jmax, parameters = beta_true, censor_strength = censor) # 3, 500, 0.7
ind = which(data$Ts <=  data$censoringtime)
t_quan = 0.5
t = quantile(data$Ts[ind], probs = t_quan)

#call for parallel, notice that the number of cores is dependent on your own laptop
stopCluster(cl)
num_cores <- detectCores()
cl <- makeCluster(num_cores-7)
registerDoParallel(cl)
packages_to_export <- c('timeROC', 'mvtnorm', 'Survivalml', "dplyr", "midasml", "survival", "MLmetrics", "pROC", "PRROC", 'survivalROC', 'dotCall64', 'glmnet', 'rlang', 'pracma')

# repetitions of the simulation
it = 500
for (n in c(100, 200, 500)) {
result = foreach(k = 1:it, .packages = packages_to_export, .export = export_env) %dopar% {
      set.seed(k*t)
      # data prepare
      data = generateData_AR( s = s, n = n, numhv = numhv, numtrue = numtrue, degree = degree, jmax = jmax, parameters = beta_true, censor_strength = censor) # 3, 500, 0.7
      data_in = testRandomLogitDataset( data, t = t )
      data_in = KM_estimate(data_in)

      #dataset for LASSO-UMIDAS
      X_orginal = data_in$X_orginal[ , 1:(1+jmax*numhv) ]
      y = 1*(data_in$time <= t)
      # IPW weights
      w_in = 1*( pmin(data_in$Ts, t) <= data_in$censoringtime)/data_in$Gp

      #LASSO-UMIDAS
      fit_cv_LASSO = survival_sparsegl(as.matrix(X_orginal[,-1]), y, group = seq(jmax*numhv), weight = w_in, lambda = 0,  asparse = 1, intercept_zero = 0, standardize = TRUE)
      est_cof = fit_cv_LASSO$beta
      results = list( est_cof = est_cof)
      return(results)
    }

b_truec = c(beta_true[1] + log((t - s)), beta_true[-1]  +  log((t - s)))
b_truec = c(b_truec[1:(1)], w1*b_truec[1+1])
est_all = 0
for (i in seq(it)) {
  est_all = est_all + result[[i]]$est_cof
}
est = est_all/it
mse = mean((est- b_truec[-1])^2)
print(mse)
}


