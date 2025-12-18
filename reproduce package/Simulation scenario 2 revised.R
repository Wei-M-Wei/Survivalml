#============================================================
# Corporate Survival Analysis with Machine Learning Methods
# Import functions for simulation and empirical application
#============================================================


#============================================================
# Reproduce simulation
#============================================================



rm(list = ls())

source('import functions for empirical application.R')



#load necessary packages
#library(midasml)
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

# data generating process, scenario 2
generateData_AR <-
  function(s, n, numhv, numtrue, degree, jmax, parameters = NA, censor_strength){
    TN = n
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
    rho = 0.9
    corr_matrix <- diag(1, p)

    for (v in 1:(p)) {
      for (u in 1:(p)) {
        corr_matrix[v, u] <- (phi)^abs(v - u)
      }
    }

    for ( t in 1:(jmax)){
      if (t == 1){
        Xdd[,t,(1:p)] <- MASS::mvrnorm(TN, mu = rep(0, p), corr_matrix) #scenario 2
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
    dataset = cbind(rep(1, n),  Xdw)

    # Inverse logit function
    logit_inv <- function(p) log( p / (1 - p))  # Inverse logit function
    uniform_rand <- runif(n)

    # X^T%*%l
    div = cbind(rep(1, n) ,  Xdw_true) %*% c(1, rep(1, length(parameters) - 1))

    dataset %>%
      as_tibble() %>%
      mutate( Ts = s + exp( (logit_inv(uniform_rand) - (cbind(rep(1,n), Xdw_true) %*% parameters)) / (div)), # equation (10) in the paper
              censoringtime = s + rexp(n, censor_strength),
              time = pmin(Ts, censoringtime),
              status = as.numeric(Ts <= censoringtime),
              X_orginal = X_orginal)
  }


#check the censoring rate
test_censoring_AR = function(s, n, censor_strength){
  res = 0
  for (i in seq(100)){
    data =generateData_AR( s = s, n = n, numhv= numtrue, numtrue=numtrue, degree=degree, jmax=jmax, parameters = beta_true, censor_strength = censor_strength)
    res = res + length(which(data$status==0))/n
  }
  return(res/100)
}


# initial setting
s = 6

#quratker/year frequency
high_frequency = 4

# degree of the polynomials
degree = 2
jmax <- s * high_frequency

# Number of covariates K
numhv = 50

# Number of relevant covariates
numtrue = 2

# Gegenbauer polynomials W, L = 3 <=> degree is set to 2
w <- gb(degree = degree, alpha = -1/2, a = 0, b = 1, jmax = jmax)/jmax

# Underlying weight function
w1 <- dbeta((1:jmax)/(jmax), shape1 = 1, shape2 = 3)
w2 <- dbeta((1:jmax)/(jmax), shape1 = 2, shape2 = 3)

#true parameters including the intercept, see Example 5.1 in the paper
beta_true = c(1, 1, -1)

#group structure
gindex = NULL
for (z in seq(numhv)){
  gindex <- c(gindex, rep(z, times = degree + 1))
}

# starting point of the gradient decent algorithm
intercept_zero = 0

# cross-validation
alpha = c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1)

#censoring strength gamma, which can generate approximately 81% censoring
censor = 1.7
test_censoring_AR( s = s, n = 1000, censor_strength = censor)

# Generate the prediction horizon t
t_quantile = 0
for (k in 1:100) {
  set.seed(k)
  print(k)
  data = generateData_AR( s = s, n = 1200, numhv = numhv, numtrue = numtrue, degree = degree, jmax = jmax, parameters = beta_true, censor_strength = censor) 
  ind = which(data$Ts <=  data$censoringtime)
  t_quantile = t_quantile + c(quantile(data$Ts[ind], probs = 0.1), quantile(data$Ts[ind], probs = 0.3), quantile(data$Ts[ind], probs = 0.5))
}


#
t_quantile = t_quantile/100 # c(quantile(data$Ts[ind], probs = 0.1), quantile(data$Ts[ind], probs = 0.3), quantile(data$Ts[ind], probs = 0.5))
censor_strength = c(censor)

#call for parallel
stopCluster(cl)
num_cores <- detectCores()
cl <- makeCluster(num_cores-22)
registerDoParallel(cl)
packages_to_export <- c('timeROC', 'mvtnorm', 'Survivalml', "dplyr", "survival", "MLmetrics", "pROC", "PRROC", 'survivalROC', 'dotCall64', 'glmnet', 'rlang', 'pracma')

# repetitions of the simulation
it = 100

#scenario 2 N = 800
for (censor in censor_strength){
  n = 800
  i=0
  pa = matrix(0, nrow = length(t_quantile), ncol = 43)
  for (t in t_quantile){
    i = i + 1
    it = it

    result = foreach(k = 1:it, .errorhandling= 'remove', .packages = packages_to_export, .export = export_env) %dopar% {
      set.seed(s*k)

      # data prepare
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

      #5-fold cross validation
      nfold = 5
      index_fold = which(1*(train_dataset$time <= t) * 1*(train_dataset$status == 1) == 1)
      if (length(index_fold) < 5){
        foldid = NULL
      }else{
        foldid = c( form_folds(nrow(train_dataset[index_fold,]), nfold), form_folds(nrow(train_dataset[-index_fold,]), nfold) )
      }

      #dataset for LASSO-UMIDAS
      X_orginal_train = train_dataset$X_orginal[ , 1:(1+jmax*numhv) ]
      X_orginal_test = test_dataset$X_orginal[ , 1:(1+jmax*numhv) ]

      #dataset for sg-LASSO-MIDAS
      X_train = train_dataset[ , 1:(numhv*(degree+1)+1) ]
      y_train = 1*(train_dataset$time <= t)
      X_test = as.matrix( test_dataset[, 1:(numhv*(degree+1)+1)] )
      y_test = 1*( test_dataset$time <= t )

      # IPW weights
      w_train = 1*( pmin(train_dataset$Ts, t) <= train_dataset$censoringtime)/train_dataset$Gp

      # sg-LASSO-MIDAS, cross-validation for AUC
      fit_cv_MIDAS = alpha_cv_sparsegl(as.matrix(X_train[,-1]), y_train, group = gindex, nlambda = 50, weight = w_train, alpha = alpha, foldid = foldid, nfolds = nfold,  pred.loss = 'censor', intercept_zero = intercept_zero, standardize = TRUE, AUC = TRUE, data = train_dataset, t = t)
      test_predictions_MIDAS = unlist( as.list(plogis(cbind(1, X_test[,-1]) %*% fit_cv_MIDAS$coff_AUC)))
      est_MIDAS = fit_cv_MIDAS$coff_AUC

      #LASSO-MIDAS
      if (fit_cv_MIDAS$alpha_out == 1){
        fit_cv_MIDAS_LASSO = fit_cv_MIDAS
        test_predictions_MIDAS_LASSO = test_predictions_MIDAS
        est_MIDAS_LASSO = est_MIDAS
      } else {
        fit_cv_MIDAS_LASSO = fit_cv_MIDAS$fit_MIDAS_LASSO
        cesnor_min_index = which(fit_cv_MIDAS_LASSO$AUC_censor == max(fit_cv_MIDAS_LASSO$AUC_censor))[1]
        est_MIDAS_LASSO = unlist(c(fit_cv_MIDAS_LASSO$survival_sparsegl.fit$b0[,cesnor_min_index], fit_cv_MIDAS_LASSO$survival_sparsegl.fit$beta[,cesnor_min_index]))
        test_predictions_MIDAS_LASSO = unlist( as.list(plogis(cbind(1, X_test[,-1]) %*% est_MIDAS_LASSO)))
      }

      #LASSO-UMIDAS
      fit_cv_LASSO = cv.survival_sparsegl(as.matrix(X_orginal_train[,-1]), y_train, group = seq(jmax*numhv), nlambda = 50, weight = w_train, asparse = 1, foldid = foldid, nfolds = nfold, pred.loss = 'censor', intercept_zero = intercept_zero, standardize = TRUE, AUC = TRUE, data = train_dataset, t = t, maxit = 30000)
      cesnor_min_index_LASSO = which(fit_cv_LASSO$AUC_censor == max(fit_cv_LASSO$AUC_censor))[1]
      est_LASSO = unlist(c(fit_cv_LASSO$survival_sparsegl$b0[,cesnor_min_index_LASSO], fit_cv_LASSO$survival_sparsegl$beta[,cesnor_min_index_LASSO]))
      test_predictions_LASSO = unlist( as.list(plogis(cbind(1, X_orginal_test[,-1]) %*% est_LASSO)))

      # oracle AUC
      b_truec = c(beta_true[1] + log((t - s)), beta_true[-1]  +  log((t - s)))
      b_truec = c(b_truec[1:(1)], w1*b_truec[1+1], w2*b_truec[1+1+1])
      predictions_true = unlist(as.list(plogis(X_orginal_test[, 1: (1+jmax*numtrue)] %*% b_truec)))
      AUC_true_N_t = survivalROC(Stime=test_dataset$time, status= test_dataset$status, marker = predictions_true, predict.time = t, span = 0.25*nrow(data)^(-1/2))
      AUC_true_N = AUC_true_N_t$AUC

      # prediction for 3 approaches
      AUC_MIDAS_N = ROC_censor_N(data = test_dataset, prediction = test_predictions_MIDAS, t = t )$AUC
      AUC_MIDAS_LASSO_N = ROC_censor_N(data = test_dataset, prediction = test_predictions_MIDAS_LASSO, t = t )$AUC
      AUC_LASSO_N = ROC_censor_N(data = test_dataset, prediction = test_predictions_LASSO, t = t )$AUC

      balance = length(data$Ts[data$Ts <= t])/n
      censoring_proportions = length(which(data$status==0))/n

      results = list( AUC_true_N = round( mean(AUC_true_N), 3 ), AUC_LASSO_N = round( mean(AUC_LASSO_N), 3 ),
                      AUC_MIDAS_N = round(mean(AUC_MIDAS_N), 3 ), AUC_MIDAS_LASSO_N = round(mean(AUC_MIDAS_LASSO_N), 3 ),
                      censoring_proportions = censoring_proportions, t = t, balance = round( mean(balance), 3 ),
                      fit_LASSO = est_LASSO, fit_MIDAS = est_MIDAS, fit_MIDAS_LASSO = est_MIDAS_LASSO)
      return(results)
    }

    res = paral_independent(result, length(result))
    pa[i, ] = c(res$censoring_proportions, res$balance,
                'AUC',
                res$AUC_true_N, res$AUC_LASSO_N, res$AUC_MIDAS_LASSO_N, res$AUC_MIDAS_N,
                res$AUC_true_N_sd, res$AUC_LASSO_N_sd, res$AUC_MIDAS_LASSO_N_sd, res$AUC_MIDAS_N_sd,
                'w1',
                round(res$w1_es_MSE_LASSO,3), round(res$w1_es_MSE_MIDAS_LASSO,3), round(res$w1_es_MSE,3),
                round(res$w1_es_sd_LASSO,3), round(res$w1_es_sd_MIDAS_LASSO,3), round(res$w1_es_sd,3),
                'w2',
                round(res$w2_es_MSE_LASSO,3), round(res$w2_es_MSE_MIDAS_LASSO,3), round(res$w2_es_MSE,3),
                round(res$w2_es_sd_LASSO,3), round(res$w2_es_sd_MIDAS_LASSO,3), round(res$w2_es_sd,3),
                'LASSO metric',
                round(res$metrics_LASSO$TPR, 3), round(res$metrics_LASSO$FPR, 3), round(res$metrics_LASSO$Accuracy, 3), 
                round(res$metrics_LASSO$Precision, 3), round(res$metrics_LASSO$F1, 3),
                'MIDAS_LASSO metric', 
                round(res$metrics_MIDAS_LASSO$TPR, 3), round(res$metrics_MIDAS_LASSO$FPR, 3), round(res$metrics_MIDAS_LASSO$Accuracy, 3), 
                round(res$metrics_MIDAS_LASSO$Precision, 3), round(res$metrics_MIDAS_LASSO$F1, 3),
                'MIDAS metric', 
                round(res$metrics_MIDAS$TPR, 3), round(res$metrics_MIDAS$FPR, 3), round(res$metrics_MIDAS$Accuracy, 3), 
                round(res$metrics_MIDAS$Precision, 3), round(res$metrics_MIDAS$F1, 3))
    print(pa[i, ])
  }
  table_name <- paste0("MIDASAUC_scenario2_new_revised", censor, 800, s,  ".csv")
  write.table(pa, table_name, sep = ";", row.names = FALSE, quote = FALSE)
}

#scenario 2 N = 1200
for (censor in censor_strength){
  n = 1200
  i=1
  pa = matrix(0, nrow = length(t_quantile), ncol = 43)
  for (t in t_quantile){
    i = i + 1
    it = it

    result = foreach(k = 1:it, .errorhandling= 'remove', .packages = packages_to_export, .export = export_env) %dopar% {
      set.seed(s*k)

      # data prepare
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

      #5-fold cross validation
      nfold = 5
      index_fold = which(1*(train_dataset$time <= t) * 1*(train_dataset$status == 1) == 1)
      if (length(index_fold) < 5){
        foldid = NULL
      }else{
        foldid = c( form_folds(nrow(train_dataset[index_fold,]), nfold), form_folds(nrow(train_dataset[-index_fold,]), nfold) )
      }

      #dataset for LASSO-UMIDAS
      X_orginal_train = train_dataset$X_orginal[ , 1:(1+jmax*numhv) ]
      X_orginal_test = test_dataset$X_orginal[ , 1:(1+jmax*numhv) ]

      #dataset for sg-LASSO-MIDAS
      X_train = train_dataset[ , 1:(numhv*(degree+1)+1) ]
      y_train = 1*(train_dataset$time <= t)
      X_test = as.matrix( test_dataset[, 1:(numhv*(degree+1)+1)] )
      y_test = 1*( test_dataset$time <= t )

      # IPW weights
      w_train = 1*( pmin(train_dataset$Ts, t) <= train_dataset$censoringtime)/train_dataset$Gp

      # sg-LASSO-MIDAS, cross-validation for AUC
      fit_cv_MIDAS = alpha_cv_sparsegl(as.matrix(X_train[,-1]), y_train, group = gindex, nlambda = 50, weight = w_train, alpha = alpha, foldid = foldid, nfolds = nfold,  pred.loss = 'censor', intercept_zero = intercept_zero, standardize = TRUE, AUC = TRUE, data = train_dataset, t = t)
      test_predictions_MIDAS = unlist( as.list(plogis(cbind(1, X_test[,-1]) %*% fit_cv_MIDAS$coff_AUC)))
      est_MIDAS = fit_cv_MIDAS$coff_AUC

      #LASSO-MIDAS
      if (fit_cv_MIDAS$alpha_out == 1){
        fit_cv_MIDAS_LASSO = fit_cv_MIDAS
        test_predictions_MIDAS_LASSO = test_predictions_MIDAS
        est_MIDAS_LASSO = est_MIDAS
      } else {
        fit_cv_MIDAS_LASSO = fit_cv_MIDAS$fit_MIDAS_LASSO
        cesnor_min_index = which(fit_cv_MIDAS_LASSO$AUC_censor == max(fit_cv_MIDAS_LASSO$AUC_censor))[1]
        est_MIDAS_LASSO = unlist(c(fit_cv_MIDAS_LASSO$survival_sparsegl.fit$b0[,cesnor_min_index], fit_cv_MIDAS_LASSO$survival_sparsegl.fit$beta[,cesnor_min_index]))
        test_predictions_MIDAS_LASSO = unlist( as.list(plogis(cbind(1, X_test[,-1]) %*% est_MIDAS_LASSO)))
      }

      #LASSO-UMIDAS
      fit_cv_LASSO = cv.survival_sparsegl(as.matrix(X_orginal_train[,-1]), y_train, group = seq(jmax*numhv), nlambda = 50, weight = w_train, asparse = 1, foldid = foldid, nfolds = nfold, pred.loss = 'censor', intercept_zero = intercept_zero, standardize = TRUE, AUC = TRUE, data = train_dataset, t = t, maxit = 30000)
      cesnor_min_index_LASSO = which(fit_cv_LASSO$AUC_censor == max(fit_cv_LASSO$AUC_censor))[1]
      est_LASSO = unlist(c(fit_cv_LASSO$survival_sparsegl$b0[,cesnor_min_index_LASSO], fit_cv_LASSO$survival_sparsegl$beta[,cesnor_min_index_LASSO]))
      test_predictions_LASSO = unlist( as.list(plogis(cbind(1, X_orginal_test[,-1]) %*% est_LASSO)))

      # oracle AUC
      b_truec = c(beta_true[1] + log((t - s)), beta_true[-1]  +  log((t - s)))
      b_truec = c(b_truec[1:(1)], w1*b_truec[1+1], w2*b_truec[1+1+1])
      predictions_true = unlist(as.list(plogis(X_orginal_test[, 1: (1+jmax*numtrue)] %*% b_truec)))
      AUC_true_N_t = survivalROC(Stime=test_dataset$time, status= test_dataset$status, marker = predictions_true, predict.time = t, span = 0.25*nrow(data)^(-1/2))
      AUC_true_N = AUC_true_N_t$AUC

      # prediction for 3 approaches
      AUC_MIDAS_N = ROC_censor_N(data = test_dataset, prediction = test_predictions_MIDAS, t = t )$AUC
      AUC_MIDAS_LASSO_N = ROC_censor_N(data = test_dataset, prediction = test_predictions_MIDAS_LASSO, t = t )$AUC
      AUC_LASSO_N = ROC_censor_N(data = test_dataset, prediction = test_predictions_LASSO, t = t )$AUC

      balance = length(data$Ts[data$Ts <= t])/n
      censoring_proportions = length(which(data$status==0))/n

      results = list( AUC_true_N = round( mean(AUC_true_N), 3 ), AUC_LASSO_N = round( mean(AUC_LASSO_N), 3 ),
                      AUC_MIDAS_N = round(mean(AUC_MIDAS_N), 3 ), AUC_MIDAS_LASSO_N = round(mean(AUC_MIDAS_LASSO_N), 3 ),
                      censoring_proportions = censoring_proportions, t = t, balance = round( mean(balance), 3 ),
                      fit_LASSO = est_LASSO, fit_MIDAS = est_MIDAS, fit_MIDAS_LASSO = est_MIDAS_LASSO)
      return(results)
    }

    res = paral_independent(result, length(result))
    pa[i, ] = c(res$censoring_proportions, res$balance,
                'AUC',
                res$AUC_true_N, res$AUC_LASSO_N, res$AUC_MIDAS_LASSO_N, res$AUC_MIDAS_N,
                res$AUC_true_N_sd, res$AUC_LASSO_N_sd, res$AUC_MIDAS_LASSO_N_sd, res$AUC_MIDAS_N_sd,
                'w1',
                round(res$w1_es_MSE_LASSO,3), round(res$w1_es_MSE_MIDAS_LASSO,3), round(res$w1_es_MSE,3),
                round(res$w1_es_sd_LASSO,3), round(res$w1_es_sd_MIDAS_LASSO,3), round(res$w1_es_sd,3),
                'w2',
                round(res$w2_es_MSE_LASSO,3), round(res$w2_es_MSE_MIDAS_LASSO,3), round(res$w2_es_MSE,3),
                round(res$w2_es_sd_LASSO,3), round(res$w2_es_sd_MIDAS_LASSO,3), round(res$w2_es_sd,3),
                'LASSO metric',
                round(res$metrics_LASSO$TPR, 3), round(res$metrics_LASSO$FPR, 3), round(res$metrics_LASSO$Accuracy, 3), 
                round(res$metrics_LASSO$Precision, 3), round(res$metrics_LASSO$F1, 3),
                'MIDAS_LASSO metric', 
                round(res$metrics_MIDAS_LASSO$TPR, 3), round(res$metrics_MIDAS_LASSO$FPR, 3), round(res$metrics_MIDAS_LASSO$Accuracy, 3), 
                round(res$metrics_MIDAS_LASSO$Precision, 3), round(res$metrics_MIDAS_LASSO$F1, 3),
                'MIDAS metric', 
                round(res$metrics_MIDAS$TPR, 3), round(res$metrics_MIDAS$FPR, 3), round(res$metrics_MIDAS$Accuracy, 3), 
                round(res$metrics_MIDAS$Precision, 3), round(res$metrics_MIDAS$F1, 3))
    print(pa[i, ])
  }
  table_name <- paste0("MIDASAUC_scenario2_new_revised", censor, 1200, s,  ".csv")
  write.table(pa, table_name, sep = ";", row.names = FALSE, quote = FALSE)
}
