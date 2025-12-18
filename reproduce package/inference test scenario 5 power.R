#============================================================
# Corporate Survival Analysis with Machine Learning Methods
# Import functions for simulation and empirical application
#============================================================


#============================================================
# Reproduce simulation, scenario 5 Table 7
#============================================================


# clear the environment and import necessary functions
rm(list = ls())

source('import functions for empirical application.R')



#load necessary packages
library(dplyr)
#library(RSpectra)
library(pROC)
library(openxlsx)
library(ggplot2)
library(xtable)
library(caret)
library(survival)
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
library(glmnet)

# data generating process, scenario 5
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
    ##########################################################
    
    # level of cross-sectional dependence across all variables
    phi = 0.9
    
    # degree of time series dependence among its lag
    rho = 0.9
    corr_matrix <- diag(1, p)
    
    for (v in 1:(p)) {
      for (u in 1:(p)) {
        corr_matrix[v, u] <- (phi)^abs(v - u)
      }
    }
    
    for (tt in 1:(jmax)){
      if (tt == 1){
        Xdd[,tt,(1:p)] <- MASS::mvrnorm(TN, mu = rep(0, p), corr_matrix) #scenario 2
      }
      else{
        Xdd[,tt,(1:p)] <- rho * Xdd[,tt-1,(1:p)] + MASS::mvrnorm(TN, mu = rep(0, p), corr_matrix*(1 - rho^2))
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
    
    for (z in seq(p)){
      z_idx <- idx + (z-1)*(degree+1)
      Xdw[,z_idx] <- Xdd[,,z]%*%w
      if (z == 1){
        Xdw_true[,z] <- Xdd[,,z]%*%w1
      }
      if (z == 2){
        Xdw_true[,z] <- Xdd[,,z]%*%w2
      }
    }
    
    #
    dataset = cbind(rep(1, n),  Xdw)
    
    # Inverse logit function
    logit_inv <- function(p) log( p / (1 - p))  # Inverse logit function
    uniform_rand <- runif(n)
    
    # X^T%*%l
    div = cbind(rep(1, n) ,  Xdw_true) %*% c(1, rep(1, length(parameters) - 1))
    
    # Calculate survival and censoring times
    Ts = s + exp( (logit_inv(uniform_rand) - (cbind(rep(1,n), Xdw_true) %*% parameters)) / (div))
    censoringtime = s + rexp(n, censor_strength)
    time_obs = pmin(Ts, censoringtime)
    status = as.numeric(Ts <= censoringtime)
    
    # Create result data frame
    result_df = as.data.frame(dataset)
    result_df$Ts = as.vector(Ts)
    result_df$censoringtime = censoringtime
    result_df$time = as.vector(time_obs)
    result_df$status = status
    result_df$X_orginal = X_orginal
    
    return(result_df)
  }

# grab the prediction metric for scenario 5
inf_result = function(result,it){
  se = NULL
  cover = NULL
  res_fit_MIDAS = NULL
  res_fit_MIDAS_LASSO = NULL
  res_fit_LASSO = NULL
  reject_MIDAS = NULL
  reject_LASSO = NULL
  reject_MIDAS_LASSO = NULL
  test_stat_MIDAS = NULL
  test_stat_MIDAS_LASSO = NULL
  test_stat_LASSO = NULL
  balance = NULL
  censoring_proportions = NULL
  it = length(result)
  for (i in seq(it)){
    res_fit_MIDAS= rbind(res_fit_MIDAS , result[[i]]$ref_MIDAS$est)
    res_fit_MIDAS_LASSO= rbind(res_fit_MIDAS_LASSO , result[[i]]$ref_MIDAS_LASSO$est)
    res_fit_LASSO= rbind(res_fit_LASSO , result[[i]]$ref_LASSO$est)
    reject_MIDAS = c(reject_MIDAS , as.numeric( result[[i]]$reject_MIDAS ))
    reject_MIDAS_LASSO = c(reject_MIDAS_LASSO , as.numeric( result[[i]]$reject_MIDAS_LASSO ))
    reject_LASSO = c(reject_LASSO , as.numeric( result[[i]]$reject_LASSO ))
    test_stat_MIDAS = c(test_stat_MIDAS, as.numeric( result[[i]]$test_stat_MIDAS ))
    test_stat_MIDAS_LASSO = c(test_stat_MIDAS_LASSO, as.numeric( result[[i]]$test_stat_MIDAS_LASSO ))
    test_stat_LASSO = c(test_stat_LASSO, as.numeric( result[[i]]$test_stat_LASSO ))
    se = rbind(se, result[[i]]$ref_MIDAS$se)
    cover = rbind(cover, result[[i]]$cover)
    censoring_proportions = rbind(censoring_proportions,  as.numeric( result[[i]]$censoring_proportions ))
    balance  = rbind(balance , as.numeric( result[[i]]$balance ))
  }
  se = colMeans(se)
  sd = apply(res_fit_MIDAS, 2, sd)
  cover = apply(cover, 2, mean)[(1+1):(1+1+degree)]
  ratio = (se/sd)[(1+1):(1+1+degree)]
  res_fit_av_MIDAS = colMeans(res_fit_MIDAS)
  res_fit_av_LASSO = 0#colMeans(res_fit_LASSO)
  res_fit_av_MIDAS_LASSO = colMeans(res_fit_MIDAS_LASSO)
  res_fit_var_LASSO = 0#mean(apply(res_fit_LASSO, 2, var))
  b_true = c(beta_true[1] + log((t - s)), beta_true[-1]  +  log((t - s)))
  
  # parameter estimation evaluation for w1
  w1_tr = w1 * b_true[1+1]
  w1_es = w%*%as.numeric(res_fit_av_MIDAS[(1+1):(1+1+degree)])
  w1_es_MIDAS_LASSO = w%*%as.numeric(res_fit_av_MIDAS_LASSO[(1+1):(1+1+degree)])
  w1_es_LASSO = 0 #as.numeric(res_fit_av_LASSO[(1+1):(1+1+jmax-1)])
  w1_es_sd = round(mean(apply(as.matrix(res_fit_MIDAS[, (1  + 1):(1  + 1 + degree)])%*%t(w), 2, sd)), 5)
  # here we use MSE = bias^2 + variance, which is close to direct calculation of MASE when the simulation repetation goes large enough
  # w1_es_MSE = mean(colMeans( (w1_tr - w %*% t(res_fit_MIDAS[ , (1+1):(1+1+degree)]) )^2 ))
  w1_es_MSE = round(mean(as.numeric((w1_tr - w1_es)^2) + apply(as.matrix(res_fit_MIDAS[, (1  + 1):(1  + 1 + degree)])%*%t(w), 2, var)), 5)
  w1_es_sd_LASSO = 0#round(mean(apply(as.matrix(res_fit_LASSO[, (1  + 1):(1  + 1 + jmax - 1)]), 2, sd)), 5)
  # w1_es_MSE_LASSO = mean(colMeans( (w1_tr -  t(res_fit_LASSO[, (1  + 1):(1  + 1 + jmax - 1)]) )^2 ))
  w1_es_MSE_LASSO = 0#round(mean(as.numeric((w1_tr - w1_es_LASSO)^2) + apply(as.matrix(res_fit_LASSO[, (1  + 1):(1  + 1 + jmax - 1)]), 2, var)), 5)
  w1_es_sd_MIDAS_LASSO = round(mean(apply(as.matrix(res_fit_MIDAS_LASSO[, (1  + 1):(1  + 1 + degree)])%*%t(w), 2, sd)), 5)
  # w1_es_MSE_MIDAS_LASSO = mean(colMeans( (w1_tr - w %*% t(res_fit_MIDAS_LASSO[ , (1+1):(1+1+degree)]) )^2 ))
  w1_es_MSE_MIDAS_LASSO = round(mean(as.numeric((w1_tr - w1_es_MIDAS_LASSO)^2) + apply(as.matrix(res_fit_MIDAS_LASSO[, (1  + 1):(1  + 1 + degree)])%*%t(w), 2, var)), 5)
  
  # parameter estimation evaluation for w1
  w2_tr = w2 * b_true[1+1+1]
  w2_es = w%*%as.numeric(res_fit_av_MIDAS[(1+1+degree+1):(1+1+degree+1+degree)])
  w2_es_LASSO = 0#as.numeric(res_fit_av_LASSO[(1+1+jmax):(1+1+jmax*2-1)])
  w2_es_MIDAS_LASSO = w%*%as.numeric(res_fit_av_MIDAS_LASSO[(1+1+degree+1):(1+1+degree+1+degree)])
  w2_es_sd = round(mean(apply(as.matrix(res_fit_MIDAS[, (1+1+degree+1):(1+1+degree+1+degree)])%*%t(w), 2, sd)), 5)
  w2_es_sd_LASSO = 0#round(mean(apply(as.matrix(res_fit_LASSO[, (1+1+jmax):(1+1+jmax*2-1)]), 2, sd)), 5)
  w2_es_sd_MIDAS_LASSO = round(mean(apply(as.matrix(res_fit_MIDAS_LASSO[, (1+1+degree+1):(1+1+degree+1+degree)])%*%t(w), 2, sd)), 5)
  # w2_es_MSE = mean(colMeans( (w2_tr - w %*% t(res_fit_MIDAS[, (1+1+degree+1):(1+1+degree+1+degree)]) )^2 ))
  # w2_es_MSE_LASSO = mean(colMeans( (w2_tr -  t(res_fit_LASSO[, (1+1+jmax):(1+1+jmax*2-1)]) )^2 ))
  # w2_es_MSE_MIDAS_LASSO = mean(colMeans( (w2_tr - w %*% t(res_fit_MIDAS_LASSO[, (1+1+degree+1):(1+1+degree+1+degree)] ))^2 ))
  w2_es_MSE = round(mean( as.numeric((w2_tr - w2_es)^2) + apply(as.matrix(res_fit_MIDAS[, (1+1+degree+1):(1+1+degree+1+degree)])%*%t(w), 2, var)), 5)
  w2_es_MSE_LASSO = 0#round(mean(as.numeric((w2_tr - w2_es_LASSO)^2) + apply(as.matrix(res_fit_LASSO[, (1+1+jmax):(1+1+jmax*2-1)]), 2, var)), 5)
  w2_es_MSE_MIDAS_LASSO = round(mean(as.numeric((w2_tr - w2_es_MIDAS_LASSO)^2) + apply(as.matrix(res_fit_MIDAS_LASSO[, (1+1+degree+1):(1+1+degree+1+degree)])%*%t(w), 2, var)), 5)
  
  # return the result
  res = list(reject_MIDAS = mean(reject_MIDAS), reject_MIDAS_LASSO = mean(reject_MIDAS_LASSO), reject_LASSO = mean(reject_LASSO),
             test_stat_MIDAS = test_stat_MIDAS, test_stat_MIDAS_LASSO = test_stat_MIDAS_LASSO, test_stat_LASSO = test_stat_LASSO,
             w1_es_MSE = w1_es_MSE, w1_es_sd = w1_es_sd,
             w2_es_MSE = w2_es_MSE, w2_es_sd = w2_es_sd,
             w1_es_MSE_LASSO = w1_es_MSE_LASSO, w1_es_sd_LASSO = w1_es_sd_LASSO,
             w2_es_MSE_LASSO = w2_es_MSE_LASSO, w2_es_sd_LASSO = w2_es_sd_LASSO,
             w1_es_MSE_MIDAS_LASSO = w1_es_MSE_MIDAS_LASSO, w1_es_sd_MIDAS_LASSO = w1_es_sd_MIDAS_LASSO,
             w2_es_MSE_MIDAS_LASSO = w2_es_MSE_MIDAS_LASSO, w2_es_sd_MIDAS_LASSO = w2_es_sd_MIDAS_LASSO,
             se = se, sd =sd, ratio = ratio, cover = cover,
             censoring_proportions = colMeans(censoring_proportions), balance = colMeans(balance))
  return(res)
}

# calculate the censoring rate of a simulation dataset (scenario 5)
test_censoring_AR = function(s, n, censor_strength){
  res = 0
  for (i in seq(100)){
    data =generateData_AR( s = s, n = n, numhv= numtrue, numtrue=numtrue, degree=degree, jmax=jmax, parameters = beta_true, censor_strength = censor_strength) # 3, 500, 0.7
    res = res + length(which(data$status==0))/n
  }
  return(res/100)
}


#============================================================
# main part
#============================================================

# initial setting
s = 6

#quarter/year frequency
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
w1 <- dbeta((1:jmax)/(jmax), shape1 = 1, shape2 = 3)*0.1
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
alpha = c(0, 0.5, 1)

#censoring strength gamma, which can generate approximately 81% censoring
censor = 0.81 
test_censoring_AR( s = s, n = 1000, censor_strength = censor)

# Generate the prediction horizon t
t = 0
for (k in 1:100) {
  set.seed(k)
  print(k)
  data = generateData_AR( s = s, n = 1200, numhv = numhv, numtrue = numtrue, degree = degree, jmax = jmax, parameters = beta_true, censor_strength = censor) # 3, 500, 0.7
  ind = which(data$Ts <=  data$censoringtime)
  t_quantile = c(quantile(data$Ts[ind], probs = 0.1), quantile(data$Ts[ind], probs = 0.3), quantile(data$Ts[ind], probs = 0.5))
  censor_strength = c(censor)
  t = t + quantile(data$Ts[ind], probs = 0.3)
}

t = t/100
t
b_truec = c(beta_true[1] + log((t - s)), beta_true[-1]  +  log((t - s)))
b_truec = c(b_truec[1:(1)], w1*b_truec[1+1], w2*b_truec[1+1+1])

#call for parallel, notice that the number of cores is dependent on your own laptop
if (exists("cl")) { try(stopCluster(cl), silent = TRUE) }
num_cores <- detectCores()
cl <- makeCluster(max(1, num_cores - 22))
registerDoParallel(cl)
packages_to_export <- c('MASS', 'timeROC', 'mvtnorm', 'Survivalml', "dplyr",  "survival", "MLmetrics", "pROC", "PRROC", 'survivalROC', 'dotCall64', 'glmnet', 'rlang', 'pracma')
vars_to_export <- ls(globalenv())

# repetitions of the simulation
it = 500
sig_level = 0.05
v <- qnorm(sig_level/2, lower.tail=F)

#scenario 5 N = 800
i=0
pa = matrix(0, nrow = length(t_quantile), ncol = 20)
i = i + 1
it = it
result = foreach(k = 1:it, .errorhandling= 'remove', .packages = packages_to_export, .export = vars_to_export) %dopar% {
  n = 1200
  set.seed(k)
  
  # data prepare
  data = generateData_AR( s = s, n = n, numhv = numhv, numtrue = numtrue, degree = degree, jmax = jmax, parameters = beta_true, censor_strength = censor) # 3, 500, 0.7
  dataset_p = data[ which( (1*( data$time <= t ))*(1*( data$status == 1 )) == 1), ]
  dataset_n = data[ which( (1*( data$time <= t ))*(1*( data$status == 1 )) == 0), ]
  dataset = rbind(dataset_p, dataset_n )
  dataset = testRandomLogitDataset(dataset, t = t )
  dataset = KM_estimate(dataset)
  
  #5-fold cross validation
  nfold = 5
  index_fold = which(1*(dataset$time <= t) * 1*(dataset$status == 1) == 1)
  length(index_fold)
  if (length(index_fold) < 5){
    foldid = NULL
  }else{
    foldid = c( form_folds(nrow(dataset[index_fold,]), nfold), form_folds(nrow(dataset[-index_fold,]), nfold) )
  }
  
  #dataset for LASSO-UMIDAS
  X_orginal = dataset$X_orginal[ , 1:(1+jmax*numhv) ]
  
  #dataset for sg-LASSO-MIDAS
  X = dataset[ , 1:(numhv*(degree+1)+1) ]
  y = 1*(dataset$time <= t)
  
  # IPW weights
  w_train = 1*( pmin(dataset$Ts, t) <= dataset$censoringtime)/dataset$Gp
  IF = IF_estimate(dataset)
  w_km = 1*( pmin(dataset$Ts, t) <= dataset$censoringtime)/dataset$Gp^2 * y
  V_weighted <- IF * matrix(w_km, n, n)
  var_km <- (V_weighted %*% as.matrix(X))/(n)
  
  # sg-LASSO-MIDAS, cross-validation for AUC
  fit_cv_MIDAS = alpha_cv_sparsegl(as.matrix(X[,-1]), y, group = gindex, nlambda = 50, weight = w_train, alpha = alpha, foldid = foldid, nfolds = nfold,  pred.loss = 'censor', intercept_zero = intercept_zero, standardize = TRUE, AUC = FALSE, data = dataset, t = t, maxit = 30000)
  est_MIDAS = as.vector( fit_cv_MIDAS$coff)
  ref_MIDAS = ORIG_DS_inf(x=as.matrix(X[,-1]), y = y, family="binomial", lasso_est = est_MIDAS, weight = w_train, interest_number = degree + 1, foldid = foldid, variance_km = var_km)
  ci_upper_orig_ds <- ref_MIDAS$est + v*ref_MIDAS$se
  ci_lower_orig_ds <- ref_MIDAS$est - v*ref_MIDAS$se
  cover = 1*(0 <= ci_upper_orig_ds) * 1*(0>=ci_lower_orig_ds)
  test_stat_MIDAS = n * t(ref_MIDAS$est[(1+1):(1+1+degree)])%*% solve(ref_MIDAS$sigma[(1+1):(1+1+degree), (1+1):(1+1+degree)])%*% ref_MIDAS$est[(1+1):(1+1+degree)]
  reject_MIDAS = 1*(test_stat_MIDAS > qchisq(p = 1 - sig_level, df = degree+1))
  
  #LASSO-MIDAS
  if (fit_cv_MIDAS$alpha_out == 1){
    fit_cv_MIDAS_LASSO = fit_cv_MIDAS
    est_MIDAS_LASSO = est_MIDAS
  } else {
    fit_cv_MIDAS_LASSO = fit_cv_MIDAS$fit_MIDAS_LASSO
    est_MIDAS_LASSO = as.vector( coef(fit_cv_MIDAS_LASSO, s=c('lambda.min') ) )
  }
  
  balance = length(data$Ts[data$Ts <= t])/n
  censoring_proportions = length(which(data$status==0))/n
  
  results = list( censoring_proportions = censoring_proportions, t = t, balance = round( mean(balance), 3 ),
                  fit_LASSO = est_MIDAS, fit_MIDAS = est_MIDAS, fit_MIDAS_LASSO = est_MIDAS_LASSO,
                  ref_MIDAS = ref_MIDAS, ref_MIDAS_LASSO = ref_MIDAS, ref_LASSO = ref_MIDAS,
                  reject_MIDAS =  reject_MIDAS,  reject_MIDAS_LASSO =  reject_MIDAS,  reject_LASSO =  reject_MIDAS, cover = cover,
                  test_stat_MIDAS = test_stat_MIDAS, test_stat_MIDAS_LASSO = test_stat_MIDAS, test_stat_LASSO = test_stat_MIDAS)
  return(results)
}

res = inf_result(result, length(result))
pa[i, ] = c(res$censoring_proportions, res$balance,
            'p value',
            round(res$reject_LASSO,3), round(res$reject_MIDAS_LASSO,3), round(res$reject_MIDAS,3),
            'w1',
            round(res$w1_es_MSE_LASSO,3), round(res$w1_es_MSE_MIDAS_LASSO,3), round(res$w1_es_MSE,3),
            round(res$w1_es_sd_LASSO,3), round(res$w1_es_sd_MIDAS_LASSO,3), round(res$w1_es_sd,3),
            'w2',
            round(res$ratio,3), round(res$cover,3))
print(pa[i, ])
table_name <- paste0("MIDASAUC_scenario5_inference_power", censor, 1200, s,  ".csv")
write.table(pa, table_name, sep = ";", row.names = FALSE, quote = FALSE)


