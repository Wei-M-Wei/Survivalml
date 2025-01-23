#============================================================
# Corporate Survival Analysis with Machine Learning Methods
# Import functions for simulation and empirical application
#============================================================


#============================================================
# Reproduce empirical application when s = 6 years and t = 8,8.5,9 years
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


#empirical settings
s = 6
lag_use_year = 6
quarter = 4
lags = lag_use_year * quarter - 1 ############we don't have the information of the last lag.
jmax <- lags



# read the data and check basic information
data_financial = read.csv2("period_final_new6.csv")  # read function from package readxl
n = dim(data_financial)[1]
p_fin = dim(data_financial)[2] - 5
for (col in names(data_financial)[1:p_fin]) {
  if (grepl("%", col)) {
    data_financial[[col]] <- data_financial[[col]] * 0.01
  }
}
data = data_financial
which('true' == is.na(data))# if null values exist
p = dim(data)[2] - 5 #'start_day, banktuptcy_day, T, censoring_status, C' variables in the last 5 columns
n = dim(data)[1]



# Define survival time, censoring time to calculate KM weights
end_observation = '2020-12-31'
data$start_day = sapply(data$start_day, convert_to_ymd)
data$censoringtime = as.numeric( as.Date(end_observation) - as.Date(data$start_day) )/365
data$Ts = data$survival_time
data$time = pmin(data$Ts, data$censoringtime)
data$status = data$status
censoring_rate = sum(data$status)/dim(data)[1]


# MIDAS weights
intercept_zero = 0
degree = 2
w_fin <- gb(degree = degree, alpha = -1/2, a = 0, b = 1, jmax = jmax)/jmax


# tuning parameter grid search
alpha = c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1)

#bootstrap size
bootstrap_number = 1000

# parallel settings: we need 10 cores
num_cores <- detectCores()
cl <- makeCluster(num_cores-2)
registerDoParallel(cl)
packages_to_export <- c('timeROC', 'Survivalml', "dplyr", "midasml", "survival", "MLmetrics", "pROC", "PRROC", 'survivalROC', 'dotCall64', 'glmnet', 'rlang', 'pracma')
##############################################################


# prediction horizons t = 8, 8.5, 9 years, Estimator(5) in this paper

for (t in c( 9)) {
  # iteration number
  it = 10
  result_para = foreach(k = 1:it, .errorhandling= 'remove', .packages = packages_to_export, .export = export_env) %dopar% {
    set.seed( s * k)
    dataset_p = data[ which( ( 1*(data$time <= t) * 1*(data$status==1)) == 1), ]
    dataset_n = data[ which( (  1*(data$time <= t) * 1*(data$status==1)) == 0), ]
    indices_p <- sample( nrow( dataset_p ), nrow( dataset_p ) * 0.8 )
    indices_n = sample( nrow( dataset_n ), nrow( dataset_n ) * 0.8 )
    train_dataset_p <- dataset_p[indices_p,]
    train_dataset_n <- dataset_n[ indices_n, ]
    train_dataset_in = rbind( train_dataset_p, train_dataset_n )
    test_dataset_p <- dataset_p[ -indices_p, ]
    test_dataset_n <- dataset_n[ -indices_n, ]
    test_dataset_in = rbind(test_dataset_p, test_dataset_n)

    ##########################
    train_dataset = testRandomLogitDataset( train_dataset_in, t = t )
    train_dataset = KM_estimate(train_dataset)
    test_dataset = testRandomLogitDataset( test_dataset_in, t = t )
    test_dataset = KM_estimate(test_dataset)

    #################################Apply MIDAS

    # extract the latest lag of each covariate
    X_train = train_dataset[,seq(lags, dim(train_dataset)[2],lags)]
    y_train = 1*(train_dataset$time <= t)
    X_test = test_dataset[,seq(lags, dim(test_dataset)[2],lags)]
    X_train = apply(X_train, 2, as.numeric)
    X_test = apply(X_test, 2, as.numeric)
    w_train = (1 - 1*( train_dataset$time <= t )*(1 - train_dataset$status))/train_dataset$Gp
    w_test = (1 - 1*( test_dataset$time <= t )*(1 - test_dataset$status))/test_dataset$Gp


    fit_bench = survival_sparsegl(X_train, y_train, group = seq(dim(X_train)[2]), lambda = c(0), weight = w_train, asparse = 1, intercept_zero = intercept_zero, standardize = TRUE, eps = 1e-8, maxit = 30000)

    est_LASSO_fin = unlist(c(as.numeric(fit_bench$b0), as.numeric(fit_bench$beta)))
    test_predictions_MIDAS =unlist( as.list(plogis( cbind(rep(1, nrow(X_test)), X_test) %*% est_LASSO_fin)) )


    AUC_MIDAS = ROC_censor_N( data = test_dataset, prediction = test_predictions_MIDAS, t = t )$AUC
    AUC_MIDAS_bootstrap =  ROC_N_bootstrap(data = test_dataset, prediction = test_predictions_MIDAS, t = t, sim_number = bootstrap_number )
    print(AUC_MIDAS)
    results = list(  AUC_MIDAS = round(mean(AUC_MIDAS), 3 ),  AUC_MIDAS_bootstrap = AUC_MIDAS_bootstrap)
    return(results)
  }
  result = result_para
  AUC_MIDAS = NULL
  AUC_MIDAS_bootstrap = NULL
  result = result[!sapply(result, is.null)]
  it = length(result)
  print(it)
  for (i in seq(it)){
    AUC_MIDAS = rbind(AUC_MIDAS , as.numeric( result[[i]]$AUC_MIDAS ))
    AUC_MIDAS_bootstrap = rbind(AUC_MIDAS_bootstrap , as.numeric( result[[i]]$AUC_MIDAS_bootstrap ))
  }
  AUC_MIDAS_bootstrap_av = colMeans(AUC_MIDAS_bootstrap)
  res = list(AUC_MIDAS = round(colMeans(AUC_MIDAS),3), AUC_MIDAS_VAR = round(apply(AUC_MIDAS, 2, var),3), AUC_MIDAS_bootstrap_av = round(AUC_MIDAS_bootstrap_av, 3 ))


  # save the table
  table_name1 <- paste0("bench_MIDAS_AUC", s, "_AUC_cvm", t, ".csv")
  write.csv(rbind( c(res$AUC_MIDAS, res$AUC_MIDAS_VAR) ,
                   c( quantile(res$AUC_MIDAS_bootstrap_av, 0.025), quantile(res$AUC_MIDAS_bootstrap_av, 0.975)
                   )) , table_name1)

}

########################################
stopCluster()

