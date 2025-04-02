rm(list = ls())
###############################################


KM_estimate <-
  function(data){
    data = data[order(data$m_time_t),]
    survival_data <- Surv(data$T_hat, rep(1,nrow(data)) - data$status)
    km_fit <- survfit(survival_data ~ 1)
    data$Gp = summary(km_fit, times = data$m_time_t)$surv
    set.seed(123) ## s = 6 and seed 1
    data = data[sample(nrow(data)),]
    return(data)
  }


######################################
s = 6
lag_use_year = 6
quarter = 4
lags = lag_use_year * quarter - 1
jmax <- lags
##########################read the data and check basic information
data_financial = read.csv2("period_final_new6.csv")  # read function from package readxl
#data_financial = X_final
dim(data_financial)
sum(data_financial$status)
n = dim(data_financial)[1]
p_fin = dim(data_financial)[2] - 5
p_fin/lags
min(data_financial$censoringtime)
for (col in names(data_financial)[1:p_fin]) {
  if (grepl("%", col)) {
    data_financial[[col]] <- data_financial[[col]] * 0.01
  }
}

###########check the median of survival time and censoring time
data = data_financial

end_observation = '2020-12-31'
data$start_day = sapply(data$start_day, convert_to_ymd)
data$censoringtime = as.numeric( as.Date(end_observation) - as.Date(data$start_day) )/365
data$Ts = data$survival_time
data$time = pmin(data$Ts, data$censoringtime)
data$status = data$status
censoring_rate = sum(data$status)/dim(data)[1]


cut_point = as.Date('2016-12-31')

test = data[which( as.numeric( cut_point - as.Date(data$start_day) )/365 < s),]
train = data[-which( as.numeric( cut_point - as.Date(data$start_day) )/365 < s),]

train$censoringtime = as.numeric( as.Date(cut_point) - as.Date(train$start_day) )/365
train$time = pmin(train$Ts, train$censoringtime)
ind_status = which(train$Ts <= train$censoringtime)
train$status[ind_status] = 1
train$status[-ind_status] = 0

censor_firm = which(train$status == 0)
min(train$censoringtime)
sort(train[censor_firm,]$survival_time)

dim(train)
dim(test)

#bootstrap size
bootstrap_number = 1000

###################

for (t in c(8, 8.5, 9)) {
  set.seed(s)
  ############################################allocate training dataset and test dataset
  train_dataset_in = train
  test_dataset_in = test

  ###########################################calculate the KM weights for each observation
  train_dataset = testRandomLogitDataset( train_dataset_in, t = t )
  train_dataset = KM_estimate(train_dataset)
  test_dataset = testRandomLogitDataset( test_dataset_in, t = t )
  test_dataset = KM_estimate(test_dataset)


  #we have 14 variables which are related to T,C and etc. Then delete them
  X_train = train_dataset[,seq(lags, dim(train_dataset)[2],lags)]
  y_train = 1*(train_dataset$time <= t)
  X_test = test_dataset[,seq(lags, dim(test_dataset)[2],lags)]
  X_train = apply(X_train, 2, as.numeric)
  X_test = apply(X_test, 2, as.numeric)
  w_train = (1 - 1*( train_dataset$time <= t )*(1 - train_dataset$status))/train_dataset$Gp
  w_test = (1 - 1*( test_dataset$time <= t )*(1 - test_dataset$status))/test_dataset$Gp


  fit_bench = survival_sparsegl(X_train, y_train, group = seq(dim(X_train)[2]), lambda = c(0), weight = w_train, asparse = 1, intercept_zero = intercept_zero, standardize = TRUE, eps = 20e-7, maxit = 100000)

  est_LASSO_fin = unlist(c(as.numeric(fit_bench$b0), as.numeric(fit_bench$beta)))
  test_predictions_MIDAS =unlist( as.list(plogis( cbind(rep(1, nrow(X_test)), X_test) %*% est_LASSO_fin)) )


  AUC_MIDAS = ROC_censor_N( data = test_dataset, prediction = test_predictions_MIDAS, t = t )$AUC
  AUC_MIDAS_bootstrap =  ROC_N_bootstrap(data = test_dataset, prediction = test_predictions_MIDAS, t = t, sim_number = bootstrap_number )
  print(AUC_MIDAS)
  results = list(  AUC_MIDAS = round(mean(AUC_MIDAS), 3 ),  AUC_MIDAS_bootstrap = AUC_MIDAS_bootstrap)


  result = results

  AUC_MIDAS = NULL
  AUC_MIDAS_bootstrap = NULL

  result = result[!sapply(result, is.null)]

  for (i in seq(1)){
    AUC_MIDAS = rbind(AUC_MIDAS , as.numeric( result$AUC_MIDAS ))
    AUC_MIDAS_bootstrap = rbind(AUC_MIDAS_bootstrap , as.numeric( result$AUC_MIDAS_bootstrap ))
  }

  AUC_MIDAS_bootstrap_av = colMeans(AUC_MIDAS_bootstrap)
  res = list(AUC_MIDAS = round(colMeans(AUC_MIDAS),3), AUC_MIDAS_VAR = round(apply(AUC_MIDAS, 2, var),3), AUC_MIDAS_bootstrap_av = round(AUC_MIDAS_bootstrap_av, 3 ))


  # save the table
  table_name1 <- paste0("application2_bench_MIDAS_AUC", s, "_AUC_cvm", t, ".csv")
  write.csv(rbind( c(res$AUC_MIDAS, res$AUC_MIDAS_VAR) ,
                   c( quantile(res$AUC_MIDAS_bootstrap_av, 0.025), quantile(res$AUC_MIDAS_bootstrap_av, 0.975)
                   )) , table_name1)

}
