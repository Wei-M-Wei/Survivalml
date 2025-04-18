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
dim(data_financial)
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
# Legendre weights
intercept_zero = 0
degree = 2
degree_macro = 2
w_fin <- gb(degree = degree, alpha = 1/2, a = 0, b = 1, jmax = jmax)/jmax #lb(degree = degree, a = 0, b = 1, jmax = jmax)/jmax # gb(degree = degree, alpha = alpha, a = 0, b = 1, jmax = jmax)
w_macro = gb(degree = degree_macro, alpha = 1/2, a = 0, b = 1, jmax = jmax)/jmax #lb(degree = degree, a = 0, b = 1, jmax = jmax)/jmax # gb(degree = degree, alpha = alpha, a = 0, b = 1, jmax = jmax)
# tuning parameter grid search
alpha = c(0.8, 0.85, 0.9, 0.95, 1)

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

  # stratified cross validation, see Section 6.2
  nfold = 5
  index_fold = which(1*(train_dataset$time <= t) * 1*(train_dataset$status == 1) == 1)
  if (length(index_fold) <5){
    foldid = NULL
  }else{
    foldid = c( form_folds(nrow(train_dataset[index_fold,]), nfold), form_folds(nrow(train_dataset[-index_fold,]), nfold) )
  }

  #we have 14 variables which are related to T,C and etc. Then delete them
  X_train = train_dataset[ , 1: (dim(train_dataset)[2]-14) ]
  y_train = 1*(train_dataset$time <= t)
  X_test = as.matrix( test_dataset[ , 1: (dim(train_dataset)[2]-14) ] )
  y_test = 1*( test_dataset$time <= t )
  X_train = apply(X_train, 2, as.numeric)
  X_test = apply(X_test, 2, as.numeric)
  w_train = (1 - 1*( train_dataset$time <= t )*(1 - train_dataset$status))/train_dataset$Gp
  w_test = (1 - 1*( test_dataset$time <= t )*(1 - test_dataset$status))/test_dataset$Gp

  # MIDAS setting
  idx <- 1:jmax
  gindex <- NULL
  Xdw = NULL
  for (z in seq(dim(X_train)[2]/jmax)){
    z_idx <- (1 + (z - 1) * jmax) : (z * jmax)
    Xdw <- cbind(Xdw, X_train[,z_idx] %*% w_fin)
    gindex <- c(gindex, rep(z, times = degree + 1))
  }

  #training dataset which are used to in the SGL-MIDAS
  X_in = as.matrix(cbind(rep(1, nrow(X_train)), Xdw))
  fit_cv_MIDAS = alpha_cv_sparsegl( X_in[,-1], y_train, group = gindex, nlambda = 200, weight = w_train, alpha = alpha, nfolds = 5, foldid = foldid, pred.loss = 'censor', intercept_zero = intercept_zero, standardize = TRUE, AUC = TRUE, data = train_dataset, t = t)
  est_MIDAS = fit_cv_MIDAS$coff
  est_MIDAS_AUC = fit_cv_MIDAS$coff_AUC

  #test dataset
  X_t = NULL
  for (z in seq(dim(X_test)[2]/jmax)){
    z_idx <- (1 + (z - 1)*jmax) : (z * jmax)
    X_t <- cbind(X_t, X_test[,z_idx]%*%w_fin)
  }
  X_te = as.matrix(cbind(rep(1, nrow(X_test)),  X_t))

  #time dependent AUC and bootstrap
  test_predictions_MIDAS = unlist( as.list(plogis(X_te %*% est_MIDAS)) )
  test_predictions_MIDAS_AUC = unlist( as.list(plogis(X_te %*% est_MIDAS_AUC)) )

  AUC_MIDAS = ROC_censor_N( data = test_dataset, prediction = test_predictions_MIDAS, t = t )$AUC
  AUC_MIDAS_AUC = ROC_censor_N( data = test_dataset, prediction = test_predictions_MIDAS_AUC, t = t )$AUC

  AUC_MIDAS_bootstrap =  ROC_N_bootstrap(data = test_dataset, prediction = test_predictions_MIDAS, t = t, sim_number = bootstrap_number )
  AUC_MIDAS_AUC_bootstrap =  ROC_N_bootstrap(data = test_dataset, prediction = test_predictions_MIDAS_AUC, t = t, sim_number = bootstrap_number )


  #LASSO-MIDAS, cross-validation for log-likelihood score
  if (fit_cv_MIDAS$alpha_out == 1){
    fit_cv_MIDAS_LASSO = fit_cv_MIDAS
    test_predictions_MIDAS_LASSO = test_predictions_MIDAS
    est_MIDAS_LASSO = est_MIDAS
  } else {
    fit_cv_MIDAS_LASSO = fit_cv_MIDAS$fit_MIDAS_LASSO
    cesnor_min_index = which(fit_cv_MIDAS_LASSO$cvm == min(fit_cv_MIDAS_LASSO$cvm))[1]
    est_MIDAS_LASSO = unlist(c(fit_cv_MIDAS_LASSO$survival_sparsegl.fit$b0[,cesnor_min_index], fit_cv_MIDAS_LASSO$survival_sparsegl.fit$beta[,cesnor_min_index]))
    test_predictions_MIDAS_LASSO = unlist( as.list(plogis(X_te %*% est_MIDAS_LASSO )) ) #predict(fit_cv_MIDAS_LASSO, X_te[,-1], s = 'lambda.min', type = 'response')
  }
  AUC_MIDAS_LASSO = ROC_censor_N(data = test_dataset, prediction = test_predictions_MIDAS_LASSO, t = t )$AUC

  #LASSO-MIDAS, cross-validation for AUC
  if (fit_cv_MIDAS$alpha_AUC == 1){
    fit_cv_MIDAS_LASSO = fit_cv_MIDAS
    test_predictions_MIDAS_LASSO_AUC = test_predictions_MIDAS_AUC
    est_MIDAS_LASSO_AUC = est_MIDAS_AUC
  } else {
    fit_cv_MIDAS_LASSO = fit_cv_MIDAS$fit_MIDAS_LASSO
    cesnor_min_index = which(fit_cv_MIDAS_LASSO$AUC_censor == max(fit_cv_MIDAS_LASSO$AUC_censor))[1]
    est_MIDAS_LASSO_AUC = unlist(c(fit_cv_MIDAS_LASSO$survival_sparsegl.fit$b0[,cesnor_min_index], fit_cv_MIDAS_LASSO$survival_sparsegl.fit$beta[,cesnor_min_index]))
    test_predictions_MIDAS_LASSO_AUC = unlist( as.list(plogis(X_te %*% est_MIDAS_LASSO_AUC )) ) #predict(fit_cv_MIDAS_LASSO, X_te[,-1], s = 'lambda.min', type = 'response')
  }
  AUC_MIDAS_LASSO_AUC = ROC_censor_N(data = test_dataset, prediction = test_predictions_MIDAS_LASSO_AUC, t = t )$AUC

  # bootstrap for LASSO-MIDAS
  AUC_MIDAS_LASSO_bootstrap =  ROC_N_bootstrap(data = test_dataset, prediction = test_predictions_MIDAS_LASSO, t = t, sim_number = bootstrap_number )
  AUC_MIDAS_LASSO_AUC_bootstrap =  ROC_N_bootstrap(data = test_dataset, prediction = test_predictions_MIDAS_LASSO_AUC, t = t, sim_number = bootstrap_number )

  #LASSO-UMIDAS
  fit_cv_LASSO = cv.survival_sparsegl(X_train, y_train, group = seq(dim(X_train)[2]), nlambda = 200, weight = w_train, asparse = 1, foldid = foldid, nfolds = 5, pred.loss = 'censor', intercept_zero = intercept_zero, standardize = TRUE, maxit = 30000, AUC = TRUE, data = train_dataset, t = t)
  nan_indices <- sapply(fit_cv_LASSO$cvm, is.nan)
  fit_cv_LASSO$cvm[nan_indices] <- 1000
  cesnor_min_index_LASSO = which(fit_cv_LASSO$cvm== min(fit_cv_LASSO$cvm))[1]
  est_LASSO = unlist(c(fit_cv_LASSO$survival_sparsegl.fit$b0[,cesnor_min_index_LASSO], fit_cv_LASSO$survival_sparsegl.fit$beta[,cesnor_min_index_LASSO]))

  nan_indices <- sapply(fit_cv_LASSO$AUC_censor, is.nan)
  fit_cv_LASSO$AUC_censor[nan_indices] <- 0
  cesnor_min_index_LASSO = which(fit_cv_LASSO$AUC_censor== max(fit_cv_LASSO$AUC_censor))[1]
  est_LASSO_AUC = unlist(c(fit_cv_LASSO$survival_sparsegl.fit$b0[,cesnor_min_index_LASSO], fit_cv_LASSO$survival_sparsegl.fit$beta[,cesnor_min_index_LASSO]))

  test_predictions_LASSO = unlist( as.list(plogis( cbind(rep(1, nrow(X_test)), X_test) %*% est_LASSO )) ) # predict(fit_cv_LASSO, X_test, s = 'lambda.min', type = 'response') # unlist( as.list(plogis(X_test %*% est_LASSO )) )
  test_predictions_LASSO_AUC = unlist( as.list(plogis( cbind(rep(1, nrow(X_test)), X_test) %*% est_LASSO_AUC )) ) # predict(fit_cv_LASSO, X_test, s = 'lambda.min', type = 'response') # unlist( as.list(plogis(X_test %*% est_LASSO )) )

  AUC_LASSO = ROC_censor_N(data = test_dataset, prediction = test_predictions_LASSO, t = t )$AUC
  AUC_LASSO_AUC = ROC_censor_N(data = test_dataset, prediction = test_predictions_LASSO_AUC, t = t )$AUC

  AUC_LASSO_bootstrap =  ROC_N_bootstrap(data = test_dataset, prediction = test_predictions_LASSO, t = t, sim_number = bootstrap_number )
  AUC_LASSO_AUC_bootstrap =  ROC_N_bootstrap(data = test_dataset, prediction = test_predictions_LASSO_AUC, t = t, sim_number = bootstrap_number )

  #return results
  results = list(
    AUC_LASSO = round( AUC_LASSO, 3 ), AUC_MIDAS = round(AUC_MIDAS, 3 ), AUC_MIDAS_LASSO = round(AUC_MIDAS_LASSO, 3 ),
    AUC_LASSO_bootstrap = AUC_LASSO_bootstrap, AUC_MIDAS_bootstrap = AUC_MIDAS_bootstrap, AUC_MIDAS_LASSO_bootstrap = AUC_MIDAS_LASSO_bootstrap,
    fit_LASSO = est_LASSO, fit_MIDAS = est_MIDAS, fit_MIDAS_LASSO = est_MIDAS_LASSO,

    AUC_LASSO_AUC = round( AUC_LASSO_AUC, 3 ), AUC_MIDAS_AUC = round(AUC_MIDAS_AUC, 3 ), AUC_MIDAS_LASSO_AUC = round(AUC_MIDAS_LASSO_AUC, 3 ),
    AUC_LASSO_AUC_bootstrap = AUC_LASSO_AUC_bootstrap, AUC_MIDAS_AUC_bootstrap = AUC_MIDAS_AUC_bootstrap, AUC_MIDAS_LASSO_AUC_bootstrap = AUC_MIDAS_LASSO_AUC_bootstrap,
    fit_LASSO_AUC = est_LASSO_AUC, fit_MIDAS_AUC = est_MIDAS_AUC, fit_MIDAS_LASSO_AUC = est_MIDAS_LASSO_AUC, k = k
  )

  ###################################################################################
  result = results

  res_fit_MIDAS = NULL
  res_fit_MIDAS_LASSO = NULL
  res_fit_LASSO = NULL
  res_fit_MIDAS_AUC = NULL
  res_fit_MIDAS_LASSO_AUC = NULL
  res_fit_LASSO_AUC = NULL

  AUC_LASSO = NULL
  AUC_MIDAS = NULL
  AUC_MIDAS_LASSO = NULL
  AUC_LASSO_AUC = NULL
  AUC_MIDAS_AUC = NULL
  AUC_MIDAS_LASSO_AUC = NULL

  AUC_LASSO_bootstrap = NULL
  AUC_MIDAS_bootstrap = NULL
  AUC_MIDAS_LASSO_bootstrap = NULL
  AUC_LASSO_AUC_bootstrap = NULL
  AUC_MIDAS_AUC_bootstrap = NULL
  AUC_MIDAS_LASSO_AUC_bootstrap = NULL


  for (i in seq(1)){

    ######################################################
    res_fit_MIDAS= rbind(res_fit_MIDAS , result$fit_MIDAS)
    res_fit_MIDAS_LASSO= rbind(res_fit_MIDAS_LASSO , result$fit_MIDAS_LASSO)
    res_fit_LASSO= rbind(res_fit_LASSO , result$fit_LASSO)
    res_fit_MIDAS_AUC= rbind(res_fit_MIDAS_AUC , result$fit_MIDAS_AUC)
    res_fit_MIDAS_LASSO_AUC= rbind(res_fit_MIDAS_LASSO_AUC , result$fit_MIDAS_LASSO_AUC)
    res_fit_LASSO_AUC= rbind(res_fit_LASSO_AUC , result$fit_LASSO_AUC)


    AUC_MIDAS = rbind(AUC_MIDAS , as.numeric( result$AUC_MIDAS ))
    AUC_MIDAS_LASSO = rbind(AUC_MIDAS_LASSO , as.numeric( result$AUC_MIDAS_LASSO ))
    AUC_LASSO = rbind(AUC_LASSO , as.numeric( result$AUC_LASSO ))
    AUC_MIDAS_AUC = rbind(AUC_MIDAS_AUC , as.numeric( result$AUC_MIDAS_AUC ))
    AUC_MIDAS_LASSO_AUC = rbind(AUC_MIDAS_LASSO_AUC , as.numeric( result$AUC_MIDAS_LASSO_AUC ))
    AUC_LASSO_AUC = rbind(AUC_LASSO_AUC , as.numeric( result$AUC_LASSO_AUC ))

    AUC_MIDAS_bootstrap = rbind(AUC_MIDAS_bootstrap , as.numeric( result$AUC_MIDAS_bootstrap ))
    AUC_MIDAS_LASSO_bootstrap = rbind(AUC_MIDAS_LASSO_bootstrap , as.numeric( result$AUC_MIDAS_LASSO_bootstrap ))
    AUC_LASSO_bootstrap = rbind(AUC_LASSO_bootstrap , as.numeric( result$AUC_LASSO_bootstrap ))
    AUC_MIDAS_AUC_bootstrap = rbind(AUC_MIDAS_AUC_bootstrap , as.numeric( result$AUC_MIDAS_AUC_bootstrap ))
    AUC_MIDAS_LASSO_AUC_bootstrap = rbind(AUC_MIDAS_LASSO_AUC_bootstrap , as.numeric( result$AUC_MIDAS_LASSO_AUC_bootstrap ))
    AUC_LASSO_AUC_bootstrap = rbind(AUC_LASSO_AUC_bootstrap , as.numeric( result$AUC_LASSO_AUC_bootstrap ))
    ######################################################################

  }

  ###########################################
  res_fit_av_MIDAS = colMeans(res_fit_MIDAS)
  res_fit_av_LASSO = colMeans(res_fit_LASSO)
  res_fit_av_MIDAS_LASSO = colMeans(res_fit_MIDAS_LASSO)
  res_fit_av_MIDAS_AUC = colMeans(res_fit_MIDAS_AUC)
  res_fit_av_LASSO_AUC = colMeans(res_fit_LASSO_AUC)
  res_fit_av_MIDAS_LASSO_AUC = colMeans(res_fit_MIDAS_LASSO_AUC)

  AUC_LASSO_bootstrap_av = colMeans(AUC_LASSO_bootstrap)
  AUC_MIDAS_bootstrap_av = colMeans(AUC_MIDAS_bootstrap)
  AUC_MIDAS_LASSO_bootstrap_av = colMeans(AUC_MIDAS_LASSO_bootstrap)
  AUC_LASSO_AUC_bootstrap_av = colMeans(AUC_LASSO_AUC_bootstrap)
  AUC_MIDAS_AUC_bootstrap_av = colMeans(AUC_MIDAS_AUC_bootstrap)
  AUC_MIDAS_LASSO_AUC_bootstrap_av = colMeans(AUC_MIDAS_LASSO_AUC_bootstrap)
  ###########################################

  #Final results for s = 6 years and t = 8,8.5,9 years
  res = list(
    AUC_MIDAS = round(colMeans(AUC_MIDAS),3), AUC_MIDAS_VAR = round(apply(AUC_MIDAS, 2, var),3),
    AUC_LASSO = round(colMeans(AUC_LASSO),3), AUC_LASSO_VAR = round(apply(AUC_LASSO, 2, var),3),
    AUC_MIDAS_LASSO = round(colMeans(AUC_MIDAS_LASSO),3), AUC_MIDAS_LASSO_VAR = round(apply(AUC_MIDAS_LASSO, 2, var),3),

    AUC_MIDAS_AUC = round(colMeans(AUC_MIDAS_AUC),3), AUC_MIDAS_AUC_VAR = round(apply(AUC_MIDAS_AUC, 2, var),3),
    AUC_LASSO_AUC = round(colMeans(AUC_LASSO_AUC),3), AUC_LASSO_AUC_VARC = round(apply(AUC_LASSO_AUC, 2, var),3),
    AUC_MIDAS_LASSO_AUC = round(colMeans(AUC_MIDAS_LASSO_AUC),3), AUC_MIDAS_LASSO_AUC_VAR = round(apply(AUC_MIDAS_LASSO_AUC, 2, var),3),

    AUC_LASSO_bootstrap_av = round( AUC_LASSO_bootstrap_av, 3 ), AUC_MIDAS_bootstrap_av = round(AUC_MIDAS_bootstrap_av, 3 ), AUC_MIDAS_LASSO_bootstrap_av = round(AUC_MIDAS_LASSO_bootstrap_av, 3 ),
    AUC_LASSO_AUC_bootstrap_av = round( AUC_LASSO_AUC_bootstrap_av, 3 ), AUC_MIDAS_AUC_bootstrap_av = round(AUC_MIDAS_AUC_bootstrap_av, 3 ), AUC_MIDAS_LASSO_AUC_bootstrap_av = round(AUC_MIDAS_LASSO_AUC_bootstrap_av, 3 ),

    res_fit_av_MIDAS = res_fit_av_MIDAS, res_fit_av_LASSO = res_fit_av_LASSO, res_fit_av_MIDAS_LASSO = res_fit_av_MIDAS_LASSO,
    res_fit_av_MIDAS_AUC = res_fit_av_MIDAS_AUC, res_fit_av_LASSO_AUC = res_fit_av_LASSO_AUC, res_fit_av_MIDAS_LASSO_AUC = res_fit_av_MIDAS_LASSO_AUC
  )

  # save the table
  table_name1 <- paste0("Application_2_MIDAS_AUC", s, "_AUC_cvm", t, ".csv")
  table_name2 <- paste0("Application_2_MIDAS_AUC", s, "_AUC", t, ".csv")
  table_name3 <- paste0("Application_2_MIDAS_bootstrap", s, "_AUC_cvm", t, ".csv")
  table_name4 <- paste0("Application_2_MIDAS_bootstrap", s, "_AUC", t, ".csv")
  table_name5 <- paste0("Application_2_MIDAS", s, "_cof", t, ".csv")
  table_name6 <- paste0("Application_2_MIDAS", s, "_cof_AUC", t, ".csv")
  table_name7 <- paste0("Application_2_MIDAS_LASSO", s, "_cof", t, ".csv")
  table_name8 <- paste0("Application_2_MIDAS_LASSO", s, "_cof_AUC", t, ".csv")
  table_name9 <- paste0("Application_2_LASSO", s, "_cof", t, ".csv")
  table_name10 <- paste0("Application_2_LASSO", s, "_cof_AUC", t, ".csv")

  write.csv(rbind( c(res$AUC_MIDAS, res$AUC_MIDAS_VAR, res$AUC_LASSO, res$AUC_LASSO_VAR, res$AUC_MIDAS_LASSO, res$AUC_MIDAS_LASSO_VAR) ,
                   c( quantile(res$AUC_MIDAS_bootstrap_av, 0.025), quantile(res$AUC_MIDAS_bootstrap_av, 0.975),
                      quantile(res$AUC_LASSO_bootstrap_av, 0.025), quantile(res$AUC_LASSO_bootstrap_av, 0.975),
                      quantile(res$AUC_MIDAS_LASSO_bootstrap_av, 0.025), quantile(res$AUC_MIDAS_LASSO_bootstrap_av, 0.975)) ) , table_name1)

  write.csv(rbind( c(res$AUC_MIDAS_AUC, res$AUC_MIDAS_AUC_VAR, res$AUC_LASSO_AUC, res$AUC_LASSO_AUC_VAR, res$AUC_MIDAS_LASSO_AUC, res$AUC_MIDAS_LASSO_AUC_VAR) ,
                   c( quantile(res$AUC_MIDAS_AUC_bootstrap_av, 0.025), quantile(res$AUC_MIDAS_AUC_bootstrap_av, 0.975),
                      quantile(res$AUC_LASSO_AUC_bootstrap_av, 0.025), quantile(res$AUC_LASSO_AUC_bootstrap_av, 0.975),
                      quantile(res$AUC_MIDAS_LASSO_AUC_bootstrap_av, 0.025), quantile(res$AUC_MIDAS_LASSO_AUC_bootstrap_av, 0.975)) ) , table_name2)

  write.csv(rbind( res$AUC_MIDAS_bootstrap_av, res$AUC_LASSO_bootstrap_av,  res$AUC_MIDAS_LASSO_bootstrap_av ), table_name3)

  write.csv(rbind( res$AUC_MIDAS_AUC_bootstrap_av, res$AUC_LASSO_AUC_bootstrap_av,  res$AUC_MIDAS_LASSO_AUC_bootstrap_av ), table_name4)

  write.csv(res$res_fit_av_MIDAS, table_name5)

  write.csv(res$res_fit_av_MIDAS_AUC, table_name6)

  write.csv(res$res_fit_av_MIDAS_LASSO, table_name7)

  write.csv(res$res_fit_av_MIDAS_LASSO_AUC, table_name8)

  write.csv(res$res_fit_av_LASSO, table_name9)

  write.csv(res$res_fit_av_LASSO_AUC, table_name10)

}

########################################
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

  # stratified cross validation, see Section 6.2
  index_train = which(train_dataset$censoringtime >= t)
  nfold = 5
  train_dataset = train_dataset[index_train, ]
  index_fold = which(1*(train_dataset$time <= t) * 1*(train_dataset$status == 1) == 1)
  if (length(index_fold) <5){
    foldid = NULL
  }else{
    foldid = c( form_folds(nrow(train_dataset[index_fold,]), nfold), form_folds(nrow(train_dataset[-index_fold,]), nfold) )
  }

  #we have 14 variables which are related to T,C and etc. Then delete them
  X_train = train_dataset[ , 1: (dim(train_dataset)[2]-14) ]
  y_train = 1*(train_dataset$time <= t)
  X_test = as.matrix( test_dataset[ , 1: (dim(train_dataset)[2]-14) ] )
  y_test = 1*( test_dataset$time <= t )
  X_train = apply(X_train, 2, as.numeric)
  X_test = apply(X_test, 2, as.numeric)
  w_train = (1 - 1*( train_dataset$time <= t )*(1 - train_dataset$status))/train_dataset$Gp
  w_test = (1 - 1*( test_dataset$time <= t )*(1 - test_dataset$status))/test_dataset$Gp

  # MIDAS setting
  idx <- 1:jmax
  gindex <- NULL
  Xdw = NULL
  for (z in seq(dim(X_train)[2]/jmax)){
    z_idx <- (1 + (z - 1) * jmax) : (z * jmax)
    Xdw <- cbind(Xdw, X_train[,z_idx] %*% w_fin)
    gindex <- c(gindex, rep(z, times = degree + 1))
  }

  #training dataset which are used to in the SGL-MIDAS
  X_in = as.matrix(cbind(rep(1, nrow(X_train)), Xdw))
  fit_cv_MIDAS = alpha_cv_sparsegl( X_in[,-1], y_train, group = gindex, nlambda = 200, weight = w_train, alpha = alpha, nfolds = 5, foldid = foldid, pred.loss = 'censor', intercept_zero = intercept_zero, standardize = TRUE, AUC = TRUE, data = train_dataset, t = t)
  est_MIDAS = fit_cv_MIDAS$coff
  est_MIDAS_AUC = fit_cv_MIDAS$coff_AUC

  #test dataset
  X_t = NULL
  for (z in seq(dim(X_test)[2]/jmax)){
    z_idx <- (1 + (z - 1)*jmax) : (z * jmax)
    X_t <- cbind(X_t, X_test[,z_idx]%*%w_fin)
  }
  X_te = as.matrix(cbind(rep(1, nrow(X_test)),  X_t))

  #time dependent AUC and bootstrap
  test_predictions_MIDAS = unlist( as.list(plogis(X_te %*% est_MIDAS)) )
  test_predictions_MIDAS_AUC = unlist( as.list(plogis(X_te %*% est_MIDAS_AUC)) )

  AUC_MIDAS = ROC_censor_N( data = test_dataset, prediction = test_predictions_MIDAS, t = t )$AUC
  AUC_MIDAS_AUC = ROC_censor_N( data = test_dataset, prediction = test_predictions_MIDAS_AUC, t = t )$AUC

  AUC_MIDAS_bootstrap =  ROC_N_bootstrap(data = test_dataset, prediction = test_predictions_MIDAS, t = t, sim_number = bootstrap_number )
  AUC_MIDAS_AUC_bootstrap =  ROC_N_bootstrap(data = test_dataset, prediction = test_predictions_MIDAS_AUC, t = t, sim_number = bootstrap_number )


  #LASSO-MIDAS, cross-validation for log-likelihood score
  if (fit_cv_MIDAS$alpha_out == 1){
    fit_cv_MIDAS_LASSO = fit_cv_MIDAS
    test_predictions_MIDAS_LASSO = test_predictions_MIDAS
    est_MIDAS_LASSO = est_MIDAS
  } else {
    fit_cv_MIDAS_LASSO = fit_cv_MIDAS$fit_MIDAS_LASSO
    cesnor_min_index = which(fit_cv_MIDAS_LASSO$cvm == min(fit_cv_MIDAS_LASSO$cvm))[1]
    est_MIDAS_LASSO = unlist(c(fit_cv_MIDAS_LASSO$survival_sparsegl.fit$b0[,cesnor_min_index], fit_cv_MIDAS_LASSO$survival_sparsegl.fit$beta[,cesnor_min_index]))
    test_predictions_MIDAS_LASSO = unlist( as.list(plogis(X_te %*% est_MIDAS_LASSO )) ) #predict(fit_cv_MIDAS_LASSO, X_te[,-1], s = 'lambda.min', type = 'response')
  }
  AUC_MIDAS_LASSO = ROC_censor_N(data = test_dataset, prediction = test_predictions_MIDAS_LASSO, t = t )$AUC

  #LASSO-MIDAS, cross-validation for AUC
  if (fit_cv_MIDAS$alpha_AUC == 1){
    fit_cv_MIDAS_LASSO = fit_cv_MIDAS
    test_predictions_MIDAS_LASSO_AUC = test_predictions_MIDAS_AUC
    est_MIDAS_LASSO_AUC = est_MIDAS_AUC
  } else {
    fit_cv_MIDAS_LASSO = fit_cv_MIDAS$fit_MIDAS_LASSO
    cesnor_min_index = which(fit_cv_MIDAS_LASSO$AUC_censor == max(fit_cv_MIDAS_LASSO$AUC_censor))[1]
    est_MIDAS_LASSO_AUC = unlist(c(fit_cv_MIDAS_LASSO$survival_sparsegl.fit$b0[,cesnor_min_index], fit_cv_MIDAS_LASSO$survival_sparsegl.fit$beta[,cesnor_min_index]))
    test_predictions_MIDAS_LASSO_AUC = unlist( as.list(plogis(X_te %*% est_MIDAS_LASSO_AUC )) ) #predict(fit_cv_MIDAS_LASSO, X_te[,-1], s = 'lambda.min', type = 'response')
  }
  AUC_MIDAS_LASSO_AUC = ROC_censor_N(data = test_dataset, prediction = test_predictions_MIDAS_LASSO_AUC, t = t )$AUC

  # bootstrap for LASSO-MIDAS
  AUC_MIDAS_LASSO_bootstrap =  ROC_N_bootstrap(data = test_dataset, prediction = test_predictions_MIDAS_LASSO, t = t, sim_number = bootstrap_number )
  AUC_MIDAS_LASSO_AUC_bootstrap =  ROC_N_bootstrap(data = test_dataset, prediction = test_predictions_MIDAS_LASSO_AUC, t = t, sim_number = bootstrap_number )

  #LASSO-UMIDAS
  fit_cv_LASSO = cv.survival_sparsegl(X_train, y_train, group = seq(dim(X_train)[2]), nlambda = 200, weight = w_train, asparse = 1, foldid = foldid, nfolds = 5, pred.loss = 'censor', intercept_zero = intercept_zero, standardize = TRUE, maxit = 30000, AUC = TRUE, data = train_dataset, t = t)
  nan_indices <- sapply(fit_cv_LASSO$cvm, is.nan)
  fit_cv_LASSO$cvm[nan_indices] <- 1000
  cesnor_min_index_LASSO = which(fit_cv_LASSO$cvm== min(fit_cv_LASSO$cvm))[1]
  est_LASSO = unlist(c(fit_cv_LASSO$survival_sparsegl.fit$b0[,cesnor_min_index_LASSO], fit_cv_LASSO$survival_sparsegl.fit$beta[,cesnor_min_index_LASSO]))

  nan_indices <- sapply(fit_cv_LASSO$AUC_censor, is.nan)
  fit_cv_LASSO$AUC_censor[nan_indices] <- 0
  cesnor_min_index_LASSO = which(fit_cv_LASSO$AUC_censor== max(fit_cv_LASSO$AUC_censor))[1]
  est_LASSO_AUC = unlist(c(fit_cv_LASSO$survival_sparsegl.fit$b0[,cesnor_min_index_LASSO], fit_cv_LASSO$survival_sparsegl.fit$beta[,cesnor_min_index_LASSO]))

  test_predictions_LASSO = unlist( as.list(plogis( cbind(rep(1, nrow(X_test)), X_test) %*% est_LASSO )) ) # predict(fit_cv_LASSO, X_test, s = 'lambda.min', type = 'response') # unlist( as.list(plogis(X_test %*% est_LASSO )) )
  test_predictions_LASSO_AUC = unlist( as.list(plogis( cbind(rep(1, nrow(X_test)), X_test) %*% est_LASSO_AUC )) ) # predict(fit_cv_LASSO, X_test, s = 'lambda.min', type = 'response') # unlist( as.list(plogis(X_test %*% est_LASSO )) )

  AUC_LASSO = ROC_censor_N(data = test_dataset, prediction = test_predictions_LASSO, t = t )$AUC
  AUC_LASSO_AUC = ROC_censor_N(data = test_dataset, prediction = test_predictions_LASSO_AUC, t = t )$AUC

  AUC_LASSO_bootstrap =  ROC_N_bootstrap(data = test_dataset, prediction = test_predictions_LASSO, t = t, sim_number = bootstrap_number )
  AUC_LASSO_AUC_bootstrap =  ROC_N_bootstrap(data = test_dataset, prediction = test_predictions_LASSO_AUC, t = t, sim_number = bootstrap_number )

  #return results
  results = list(
    AUC_LASSO = round( AUC_LASSO, 3 ), AUC_MIDAS = round(AUC_MIDAS, 3 ), AUC_MIDAS_LASSO = round(AUC_MIDAS_LASSO, 3 ),
    AUC_LASSO_bootstrap = AUC_LASSO_bootstrap, AUC_MIDAS_bootstrap = AUC_MIDAS_bootstrap, AUC_MIDAS_LASSO_bootstrap = AUC_MIDAS_LASSO_bootstrap,
    fit_LASSO = est_LASSO, fit_MIDAS = est_MIDAS, fit_MIDAS_LASSO = est_MIDAS_LASSO,

    AUC_LASSO_AUC = round( AUC_LASSO_AUC, 3 ), AUC_MIDAS_AUC = round(AUC_MIDAS_AUC, 3 ), AUC_MIDAS_LASSO_AUC = round(AUC_MIDAS_LASSO_AUC, 3 ),
    AUC_LASSO_AUC_bootstrap = AUC_LASSO_AUC_bootstrap, AUC_MIDAS_AUC_bootstrap = AUC_MIDAS_AUC_bootstrap, AUC_MIDAS_LASSO_AUC_bootstrap = AUC_MIDAS_LASSO_AUC_bootstrap,
    fit_LASSO_AUC = est_LASSO_AUC, fit_MIDAS_AUC = est_MIDAS_AUC, fit_MIDAS_LASSO_AUC = est_MIDAS_LASSO_AUC, k = k
  )

  ###################################################################################
  result = results

  res_fit_MIDAS = NULL
  res_fit_MIDAS_LASSO = NULL
  res_fit_LASSO = NULL
  res_fit_MIDAS_AUC = NULL
  res_fit_MIDAS_LASSO_AUC = NULL
  res_fit_LASSO_AUC = NULL

  AUC_LASSO = NULL
  AUC_MIDAS = NULL
  AUC_MIDAS_LASSO = NULL
  AUC_LASSO_AUC = NULL
  AUC_MIDAS_AUC = NULL
  AUC_MIDAS_LASSO_AUC = NULL

  AUC_LASSO_bootstrap = NULL
  AUC_MIDAS_bootstrap = NULL
  AUC_MIDAS_LASSO_bootstrap = NULL
  AUC_LASSO_AUC_bootstrap = NULL
  AUC_MIDAS_AUC_bootstrap = NULL
  AUC_MIDAS_LASSO_AUC_bootstrap = NULL

  for (i in seq(1)){

    ######################################################
    res_fit_MIDAS= rbind(res_fit_MIDAS , result$fit_MIDAS)
    res_fit_MIDAS_LASSO= rbind(res_fit_MIDAS_LASSO , result$fit_MIDAS_LASSO)
    res_fit_LASSO= rbind(res_fit_LASSO , result$fit_LASSO)
    res_fit_MIDAS_AUC= rbind(res_fit_MIDAS_AUC , result$fit_MIDAS_AUC)
    res_fit_MIDAS_LASSO_AUC= rbind(res_fit_MIDAS_LASSO_AUC , result$fit_MIDAS_LASSO_AUC)
    res_fit_LASSO_AUC= rbind(res_fit_LASSO_AUC , result$fit_LASSO_AUC)


    AUC_MIDAS = rbind(AUC_MIDAS , as.numeric( result$AUC_MIDAS ))
    AUC_MIDAS_LASSO = rbind(AUC_MIDAS_LASSO , as.numeric( result$AUC_MIDAS_LASSO ))
    AUC_LASSO = rbind(AUC_LASSO , as.numeric( result$AUC_LASSO ))
    AUC_MIDAS_AUC = rbind(AUC_MIDAS_AUC , as.numeric( result$AUC_MIDAS_AUC ))
    AUC_MIDAS_LASSO_AUC = rbind(AUC_MIDAS_LASSO_AUC , as.numeric( result$AUC_MIDAS_LASSO_AUC ))
    AUC_LASSO_AUC = rbind(AUC_LASSO_AUC , as.numeric( result$AUC_LASSO_AUC ))

    AUC_MIDAS_bootstrap = rbind(AUC_MIDAS_bootstrap , as.numeric( result$AUC_MIDAS_bootstrap ))
    AUC_MIDAS_LASSO_bootstrap = rbind(AUC_MIDAS_LASSO_bootstrap , as.numeric( result$AUC_MIDAS_LASSO_bootstrap ))
    AUC_LASSO_bootstrap = rbind(AUC_LASSO_bootstrap , as.numeric( result$AUC_LASSO_bootstrap ))
    AUC_MIDAS_AUC_bootstrap = rbind(AUC_MIDAS_AUC_bootstrap , as.numeric( result$AUC_MIDAS_AUC_bootstrap ))
    AUC_MIDAS_LASSO_AUC_bootstrap = rbind(AUC_MIDAS_LASSO_AUC_bootstrap , as.numeric( result$AUC_MIDAS_LASSO_AUC_bootstrap ))
    AUC_LASSO_AUC_bootstrap = rbind(AUC_LASSO_AUC_bootstrap , as.numeric( result$AUC_LASSO_AUC_bootstrap ))
    ######################################################################

  }

  ###########################################
  res_fit_av_MIDAS = colMeans(res_fit_MIDAS)
  res_fit_av_LASSO = colMeans(res_fit_LASSO)
  res_fit_av_MIDAS_LASSO = colMeans(res_fit_MIDAS_LASSO)
  res_fit_av_MIDAS_AUC = colMeans(res_fit_MIDAS_AUC)
  res_fit_av_LASSO_AUC = colMeans(res_fit_LASSO_AUC)
  res_fit_av_MIDAS_LASSO_AUC = colMeans(res_fit_MIDAS_LASSO_AUC)

  AUC_LASSO_bootstrap_av = colMeans(AUC_LASSO_bootstrap)
  AUC_MIDAS_bootstrap_av = colMeans(AUC_MIDAS_bootstrap)
  AUC_MIDAS_LASSO_bootstrap_av = colMeans(AUC_MIDAS_LASSO_bootstrap)
  AUC_LASSO_AUC_bootstrap_av = colMeans(AUC_LASSO_AUC_bootstrap)
  AUC_MIDAS_AUC_bootstrap_av = colMeans(AUC_MIDAS_AUC_bootstrap)
  AUC_MIDAS_LASSO_AUC_bootstrap_av = colMeans(AUC_MIDAS_LASSO_AUC_bootstrap)
  ###########################################

  #Final results for s = 6 years and t = 8,8.5,9 years
  res = list(
    AUC_MIDAS = round(colMeans(AUC_MIDAS),3), AUC_MIDAS_VAR = round(apply(AUC_MIDAS, 2, var),3),
    AUC_LASSO = round(colMeans(AUC_LASSO),3), AUC_LASSO_VAR = round(apply(AUC_LASSO, 2, var),3),
    AUC_MIDAS_LASSO = round(colMeans(AUC_MIDAS_LASSO),3), AUC_MIDAS_LASSO_VAR = round(apply(AUC_MIDAS_LASSO, 2, var),3),

    AUC_MIDAS_AUC = round(colMeans(AUC_MIDAS_AUC),3), AUC_MIDAS_AUC_VAR = round(apply(AUC_MIDAS_AUC, 2, var),3),
    AUC_LASSO_AUC = round(colMeans(AUC_LASSO_AUC),3), AUC_LASSO_AUC_VARC = round(apply(AUC_LASSO_AUC, 2, var),3),
    AUC_MIDAS_LASSO_AUC = round(colMeans(AUC_MIDAS_LASSO_AUC),3), AUC_MIDAS_LASSO_AUC_VAR = round(apply(AUC_MIDAS_LASSO_AUC, 2, var),3),

    AUC_LASSO_bootstrap_av = round( AUC_LASSO_bootstrap_av, 3 ), AUC_MIDAS_bootstrap_av = round(AUC_MIDAS_bootstrap_av, 3 ), AUC_MIDAS_LASSO_bootstrap_av = round(AUC_MIDAS_LASSO_bootstrap_av, 3 ),
    AUC_LASSO_AUC_bootstrap_av = round( AUC_LASSO_AUC_bootstrap_av, 3 ), AUC_MIDAS_AUC_bootstrap_av = round(AUC_MIDAS_AUC_bootstrap_av, 3 ), AUC_MIDAS_LASSO_AUC_bootstrap_av = round(AUC_MIDAS_LASSO_AUC_bootstrap_av, 3 ),

    res_fit_av_MIDAS = res_fit_av_MIDAS, res_fit_av_LASSO = res_fit_av_LASSO, res_fit_av_MIDAS_LASSO = res_fit_av_MIDAS_LASSO,
    res_fit_av_MIDAS_AUC = res_fit_av_MIDAS_AUC, res_fit_av_LASSO_AUC = res_fit_av_LASSO_AUC, res_fit_av_MIDAS_LASSO_AUC = res_fit_av_MIDAS_LASSO_AUC
  )

  # save the table
  table_name1 <- paste0("Application_2_comparison_MIDAS_AUC", s, "_AUC_cvm", t, ".csv")
  table_name2 <- paste0("Application_2_comparison_MIDAS_AUC", s, "_AUC", t, ".csv")
  table_name3 <- paste0("Application_2_comparison_MIDAS_bootstrap", s, "_AUC_cvm", t, ".csv")
  table_name4 <- paste0("Application_2_comparison_MIDAS_bootstrap", s, "_AUC", t, ".csv")
  table_name5 <- paste0("Application_2_comparison_MIDAS", s, "_cof", t, ".csv")
  table_name6 <- paste0("Application_2_comparison_MIDAS", s, "_cof_AUC", t, ".csv")
  table_name7 <- paste0("Application_2_comparison_MIDAS_LASSO", s, "_cof", t, ".csv")
  table_name8 <- paste0("Application_2_comparison_MIDAS_LASSO", s, "_cof_AUC", t, ".csv")
  table_name9 <- paste0("Application_2_comparison_LASSO", s, "_cof", t, ".csv")
  table_name10 <- paste0("Application_2_comparison_LASSO", s, "_cof_AUC", t, ".csv")

  write.csv(rbind( c(res$AUC_MIDAS, res$AUC_MIDAS_VAR, res$AUC_LASSO, res$AUC_LASSO_VAR, res$AUC_MIDAS_LASSO, res$AUC_MIDAS_LASSO_VAR) ,
                   c( quantile(res$AUC_MIDAS_bootstrap_av, 0.025), quantile(res$AUC_MIDAS_bootstrap_av, 0.975),
                      quantile(res$AUC_LASSO_bootstrap_av, 0.025), quantile(res$AUC_LASSO_bootstrap_av, 0.975),
                      quantile(res$AUC_MIDAS_LASSO_bootstrap_av, 0.025), quantile(res$AUC_MIDAS_LASSO_bootstrap_av, 0.975)) ) , table_name1)

  write.csv(rbind( c(res$AUC_MIDAS_AUC, res$AUC_MIDAS_AUC_VAR, res$AUC_LASSO_AUC, res$AUC_LASSO_AUC_VAR, res$AUC_MIDAS_LASSO_AUC, res$AUC_MIDAS_LASSO_AUC_VAR) ,
                   c( quantile(res$AUC_MIDAS_AUC_bootstrap_av, 0.025), quantile(res$AUC_MIDAS_AUC_bootstrap_av, 0.975),
                      quantile(res$AUC_LASSO_AUC_bootstrap_av, 0.025), quantile(res$AUC_LASSO_AUC_bootstrap_av, 0.975),
                      quantile(res$AUC_MIDAS_LASSO_AUC_bootstrap_av, 0.025), quantile(res$AUC_MIDAS_LASSO_AUC_bootstrap_av, 0.975)) ) , table_name2)

  write.csv(rbind( res$AUC_MIDAS_bootstrap_av, res$AUC_LASSO_bootstrap_av,  res$AUC_MIDAS_LASSO_bootstrap_av ), table_name3)

  write.csv(rbind( res$AUC_MIDAS_AUC_bootstrap_av, res$AUC_LASSO_AUC_bootstrap_av,  res$AUC_MIDAS_LASSO_AUC_bootstrap_av ), table_name4)

  write.csv(res$res_fit_av_MIDAS, table_name5)

  write.csv(res$res_fit_av_MIDAS_AUC, table_name6)

  write.csv(res$res_fit_av_MIDAS_LASSO, table_name7)

  write.csv(res$res_fit_av_MIDAS_LASSO_AUC, table_name8)

  write.csv(res$res_fit_av_LASSO, table_name9)

  write.csv(res$res_fit_av_LASSO_AUC, table_name10)

}



