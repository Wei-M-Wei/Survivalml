#============================================================
# Corporate Survival Analysis with Machine Learning Methods
# Import functions for simulation and empirical application
#============================================================



# cross-validation for sparse group lasso
alpha_cv_sparsegl = function(X, y, group, nlambda = 50, lambda = NULL, weight , alpha , foldid = NULL , nfolds, pred.loss = 'censor', intercept_zero = 0, standardize = TRUE, AUC = FALSE, data = NULL, t = 0, maxit = 30000){
  cv_error = matrix(nlambda, length(alpha), nlambda )
  AUC_c = matrix(0, length(alpha), nlambda )
  coff = matrix(0,length(alpha), dim(X)[2] + 1)
  coff_AUC = matrix(0,length(alpha), dim(X)[2] + 1)
  lamb = rep(0, nlambda)
  lamb_AUC = rep(0, nlambda)
  alpha_out = rep(0, length(alpha))
  alpha_AUC_out = rep(0, length(alpha))
  for (alp in seq(length(alpha))){
    print(alp)
    if (alpha[alp] == 1){
      fit = cv.survival_sparsegl( X, y, group = seq(dim(X)[2]), nlambda = nlambda, lambda = lambda,  weight = weight, asparse= alpha[alp], foldid = foldid, nfolds = nfolds, pred.loss = 'censor', intercept_zero = intercept_zero, standardize = standardize, AUC = AUC, data = data, t = t, maxit = maxit)
    }
    else{
      fit = cv.survival_sparsegl(X, y, group = group, nlambda = nlambda, lambda = lambda, weight = weight, asparse= alpha[alp], foldid = foldid, nfolds = nfolds, pred.loss = 'censor', intercept_zero = intercept_zero, standardize = standardize, AUC = AUC, data = data, t = t, maxit = maxit)
    }
    est = coef.cv.survival_sparsegl(fit, s = 'lambda.min')
    fit$cvm = na.omit(fit$cvm)
    print(fit$lambda)
    print(fit$cvm)
    AUC_cal = na.omit(fit$AUC_censor)
    print(AUC_cal)
    AUC_max_index = which(AUC_cal == max(AUC_cal))[1]
    cesnor_min_index = which(fit$cvm == min(fit$cvm))[1]
    est = unlist(c(fit$survival_sparsegl.fit$b0[,cesnor_min_index], fit$survival_sparsegl.fit$beta[,cesnor_min_index]))
    est_AUC = unlist(c(fit$survival_sparsegl.fit$b0[,AUC_max_index], fit$survival_sparsegl.fit$beta[,AUC_max_index]))
    cv_error[alp, 1:length(fit$cvm)] = fit$cvm
    AUC_c[alp, 1:length(AUC_cal)] = AUC_cal
    coff[alp, ] = est
    coff_AUC[alp, ] = est_AUC
    lamb[alp] = fit$lambda.min
    lamb_AUC[alp] = fit$lambda[AUC_max_index]
    alpha_out[alp] = alp
  }
  min_index <- which(cv_error == min(cv_error), arr.ind = TRUE)
  max_index = which(AUC_c == max(AUC_c), arr.ind = TRUE)
  res = list( fit_MIDAS_LASSO = fit, coff = coff[min_index[1], ], coff_AUC = coff_AUC[max_index[1], ], alpha_out = alpha[min_index[1]], alpha_AUC = alpha[max_index[1]], lambda_out = lamb[min_index[1]], lambda_AUC = lamb_AUC[max_index[1]], cv_error = cv_error, AUC = AUC_c)
  return(res)
}


lb <- function(degree,a=0,b=1,jmax=NULL,X=NULL){
  if (!is.null(jmax)){
    X <- seq(0,1,length.out=jmax)
  }
  if (is.null(X)){
    stop("X is not provided. Either set X or set jmax.")
  }
  n <- length(X)
  P <- matrix(1,nrow=n,ncol=degree+2)
  Psi <- matrix(1,nrow=n,ncol=degree+1)/ sqrt(b-a)
  P[, 2] <-  2*X/(b-a) - ((b+a)/(b-a))
  if(degree>0){
    for (i in 1:degree){
      P[, i+2]   <- ((2*i+1)/(i+1)) * P[, 2]*P[, i+1] - i/(i+1) * P[, i]
      Psi[, i+1] <- sqrt((2*i + 1) / (b-a)) %*% P[, i+1]
    }
  }
  col.name <- NULL
  for (i in 0:degree)
    col.name <- c(col.name, paste0("poly-degree-",i))
  
  row.name <- NULL
  for (i in 1:n)
    row.name <- c(row.name, paste0("lag-",n))
  
  colnames(Psi) <- col.name
  rownames(Psi) <- row.name
  class(Psi) <- "midasml"
  Psi
  return(Psi)
}


gb <- function(degree,alpha,a=0,b=1,jmax=NULL,X=NULL){
  if (!is.null(jmax)){
    X <- seq(0,1,length.out=jmax)
  }
  if (is.null(X)){
    stop("X is not provided. Either set X or set jmax.")
  }
  n <- length(X)
  P <- matrix(1,nrow=n,ncol=degree+2)
  Psi <- matrix(1,nrow=n,ncol=degree+1)/ sqrt(b-a)
  P[, 2] <-  2*X/(b-a) - ((b+a)/(b-a))
  if(degree>0){
    for (i in 1:degree){
      d <- (2*i + 2*alpha + 1)*(2*i + 2*alpha + 2)/(2*(i+1)*(i + 2*alpha + 1))
      c <- (alpha + i)^2*(2*i + 2*alpha + 2)/( (i+1)*(i + 2*alpha + 1)*(2*i + 2*alpha))
      P[, i+2]   <- d * P[, 2]*P[, i+1] - c * P[, i]
      Psi[, i+1] <- sqrt((2*i + 1) / (b-a)) %*% P[, i+1]
    }
  }
  col.name <- NULL
  for (i in 0:degree)
    col.name <- c(col.name, paste0("poly-degree-",i))
  
  row.name <- NULL
  for (i in 1:n)
    row.name <- c(row.name, paste0("lag-",n))
  
  colnames(Psi) <- col.name
  rownames(Psi) <- row.name
  class(Psi) <- "midasml"
  Psi
  return(Psi)
}

# calculate the Kaplan-meier estimator for censoring time C
testRandomLogitDataset <-
  function(data,t){
    data %>%
      as_tibble() %>%
      mutate(
        T_hat = pmin(data$Ts, data$censoringtime),
        dead =  runif(nrow(data)),
        re=runif(nrow(data)),
        m_time_t = pmin(t,T_hat),
        cp = runif(nrow(data)),
        Gp = runif(nrow(data)),
        pr = runif(nrow(data)))
  }

KM_estimate <-
  function(data){
    data = data[order(data$m_time_t),]
    survival_data <- Surv(data$T_hat, rep(1,nrow(data)) - data$status)
    km_fit <- survfit(survival_data ~ 1)
    data$Gp = summary(km_fit, times = data$m_time_t)$surv
    return(data)
  }

IF_estimate <-
  function(data){
    data = data[order(data$m_time_t),]
    survival_data <- Surv(data$T_hat, rep(1,nrow(data)) - data$status)
    km_fit <- survfit(survival_data ~ 1, influence = TRUE)
    times_fit <- km_fit$time
    IF_mat <- km_fit$influence.surv*dim(data)[1]
    idx <- findInterval(data$m_time_t, times_fit)
    idx[idx == 0] <- 1
    IF_est <- IF_mat[, idx, drop = FALSE]
    return(IF_est)
  }


# cut folds for cross-validation
form_folds <- function(n_obs, nfold){
  folds <- cut(seq(1,n_obs),breaks=nfold,labels=FALSE)
  folds_final <- sample(folds, replace = FALSE)
  return(folds_final)
}


# time-dependent AUC estimator, from package 'survivalROC'
ROC_censor_N = function(data, prediction, t, cut.values = NULL ){
  AUC_N = survivalROC(Stime=data$time, status= data$status, marker = prediction, predict.time = t, span = 0.25*nrow(data)^(-1/2))
  plot(AUC_N$FP, AUC_N$TP, type="l", xlim=c(0,1), ylim=c(0,1), xlab=paste( "FP", "\n", "AUC = ",round(AUC_N$AUC,3)), ylab="TP",main="Method = NNE ")
  return(AUC_N)
}

# Bootstrap for the time-dependent AUC estimator
ROC_N_bootstrap = function(data, prediction, t, sim_number){
  set.seed(123)
  AUC_bootstrap = c()
  data_to_use = data
  data_to_use$predictions = prediction
  for (b in seq(sim_number)){
    bootstrap_sample =  data_to_use[sample(nrow(data_to_use), size = length(data), replace = TRUE), ]
    AUC_N = survivalROC(Stime=bootstrap_sample$time, status= bootstrap_sample$status, marker = bootstrap_sample$predictions, predict.time = t, span = 0.25*nrow(data)^(-1/2))
    AUC_bootstrap = cbind(AUC_bootstrap, AUC_N$AUC)
  }
  return(AUC_bootstrap)
}


# Transfer the format of time
convert_to_ymd <- function(date_str) {
  # Check if the year is at the beginning (yyyy-mm-dd)
  if (grepl("^[0-9]{4}", date_str)) {
    # Year is at the beginning, so the format is yyyy-mm-dd
    date_parsed <- as.Date(date_str, format = "%Y-%m-%d")
  }
  # Check if the format is dd/mm/yyyy (e.g. 30/04/2015)
  else if (grepl("[0-9]{4}$", date_str) && grepl("/", date_str)) {
    # Check if it's mm/dd/yyyy or dd/mm/yyyy based on the middle value
    parts <- strsplit(date_str, "/")[[1]]
    if (as.numeric(parts[1]) > 12) {
      # Day comes first: dd/mm/yyyy
      date_parsed <- as.Date(date_str, format = "%d/%m/%Y")
    } else {
      # Month comes first: mm/dd/yyyy
      date_parsed <- as.Date(date_str, format = "%m/%d/%Y")
    }
  }
  else {
    stop("Unknown date format")  # If neither format is matched
  }
  
  # Return the date in yyyy-mm-dd format with leading zeros for day and month
  return(format(date_parsed, "%Y-%m-%d"))
}


# Merge the macro data with the financial data, 'process.xls' is the macro data
macro_to_be_merge = function(data_financial, lag_use_year){
  
  ########################################initial settings and read the macro data
  quarter = 4
  macro_lags = lag_use_year * quarter
  data_in = read_xls('process.xls')
  macro = c()
  for ( i in seq(2,dim(data_in)[2])){
    macro = cbind(macro,t(data_in[,i]))
  }
  data_macro = as.data.frame(matrix(rep(macro, n), nrow = n, byrow = TRUE))
  
  data_macro$start_day = sapply(data_financial$start_day, convert_to_ymd)
  n = dim(data_financial)[1]
  
  ########################start year and the end year of the covariates
  end_in = 2023
  sat_in = 1985
  
  
  ###################Total number of macro covariates is 98
  macro_num = 98
  index_macro = seq( 1, 1  + ((2023 - 1985) + 1)*4*macro_num - 1, 1 )
  total_macro = length(index_macro)
  
  
  ######################################################merge the macro data
  lags = lag_use_year * 4
  X_g = data.frame(matrix(ncol = 0, nrow = 0))
  data_pre = data_macro
  for (i in seq(n)) {
    nex_year = next_quarter(year(data_pre$start_day[i]), month(data_pre$start_day[i]))$year
    nex_qu = next_quarter(year(data_pre$start_day[i]), month(data_pre$start_day[i]))$quarter
    start_index = (nex_year - sat_in ) * 4 + nex_qu + 4 * (s - lag_use_year)
    if (i == 1){
      X_g = data_extract(data_pre[i, index_macro], start_index, start_index + lags - 1, macro_num)
    }
    else{
      X_g = rbind(as.matrix(X_g), as.matrix(data_extract(data_pre[i, index_macro], start_index, start_index + lags - 1, macro_num))
      )
    }
  }
  X = data.frame(X_g)
  
  
  ###############################
  # delete covariates which have missing values
  empty_columns <- unique(which(colSums(is.na(X_g)) > 0))
  ind_remove = NULL
  for (k in empty_columns) {
    i_index = k / (macro_lags)
    
    if (i_index  %% 1 == 0) {
      ind_remove = c(ind_remove, i_index)
    } else {
      ind_remove = c(ind_remove, floor(i_index) + 1)
    }
  }
  remove <- unique(ind_remove)
  remove_all = NULL
  for (h in remove) {
    id = seq((h-1)*macro_lags + 1, h*macro_lags)
    remove_all = c(remove_all, id)
  }
  X_macro = X_g[, -remove_all]
  
  #Remove the last lag
  lag_of_covariates = (lag_use_year-0.25)
  interval <- lag_of_covariates*quarter
  
  # Get the total number of columns in the data frame
  n_cols <- ncol(X_macro)
  
  # Create a vector to store the selected columns
  selected_columns <- c()
  
  # Generate the column indices with the specified interval
  for (i in seq(1, n_cols, by = interval + 1)) {
    selected_columns <- c(selected_columns, i:(min(i + interval - 1, n_cols)))
  }
  X_macro = X_macro[, selected_columns]
  
  return(X_macro)
  
}

# data preprocess functions
next_quarter <- function(year, month) {
  
  current_quarter <- (month - 1) %/% 3 + 1
  
  if (current_quarter == 4) {
    next_year <- year + 1
    next_quarter <- 1
  } else {
    next_year <- year
    next_quarter <- current_quarter + 1
  }
  
  return(list(year = next_year, quarter = next_quarter))
}

add_column <- function(df, col_name, col_data) {
  if (ncol(df) == 0) {
    # If the data frame is empty, add the new column with column name and data
    df[[col_name]] <- col_data
  } else {
    # If the data frame already has columns, bind the new column to it
    df <- cbind(df, col_data)
    colnames(df)[ncol(df)] <- col_name
  }
  return(df)
}

data_extract = function(data, start, end, total){
  data_e = data.frame(matrix(ncol = 0, nrow = 0))
  end_in = 2023
  sat_in = 1985
  for (i in seq(total)){
    data_ex = data[, (start + (end_in - sat_in + 1)*4*(i-1)) : (end + (end_in - sat_in + 1)*4*(i-1)) ]
    if (i == 1){
      data_e = data_ex
    }
    else{
      data_e = cbind(data_e, data_ex)
    }
  }
  return(data_e)
}

reverse_matrix <- function(mat, p, group_size) {
  indices <- seq_len(p)
  
  # Split the indices into groups of 'group_size'
  groups <- split(indices, ceiling(indices / group_size))
  
  # Reverse each group and combine the indices
  new_order <- unlist(lapply(groups, rev))
  
  # Subset the matrix with the reordered indices
  return(cbind(mat[, new_order], mat[,(p+1):ncol(mat)]))
}

# Simulation section
# grab the prediction metric for simulation
paral_independent = function(result,it){
  res_fit_MIDAS = NULL
  res_fit_MIDAS_LASSO = NULL
  res_fit_LASSO = NULL
  AUC_true = NULL
  AUC_LASSO = NULL
  AUC_MIDAS = NULL
  AUC_MIDAS_LASSO = NULL
  AUC_true_N = NULL
  AUC_LASSO_N = NULL
  AUC_MIDAS_N = NULL
  AUC_MIDAS_LASSO_N = NULL
  PRAUC_true_N = NULL
  PRAUC_LASSO_N = NULL
  PRAUC_MIDAS_N = NULL
  PRAUC_MIDAS_LASSO_N = NULL
  balance = NULL
  censoring_proportions = NULL
  it = length(result)
  for (i in seq(it)){
    res_fit_MIDAS= rbind(res_fit_MIDAS , result[[i]]$fit_MIDAS)
    res_fit_MIDAS_LASSO= rbind(res_fit_MIDAS_LASSO , result[[i]]$fit_MIDAS_LASSO)
    res_fit_LASSO= rbind(res_fit_LASSO , result[[i]]$fit_LASSO)
    AUC_true = rbind(AUC_true , as.numeric( result[[i]]$AUC_true ))
    AUC_MIDAS = rbind(AUC_MIDAS , as.numeric( result[[i]]$AUC_MIDAS ))
    AUC_MIDAS_LASSO = rbind(AUC_MIDAS_LASSO , as.numeric( result[[i]]$AUC_MIDAS_LASSO ))
    AUC_LASSO = rbind(AUC_LASSO , as.numeric( result[[i]]$AUC_LASSO ))
    AUC_true_N = rbind(AUC_true_N , as.numeric( result[[i]]$AUC_true_N ))
    AUC_MIDAS_N = rbind(AUC_MIDAS_N , as.numeric( result[[i]]$AUC_MIDAS_N ))
    AUC_MIDAS_LASSO_N = rbind(AUC_MIDAS_LASSO_N , as.numeric( result[[i]]$AUC_MIDAS_LASSO_N ))
    AUC_LASSO_N = rbind(AUC_LASSO_N , as.numeric( result[[i]]$AUC_LASSO_N ))
    censoring_proportions = rbind(censoring_proportions,  as.numeric( result[[i]]$censoring_proportions ))
    balance  = rbind(balance , as.numeric( result[[i]]$balance ))
    PRAUC_true_N = rbind(PRAUC_true_N , as.numeric( result[[i]]$PRAUC_true_N ))
    PRAUC_MIDAS_N = rbind(PRAUC_MIDAS_N , as.numeric( result[[i]]$PRAUC_MIDAS_N ))
    PRAUC_MIDAS_LASSO_N = rbind(PRAUC_MIDAS_LASSO_N , as.numeric( result[[i]]$PRAUC_MIDAS_LASSO_N ))
    PRAUC_LASSO_N = rbind(PRAUC_LASSO_N , as.numeric( result[[i]]$PRAUC_LASSO_N ))
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
             AUC_true = round(colMeans(AUC_true),3), AUC_true_VAR = round(apply(AUC_true, 2, var),3),  AUC_MIDAS = round(colMeans(AUC_MIDAS),3), AUC_MIDAS_VAR = round(apply(AUC_MIDAS, 2, var),3), AUC_LASSO = round(colMeans(AUC_LASSO),3), AUC_LASSO_VAR = round(apply(AUC_LASSO, 2, var),3), AUC_MIDAS_LASSO = round(colMeans(AUC_MIDAS_LASSO),3), AUC_MIDAS_LASSO_VAR = round(apply(AUC_MIDAS_LASSO, 2, var),3),
             PRAUC_true_N = round(colMeans(PRAUC_true_N),3), PRAUC_true_N_VAR = round(apply(PRAUC_true_N, 2, var),3),  PRAUC_MIDAS_N = round(colMeans(PRAUC_MIDAS_N),3), PRAUC_MIDAS_N_VAR = round(apply(PRAUC_MIDAS_N, 2, var),3), PRAUC_LASSO_N = round(colMeans(PRAUC_LASSO_N),3), PRAUC_LASSO_N_VAR = round(apply(PRAUC_LASSO_N, 2, var),3), PRAUC_MIDAS_LASSO_N = round(colMeans(PRAUC_MIDAS_LASSO_N),3), PRAUC_MIDAS_LASSO_N_VAR = round(apply(PRAUC_MIDAS_LASSO_N, 2, var),3),
             w1_es_MSE = w1_es_MSE, w1_es_var = w1_es_var,
             w2_es_MSE = w2_es_MSE, w2_es_var = w2_es_var,
             w1_es_MSE_LASSO = w1_es_MSE_LASSO, w1_es_var_LASSO = w1_es_var_LASSO,
             w2_es_MSE_LASSO = w2_es_MSE_LASSO, w2_es_var_LASSO = w2_es_var_LASSO,
             w1_es_MSE_MIDAS_LASSO = w1_es_MSE_MIDAS_LASSO, w1_es_var_MIDAS_LASSO = w1_es_var_MIDAS_LASSO,
             w2_es_MSE_MIDAS_LASSO = w2_es_MSE_MIDAS_LASSO, w2_es_var_MIDAS_LASSO = w2_es_var_MIDAS_LASSO,
             censoring_proportions = colMeans(censoring_proportions), balance = colMeans(balance))
  return(res)
}

#check the censoring rate
test_censoring_AR = function(s = s, n, censor_strength){
  res = 0
  for (i in seq(100)){
    data =generateData_AR( s = s, numberOfObservations = n, numhv= numtrue, numtrue=numtrue, degree=degree, jmax=jmax, parameters = beta_true, censor_strength = censor_strength) # 3, 500, 0.7
    res = res + length(which(data$status==0))/n
  }
  return(res/100)
}



###########################################################################
# De-biased sg-LASSO estimator for GLM with observation-level influence function adjustment

ORIG_DS_inf <- function(x, y, family, lasso_est, weight, nfold=5, n_lambda=50,
                        lambda_ratio=0.005, interest_number = NULL, foldid, variance_km) {

  nn <- length(y)
  pp <- ncol(x)
  X <- cbind(rep(1, nrow(x)), x)

  # Validate input dimensions
  if(ncol(X) != length(lasso_est)) {
    stop("The length of lasso_est is incompatible with the covariate matrix.")
  }

  # Compute gradient and Hessian of negative log-likelihood
  if(family == "binomial") {
    mu <- as.vector(exp(X%*%lasso_est)/(1+exp(X%*%lasso_est)))
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(weight*y - mu))/nn
    neg_ddloglik_glmnet <- t(X)%*%diag(mu*(1-mu))%*%X/nn
    C_glmnet <- sqrt(diag(mu*(1-mu))/nn)%*%X
  } else if(family == "poisson") {
    mu <- as.vector(exp(X%*%lasso_est))
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(weight*y - mu))/nn
    neg_ddloglik_glmnet <- t(X)%*%diag(mu)%*%X/nn
    C_glmnet <- sqrt(diag(mu)/nn)%*%X
  } else {
    stop("Input family is not supported.")
  }

  # Determine which parameters to de-bias
  n_debias <- if(is.null(interest_number)) pp+1 else interest_number+1

  # Nodewise LASSO for Riesz representor estimation
  theta_glmnet <- diag(pp+1)
  tau_glmnet <- rep(NA, pp+1)

  for(j in 1:n_debias) { # nodewise lasso for parameters of interest
    current_x <- sqrt(nn)*C_glmnet[, -j]
    current_y <- sqrt(nn)*as.vector(C_glmnet[, j])

    # Cross-validated LASSO
    gamma_j_glmnet <- cv.glmnet(x = current_x, y = current_y,
                                family = "gaussian", alpha = 1,
                                standardize = TRUE, intercept = FALSE,
                                nfolds = nfold, nlambda = n_lambda, foldid = foldid)
    gamma_j_glmnet <- as.vector(coef(gamma_j_glmnet, s = "lambda.min"))[-1]

    # Update Riesz representor
    theta_glmnet[j, -j] <- (-1)*t(gamma_j_glmnet)
    tau_glmnet[j] <- as.numeric(neg_ddloglik_glmnet[j, j] -
                                  neg_ddloglik_glmnet[j, -j]%*%gamma_j_glmnet)
  }

  # Scale Riesz representor by diagonal elements (only for de-biased parameters)
  theta_glmnet[1:n_debias, ] <- diag(1/tau_glmnet[1:n_debias])%*%theta_glmnet[1:n_debias, ]

  # De-biased estimator
  b_hat_nw <- as.vector(lasso_est - theta_glmnet%*%neg_dloglik_glmnet)

  # Compute score contribution for each observation: sigma_i = X_i * (-weight*y + mu) - variance_km
  # Shape: N x (pp+1)
  sigma <- as.matrix(X) * matrix((-weight * y + mu), nrow = nrow(X), ncol = ncol(X)) - variance_km

  # Compute covariance matrix: (1/N) * sum_i(sigma_i * sigma_i^T)
  # Using t(sigma) %*% sigma / nn to efficiently compute the outer product sum
  cov_matrix <- t(sigma) %*% sigma / nn

  # Compute sandwich variance estimate
  var_matrix <- theta_glmnet %*% cov_matrix %*% t(theta_glmnet)

  # Standard errors for all parameters
  se_nw <- sqrt(diag(var_matrix) / nn)

  # P-values (two-sided test)
  pval_nw <- 2*pnorm(abs(b_hat_nw/se_nw), lower.tail = FALSE)

  return(list(est = b_hat_nw,
              se = se_nw,
              pvalue = pval_nw,
              theta = theta_glmnet,
              sigma = var_matrix))
}


export_env <- ls(globalenv())


