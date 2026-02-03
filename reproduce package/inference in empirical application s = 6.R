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
library(glmnet)


#empirical settings
s = 6
lag_use_year = 6
quarter = 4
lags = lag_use_year * quarter - 1 ############we don't have the information of the latest lag.
jmax <- lags



# read the data and check basic information
#### notice that if you want to replicate, you have to import the sub-dataset which is extracted by the algorithm in te paper.
data_financial = read.csv2("period_final_new6.csv")  # read function from package readxl, the name of the dataset file is 'period_final_new6.csv'. Change the name as you want.
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
K = p/jmax
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

# The polynomial is symmteric
w_fin <- gb(degree = degree, alpha = -1/2, a = 0, b = 1, jmax = jmax)/jmax
t(w_fin) %*% w_fin

# tuning parameter grid search
alpha = c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1)
sig_level = 0.05
v <- qnorm(sig_level/2, lower.tail=F)

# prediction horizons t = 8, 8.5, 9 years, Estimator(5) in this paper
rej_vector = c()
for (t in c( 8.5, 9)) {
    dataset = testRandomLogitDataset( data, t = t )
    dataset = KM_estimate(dataset)
    # stratified cross validation, see Section 6.2
    index_train = which(dataset$censoringtime >= t)
    nfold = 5
    index_fold = which(1*(dataset$time <= t) * 1*(dataset$status == 1) == 1)
    if (length(index_fold) <5){
      foldid = NULL
    }else{
      foldid = c( form_folds(nrow(dataset[index_fold,]), nfold), form_folds(nrow(dataset[-index_fold,]), nfold) )
    }

    #we have 14 variables which are related to T,C and etc. Then delete them
    X = dataset[ , 1: (dim(dataset)[2]-14) ]
    X = apply(X, 2, as.numeric)
    y = 1*(dataset$time <= t)
    w_train = 1*( pmin(dataset$Ts, t) <= dataset$censoringtime)/dataset$Gp
    IF = IF_estimate(dataset)
    w_km = 1*( pmin(dataset$Ts, t) <= dataset$censoringtime)/dataset$Gp^2 * y
    V_weighted <- IF * matrix(rep(w_km, each = n), n, n)
    # MIDAS setting
    idx <- 1:jmax
    gindex <- NULL
    Xdw = NULL
    for (z in seq(dim(X)[2]/jmax)){
      z_idx <- (1 + (z - 1) * jmax) : (z * jmax)
      Xdw <- cbind(Xdw, X[,z_idx] %*% w_fin)
      gindex <- c(gindex, rep(z, times = degree + 1))
    }
    # dataset which are used to in the SGL-MIDAS
    X_in = as.matrix(cbind(rep(1, nrow(X)), Xdw))
    var_km <- (V_weighted %*% as.matrix(X_in))/(n)
    fit_cv_MIDAS = alpha_cv_sparsegl( X_in[,-1], y, group = gindex, nlambda = 50, weight = w_train, alpha = alpha, nfolds = 5, foldid = foldid, pred.loss = 'censor', intercept_zero = intercept_zero, standardize = TRUE, AUC = TRUE, data = dataset, t = t)
    est_MIDAS = as.vector( fit_cv_MIDAS$coff)
    ref_MIDAS = ORIG_DS_inf(x=as.matrix(X_in[,-1]), y = y, family="binomial", lasso_est = est_MIDAS, weight = w_train, foldid = foldid, variance_km = var_km)
    ci_upper_orig_ds <- ref_MIDAS$est + v*ref_MIDAS$se
    ci_lower_orig_ds <- ref_MIDAS$est - v*ref_MIDAS$se
    cover = 1*(0 <= ci_upper_orig_ds) * 1*(0>=ci_lower_orig_ds)
    reject_MIDAS = NULL
    for (k in seq(K)){
    test_stat_MIDAS = n * t(ref_MIDAS$est[((k - 1) * (degree + 1) + 1 + 1):(k * (degree + 1) + 1)])%*% solve(ref_MIDAS$sigma[(((k - 1) * (degree + 1) + 1 + 1):(k * (degree + 1) + 1)),(((k - 1) * (degree + 1) + 1 + 1):(k * (degree + 1) + 1))])%*% ref_MIDAS$est[((k - 1) * (degree + 1) + 1 + 1):(k * (degree + 1) + 1)]
    reject_MIDAS[k] = 1*(test_stat_MIDAS > qchisq(p = 1 - sig_level, df = degree+1))
    }
    rej_vector = rbind(rej_vector, reject_MIDAS)
  }

table_name <- paste0("inference result for 6 years",  ".csv")
write.table(rej_vector, table_name, sep = ";", row.names = FALSE, quote = FALSE)

library(reshape2)
library(ggplot2)

rej_vector = read.csv('inference result for 6 years.csv')
### ----------------------------------------------------------
### 1. Column labels (time horizons)
### ----------------------------------------------------------
col_labels <- c("8 years", "8.5 years", "9 years")

### ----------------------------------------------------------
### 2. Variable names (32 variables)
### ----------------------------------------------------------
names_pic = c("Accounts receivable turnover ratio",
              "Current assets turnover ratio"
              , "Fixed assets turnover rate"
              ,"Total assets turnover ratio"
              , "Current ratio"
              ,"Quick ratio"
              ,"Equity ratio"
              ,"Total tangible assetss / Total liabilities"
              ,"Tangible assetss / Net debt"
              ,"Earnings before interest, tax, depreciation, and amortization / Total liabilities"
              ,"Cash flow debt ratio"

              ,"Return on total assetss ROA"
              ,"Net profit on assetss"
              ,"Net profit / Total operating income"
              ,"Earnings before interest and taxes / Total operating income"

              ,"X1 Working Capital / Total Assets"
              ,"X2 Retained Earnings / Total Assets"
              ,"X3 Earnings Before Interest and Taxes / Total Assets"
              ,"X4 Market Value of Equity / Book Value of Total Liabilities"
              ,"X5 Sales / Total Assets"

              ,"Total shareholders equity / Total liabilities"
              ,"Debt ratio"
              ,"Interest bearing debt ratio"
              ,"Equity multiplier"
              ,"Current assetss / Total assetss"
              ,"Current liabilities / Total liabilities"

              ,"Net cash flow from operating activities per share"
              ,"Operating income per share"
              ,"Profit before tax per share"
              ,"Net assetss per share BPS"

              ,"Net cash flow from operating activities / Operating income"
              ,"Net operating Cash flow / Operating income")

### ----------------------------------------------------------
### 3. Your three indicator vectors (length 32 each)
### ----------------------------------------------------------
v1 <- as.numeric(strsplit(rej_vector[1,], ";")[[1]])# c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0)
v2 <- as.numeric(strsplit(rej_vector[2,], ";")[[1]])# c(0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0)
v3 <- as.numeric(strsplit(rej_vector[3,], ";")[[1]])# c(0,0,0,0,1,0,0,1,0,1,0,0,1,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0)

### ----------------------------------------------------------
### 4. Create matrix for the heatmap
### ----------------------------------------------------------
cof_analysis <- cbind(v1, v2, v3)
colnames(cof_analysis) <- col_labels
rownames(cof_analysis) <- names_pic

### ----------------------------------------------------------
### 5. Convert to long format for ggplot
### ----------------------------------------------------------
melted_data <- melt(cof_analysis)
colnames(melted_data) <- c("Variables", "Time", "Value")

### ----------------------------------------------------------
### 6. Draw heatmap (1 = red, 0 = white)
### ----------------------------------------------------------
ggplot(melted_data, aes(Time, Variables, fill = factor(Value))) +  # Use factor for discrete values
  geom_tile() +  # Add border to tiles
  scale_fill_manual(values = c("0" = "white", "1" = "blue")) +  # Custom colors
  labs(
    x = "Prediction horizon",
    y = "Variable") +
  theme_minimal() +
  theme(legend.position = "none")
