#============================================================
# Corporate Survival Analysis with Machine Learning Methods
# Import functions for simulation and empirical application
#============================================================


#============================================================
# Reproduce pair-wise difference test
#============================================================




###################KM
data_test = read.csv2("Application_2_MIDAS_bootstrap10_AUC13.csv")
dim(data_test)
data_test[1,1]

MIDAS <- as.numeric(strsplit(data_test[1,1], ",")[[1]])[-1]
LASSO = as.numeric(strsplit(data_test[2,1], ",")[[1]])[-1]
MIDAS_LASSO = as.numeric(strsplit(data_test[3,1], ",")[[1]])[-1]


1-max(length(which(MIDAS > LASSO))/1000)
quantile( (MIDAS - LASSO)/sd(MIDAS - LASSO), 0.1)


####################Comparison

data_test_2 = read.csv2("Application_2_comparison_MIDAS_bootstrap10_AUC13.csv")
dim(data_test_2)

MIDAS_2 <- as.numeric(strsplit(data_test_2[1,1], ",")[[1]])[-1]
LASSO_2 = as.numeric(strsplit(data_test_2[2,1], ",")[[1]])[-1]
MIDAS_LASSO_2 = as.numeric(strsplit(data_test_2[3,1], ",")[[1]])[-1]

1 - max(length(which(MIDAS > MIDAS_2))/1000)
quantile( (MIDAS - MIDAS_2)/sd(MIDAS - MIDAS_2), 0.05)


