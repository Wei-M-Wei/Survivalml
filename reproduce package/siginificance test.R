#============================================================
# Corporate Survival Analysis with Machine Learning Methods
# Import functions for simulation and empirical application
#============================================================


#============================================================
# Reproduce pair-wise difference test
#============================================================


# s = 6 years, t = 8, 8.5, 9 years

#sg-LASSO-M v.s. LASSOU

#
data_test = read.csv2("MIDAS_bootstrap6_AUC8.csv")
dim(data_test)
data_test[1,1]
MIDAS <- as.numeric(strsplit(data_test[1,1], ",")[[1]])[-1]
LASSO = as.numeric(strsplit(data_test[2,1], ",")[[1]])[-1]
MIDAS_LASSO = as.numeric(strsplit(data_test[3,1], ",")[[1]])[-1]
n = length(MIDAS)
p_value = 1-max(length(which(MIDAS > LASSO))/n)

#
data_test = read.csv2("MIDAS_bootstrap6_AUC8.5.csv")
dim(data_test)
data_test[1,1]
MIDAS <- as.numeric(strsplit(data_test[1,1], ",")[[1]])[-1]
LASSO = as.numeric(strsplit(data_test[2,1], ",")[[1]])[-1]
MIDAS_LASSO = as.numeric(strsplit(data_test[3,1], ",")[[1]])[-1]
n = length(MIDAS)
p_value = 1-max(length(which(MIDAS > LASSO))/n)


#
data_test = read.csv2("MIDAS_bootstrap6_AUC9.csv")
dim(data_test)
data_test[1,1]
MIDAS <- as.numeric(strsplit(data_test[1,1], ",")[[1]])[-1]
LASSO = as.numeric(strsplit(data_test[2,1], ",")[[1]])[-1]
MIDAS_LASSO = as.numeric(strsplit(data_test[3,1], ",")[[1]])[-1]
n = length(MIDAS)
p_value = 1-max(length(which(MIDAS > LASSO))/n)


# s = 10 years, t = 13, 13.5, 14 years

#sg-LASSO-M v.s. LASSOU
#
data_test = read.csv2("MIDAS_bootstrap10_AUC13.csv")
dim(data_test)
data_test[1,1]
MIDAS <- as.numeric(strsplit(data_test[1,1], ",")[[1]])[-1]
LASSO = as.numeric(strsplit(data_test[2,1], ",")[[1]])[-1]
MIDAS_LASSO = as.numeric(strsplit(data_test[3,1], ",")[[1]])[-1]
n = length(MIDAS)
p_value = 1-max(length(which(MIDAS > LASSO))/n)


#
data_test = read.csv2("MIDAS_bootstrap10_AUC13.5.csv")
dim(data_test)
data_test[1,1]
MIDAS <- as.numeric(strsplit(data_test[1,1], ",")[[1]])[-1]
LASSO = as.numeric(strsplit(data_test[2,1], ",")[[1]])[-1]
MIDAS_LASSO = as.numeric(strsplit(data_test[3,1], ",")[[1]])[-1]
n = length(MIDAS)
p_value = 1-max(length(which(MIDAS > LASSO))/n)


#
data_test = read.csv2("MIDAS_bootstrap10_AUC14.csv")
dim(data_test)
data_test[1,1]
MIDAS <- as.numeric(strsplit(data_test[1,1], ",")[[1]])[-1]
LASSO = as.numeric(strsplit(data_test[2,1], ",")[[1]])[-1]
MIDAS_LASSO = as.numeric(strsplit(data_test[3,1], ",")[[1]])[-1]
n = length(MIDAS)
p_value = 1-max(length(which(MIDAS > LASSO))/n)




# s = 6 years, t = 8, 8.5, 9 years

#sg-LASSO-M v.s. data without censored firms

#
data_test_2 = read.csv2("comparison_MIDAS_bootstrap6_AUC8.csv")
dim(data_test_2)
MIDAS_2 <- as.numeric(strsplit(data_test_2[1,1], ",")[[1]])[-1]
LASSO_2 = as.numeric(strsplit(data_test_2[2,1], ",")[[1]])[-1]
MIDAS_LASSO_2 = as.numeric(strsplit(data_test_2[3,1], ",")[[1]])[-1]
n = length(MIDAS_2)
P_value = 1 - max(length(which(MIDAS > MIDAS_2))/1000)

#
data_test_2 = read.csv2("comparison_MIDAS_bootstrap6_AUC8.5.csv")
dim(data_test_2)
MIDAS_2 <- as.numeric(strsplit(data_test_2[1,1], ",")[[1]])[-1]
LASSO_2 = as.numeric(strsplit(data_test_2[2,1], ",")[[1]])[-1]
MIDAS_LASSO_2 = as.numeric(strsplit(data_test_2[3,1], ",")[[1]])[-1]
n = length(MIDAS_2)
P_value = 1 - max(length(which(MIDAS > MIDAS_2))/1000)

#
data_test_2 = read.csv2("comparison_MIDAS_bootstrap6_AUC9.csv")
dim(data_test_2)
MIDAS_2 <- as.numeric(strsplit(data_test_2[1,1], ",")[[1]])[-1]
LASSO_2 = as.numeric(strsplit(data_test_2[2,1], ",")[[1]])[-1]
MIDAS_LASSO_2 = as.numeric(strsplit(data_test_2[3,1], ",")[[1]])[-1]
n = length(MIDAS_2)
P_value = 1 - max(length(which(MIDAS > MIDAS_2))/1000)

# s = 10 years, t = 13, 13.5, 14 years

#sg-LASSO-M v.s. data without censored firms
#
data_test_2 = read.csv2("comparison_MIDAS_bootstrap10_AUC13.csv")
dim(data_test_2)
MIDAS_2 <- as.numeric(strsplit(data_test_2[1,1], ",")[[1]])[-1]
LASSO_2 = as.numeric(strsplit(data_test_2[2,1], ",")[[1]])[-1]
MIDAS_LASSO_2 = as.numeric(strsplit(data_test_2[3,1], ",")[[1]])[-1]
n = length(MIDAS_2)
P_value = 1 - max(length(which(MIDAS > MIDAS_2))/1000)

#
data_test_2 = read.csv2("comparison_MIDAS_bootstrap10_AUC13.5.csv")
dim(data_test_2)
MIDAS_2 <- as.numeric(strsplit(data_test_2[1,1], ",")[[1]])[-1]
LASSO_2 = as.numeric(strsplit(data_test_2[2,1], ",")[[1]])[-1]
MIDAS_LASSO_2 = as.numeric(strsplit(data_test_2[3,1], ",")[[1]])[-1]
n = length(MIDAS_2)
P_value = 1 - max(length(which(MIDAS > MIDAS_2))/1000)

#
data_test_2 = read.csv2("comparison_MIDAS_bootstrap10_AUC14.csv")
dim(data_test_2)
MIDAS_2 <- as.numeric(strsplit(data_test_2[1,1], ",")[[1]])[-1]
LASSO_2 = as.numeric(strsplit(data_test_2[2,1], ",")[[1]])[-1]
MIDAS_LASSO_2 = as.numeric(strsplit(data_test_2[3,1], ",")[[1]])[-1]
n = length(MIDAS_2)
P_value = 1 - max(length(which(MIDAS > MIDAS_2))/1000)


# Additional application

# s = 6 years, t = 8, 8.5, 9 years

#sg-LASSO-M v.s. LASSOU

#
data_test = read.csv2("Application_2_MIDAS_bootstrap6_AUC8.csv")
dim(data_test)
data_test[1,1]
MIDAS <- as.numeric(strsplit(data_test[1,1], ",")[[1]])[-1]
LASSO = as.numeric(strsplit(data_test[2,1], ",")[[1]])[-1]
MIDAS_LASSO = as.numeric(strsplit(data_test[3,1], ",")[[1]])[-1]
n = length(MIDAS)
p_value = 1-max(length(which(MIDAS > LASSO))/n)

#
data_test = read.csv2("Application_2_MIDAS_bootstrap6_AUC8.5.csv")
dim(data_test)
data_test[1,1]
MIDAS <- as.numeric(strsplit(data_test[1,1], ",")[[1]])[-1]
LASSO = as.numeric(strsplit(data_test[2,1], ",")[[1]])[-1]
MIDAS_LASSO = as.numeric(strsplit(data_test[3,1], ",")[[1]])[-1]
n = length(MIDAS)
p_value = 1-max(length(which(MIDAS > LASSO))/n)


#
data_test = read.csv2("Application_2_MIDAS_bootstrap6_AUC9.csv")
dim(data_test)
data_test[1,1]
MIDAS <- as.numeric(strsplit(data_test[1,1], ",")[[1]])[-1]
LASSO = as.numeric(strsplit(data_test[2,1], ",")[[1]])[-1]
MIDAS_LASSO = as.numeric(strsplit(data_test[3,1], ",")[[1]])[-1]
n = length(MIDAS)
p_value = 1-max(length(which(MIDAS > LASSO))/n)


# s = 10 years, t = 13, 13.5, 14 years

#sg-LASSO-M v.s. LASSOU
#
data_test = read.csv2("Application_2_MIDAS_bootstrap10_AUC13.csv")
dim(data_test)
data_test[1,1]
MIDAS <- as.numeric(strsplit(data_test[1,1], ",")[[1]])[-1]
LASSO = as.numeric(strsplit(data_test[2,1], ",")[[1]])[-1]
MIDAS_LASSO = as.numeric(strsplit(data_test[3,1], ",")[[1]])[-1]
n = length(MIDAS)
p_value = 1-max(length(which(MIDAS > LASSO))/n)


#
data_test = read.csv2("Application_2_MIDAS_bootstrap10_AUC13.5.csv")
dim(data_test)
data_test[1,1]
MIDAS <- as.numeric(strsplit(data_test[1,1], ",")[[1]])[-1]
LASSO = as.numeric(strsplit(data_test[2,1], ",")[[1]])[-1]
MIDAS_LASSO = as.numeric(strsplit(data_test[3,1], ",")[[1]])[-1]
n = length(MIDAS)
p_value = 1-max(length(which(MIDAS > LASSO))/n)


#
data_test = read.csv2("Application_2_MIDAS_bootstrap10_AUC14.csv")
dim(data_test)
data_test[1,1]
MIDAS <- as.numeric(strsplit(data_test[1,1], ",")[[1]])[-1]
LASSO = as.numeric(strsplit(data_test[2,1], ",")[[1]])[-1]
MIDAS_LASSO = as.numeric(strsplit(data_test[3,1], ",")[[1]])[-1]
n = length(MIDAS)
p_value = 1-max(length(which(MIDAS > LASSO))/n)




# s = 6 years, t = 8, 8.5, 9 years

#sg-LASSO-M v.s. data without censored firms

#
data_test_2 = read.csv2("Application_2_comparison_MIDAS_bootstrap6_AUC8.csv")
dim(data_test_2)
MIDAS_2 <- as.numeric(strsplit(data_test_2[1,1], ",")[[1]])[-1]
LASSO_2 = as.numeric(strsplit(data_test_2[2,1], ",")[[1]])[-1]
MIDAS_LASSO_2 = as.numeric(strsplit(data_test_2[3,1], ",")[[1]])[-1]
n = length(MIDAS_2)
P_value = 1 - max(length(which(MIDAS > MIDAS_2))/1000)

#
data_test_2 = read.csv2("Application_2_comparison_MIDAS_bootstrap6_AUC8.5.csv")
dim(data_test_2)
MIDAS_2 <- as.numeric(strsplit(data_test_2[1,1], ",")[[1]])[-1]
LASSO_2 = as.numeric(strsplit(data_test_2[2,1], ",")[[1]])[-1]
MIDAS_LASSO_2 = as.numeric(strsplit(data_test_2[3,1], ",")[[1]])[-1]
n = length(MIDAS_2)
P_value = 1 - max(length(which(MIDAS > MIDAS_2))/1000)

#
data_test_2 = read.csv2("Application_2_comparison_MIDAS_bootstrap6_AUC9.csv")
dim(data_test_2)
MIDAS_2 <- as.numeric(strsplit(data_test_2[1,1], ",")[[1]])[-1]
LASSO_2 = as.numeric(strsplit(data_test_2[2,1], ",")[[1]])[-1]
MIDAS_LASSO_2 = as.numeric(strsplit(data_test_2[3,1], ",")[[1]])[-1]
n = length(MIDAS_2)
P_value = 1 - max(length(which(MIDAS > MIDAS_2))/1000)

# s = 10 years, t = 13, 13.5, 14 years

#sg-LASSO-M v.s. data without censored firms
#
data_test_2 = read.csv2("Application_2_comparison_MIDAS_bootstrap10_AUC13.csv")
dim(data_test_2)
MIDAS_2 <- as.numeric(strsplit(data_test_2[1,1], ",")[[1]])[-1]
LASSO_2 = as.numeric(strsplit(data_test_2[2,1], ",")[[1]])[-1]
MIDAS_LASSO_2 = as.numeric(strsplit(data_test_2[3,1], ",")[[1]])[-1]
n = length(MIDAS_2)
P_value = 1 - max(length(which(MIDAS > MIDAS_2))/1000)

#
data_test_2 = read.csv2("Application_2_comparison_MIDAS_bootstrap10_AUC13.5.csv")
dim(data_test_2)
MIDAS_2 <- as.numeric(strsplit(data_test_2[1,1], ",")[[1]])[-1]
LASSO_2 = as.numeric(strsplit(data_test_2[2,1], ",")[[1]])[-1]
MIDAS_LASSO_2 = as.numeric(strsplit(data_test_2[3,1], ",")[[1]])[-1]
n = length(MIDAS_2)
P_value = 1 - max(length(which(MIDAS > MIDAS_2))/1000)

#
data_test_2 = read.csv2("Application_2_comparison_MIDAS_bootstrap10_AUC14.csv")
dim(data_test_2)
MIDAS_2 <- as.numeric(strsplit(data_test_2[1,1], ",")[[1]])[-1]
LASSO_2 = as.numeric(strsplit(data_test_2[2,1], ",")[[1]])[-1]
MIDAS_LASSO_2 = as.numeric(strsplit(data_test_2[3,1], ",")[[1]])[-1]
n = length(MIDAS_2)
P_value = 1 - max(length(which(MIDAS > MIDAS_2))/1000)
