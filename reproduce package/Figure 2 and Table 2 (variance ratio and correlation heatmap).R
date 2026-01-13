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
library(mvtnorm)

po = c()
cor1 = 0
cor2 = 0
var1 = 0
var2 = 0

# scenario 1
# data generating process, scenario 1
generateData_AR <-
  function(s , numberOfObservations, numhv, numtrue, degree, jmax, parameters = NA, censor_strength){
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

    rho = 0.1
    # degree of time series dependence among its lag
    corr_matrix <- diag(1, p)

    for (v in 1:(p)) {
      for (u in 1:(p)) {
        corr_matrix[v, u] <- (0.1)^abs(v - u)
      }
    }

    for ( t in 1:(jmax)){
      if (t == 1){
        Xdd[,t,(1:p)] <- MASS::mvrnorm(TN, mu = rep(0, p), corr_matrix) #scenario 1
      }
      else{
        Xdd[,t,(1:p)] <- rho * Xdd[,t-1,(1:p)] + MASS::mvrnorm(TN, mu = rep(0, p), corr_matrix*(1 - rho^2))
      }
    }

    # transform to its absolute value, as described in the paper
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

# calculate the censoring rate of a simulation dataset
test_censoring_AR = function(s , n, censor_strength){
  res = 0
  for (i in seq(100)){
    data =generateData_AR( s = s, numberOfObservations = n, numhv= numtrue, numtrue=numtrue, degree=degree, jmax=jmax, parameters = beta_true, censor_strength = censor_strength) # 3, 500, 0.7
    res = res + length(which(data$status==0))/n
  }
  return(res/100)
}
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
alpha = c(1)

#censoring strength gamma, which can generate approximately 81% censoring
censor = 1.68 # 3.9
test_censoring_AR( s = s, n = 1200, censor_strength = censor)

for (j in seq(100)) {
  n = 800
  set.seed(j)
  print(j)
  # data prepare
  data = generateData_AR( s = s, numberOfObservations = n, numhv = numhv, numtrue = numtrue, degree = degree, jmax = jmax, parameters = beta_true, censor_strength = censor) # 3, 500, 0.7
  dataset_p = data[ which( (1*( data$time <= t ))*(1*( data$status == 1 )) == 1), ]
  dataset_n = data[ which( (1*( data$time <= t ))*(1*( data$status == 1 )) == 0), ]
  dataset = rbind(dataset_p, dataset_n )
  dataset = testRandomLogitDataset(dataset, t = t )
  dataset = KM_estimate(dataset)

  #5-fold cross validation
  nfold = 5
  index_fold = which(1*(dataset$time <= t) * 1*(dataset$status == 1) == 1)
  #dataset for LASSO-UMIDAS
  X_orginal = dataset$X_orginal[ , 1:(1+jmax*numhv) ]
  cor1 = cor1 + cor(X_orginal[,(jmax + 2):(2*jmax+1)])
  var1 = var1 + var(X_orginal[,2])
  #dataset for sg-LASSO-MIDAS
  X = dataset[ , 1:(numhv*(degree+1)+1) ]
  cor2 = cor2 + cor(X[,2:4])
  var2 = var2 + apply(X[,2:4],2,var)
  po = c(po, as.numeric(length(index_fold)))
}
var_scenario1 = var2/500
c2_1 = cor2/100



########scenario 2

# data generating process, scenario 2
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

    rho = 0.9
    # degree of time series dependence among its lag
    corr_matrix <- diag(1, p)

    for (v in 1:(p)) {
      for (u in 1:(p)) {
        corr_matrix[v, u] <- (phi)^abs(v - u)
      }
    }

    for ( t in 1:(jmax)){
      if (t == 1){
        Xdd[,t,(1:p)] <- MASS::mvrnorm(TN, mu = rep(0, p), corr_matrix) #scenario 1
      }
      else{
        Xdd[,t,(1:p)] <- rho * Xdd[,t-1,(1:p)] + MASS::mvrnorm(TN, mu = rep(0, p), corr_matrix*(1 - rho^2))
      }
    }

    # transform to its absolute value, as described in the paper
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
alpha = c(1)

#censoring strength gamma, which can generate approximately 81% censoring
censor = 1.7 # 3.9
test_censoring_AR( s = s, n = 1000, censor_strength = censor)

po = c()
cor1 = 0
cor2 = 0
var1 = 0
var2 = 0

for (j in seq(100)) {
  set.seed(j)
  n = 800
  print(j)
  # data prepare
  data = generateData_AR( s = s, numberOfObservations = n, numhv = numhv, numtrue = numtrue, degree = degree, jmax = jmax, parameters = beta_true, censor_strength = censor) # 3, 500, 0.7
  dataset_p = data[ which( (1*( data$time <= t ))*(1*( data$status == 1 )) == 1), ]
  dataset_n = data[ which( (1*( data$time <= t ))*(1*( data$status == 1 )) == 0), ]
  dataset = rbind(dataset_p, dataset_n )
  dataset = testRandomLogitDataset(dataset, t = t )
  dataset = KM_estimate(dataset)

  #5-fold cross validation
  nfold = 5
  index_fold = which(1*(dataset$time <= t) * 1*(dataset$status == 1) == 1)
  #dataset for LASSO-UMIDAS
  X_orginal = dataset$X_orginal[ , 1:(1+jmax*numhv) ]
  cor1 = cor1 + cor(X_orginal[,(jmax + 2):(2*jmax+1)])
  var1 = var1 + var(X_orginal[,2])
  #dataset for sg-LASSO-MIDAS
  X = dataset[ , 1:(numhv*(degree+1)+1) ]
  cor2 = cor2 + cor(X[,2:4])
  var2 = var2 + apply(X[,2:4],2,var)
  po = c(po, as.numeric(length(index_fold)))
}
var_scenario2 = var2/100
c2_2 = cor2/100

# variance ratio
print(var_scenario2/var_scenario1)


# correlation plot
library(patchwork)
library(ggplot2)
library(reshape2)
# First correlation matrix
c2_1 <- round(c2_1,2)
cov_names <- c("Element 1", "Element 2", "Element 3")
colnames(c2_1) <- rownames(c2_1) <- cov_names
c2_1_melt <- melt(c2_1)

# Second correlation matrix (example)
c2_2 <- round(c2_2,2)
colnames(c2_2) <- rownames(c2_2) <- cov_names
c2_2_melt <- melt(c2_2)

# Shared fill scale
fill_scale <- scale_fill_gradient2(
  low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0,
  limits = c(-1,1), name = "Correlation"
)


compact_theme <- theme(
  plot.margin = margin(t = 2, r = 2, b = 2, l = 2), # minimal margins
  panel.spacing = unit(1, "mm")
)

# Heatmap 1
heatmap1 <- ggplot(c2_1_melt, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(color="grey70", size=0.3) +
  fill_scale +
  geom_text(aes(label=sprintf("%.3f", value)), color="black", size=3.5) +
  labs(title="Scenario 1", x=NULL, y=NULL) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, face="bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(size=9, face="bold", hjust=0.5),
    legend.position = "none"
  ) +
  coord_fixed(ratio = 1.5) +
  scale_y_discrete(expand = c(0,0)) +  # remove extra space at bottom/top
  scale_x_discrete(expand = c(0,0)) +
  compact_theme

# Heatmap 2
heatmap2 <- ggplot(c2_2_melt, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(color="grey70", size=0.3) +
  fill_scale +
  geom_text(aes(label=sprintf("%.3f", value)), color="black", size=3.5) +
  labs(title="Scenario 2", x=NULL, y=NULL) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, face="bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(size=9, face="bold", hjust=0.5),
    legend.position = "right"
  ) +
  coord_fixed(ratio = 1.5) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  compact_theme

# Combine heatmaps
final_plot <- heatmap1 + heatmap2 +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Display
final_plot
