# Step 1: Generate Simulated Data
library(glmnet)
library(Survivalml)
set.seed(123)  # For reproducibility

# Define the true parameters
beta_true <- c(-1.5, 0.8, -0.6)  # Intercept, coefficient for x1, coefficient for x2

# Simulate predictor variables (x1 and x2)
n <- 500  # Number of observations
x1 <- rnorm(n, mean = 2, sd = 1)
x2 <- rnorm(n, mean = -1, sd = 1)

# Create a design matrix with an intercept (column of 1s)
x <- cbind(1, x1, x2)

# Calculate the linear predictor
linear_predictor <- x %*% beta_true

# Generate the probabilities using the logistic function
p <- 1 / (1 + exp(-linear_predictor))

# Simulate the binary outcome y based on the probabilities
y <- rbinom(n, size = 1, prob = p)

# weight
wei = c(rep(1, 50), rep(1, 100), rep(1,n-150))

# Step 2: Use optim() to Estimate Parameters
# Negative log-likelihood function for logistic regression
neg_log_likelihood <- function(beta) {
  linear_predictor <- x %*% beta
  p <- 1 / (1 + exp(-linear_predictor))

  # Calculate the negative log-likelihood
  -sum(wei*y * log(p) + (1 - wei*y) * log(1 - p))
}

# Initial guesses for beta
initial_params <- rep(0, ncol(x))

# Use optim() to find the MLE estimates for beta
result_optim <- optim(par = initial_params, fn = neg_log_likelihood, method = "BFGS")

# Estimated coefficients from optim
beta_optim <- result_optim$par
print("Estimated coefficients using optim:")
print(beta_optim)

# Estimated coefficients from survivalml
fit2 = survival_sparsegl(x[,-1], y, group = seq(2), nlambda = 10, lambda = c(0), weight = wei, asparse = 1, intercept_zero = -1, standardize = TRUE, eps = 1e-8, maxit = 10000000,intercept = TRUE)
fit2$beta

# compare fit2$beta and beta_optim, we find the difference is very little
fit2$beta - beta_optim[2:3]


# comparison between glmenet and survivalml, we set alpha = 1, lambda = c(0.01, 0.08)

# Step 3: Use glm() to Estimate Parameters
# Fit the logistic regression model using glm()
model_glm <- glmnet(x[,-1], y, alpha = 1, family = 'binomial', standardize = TRUE, nlambda = 10, lambda = c(0.01,0.08))
model_glm$lambda

# Estimated coefficients from glm
beta_glm <- coef(model_glm)
print("Estimated coefficients using glm:")
print(beta_glm)


fit2 = survival_sparsegl(x[,-1], y, group = seq(2), nlambda = 10, lambda = c(0.01,0.08), weight = wei, asparse = 1, intercept_zero = -1, standardize = TRUE, eps = 1e-8, maxit = 10000000,intercept = TRUE)
fit2$beta

# compare fit2$beta and beta_glm, we find the difference is very little
fit2$beta - beta_glm[2:3,]
