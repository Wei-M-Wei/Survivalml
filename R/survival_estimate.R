#' @description
#' Fits regularization paths for sparse group-lasso penalized learning problems at a
#' sequence of regularization parameters `lambda`.
#' Note that the objective function for maximum likely hood function is
#' \deqn{MSE + \lambda penalty}
#' Users can also tweak the penalty by choosing a different penalty factor.
#'
#'
#' @param x Double. A matrix of predictors, of dimension
#'   \eqn{n \times p}{n * p}; each row
#'   is a vector of measurements and each column is a feature. Objects of class
#'   [`Matrix::sparseMatrix`] are supported.
#' @param y Double/Integer/Factor. Either a factor with two levels or
#'   a vector of integers taking 2 unique values. For a factor, the last level
#'   in alphabetical order is the target class.
#' @param weight Double vector. It represents the weighting vector of censored individuals.
#' @param interceptzero Double/Integer. The initial intercept into the optimization.
#' @param group Integer. A vector of consecutive integers describing the
#'   grouping of the coefficients (see example below).
#' @param nlambda The number of \code{lambda} values - default is 100.
#' @param lambda.factor A multiplicative factor for the minimal lambda in the
#'   `lambda` sequence, where `min(lambda) = lambda.factor * max(lambda)`.
#'   `max(lambda)` is the smallest value of `lambda` for which all coefficients
#'   are zero. The default depends on the relationship between \eqn{n}
#'   (the number of rows in the matrix of predictors) and \eqn{p}
#'   (the number of predictors). If \eqn{n \geq p}, the
#'   default is `0.0001`.  If \eqn{n < p}, the default is `0.01`.
#'   A very small value of `lambda.factor` will lead to a
#'   saturated fit. This argument has no effect if there is user-defined
#'   `lambda` sequence.
#' @param lambda A user supplied `lambda` sequence. The default, `NULL`
#'   results in an automatic computation based on `nlambda`, the smallest value
#'   of `lambda` that would give the null model (all coefficient estimates equal
#'   to zero), and `lambda.factor`. Supplying a value of `lambda` overrides
#'   this behaviour. It is likely better to supply a
#'   decreasing sequence of `lambda` values than a single (small) value. If
#'   supplied, the user-defined `lambda` sequence is automatically sorted in
#'   decreasing order.
#' @param pf_group Penalty factor on the groups, a vector of the same
#'   length as the total number of groups. Separate penalty weights can be applied
#'   to each group of \eqn{\beta}{beta's}s to allow differential shrinkage.
#'   Can be 0 for some
#'   groups, which implies no shrinkage, and results in that group always being
#'   included in the model (depending on `pf_sparse`). Default value for each
#'   entry is the square-root of the corresponding size of each group.
#'   Because this default is typical, these penalties are not rescaled.
#' @param pf_sparse Penalty factor on l1-norm, a vector the same length as the
#'   total number of columns in `x`. Each value corresponds to one predictor
#'   Can be 0 for some predictors, which
#'   implies that predictor will be receive only the group penalty.
#'   Note that these are internally rescaled so that the sum is the same as
#'   the number of predictors.
#' @param dfmax Limit the maximum number of groups in the model. Default is
#'   no limit.
#' @param pmax Limit the maximum number of groups ever to be nonzero. For
#'   example once a group enters the model, no matter how many times it exits or
#'   re-enters model through the path, it will be counted only once.
#' @param eps Convergence termination tolerance. Defaults value is `1e-8`.
#' @param maxit Maximum number of outer-loop iterations allowed at fixed lambda
#'   value. Default is `3e8`. If models do not converge, consider increasing
#'   `maxit`.
#' @param intercept Whether to include intercept in the model. Default is TRUE.
#' @param asparse The relative weight to put on the \eqn{\ell_1}-norm in
#'   sparse group lasso. Default is `0.05` (resulting in `0.95` on the
#'   \eqn{\ell_2}-norm).
#' @param standardize Logical flag for variable standardization (scaling) prior
#'   to fitting the model. Default is TRUE.
#' @param lower_bnd Lower bound for coefficient values, a vector in length of 1
#'   or of length the number of groups. Must be non-positive numbers only.
#'   Default value for each entry is `-Inf`.
#' @param upper_bnd Upper for coefficient values, a vector in length of 1
#'   or of length the number of groups. Must be non-negative numbers only.
#'   Default value for each entry is `Inf`.
#' @param offset Double vector. Optional offset (constant predictor without a
#'   corresponding coefficient). These can only be used with a
#'   [stats::family()] object.
#' @param trace_it Scalar integer. Larger values print more output during
#'   the irls loop. Typical values are `0` (no printing), `1` (some printing
#'   and a progress bar), and `2` (more detailed printing).
#'   These can only be used with a [stats::family()] object.
#'
#' @return An object with S3 class `"survival_sparsegl"`. Among the list components:
#' * `call` The call that produced this object.
#' * `b0` Intercept sequence of length `length(lambda)`.
#' * `beta` A `p` x `length(lambda)` sparse matrix of coefficients.
#' * `df` The number of features with nonzero coefficients for each value of
#'     `lambda`.
#' * `dim` Dimension of coefficient matrix.
#' * `lambda` The actual sequence of `lambda` values used.
#' * `npasses` Total number of iterations summed over all `lambda` values.
#' * `jerr` Error flag, for warnings and errors, 0 if no error.
#' * `group` A vector of consecutive integers describing the grouping of the
#'     coefficients.
#' * `nobs` The number of observations used to estimate the model.

#'
#'
#' @seealso [cv.survival_sparsegl()] and [`predict()`][predict.survival_sparsegl()], and [`coef()`][coef.survival_sparsegl()]
#'   methods for `"survival_sparsegl"` objects.
#'
#' @export
#'
#'
#' @examples
#' n <- 2000
#' p=10
#' x1 <- rnorm(n)
#' x2 <- 0.2*x1 + rnorm(n)
#' x3= matrix(rnorm(n*p), nrow = n, ncol = p)  # Define parameters for generating y
#' intercept <- 5
#' coef_x1 <- 2
#' coef_x2 <- 0
#' coef_x3 <- rep(0,p)
#' #Generate y based on logistic function
#' logistic_function <- intercept + coef_x1 * x1 + coef_x2 * x2 + x3 %*% coef_x3
#' probabilities <- plogis(logistic_function)
#' y <- rbinom(n, 1, probabilities)
#' X <- data.frame(x1, x2, x3)
#' X = data.matrix(X)
#' index = seq(p+2)
#' weight = c( rep(0.2,100), rep(1,100), rep(1,n - 200))
#' fit = survival_sparsegl(X, y, group = index, nlambda = 1, lambda = c(0), weight = weight, intercept_zero = 0, standardize = TRUE)
#' fit
survival_sparsegl <- function(
    x, y, weight, intercept_zero, group = NULL,
    nlambda = 100, lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04),
    lambda = NULL, pf_group = sqrt(bs), pf_sparse = rep(1, nvars),
    intercept = TRUE, asparse = 0.05, standardize = TRUE,
    lower_bnd = -Inf, upper_bnd = Inf,
    weights = NULL, offset = NULL, warm = NULL,
    trace_it = 0,
    dfmax = as.integer(max(group)) + 1L,
    pmax = min(dfmax * 1.2, as.integer(max(group))),
    eps = 1e-08, maxit = 100000) {
  # eps = 1e-08
  this.call <- match.call()
  if (!is.matrix(x) && !inherits(x, "sparseMatrix")) {
    cli::cli_abort("`x` must be a matrix.")
  }

  if (any(is.na(x))) cli::cli_abort("Missing values in `x` are not supported.")

  y <- drop(y)
  if (!is.null(dim(y))) cli::cli_abort("`y` must be a vector or 1-column matrix.")
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  vnames <- colnames(x)

  if (is.null(vnames)) vnames <- paste("V", seq(nvars), sep = "")

  if (length(y) != nobs) {
    cli::cli_abort("`x` has {nobs} rows while `y` has {length(y)}.")
  }

  #    group setup
  if (is.null(group)) {
    group <- 1:nvars
  } else {
    if (length(group) != nvars) {
      cli::cli_abort(c(
        "The length of `group` is {length(group)}.",
        "It must match the number of columns in `x`: {nvars}"
      ))
    }
  }

  bn <- as.integer(max(group))  # number of groups
  bs <- as.integer(as.numeric(table(group)))  # number of elements in each group

  if (!identical(as.integer(sort(unique(group))), as.integer(1:bn))) {
    cli::cli_abort("Groups must be consecutively numbered 1, 2, 3, ...")
  }

  if (asparse > 1) {
    cli::cli_abort(c(
      "`asparse` must be less than or equal to 1.",
      i = "You may want {.fn glmnet::glmnet} instead."
    ))
  }

  if (asparse < 0) {
    asparse <- 0
    cli::cli_warn("`asparse` must be in {.val [0, 1]}, running ordinary group lasso.")
  }
  if (any(pf_sparse < 0)) cli::cli_abort("`pf_sparse` must be non-negative.")
  if (any(is.infinite(pf_sparse))) {
    cli::cli_abort(
      "`pf_sparse` may not be infinite. Simply remove the column from `x`."
    )
  }
  if (any(pf_group < 0)) cli::cli_abort("`pf_group` must be non-negative.")
  if (any(is.infinite(pf_group))) {
    cli::cli_abort(c(
      "`pf_group` must be finite.",
      i = "Simply remove the group from `x`."
    ))
  }
  if (all(pf_sparse == 0)) {
    if (asparse > 0) {
      cli::cli_abort(
        "`pf_sparse` is identically 0 but `asparse` suggests some L1 penalty is desired."
      )
    } else {
      cli::cli_warn("`pf_sparse` was set to 1 because `asparse` = {.val {0}}.")
      pf_sparse = rep(1, nvars)
    }
  }

  ## Note: should add checks to see if any columns are completely unpenalized
  ## This is not currently expected.

  iy <- cumsum(bs) # last column of x in each group
  ix <- c(0, iy[-bn]) + 1 # first column of x in each group
  ix <- as.integer(ix)
  iy <- as.integer(iy)
  group <- as.integer(group)

  #parameter setup
  if (length(pf_group) != bn) {
    cli::cli_abort(
      "The length of `pf_group` must be the same as the number of groups: {.val {bn}}."
    )
  }
  if (length(pf_sparse) != nvars) {
    cli::cli_abort(
      "The length of `pf_sparse` must be equal to the number of predictors: {.val {nvars}}."
    )
  }

  pf_sparse <- pf_sparse / sum(pf_sparse) * nvars
  maxit <- as.integer(maxit)
  pf_group <- as.double(pf_group)
  pf_sparse <- as.double(pf_sparse)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)

  #lambda setup
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1) {
      cli::cli_abort("`lambda.factor` must be less than {.val {1}}.")
    }
    flmin <- as.double(lambda.factor)
    ulam <- double(1)
  } else {
    #flmin = 1 if user define lambda
    flmin <- as.double(1)
    if (any(lambda < 0)) cli::cli_abort("`lambda` must be non-negative.")
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }
  intr <- as.integer(intercept)

  ### check on upper/lower bounds
  if (any(lower_bnd > 0)) cli::cli_abort("`lower_bnd` must be non-positive.")
  if (any(upper_bnd < 0)) cli::cli_abort("`upper_bnd` must be non-negative.")
  lower_bnd[lower_bnd == -Inf] <- -9.9e30
  upper_bnd[upper_bnd == Inf] <- 9.9e30
  if (length(lower_bnd) < bn) {
    if (length(lower_bnd) == 1) lower_bnd <- rep(lower_bnd, bn)
    else cli::cli_abort("`lower_bnd` must be length {.val {1}} or length {.val {bn}}.")
  } else {
    lower_bnd <- lower_bnd[seq_len(bn)]
  }
  if (length(upper_bnd) < bn) {
    if (length(upper_bnd) == 1) upper_bnd <- rep(upper_bnd, bn)
    else cli::cli_abort("`upper_bnd` must be length {.val {1}} or length {.val {bn}}.")
  } else {
    upper_bnd <- upper_bnd[seq_len(bn)]
  }
  storage.mode(upper_bnd) <- "double"
  storage.mode(lower_bnd) <- "double"

  fit = sgl_logit_weight(
    bn, bs, ix, iy, nobs, nvars, x, y, pf_group, pf_sparse,
    dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr,
    as.double(asparse), standardize, lower_bnd, upper_bnd, weight, intercept_zero)

  # output
  if (is.null(lambda)) fit$lambda <- lamfix(fit$lambda)
  fit$call <- this.call
  fit$asparse <- asparse
  fit$nobs <- nobs
  fit$pf_group <- pf_group
  fit$pf_sparse <- pf_sparse
  class(fit) <- c(class(fit), "survival_sparsegl")
  fit
}

#' @importFrom stats glm binomial gaussian
sgl_logit_weight <- function(
    bn, bs, ix, iy, nobs, nvars, x, y, pf, pfl1,
    dfmax, pmax, nlam, flmin, ulam, eps,
    maxit, vnames, group, intr, asparse, standardize,
    lower_bnd, upper_bnd, weight, intercept_zero) {

  #y <- as.factor(y)
  lev <- levels(y)
  ntab <- table(y)
  minclass <- min(ntab)
  if (minclass <= 1)
    rlang::abort("Binomial regression: one class has 1 or 0 observations; not allowed")
  if (length(ntab) != 2)
    rlang::abort("Binomial regression: more than one class is not supported")
  if (minclass < 8)
    rlang::warn(c("Binomial regression: one class has fewer than 8",
                  "observations; dangerous ground"))
  # TODO, enable prediction with class labels if factor is passed
  if (intr == 1L && flmin < 1) b0_first <- coef(glm(y ~ 1, family = binomial()))
  #y <- 2 * (as.integer(y) - 1) - 1 # convert to -1 / 1

  #is.sparse <- FALSE
  #if (inherits(x, "sparseMatrix")) {
  #  is.sparse <- FALSE
  #  x <- as_dgCMatrix(x)
  #}
  if (standardize) {
    sm = Matrix::colMeans(x)
    sx <- apply(x, 2, sd)  #sqrt(Matrix::colSums(x^2))
    sx[sx < sqrt(.Machine$double.eps)] <- 1 # Don't divide by zero!]
    xs <- 1 / sx
    x <- sweep(x, 2, sm, FUN = "-")
    x <- x %*% Matrix::Diagonal(x = xs)
  }
  # if (is.sparse) {
  #   xidx <- as.integer(x@i + 1)
  #   xcptr <- as.integer(x@p + 1)
  #   xval <- as.double(x@x)
  #   nnz <- as.integer(utils::tail(x@p, 1))
  # }
  gamma <- 0.25 * calc_gamma(x, ix, iy, bn)
  #if (!0) {
  fit <- dotCall64::.C64(
    "log_sparse_four_weight",
    SIGNATURE = c("integer", "integer", "integer", "integer", "double",
                  "integer", "integer", "double", "double", "double",
                  "double", "integer", "integer", "integer", "double",
                  "double", "double", "integer", "integer", "integer",
                  "double", "double", "integer", "integer", "double",
                  "integer", "integer", "double", "double", "double", "double","double"),
    # Read only
    bn = bn, bs = bs, ix = ix, iy = iy, gam = gamma,
    nobs = nobs, nvars = nvars, x = as.double(x), y = as.double(y), pf = pf,
    pfl1 = pfl1,
    # Read / write
    dfmax = dfmax, pmax = pmax, nlam = nlam, flmin = flmin, ulam = ulam,
    eps = eps, maxit = maxit, intr = as.integer(intr),
    # Write only
    nalam = integer_dc(1), b0 = numeric_dc(nlam),
    beta = numeric_dc(nvars * nlam), activeGroup = integer_dc(pmax),
    nbeta = integer_dc(nlam), alam = numeric_dc(nlam), npass = integer_dc(1),
    jerr = integer_dc(1),
    # read only
    alsparse = asparse, lb = lower_bnd, ub = upper_bnd, weight = as.double(weight), intercept_zero = as.double(intercept_zero),
    INTENT = c(rep("r", 11), rep("rw", 8), rep("w", 8), rep("r", 5)))

  # output
  outlist <- getoutput(x, group, fit, maxit, pmax, nvars, vnames, eps)
  if (standardize) {
    outlist$beta <- outlist$beta * xs
    outlist$b0 <- matrix(outlist$b0 - colSums(sweep(outlist$beta, 1, sm, FUN = "*")), nrow = 1)
  }else{
    outlist$b0 <- matrix(outlist$b0, nrow = 1)
  }

  if (intr == 1L && flmin < 1) outlist$b0[1] <- b0_first
  outlist <- c(
    outlist,
    list(npasses = fit$npass,
         jerr = fit$jerr,
         group = group,
         classnames = lev)
  )
  class(outlist) <- c("logitspgl")
  outlist
}

#' @method calc_gamma
#' @export
calc_gamma <- function(x, ix, iy, bn) {
  gamma <- rep(NA, bn)
  for (g in seq_len(bn)) {
    grabcols <- ix[g]:iy[g]
    ncols <- length(grabcols)
    if (ncols > 2) gamma[g] <- RSpectra::svds(x[,grabcols], 1, 0, 0)$d^2
    else {
      if (ncols == 2) {gamma[g] <- maxeig2(x[,grabcols])}
      else {
        gamma[g] <- sum(x[,grabcols]^2)
      }
    }
  }
  return(as.double(gamma / nrow(x)))
}

#' @method maxeig2
#' @export
maxeig2 <- function(x) {
  # returns the largest squared singular value of n-by-2 matrix x
  # (the largest eigenvalue of 2-by-2 matrix mat)
  mat <- crossprod(x)
  tr <- mat[1] + mat[4]
  dt <- mat[1] * mat[4] - mat[2]^2
  return((tr + sqrt(tr^2 - 4 * dt)) / 2)
}

#######################################################
######################################################################
####################################################################
#' Cross-validation for a `survival_sparsegl` object.
#'
#' Performs k-fold cross-validation for [survival_sparsegl()].
#' This function is largely similar [glmnet::cv.glmnet()].
#'
#' The function runs [survival_sparsegl()] `nfolds + 1` times; the first to
#' get the `lambda` sequence, and then the remainder to compute the fit
#' with each of the folds omitted. The average error and standard error
#' over the folds are computed.
#'
#' @aliases cv.survival_sparsegl
#' @inheritParams survival_sparsegl
#' @param pred.loss Loss to use for cross-validation error. Valid options are:
#'  * `"default"` the same as deviance (mse for regression and deviance otherwise)
#'  * `"mse"` mean square error
#'  * `"deviance"` the default (mse for Gaussian regression, and negative
#'    log-likelihood otherwise)
#'  * `"mae"` mean absolute error, can apply to any family
#'  * `"misclass"` for classification only, misclassification error.
#' @param nfolds Number of folds - default is 10. Although `nfolds` can be
#'   as large as the sample size (leave-one-out CV), it is not recommended for
#'   large datasets. Smallest value allowable is `nfolds = 3`.
#' @param foldid An optional vector of values between 1 and `nfolds`
#'   identifying which fold each observation is in. If supplied, `nfolds` can
#'   be missing.
#' @param ... Additional arguments to [survival_sparsegl()].
#'
#' @return An object of class [cv.survival_sparsegl()] is returned, which is a
#'   list with the components describing the cross-validation error.
#'   \item{lambda}{The values of \code{lambda} used in the fits.}
#'   \item{cvm}{The mean cross-validated error - a vector of
#'     length \code{length(lambda)}.}
#'   \item{cvsd}{Estimate of standard error of \code{cvm}.}
#'   \item{cvupper}{Upper curve = \code{cvm + cvsd}.}
#'   \item{cvlower}{Lower curve = \code{cvm - cvsd}.}
#'   \item{name}{A text string indicating type of measure (for plotting
#'     purposes).}
#'   \item{nnzero}{The number of non-zero coefficients for each \code{lambda}}
#'   \item{active_grps}{The number of active groups for each \code{lambda}}
#'   \item{survival_sparsegl.fit}{A fitted [survival_sparsegl()] object for the full data.}
#'   \item{lambda.min}{The optimal value of \code{lambda} that gives
#'     minimum cross validation error \code{cvm}.}
#'   \item{lambda.1se}{The largest value of \code{lambda} such that error
#'     is within 1 standard error of the minimum.}
#'   \item{call}{The function call.}
#'
#'
#' @seealso [survival_sparsegl()], as well as [`plot()`][plot.cv.survival_sparsegl()],
#'   [`predict()`][predict.cv.survival_sparsegl()], and [`coef()`][coef.cv.survival_sparsegl()]
#'   methods for `"cv.survival_sparsegl"` objects.
#'
#' @export
#'
#'
#' @examples
#' n <- 2000
#' p=10
#' x1 <- rnorm(n)
#' x2 <- 0.2*x1 + rnorm(n)
#' x3= matrix(rnorm(n*p), nrow = n, ncol = p)  # Define parameters for generating y
#' intercept <- 5
#' coef_x1 <- 2
#' coef_x2 <- 0
#' coef_x3 <- rep(0,p)
#' #Generate y based on logistic function
#' logistic_function <- intercept + coef_x1 * x1 + coef_x2 * x2 + x3 %*% coef_x3
#' probabilities <- plogis(logistic_function)
#' y <- rbinom(n, 1, probabilities)
#' X <- data.frame(x1, x2, x3)
#' X = data.matrix(X)
#' index = seq(p+2)
#' weight = c( rep(0.2,100), rep(1,100), rep(1,n - 200))
#' fit.cv = cv.survival_sparsegl(X, y, group = index, nlambda = 1, lambda = c(0), weight = weight, intercept_zero = 0, standardize = TRUE)
#' fit.cv
#'
cv.survival_sparsegl <- function(
    x, y, group = NULL, family = c( "binomial"),
    lambda = NULL,
    pred.loss = c("default", "mse", "deviance", "mae", "misclass", 'censor', 'AUC_censor'),
    nfolds = 10, foldid = NULL, offset = NULL, weight, intercept_zero, AUC = FALSE, data = NULL, t = 0, maxit = 50000,
    ...) {

  pred.loss <- match.arg(pred.loss)

  N <- nrow(x)
  ###Fit the model once to get dimensions etc of output
  y <- drop(y)
  survival_sparsegl.object <- survival_sparsegl(x, y, group, lambda = lambda, offset = offset, weight = weight, intercept_zero = intercept_zero, maxit = maxit, ...)
  lambda <- survival_sparsegl.object$lambda
  # predict -> coef
  if (is.null(foldid)) foldid <- sample(rep(seq(nfolds), length = N))
  else nfolds <- max(foldid)
  if (nfolds < 2) {
    cli::cli_abort(
      "`nfolds` must be at least {.val {2}}. `nfolds` = {.val {10}} is recommended."
    )
  }
  outlist <- as.list(seq(nfolds))
  ###Now fit the nfold models and store them
  for (i in seq(nfolds)) {
    print(i)
    test_fold <- foldid == i
    outlist[[i]] <- survival_sparsegl(
      x = x[!test_fold, , drop = FALSE],
      y = y[!test_fold], weight = weight[!test_fold], group = group, lambda = lambda, offset = offset[!test_fold], intercept_zero = intercept_zero, maxit = maxit, ...)
  }
  ###What to do depends on the pred.loss and the model fit
  cvstuff <- cverror_survival(survival_sparsegl.object, outlist, lambda, x, y, foldid,
                     pred.loss, weight = weight, AUC = AUC, data = data, t = t)
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  AUC_censor = cvstuff$AUC_censor
  cvname <- cvstuff$name
  nz <- predict(survival_sparsegl.object, type = "nonzero")
  nnzero <- sapply(nz, length)
  active_grps <- sapply(nz, function(x) length(unique(group[x])))
  out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd,
              cvlo = cvm - cvsd, name = cvname,
              nnzero = nnzero, active_grps = active_grps, AUC_censor = AUC_censor,
              survival_sparsegl.fit = survival_sparsegl.object,
              call = match.call())
  lamin <- getmin(lambda, cvm, cvsd)
  obj <- c(out, as.list(lamin))
  class(obj) <- "cv.survival_sparsegl"
  obj
}



cverror_survival <- function(fullfit, outlist, lambda, x, y, foldid, pred.loss, weight, AUC, data, t, ...) {
  UseMethod("cverror_survival")
}

#' @export
cverror_survival.logitspgl <- function(
    fullfit, outlist, lambda, x, y, foldid, weight = weight, AUC = AUC, data = data, t = t,
    pred.loss = c("default", "mse", "deviance", "mae", "misclass", 'censor', 'AUC_censor'),
    ...) {
  typenames <- c(default = "Binomial deviance", mse = "Mean squared error",
                 deviance = "Binomial deviance", mae = "Mean absolute error",
                 misclass = "Missclassification error")
  pred.loss <- match.arg(pred.loss)
  prob_min <- 1e-05
  fmax <- log(1 / prob_min - 1)
  fmin <- -fmax
  nfolds <- max(foldid)
  predmat <- matrix(NA, length(y), length(lambda))
  nlams <- double(nfolds)
  for (i in seq_len(nfolds)) {
    test_fold <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[test_fold, , drop = FALSE], type = "response")
    nlami <- length(outlist[[i]]$lambda)
    predmat[test_fold, seq_len(nlami)] <- preds
    nlams[i] <- nlami
  }
  predmat <- pmin(pmax(predmat, fmin), fmax)
  binom_deviance <- function(m) stats::binomial()$dev.resids(y, m, 1)

  AUC_cen = matrix(0, nfolds, dim(predmat)[2])
  if (AUC == TRUE){
    library(survivalROC)
    for ( i in seq(dim(predmat)[2])){
      for (j in seq(nfolds)){
        test_fold <- foldid == j
        if (any(is.na(predmat[test_fold,i]))){
          AUC_cen[j, i] = 0
        } else{
          tr =  try({
            survivalROC::survivalROC(Stime=data$time[test_fold], status= data$status[test_fold], marker = predmat[test_fold,i], predict.time = t, span = 0.25*nrow(data[test_fold,])^(-1/2) )$AUC
          }, silent = TRUE)
          if (inherits(tr, "try-error")){
            AUC_cen[j, i] = 0}
          else{
            if (all(data$status[test_fold] == 1)){
            AUC_cen[j, i] = survivalROC::survivalROC(Stime=data$time[test_fold], status= data$status[test_fold], marker = predmat[test_fold,i], predict.time = t, span = 0)$AUC
            }else{
            AUC_cen[j, i] = survivalROC::survivalROC(Stime=data$time[test_fold], status= data$status[test_fold], marker = predmat[test_fold,i], predict.time = t, span = 0.25*nrow(data[test_fold,])^(-1/2) )$AUC
            }
            }
        }
      }
    }
  } else{
    AUC_cen = AUC_cen
  }
  AUC_censor = apply(AUC_cen, 2, mean, na.rm = TRUE)

  censor = matrix(0, length(y), dim(predmat)[2])
  for ( i in seq(dim(predmat)[2])){
    if (any(is.na(predmat[,i])) || any(is.infinite(predmat[,i])) ){
      censor[ ,i] = 1000
    } else{
      censor[ ,i] = -(weight*y*log(predmat[,i]) + (1-weight*y)*log(1-predmat[,i]))
    }
  }
  cvraw <- switch(
    pred.loss,
    mse = (y - predmat)^2,
    mae = abs(y - predmat),
    misclass = y != ifelse(predmat > 0.5, 1, 0),
    censor = censor,
    AUC_censor = AUC_censor,
    apply(predmat, 2, binom_deviance)
  )
  N <- length(y) - apply(is.na(predmat), 2, sum)
  cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  for ( j in seq(length(cvm))){
    if ( is.infinite(cvm[j]) ){
      cvm[j] = 1000
    }
  }
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, NA.RM = TRUE) /
                 (N - 1))
  list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss], AUC_censor = AUC_censor )
}

####################################################

############################################################################
############################################################################
############################################################################
#' Extract model coefficients from a `survival_sparsegl` object.
#'
#' Computes the coefficients at the requested value(s) for `lambda` from a
#' [survival_sparsegl()] object.
#'
#' `s` is the new vector of `lambda` values at which predictions are requested.
#' If `s` is not in the lambda sequence used for fitting the model, the `coef`
#' function will use linear interpolation to make predictions. The new values
#' are interpolated using a fraction of coefficients from both left and right
#' `lambda` indices.
#'
#' @param object Fitted [survival_sparsegl()] object.
#' @param s Value(s) of the penalty parameter `lambda` at which
#'  coefficients are required. Default is the entire sequence.
#' @param ... Not used.
#' @seealso [survival_sparsegl()] and [predict.survival_sparsegl()].
#'
#' @return The coefficients at the requested values for `lambda`.
#'
#' @method coef survival_sparsegl
#' @export
#' @examples
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n * p), nrow = n)
#' eps <- rnorm(n)
#' beta_star <- c(rep(5, 5), c(5, -5, 2, 0, 0), rep(-5, 5), rep(0, (p - 15)))
#' y <- X %*% beta_star + eps
#' groups <- rep(1:(p / 5), each = 5)
#' fit1 <- survival_sparsegl(X, y, group = groups)
#' coef.survival_sparsegl(fit1)
coef.survival_sparsegl <- function(object, s = NULL, ...) {
  rlang::check_dots_empty()
  b0 <- matrix(object$b0, nrow = 1)
  rownames(b0) <- "(Intercept)"
  nbeta <- rbind2(b0, object$beta)
  if (!is.null(s)) {
    vnames <- dimnames(nbeta)[[1]]
    dimnames(nbeta) <- list(NULL, NULL)
    lambda <- object$lambda
    lamlist <- lambda.interp(lambda, s)
    ls <- length(s)
    if (ls == 1) {
      nbeta = nbeta[, lamlist$left, drop = FALSE] * lamlist$frac +
        nbeta[, lamlist$right, drop = FALSE] * (1 - lamlist$frac)
    } else {
      nbeta = nbeta[, lamlist$left, drop = FALSE] %*%
        Matrix::Diagonal(ls, lamlist$frac) +
        nbeta[, lamlist$right, drop = FALSE] %*%
        Matrix::Diagonal(ls, 1 - lamlist$frac)
    }
    namess <- names(s) %||% paste0("s", seq_along(s))
    dimnames(nbeta) <- list(vnames, namess)
  }
  return(nbeta)
}




#' Make predictions from a `survival_sparsegl` object.
#'
#' Similar to other predict methods, this function produces fitted values and
#' class labels from a fitted [`survival_sparsegl`] object.
#'
#' `s` is the new vector of `lambda` values at which predictions are requested.
#' If `s` is not in the lambda sequence used for fitting the model, the `coef`
#' function will use linear interpolation to make predictions. The new values
#' are interpolated using a fraction of coefficients from both left and right
#' `lambda` indices.
#'
#'
#' @param object Fitted [survival_sparsegl()] model object.
#' @param newx Matrix of new values for `x` at which predictions are to be
#'   made. Must be a matrix. This argument is mandatory.
#' @param s Value(s) of the penalty parameter `lambda` at which
#'   predictions are required. Default is the entire sequence used to create the
#'   model.
#' @param type Type of prediction required. Type `"link"` gives the linear
#'   predictors for `"binomial"`; for `"gaussian"` models it gives the fitted
#'   values. Type `"response"` gives predictions on the scale of the response
#'   (for example, fitted probabilities for `"binomial"`); for `"gaussian"` type
#'   `"response"` is equivalent to type `"link"`. Type
#'   `"coefficients"` computes the coefficients at the requested values for
#'   `s`.
#'   Type `"class"` applies only to `"binomial"` models, and produces the
#'   class label corresponding to
#'   the maximum probability. Type `"nonzero"` returns a list of the indices
#'   of the nonzero coefficients for each value of \code{s}.
#'
#' @param ... Not used.
#' @return The object returned depends on type.
#' @seealso [survival_sparsegl()], [coef.survival_sparsegl()].
#'
#' @method predict survival_sparsegl
#' @export
predict.survival_sparsegl <- function(
    object, newx, s = NULL,
    type = c("link", "response", "coefficients", "nonzero", "class"),
    ...) {
  rlang::check_dots_empty()
  type <- match.arg(type)
  if (missing(newx)) {
    if (!match(type, c("coefficients", "nonzero"), FALSE))
      cli::cli_abort(
        "You must supply a value for `newx` when `type` == '{type}'."
      )
  }
  nbeta <- coef(object, s)
  if (type == "coefficients") return(nbeta)
  if (type == "nonzero") return(nonzeroCoef(nbeta[-1, , drop = FALSE]))
  #if (inherits(newx, "sparseMatrix")) newx <- as_dgCMatrix(newx)
  dx <- dim(newx)
  p <- object$dim[1]
  if (is.null(dx)) newx <- matrix(newx, 1, byrow = TRUE)
  if (ncol(newx) != p)
    cli::cli_abort("The number of variables in `newx` must be {p}.")
  fit <- as.matrix(cbind2(1, newx) %*% nbeta)
  fit
}


#' @export
predict.logitspgl <- function(
    object, newx, s = NULL,
    type = c("link", "response", "coefficients", "nonzero", "class"),
    ...) {
  type <- match.arg(type)
  nfit <- NextMethod("predict")
  switch(
    type,
    response = 1 / (1 + exp(-nfit)),
    class = object$classnames[ifelse(nfit > 0, 1, 0)],
    nfit
  )
}


#' @export
fitted.survival_sparsegl <- function(object, ...) {
  cli::cli_abort(c(
    "!" = "Because design matrices are typically large, these are not stored ",
    "!" = "in the estimated `survival_sparsegl` object. Use `predict()` instead, and ",
    "!" = "pass in the original data."))
}

#' @method summary survival_sparsegl
#' @export
summary.survival_sparsegl <- function(object, ...) {
  rlang::check_dots_empty()
  ns <- length(object$lambda)
  if (ns > 5) {
    xlam <- round(stats::quantile(1:ns))
    names(xlam) <- c("Max.", "3rd Qu.", "Median", "1st Qu.", "Min.")
  } else {
    xlam <- seq_len(ns)
    names(xlam) <- paste0("s", seq_len(ns))
  }
  nz <- predict(object, type = "nonzero")
  nnzero <- sapply(nz, length)
  active_grps <- sapply(nz, function(x) length(unique(object$group[x])))
  tab <- with(object, data.frame(
    lambda = lambda[xlam],
    index = xlam,
    nnzero = nnzero[xlam],
    active_grps = active_grps[xlam])
  )
  rownames(tab) <- names(xlam)
  out <- structure(list(call = object$call, table = tab),
                   class = "summary.survival_sparsegl")
  out
}

#' @method print summary.survival_sparsegl
#' @export
print.summary.survival_sparsegl <- function(
    x, digits = max(3, getOption("digits") - 3), ...) {

  rlang::check_dots_empty()
  lambda_warning <- all(x$table$nnzero == 0)

  cat("\nCall: ", deparse(x$call), "\n", fill = TRUE)


  if (lambda_warning) {
    cat("Warning: all regularization parameters resulted in empty models.\n\n")
  }

  cat("Summary of Lambda sequence:\n")
  print(x$tab, digits = digits)
  cat("\n")

}

#' @method print survival_sparsegl
#' @export
print.survival_sparsegl <- function(x, digits = min(3, getOption("digits") - 3), ...) {
  rlang::check_dots_empty()
  print(summary(x), digits = digits)
}

getmin <- function(lambda, cvm, cvsd) {
  plus = mapply(add_elements, cvm, cvsd)
  cvsd = na.omit(cvsd)
  cvm = na.omit(cvm)
  plus = na.omit(plus)
  cvmin <- min(cvm)
  idmin <- cvm <= cvmin
  lambda.min <- max(lambda[idmin])
  idmin <- match(lambda.min, lambda)
  semin <- (plus)[idmin]
  idmin <- cvm <= semin
  lambda.1se <- max(lambda[idmin])
  list(lambda.min = lambda.min, lambda.1se = lambda.1se)
}

add_elements <- function(x, y) {
  if (is.na(x) || is.na(y)) {
    # If either element is NA, return NA
    return(NA)
  } else {
    # Otherwise, return the sum of the elements
    return(x + y)
  }
}
####################################################################################
####################################################################################
####################################################################################
#' Extract coefficients from a `cv.survival_sparsegl` object.
#'
#' This function etracts coefficients from a
#' cross-validated [survival_sparsegl()] model, using the stored `"survival_sparsegl.fit"`
#' object, and the optimal value chosen for `lambda`.
#'
#' @param object Fitted [cv.survival_sparsegl()] object.
#' @param s Value(s) of the penalty parameter `lambda` at which
#'   coefficients are desired. Default is the single
#'   value `s = "lambda.1se"` stored in the CV object (corresponding to
#'   the largest value of `lambda` such that CV error estimate is within 1
#'   standard error of the minimum). Alternatively `s = "lambda.min"` can be
#'   used (corresponding to the minimum of cross validation error estimate).
#'   If `s` is numeric, it is taken as the value(s) of `lambda` to be used.
#' @param ... Not used.
#'
#'
#' @return The coefficients at the requested value(s) for `lambda`.
#' @seealso [cv.survival_sparsegl()] and [predict.cv.survival_sparsegl()].
#' @method coef cv.survival_sparsegl
#' @export
#' @examples
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n * p), nrow = n)
#' eps <- rnorm(n)
#' beta_star <- c(rep(5, 5), c(5, -5, 2, 0, 0), rep(-5, 5), rep(0, (p - 15)))
#' y <- X %*% beta_star + eps
#' groups <- rep(1:(p / 5), each = 5)
#' fit1 <- survival_sparsegl(X, y, group = groups)
#' cv_fit <- cv.survival_sparsegl(X, y, groups)
#' coef(cv_fit, s = c(0.02, 0.03))
coef.cv.survival_sparsegl <- function(object, s = c("lambda.1se", "lambda.min"), ...) {
  rlang::check_dots_empty()
  if (!(is.numeric(s) || is.character(s)))
    cli::cli_abort("Invalid form for `s`.")
  if (is.numeric(s)) lambda <- s
  else {
    s <- match.arg(s)
    lambda <- object[[s]]
  }
  coef(object$survival_sparsegl.fit, s = lambda)
}



#' Make predictions from a `cv.survival_sparsegl` object.
#'
#' This function makes predictions from a cross-validated [cv.survival_sparsegl()] object,
#' using the stored `survival_sparsegl.fit` object, and the value chosen for `lambda`.
#'
#' @inheritParams predict.survival_sparsegl
#' @param object Fitted [cv.survival_sparsegl()] object.
#' @param s Value(s) of the penalty parameter `lambda` at which
#'   coefficients are desired. Default is the single
#'   value `s = "lambda.1se"` stored in the CV object (corresponding to
#'   the largest value of `lambda` such that CV error estimate is within 1
#'   standard error of the minimum). Alternatively `s = "lambda.min"` can be
#'   used (corresponding to the minimum of cross validation error estimate).
#'   If `s` is numeric, it is taken as the value(s) of `lambda` to be used.
#'
#' @return A matrix or vector of predicted values.
#' @seealso [cv.survival_sparsegl()] and [coef.cv.survival_sparsegl()].
#'
#' @method predict cv.survival_sparsegl
#' @export
#' @examples
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n * p), nrow = n)
#' eps <- rnorm(n)
#' beta_star <- c(rep(5, 5), c(5, -5, 2, 0, 0), rep(-5, 5), rep(0, (p - 15)))
#' y <- X %*% beta_star + eps
#' groups <- rep(1:(p / 5), each = 5)
#' fit1 <- survival_sparsegl(X, y, group = groups)
#' cv_fit <- cv.survival_sparsegl(X, y, groups)
#' predict(cv_fit, newx = X[50:60, ], s = "lambda.min")
#'
predict.cv.survival_sparsegl <- function(
    object, newx,
    s = c("lambda.1se", "lambda.min"),
    type = c("link", "response", "coefficients", "nonzero", "class"), ...) {
  rlang::check_dots_empty()
  type <- match.arg(type)
  if (!(is.numeric(s) || is.character(s)))
    cli::cli_abort("Invalid form for `s`.")
  if (is.numeric(s)) lambda <- s
  else {
    s <- match.arg(s)
    lambda <- object[[s]]
  }
  predict(object$survival_sparsegl.fit, newx, s = lambda, type = type)
}

#' @method fitted cv.survival_sparsegl
#' @export
fitted.cv.survival_sparsegl <- function(object, ...) {
  cli::cli_abort(c(
    "!" = "Because design matrices are typically large, these are not stored ",
    "!" = "in the estimated `cv.survival_sparsegl` object. Use `predict()` instead, and ",
    "!" = "pass in the original data."))
}


#' @method summary cv.survival_sparsegl
#' @export
summary.cv.survival_sparsegl <- function(object, ...) {
  rlang::check_dots_empty()
  optlams <- c(object$lambda.1se, object$lambda.min)
  optlax <- c(1, match(optlams, object$lambda), length(object$lambda))
  tab <- with(object, data.frame(
    lambda = lambda[optlax],
    index = optlax,
    cvm = cvm[optlax],
    cvsd = cvsd[optlax],
    nnzero = nnzero[optlax],
    active_grps = active_grps[optlax])
  )
  rownames(tab) <- c("Max.", "lambda.1se", "lambda.min", "Min.")
  out <- structure(
    list(
      call = object$call,
      error_measure = object$name,
      table = tab
    ),
    class = "summary.cv.survival_sparsegl"
  )
  out

}

#' @method print summary.cv.survival_sparsegl
#' @export
print.summary.cv.survival_sparsegl <- function(
    x,
    digits = max(3, getOption("digits") - 3), ...
) {

  rlang::check_dots_empty()
  lambda_warning = NULL
  if (x$table$index[2] == 1) lambda_warning = "smallest"
  if (x$table$index[3] == x$table$index[4]) lambda_warning = "largest"
  cat("\nCall: ", deparse(x$call), "\n", fill = TRUE)

  cat("Error measure: ", x$error_measure, "\n\n")

  if (!is.null(lambda_warning)) {
    cat("Warning: the CV minimum occurred at the", lambda_warning,
        "lambda in the path.\n\n")
  }

  print(x$tab, digits = digits)
  cat("\n")
}

#' @method print cv.survival_sparsegl
#' @export
print.cv.survival_sparsegl <- function(x, digits = max(3, getOption("digits") - 3),
                                       ...) {

  rlang::check_dots_empty()
  print(summary(x), digits = digits)
}

#############################################################################################
#' @method ROC_censor_N
#' @export
ROC_censor_N = function(data, prediction, t){
  AUC_N = survivalROC::survivalROC(Stime=data$time, status= data$status, marker = prediction, predict.time = t, span = 0.25*length(data)^(-0.20) )
  return(AUC_N)
}

#' @method ROC_censor_KM
#' @export
ROC_censor_KM = function(data, prediction, t){
  AUC_KM = survivalROC::survivalROC(Stime=data$time, status= data$status, marker = prediction, predict.time = t, method = 'KM'  )
  return(AUC_KM)
}

as_dgCMatrix <- function(x) {
  as(as(x, "sparseMatrix"), "CsparseMatrix")
}

err <- function(n, maxit, pmax) {
  if (n == 0) msg <- ""
  if (n > 0) {
    #fatal error
    if (n < 7777)
      msg <- "Memory allocation error; contact package maintainer"
    if (n == 10000)
      msg <- "All penalty factors are <= 0"
    n <- 1
    msg <- paste("in survival_sparsegl fortran code -", msg)
  }
  if (n < 0) {
    #non fatal error
    if (n > -10000)
      msg <- paste(
        "Convergence for ", -n, "th lambda value not reached after maxit=",
        maxit, " iterations; solutions for larger lambdas returned",
        sep = "")
    if (n < -10000)
      msg <- paste(
        "Number of nonzero coefficients along the path exceeds pmax=",
        pmax, " at ", -n - 10000, "th lambda value; solutions for larger lambdas returned",
        sep = "")
    n <- -1
    msg <- paste("from gglasso fortran code -", msg)
  }
  list(n = n, msg = msg)
}

getoutput <- function(x, group, fit, maxit, pmax, nvars, vnames, eps) {
  if (fit$nalam != 0 ){
    nalam <- fit$nalam
    nbeta <- fit$nbeta[seq(nalam)]
    nbetamax <- max(nbeta)
    lam <- fit$alam[seq(nalam)]
    stepnames <- paste("s", seq(nalam) - 1, sep = "")
    errmsg <- err(fit$jerr, maxit, pmax)  ### error messages from fortran
    switch(paste(errmsg$n),
           `1` = stop(errmsg$msg, call. = FALSE),
           `-1` = print(errmsg$msg, call. = FALSE))
    dd <- c(nvars, nalam)
    if (nbetamax > 0) {
      beta <- Matrix::drop0(
        matrix(fit$beta[seq(nvars * nalam)], nvars, nalam),
        tol = eps^2)
      df <- apply(abs(beta) > 0, 2, sum) ## this is wrong, but fast
    } else {
      beta <- Matrix::Matrix(0, nvars, nalam)#, dimnames = list(vnames, stepnames))
      df <- rep(0, nalam)
    }
    b0 <- fit$b0
    if (!is.null(b0)) {
      b0 <- b0[seq(nalam)]
      #names(b0) <- stepnames
    }
    list(b0 = b0, beta = beta, df = df, dim = dd, lambda = lam)
  }else{
    nalam <- fit$nalam
    nbeta <- fit$nbeta[seq(nalam)]
    nbetamax <- max(nbeta)
    lam <- fit$alam[seq(nalam)]
    stepnames <- paste("s", seq(nalam) - 1, sep = "")
    errmsg <- err(fit$jerr, maxit, pmax)  ### error messages from fortran
    switch(paste(errmsg$n),
           `1` = stop(errmsg$msg, call. = FALSE),
           `-1` = print(errmsg$msg, call. = FALSE))
    dd <- c(nvars, nalam)
    if (nbetamax > 0) {
      beta <- Matrix::drop0(
        matrix(fit$beta[seq(nvars * nalam)], nvars, nalam,
               dimnames = list(vnames, stepnames)),
        tol = eps^2)
      df <- apply(abs(beta) > 0, 2, sum) ## this is wrong, but fast
    } else {
      beta <- Matrix::Matrix(0, nvars, nalam)#, dimnames = list(vnames, stepnames))
      df <- rep(0, nalam)
    }
    b0 <- fit$b0
    if (!is.null(b0)) {
      b0 <- b0[seq(nalam)]
      names(b0) <- stepnames
    }
    list(b0 = b0, beta = beta, df = df, dim = dd, lambda = lam)
  }
}

lambda.interp <- function(lambda, s) {
  ### lambda is the index sequence that is produced by the model
  ### s is the new vector at which evaluations are required.
  ### the value is a vector of left and right indicies, and a
  #   vector of fractions.
  ### the new values are interpolated between the two using the
  #   fraction
  ### Note: lambda decreases. you take:
  ### sfrac*left+(1-sfrac*right)
  if (length(lambda) == 1) {
    nums <- length(s)
    left <- rep(1, nums)
    right <- left
    sfrac <- rep(1, nums)
  } else {
    s[s > max(lambda)] <- max(lambda)
    s[s < min(lambda)] <- min(lambda)
    k <- length(lambda)
    sfrac <- (lambda[1] - s) / (lambda[1] - lambda[k])
    lambda <- (lambda[1] - lambda) / (lambda[1] - lambda[k])
    coord <- approx(lambda, seq(lambda), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac <- (sfrac - lambda[right]) / (lambda[left] - lambda[right])
    sfrac[left == right] <- 1
  }
  list(left = left, right = right, frac = sfrac)
}


lamfix <- function(lam) {
  llam <- log(lam)
  lam[1] <- exp(2 * llam[2] - llam[3])
  lam
}

nonzeroCoef <- function(beta) {
  nr <- nrow(beta)
  if (nr == 1) { # degenerate case
    apply(beta, 2, function(x) if (abs(x) > 0) 1 else NULL)
  } else {
    beta <- abs(beta) > 0 # this is sparse
    which <- seq(nr)
    ones <- rep(1, ncol(beta))
    nz <- as.vector((beta %*% ones) > 0)
    which <- which[nz]

    if (length(which) > 0) {
      beta <- as.matrix(beta[which, , drop = FALSE])
      nzel <- function(x, which) if (any(x)) which[x] else NULL
      which <- apply(beta, 2, nzel, which)
      if (!is.list(which)) which <- data.frame(which)
      which
    } else {
      dn <- dimnames(beta)[[2]]
      which <- vector("list", length(dn))
      names(which) <- dn
      which
    }
  }
}



###########################################################
##########################################################



