################################################################################
#' Dealing with the subgroup selection bias with confidence-distribution-based method
#' @title cfdd.debias
#'
#' @description The function is used to reduce the subgroup selection bias based on confidence-distribution method according to https://doi.org/10.48550/arXiv.2206.11868.
#'
#' @details you can use this function to estimate any quantile of the parameters of interest compare to the bootstrap-based method.
#'
#' @param thetas A numeric vector of the estimated treatment effect of the subgroups.
#' @param s A numeric vector of the estimated variance of the estimated treatment effect of the subgroups.
#' @param nk A vector which denotes the sample size of each subgroup.
#' @param sigm sigm is the estimated covariance matrix of the estimated treatment effect of the subgroups.
#' @param K The number of subgroups.
#' @param m The mth smallest parameter will be estimated.
#' @param B The bootstrap times.
#' @param cl The range of left bounds being chosen from within the estimation process.
#' @param cr The range of right bounds being chosen from within the estimation process.
#' @param M The true value of the mth smallest parameter, only for simulation use.
#' @param R Tuning parameter.
#'
#' @return
#' \item{beta.hat}{The estimated debiased treatment effect size of interest.}
#' \item{upper}{The 1 - alpha upper bound of the estimated confidence interval for the treatment effect of interest.}
#' \item{lower}{The 1 - alpha lower bound of the estimated confidence interval for the treatment effect of interest.}
#'
#' @examples
#' # Generate a random sample
#' set.seed(11); thetas <- rnorm(3); s <- runif(3); nk <- c(200,200,200); sigm <- diag(s); K <- 3; m <- 3; cl <- c(1,2,3); cr <- c(1,2,3)
#'
#' # Print the bias-reduced results
#' estim <- cfdd.debias(thetas,s,nk,sigm,K,m,B = 2000,cl,cr,M = NULL,R = 200)
#' print(estim)
#' summary(estim)
#'
#' @export
cfdd.debias <- function(thetas,s,nk,sigm,K,m,B = 2000,cl,cr,M = NULL,R = 200){
  ga = gam(B)

  temp = tuningpara_cfd(cl,cr,thetas,s,nk,sigm,K,m,ga,R,B)
  cl = temp$cl; cr = temp$cr

  temp = debias.beta_cfd(thetas,s,sigm,K,m,B,cl,cr,M)
  result = list(beta.hat = esti, upper = upper, lower = lower)
  class(result) <- "cfddclass"

  return(result)
}

################################################################################
#' Print method for subgroup selection bias-reduced inference based on confidence-distribution-based method
#'
#' @description The print method for subgroup selection bias-reduced inference objects with confidence-distribution-based method.
#'
#' @param object Class "cfddclass" object, obtained by calling \code{\link{cfdd.debias}}.
#' @param ... Other arguments.
#'
#' @return
#' \item{}{The estimated debiased beta_max results and the estimated confidence interval.}
#'
#' @examples
#' # Generate a random sample
#' set.seed(11); thetas <- rnorm(3); s <- runif(3); nk <- c(200,200,200); sigm <- diag(s); K <- 3; m <- 3; cl <- c(1,2,3); cr <- c(1,2,3)
#'
#' # Print the bias-reduced results
#' print(cfdd.debias(thetas,s,nk,sigm,K,m,B = 2000,cl,cr,M = NULL,R = 200))
#'
#' @export
print.cfddclass <- function(x, ... ){

  cat("the estimated beta_hat:\n")
  print(x$beta.hat)
  cat("the estimated confidence interval:\n")
  print(c(x$lower, x$upper))

}

################################################################################
#' Summarize method for subgroup selection bias-reduced inference based on confidence-distribution-based method
#'
#' @description The print method for subgroup selection bias-reduced inference objects with confidence-distribution-based method.
#'
#' @param object Class "cfddclass" object, obtained by calling \code{\link{cfdd.debias}}.
#' @param ... Other arguments.
#'
#' @return
#' \item{}{The estimated debiased beta_max results and the estimated confidence interval.}
#'
#' @examples
#' # Generate a random sample
#' set.seed(11); thetas <- rnorm(3); s <- runif(3); nk <- c(200,200,200); sigm <- diag(s); K <- 3; m <- 3; cl <- c(1,2,3); cr <- c(1,2,3)
#'
#' # Print the bias-reduced results
#' summary(cfdd.debias(thetas,s,nk,sigm,K,m,B = 2000,cl,cr,M = NULL,R = 200))
#'
#' @export
summary.cfddclass <- function ( object , ...) {

  temp <- cbind (object$beta.hat , object$lower, object$upper)
  temp <- round (temp, 3)
  colnames (temp) <- c(" beta.hat ", "lower", "upper")
  print (temp)
}




