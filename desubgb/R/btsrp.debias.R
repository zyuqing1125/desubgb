################################################################################
#' Dealing with the subgroup selection bias with bootstrap-based method
#' @title btsrp.debias
#'
#' @description A function used to deal with the subgroup selection bias based on bootstrap method. You can use this function to estimate the maximum of the parameters.
#'
#' @details The order in data and the order in index, which should be corresponding to each other, should both be the sequence of the subject observed.
#'
#' @param data Numeric vector that represents the outcome measures of all the independent observations.
#' @param index A list indicates which subgroup or subgroups each subject belongs to. For example, if the nth item in index is the vector c(1,2), the nth subject will be classified into group 1 and group 2 at the same time.
#' @param K The number of subgroups.
#' @param func A self-defined function that estimates the treatment effect and the standard error using the outcome measures of one subgroup.
#' @param choice a vector denotes the range of tuning parameters being chosen from.
#' @param B B is the bootstrap times.
#' @param v v is the tuning parameter for v-fold cross-validated choice of the tuning parameter.
#' @param alpha alpha denotes the level 1 - alpha lower confidence limit.
#'
#' @return
#' \item{lcl}{The estimated 1 - alpha confidence limit for the maximum treatment effect size.}
#' \item{betamr}{The bias-reduced estimation for the maximum treatment effect size. The details could be found in https://doi.org/10.1080/01621459.2020.1740096}
#'
#' @examples
#' # Generate a random sample
#' set.seed(11); data <- rnorm(300); index <- rep(list(c(1,2),c(1,3),c(2,3)),100); K <- 3; func <- function(x){mean(x)}; choice <- c(1,2,3); B <- 200; v <- 5, alpha <- 0.05
#'
#' # Print the bias-reduced results
#' estim <- btsrp.debias(data, index, K, func, choice, B, v, alpha)
#' print(estim)
#' summary(estim)
#'
#' @export
btsrp.debias <- function(data, index, K, func, choice, B, v, alpha){

  para = tuningpara_bts(data, index, K, func, choice, B, v)

  result = debias.beta_bts(trda, trin, K, func, para, B, alpha)

  class(result) = "btsrpclass"
  return(result)
}

################################################################################
#' Print method for subgroup selection bias-reduced inference with the bootstrap-based method
#'
#' @description The print method for subgroup selection bias-reduced inference objects.
#'
#' @param object Class "btsrpclass" object, obtained by calling \code{\link{btsrp.debias}}.
#' @param ... Other arguments.
#'
#' @return
#' \item{}{The lower confidence limit and the estimated debiased beta_max results.}
#'
#' @examples
#' # Generate a random sample
#' set.seed(11); data <- rnorm(300); index <- rep(list(c(1,2),c(1,3),c(2,3)),100); K <- 3; func <- function(x){mean(x)}; choice <- c(1,2,3); B <- 200; v <- 5, alpha <- 0.05
#'
#' # Print the bias-reduced results
#' print(btsrp.debias(data, index, K, func, choice, B, v, alpha))
#'
#' @export
print.btsrpclass <- function(x, ... ){

  cat("the lower confidence limit:\n")
  print(x$lcl)
  cat("the estimated debiased beta_max:\n")
  print(x$betamr)
}


################################################################################
#' Summary method for subgroup selection bias-reduced inference with the bootstrap-based method
#'
#' @description The print method for subgroup selection bias-reduced inference objects.
#'
#' @param object Class "btsrpclass" object, obtained by calling \code{\link{btsrp.debias}}.
#' @param ... Other arguments.
#'
#' @return
#' \item{}{The lower confidence limit and the estimated debiased beta_max results.}
#'
#' @examples
#' # Generate a random sample
#' set.seed(11); data <- rnorm(300); index <- rep(list(c(1,2),c(1,3),c(2,3)),100); K <- 3; func <- function(x){mean(x)}; choice <- c(1,2,3); B <- 200; v <- 5, alpha <- 0.05
#'
#' # Summarize the bias-reduced results
#' summary(btsrp.debias(data, index, K, func, choice, B, v, alpha))
#'
#' @export
summary.btsrpclass <- function ( object , ...) {

  temp <- cbind (object$lcl , object$betamr)
  temp <- round (temp, 3)
  colnames (temp) <- c(" lower.confidence.limit ", "betamr")
  print (temp)
}

















