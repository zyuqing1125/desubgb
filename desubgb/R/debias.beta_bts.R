#' An internal function of the de-bias process under the bootstrap-based method
#'
#' @description An internal function used to estimate the maximum treatment effect of the subgroups and the level 1- alpha confidence limit.
#'
#' @details The order in data and the order in index, which should be corresponding to each other, should both be the sequence of the subject observed.
#'
#' @param data Numeric vector that represents the outcome measures of all the independent observations.
#' @param index A list indicates which subgroup or subgroups each subject belongs to. For example, if the nth item in index is the vector c(1,2), the nth subject will be classified into group 1 and group 2 at the same time.
#' @param K The number of the subgroups.
#' @param func A self-defined function that estimates the treatment effect and the standard error using the outcome measures of one subgroup.
#' @param B The bootstrap times.
#' @param v The tuning parameter for v-fold cross-validated choice of tuning parameter.
#' @param alpha Denotes the level 1 - alpha lower confidence limit.

#' @return
#' \item{lcl}{The estimated 1 - alpha confidence limit for the maximum treatment effect size.}
#' \item{betamr}{The bias-reduced estimation for the maximum treatment effect size. The details could be found in https://doi.org/10.1080/01621459.2020.1740096}
#'
#' @keywords internal




debias.beta_bts <- function(data, index, K, func, r, B, alpha){
  betamm <- t <- NULL

  n = length(data)
  beta = estibeta(data, index, K, func)$beta; betam = max(beta)
  d = (1-n^(r-0.5))*(betam - beta)

  b = 1
  while(b <= B) {
    j <- sample(1:n,size = n,replace = T)
    databo <- data[j]; indexbo <- index[j]

    betabo = estibeta(databo, indexbo, K, func)$beta
    t[b] <- sqrt(n)*(max(betabo + d) - betam); betamm[b] <- max(betabo + d)
    if(is.na(t[b])) next
    b = b + 1
  }

  betamr <- betam - mean(betamm - betam)

  calpha = quantile(t, 1-alpha)
  lcl <- betam - calpha/sqrt(n)

  return(list(lcl = lcl, betamr = betamr))

}
