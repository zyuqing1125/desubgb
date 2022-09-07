#' internal function for debias process with bootstrap-based method
#'
#' This internal function is used to estimate the subgroup treatment effects with the outcome measures and indicators of which subgroup each subject belongs to.
#'
#' @details The order in data and the order in index, which should be corresponding to each other, should both be the sequence of the subject observed.
#'
#' @param data Numeric vector that represents the outcome measures of all the independent observations.
#' @param index A list indicates which subgroup or subgroups each subject belongs to. For example, if the nth item in index is the vector c(1,2), the nth subject will be classified into group 1 and group 2 at the same time.
#' @param K The number of the subgroups.
#' @param func A self-defined function that estimates the treatment effect and the standard error using the outcome measures of one subgroup.
#' @return
#' \item{beta}{A vector containing the estimated treatment effect of all subgroups.}
#' \item{sigma}{A vector containing the estimated standard error of the treatment effect of all subgroups.}
#'
#' @keywords internal



estibeta <- function(data, index, K, func){
  sigma <- beta <- NULL
  df = data.frame(
    grp = unlist(index),
    num = factor(rep(seq_along(index), lengths(index)))
  )

  for (i in 1:K) {
    where = df[df$grp == i, "num"]
    temp = func(data[where])
    beta[i] = temp$beta; sigma[i] = temp$sigma
  }

  return(list(beta = beta, sigma = sigma))
}
