#' An internal function for debias process under the bootstrap-based method
#'
#' @description The function is used to select the best tuning parameter for dealing with the subgroup selection bias based on bootstrap method.
#'
#' @details The order in data and the order in index, which should be corresponding to each other, should both be the sequence of the subject observed.
#'
#' @param data a vector represents the outcome measure of all the independent observations.
#' @param index a list indicates which subgroup or subgroups subject belongs to.
#' @param K K is the number of subgroups.
#' @param func func is a self-defined function that estimates the treatment effect and the standard error with the outcome measures of a subgroup.
#' @param choice A vector denotes the range of the tuning parameter being chosen from.
#' @param B The bootstrap times.
#' @param v The tuning parameter for v-fold cross-validated choice of tuning parameter.

#' @return The best selected parameter
#' @export
#' @keywords internal


tuningpara_bts <- function(data, index, K, func, choice, B, v){
  n = length(data); m = length(choice); bmr <- matrix(nrow = v,ncol = m); h <- matrix(nrow = v,ncol = K); hl <- matrix(nrow = m, ncol = K)

  foldid <- sample(rep(1:v,length = n))

  for (l in 1:m) {
    for (j in 1:v) {
      trda <- data[foldid != j]; reda <- data[foldid == j]; trin <- index[foldid != j]; rein <- index[foldid == j]

      bmr[j,l] <- debias.beta_bts(trda, trin, K, func, choice[l], B, 0.5)$betamr

      temp = estibeta(reda, rein, K, func)

      true_beta = temp$beta; sigma = temp$sigma

      h[j,] <- (bmr[j,l]-true_beta)^2-sigma^2
    }

    hl[l,] = colMeans(h)

  }

  mi = apply(hl, 1, min)
  ll <- which(mi==min(mi))

  return(choice[ll])
}
