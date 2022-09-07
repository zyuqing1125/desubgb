#' Internal function for debias process under the confidence-distribution-based method
#'
#' @description A function used to reduce the subgroup selection bias based on confidence-distribution method and provide a point estimate and a confidence interval.
#'
#' @details you can use this function to estimate any quantile of the parameters of interest compare to the bootstrap-based method.
#'
#' @param thetas A numeric vector which denotes the estimated treatment effect of all the subgroups.
#' @param s A numeric vector which is the estimated variance of the estimated treatment effect of each subgroup.
#' @param sigm sigm is the estimated covariance matrix of the estimated treatment effect of all the subgroups.
#' @param K The number of subgroups.
#' @param m The mth smallest parameter will be estimated.
#' @param B The bootstrap times.
#' @param cl The left bound within the estimation process.
#' @param cr The right bound within the estimation process.
#' @param M The true value of the mth smallest parameter, only for simulation use.
#'
#' @return
#' \item{esti}{The point estimate for the treatment effect of interest.}
#' \item{upper}{The 1 - alpha upper bound of the estimated confidence interval for the treatment effect of interest.}
#' \item{lower}{The 1 - alpha lower bound of the estimated confidence interval for the treatment effect of interest.}
#'
#' @keywords internal

newob2 <- function(the,s,nk){
  data=mvrnorm(n=1,the,s,tol = 2)
  return(data)
}

debias.beta_cfd <- function(thetas,s,sigm,K,m,B = 2000,cl,cr,M = NULL){


  them = sort(thetas)[m]
  epsiss = numeric(B)
  tilde.beta = numeric(B)
  w = which(thetas==them) ;d=1/4; tn = mean(s[w])^d

  for (j in 1:B) {
    epsi = newob2(thetas,sigm)
    e = epsi-sort(epsi)[m]
    wkm = numeric(K)
    temp = ((e<=cr*tn)+(e>=-cl*tn))==2
    wkm[temp] = 1
    epsiss[j] = sum(wkm*epsi)/sum(wkm)
    tilde.beta[j] = sum(wkm*thetas)/sum(wkm)
  }

  epsiss.sort = sort(epsiss)
  lower = epsiss.sort[0.025*B] # estimate the 95% confidence interval lower bound
  upper = epsiss.sort[0.975*B] # estimate the 95% confidence interval upper bound
  width = upper-lower
  esti = mean(epsiss.sort)
  esti2 = median(epsiss.sort) # estimate the mth smallest parameter
  mtbeta = median(tilde.beta)

  if (!is.null(M)){
    c2 = ifelse((M-lower>=0)&(upper-M>=0),1,0)
    return(list(esti = esti, upper = upper, lower = lower, cove = c2))
  }

  return(list(esti = esti, upper = upper, lower = lower))
}

