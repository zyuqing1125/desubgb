#' Internal function for debias process using the confidence-distribution-based method
#' @title tuningpara_cfd
#'
#' @description A function used to select the best tuning paramters for reducing the subgroup selection bias based on confidence-distribution method.
#'
#'
#' @param thetah The numeric vector which denotes the estimated treatment effect of each subgroup.
#' @param s The numeric vector which represents the estimated variance of the parameter of each subgroup.
#' @param nk nk is a vector which denotes the sample size of each subgroup.
#' @param sigm sigm is the estimated covariance matrix of the estimated treatment effect of all the subgroups.
#' @param K The number of subgroups.
#' @param m The mth smallest parameter will be estimated.
#' @param B The bootstrap times.
#' @param cl The range of left bound within the estimation process being chosen from.
#' @param cr The range of right bound within the estimation process being chosen from.
#' @param R Tuning parameter.
#'
#' @return
#' \item{cl}{The selected left bound.}
#' \item{cr}{The selected right bound.}
#'
#' @keywords internal



tuning.dist <- function(thetahh,s,nk,sigm,K,m,cl,cr,R,B){

  # thetahh:parameters of interest
  # s:the estimated variances of parameters
  # nk:the size of each subgroup
  # sigm:the estimated covariance matrix
  # K:the number of subgroups
  # m:the mth smallest parameters of interest (the biggest:m = K)

  r = mean(s)/var(thetahh); sig = sqrt(sum(s*nk)/K); d=1/4
  temp = r*(sig/sqrt(mean(s)))^0.1; tri = min(c(1,temp))
  thetah = tri*mean(thetahh)+(1-tri)*thetahh
  them = sort(thetah)[m]
  b = numeric(B)

  for (i in 1:B) {
    thetash = newob2(thetah,sigm)
    temp = sort(thetash)[m]; w = which(thetash==temp)

    for (j in 1:R) {
      epsi = newob2(thetash,sigm); tn = (s[w])^d
      e = epsi-sort(epsi)[m]
      wkm = numeric(K)
      temp = ((e<=cr*tn)+(e>=-cl*tn))==2
      wkm[temp] = 1
      epsiss = sum(wkm*epsi)/sum(wkm)
      b[i] = b[i]+ifelse(them-epsiss>=0,1,0)
    }

    b[i] = b[i]/R
  }

  l = sum((sort(b)-seq(1,B)/(B+1))^2)/B
  return(l)
}

gam = function(B){
  l = numeric(20000)
  for (i in 1:20000) {
    l[i] = sum((sort(runif(B))-(1:B)/(B+1))^2)/B
  }
  return(quantile(l,0.975))
}

tuningpara_cfd <- function(cl,cr,thetah,s,nk,sigm,K,m,gam,R,B){

  # cl:choices of c_L, cr:choices of c_R
  # return the best selected (cl,cr)

  l1 = length(cr); l2 = length(cl); d=1/4
  l = matrix(0,nrow = l2,ncol = l1)

  for (n1 in 1:l2) {
    c1=cl[n1]

    for (n2 in 1:l1) {
      c2=cr[n2]
      l[n1,n2]=tuning.dist(thetah,s,nk,sigm,K,m,c1,c2,R,B)
    }

  }
  #print(l)

  if (!any(l<=gam)) {
    mi = min(l); temp = which(l==mi,arr.ind = T)[1,]
    c1 = temp[1]; c2 = temp[2]

    return(list(cl = cl[c1],cr = cr[c2]))
  }

  temp = which(l<=gam,arr.ind = T)
  c1 = mean(cl[temp[,1]]); c2 = mean(cr[temp[,2]])

  return(list(cl = c1,cr = c2))
}
