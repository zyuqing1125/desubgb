################################################################################
#' @title desubgb: A package for reducing the subgroup selection bias.
#'
#' @description The desubgb package provides two methods for adjusting the subgroup selection bias: the bootstrap-based method and the confidence-distribution-based method.
#'
#' The bootstrap-based method is established from Guo and He (2020) (https://doi.org/10.1080/01621459.2020.1740096)
#'  and the confidence-distribution-based method comes from Inference on the Best Policies with Many Covariates (https://doi.org/10.48550/arXiv.2206.11868).
#'  You can find the coding details in the two articles.
#'
#' @author
#' Yuqing Zhou, UMich. \email{zyuqing@umich.edu}
#'
#' Waverly Wei, UC Berkeley.
#'
#' Jingshen Wang, UC Berkeley.
#'
#' @maintainer Yuqing Zhou
#'
#' @references
#' Guo and He (2020). https://doi.org/10.1080/01621459.2020.1740096
#'
#' Inference on the Best Policies with Many Covariates. https://doi.org/10.48550/arXiv.2206.11868
#'
#' @importFrom MASS mvrnorm
#' @aliases desubgb
"_PACKAGE"
