#' Multivariate normal model imputation
#'
#' This function performs multiple imputation under a multivariate normal model
#' as described by Schafer (1997), using his \code{norm} package, either with
#' or without posterior draws.
#'
#' This function imputes from a multivariate normal model with unstructured covariance
#' matrix, as described by Schafer (1997). With \code{pd=FALSE}, all imputed datasets
#' are generated conditional on the MLE of the model parameter, referred to as maximum
#' likelihood multiple imputation by von Hippel (2018).
#'
#' With \code{pd=TRUE}, regular 'proper' multiple imputation
#' is used, where each imputation is drawn from a distinct value of the model
#' parameter. Specifically, for each imputation, a single MCMC chain is run,
#' iterating for \code{steps} iterations.
#'
#' Imputed datasets can be analysed using \code{\link{withinBetween}},
#' \code{\link{scoreBased}}, or for example the
#' \href{https://cran.r-project.org/package=bootImpute}{bootImpute} package.
#'
#' @param obsData The data frame to be imputed.
#' @param M Number of imputations to generate.
#' @param pd Specify whether to use posterior draws (\code{TRUE})
#' or not (\code{FALSE}).
#' @param steps If \code{pd} is \code{TRUE}, the \code{steps} argument specifies
#' how many MCMC iterations to perform.
#' @param rseed The value to set the \code{norm} package's random number seed to,
#' using the \code{rngseed} function of \code{norm}. This function must be called at least
#' once before imputing using \code{norm}. If the user wishes to set the seed using
#' \code{rngseed} before calling \code{normImp}, set \code{rseed=NULL}.
#' @return A list of imputed datasets, or if \code{M=1}, just the imputed data frame.
#'
#' @references Schafer J.L. (1997). Analysis of incomplete multivariate data.
#' Chapman & Hall, Boca Raton, Florida, USA.
#'
#' @references von Hippel P.T. (2018) Maximum likelihood multiple imputation: faster,
#' more efficient imputation without posterior draws. \href{https://arxiv.org/abs/1210.0870v9}{arXiv:1210.0870v9}.
#'
#' @example data-raw/normExample.r
#'
#' @export
normImp <- function(obsData, M=10, pd=FALSE, steps=100, rseed) {

  if (is.data.frame(obsData)==FALSE) {
    stop("obsData argumment must be a data frame.")
  }

  if (is.null(rseed)==FALSE) {
    norm::rngseed(rseed)
  }

  imps <- vector("list", M)

  s <- norm::prelim.norm(as.matrix(obsData))
  thetahat <- norm::em.norm(s)

  if (pd==FALSE) {
    #MLMI
    for (i in 1:M) {
      imps[[i]] <- as.data.frame(norm::imp.norm(s,thetahat,as.matrix(obsData)))
    }
  } else {
    #PDMI
    for (i in 1:M) {
      theta <- norm::da.norm(s,thetahat,steps=steps)
      imps[[i]] <- as.data.frame(norm::imp.norm(s,theta,as.matrix(obsData)))
    }
  }

  if (M==1) {
    imps <- data.frame(imps[[1]])
  }
  attr(imps, "pd") <- pd

  imps
}
