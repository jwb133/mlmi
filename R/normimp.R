#' Multivariate normal model imputation
#'
#' @param obsData The data frame to be imputed.
#' @param M Number of imputations to generate.
#' @param pd Specify whether to use posterior draws (\code{TRUE})
#' or not (\code{FALSE}).
#' @param steps If \code{pd} is \code{TRUE}, the \code{steps} argument specifies
#' how many MCMC iterations to perform.
#' @return A list of imputed datasets, or if \code{M=1}, just the imputed data frame.
#' @export
normImp <- function(obsData, M=10, pd=FALSE, steps=100) {

  if (is.data.frame(obsData)==FALSE) {
    stop("obsData argumment must be a data frame.")
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
