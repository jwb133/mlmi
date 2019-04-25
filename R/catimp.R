#' Categorical data model imputation under saturated multinomial model
#'
#' This function performs multiple imputation under a saturated multinomial model
#' for categorical data, as described by Schafer (1997). If the dataset has
#' more than a small number of variables, a saturated multinomial model may not
#' be estimable, in which case the user should consider a simpler log-linear
#' model.
#'
#' The \code{rngseed} function must be called at least once before \code{catSatImp}
#' is called.
#'
#' @param obsData The data frame to be imputed. Variables must be coded taking positive
#' integer values only.
#' @param M Number of imputations to generate.
#' @param pd Specify whether to use posterior draws (\code{TRUE})
#' or not (\code{FALSE}).
#' @param steps If \code{pd} is \code{TRUE}, the \code{steps} argument specifies
#' how many MCMC iterations to perform.
#' @return A list of imputed datasets, or if \code{M=1}, just the imputed data frame.
#' @export
catSatImp <- function(obsData, M=10, pd=FALSE, steps=100) {

  if (is.data.frame(obsData)==FALSE) {
    stop("obsData argumment must be a data frame.")
  }

  imps <- vector("list", M)

  s <- cat::prelim.cat(as.matrix(obsData))
  thetahat <- cat::em.cat(s)

  if (pd==FALSE) {
    #MLMI
    for (i in 1:M) {
      imps[[i]] <- as.data.frame(cat::imp.cat(s,thetahat))
    }
  } else {
    #PDMI
    for (i in 1:M) {
      theta <- cat::da.cat(s,thetahat,steps=steps)
      imps[[i]] <- as.data.frame(cat::imp.cat(s,theta))
    }
  }

  if (M==1) {
    imps <- data.frame(imps[[1]])
  }
  attr(imps, "pd") <- pd

  imps
}
