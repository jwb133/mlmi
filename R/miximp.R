#' Imputation for a mixture of continuous and categorical variables using
#' the general location model.
#'
#' This function performs multiple imputation under a general location model
#' as described by Schafer (1997), using his \code{mix} package.
#'
#' @param obsData The data frame to be imputed. The categorical variables must be
#' in the first \code{p} columns, and they must be coded using consecutive positive
#' integers.
#' @param p The number of categorical variables in \code{obsData}.
#' @param M Number of imputations to generate.
#' @param pd Specify whether to use posterior draws (\code{TRUE})
#' or not (\code{FALSE}).
#' @param margins An optional margins argument for specifying a restricted log-linear
#' model for the categorical variables. If \code{margins} is specified, \code{design}
#' must also be specified.
#' @param design An optional argument for specifying the design matrix controlling the
#' dependence of the continuous variables on the categorical variables.  If
#' \code{design} is specified, \code{margins} must also be specified.
#' @param steps If \code{pd} is \code{TRUE}, the \code{steps} argument specifies
#' how many MCMC iterations to perform.
#' @param rseed The value to set the \code{mix} package's random number seed to,
#' using the \code{rngseed} function of \code{mix}. This function must be called at least
#' once before imputing using \code{mix}. If the user wishes to set the seed using
#' \code{rngseed} before calling \code{mixImp}, set \code{rseed=NULL}.
#' @return A list of imputed datasets, or if \code{M=1}, just the imputed data frame.
#' @export
mixImp <- function(obsData, nCat, M=10, pd=FALSE, type=1, margins=NULL, design=NULL,
                   steps=100, rseed) {

  if (is.data.frame(obsData)==FALSE) {
    stop("obsData argumment must be a data frame.")
  }

  if (is.null(rseed)==FALSE) {
    mix::rngseed(rseed)
  }

  if (((is.null(margins)==FALSE) & (is.null(design)==TRUE)) |
    ((is.null(margins)==TRUE) & (is.null(design)==FALSE))) {
    stop("Margins and design arguments must both be specified if one is specified.")
  }

  imps <- vector("list", M)

  s <- mix::prelim.mix(as.matrix(obsData), p=nCat)

  #find MLE
  if (is.null(margins)==FALSE) {
    thetahat <- mix::ecm.mix(s, margins=margins, design=design)
  } else {
    print("Imputing using unrestricted general location model")
    thetahat <- mix::em.mix(s)
  }

  if (pd==FALSE) {
    #MLMI
    for (i in 1:M) {
      imps[[i]] <- as.data.frame(mix::imp.mix(s,thetahat,as.matrix(obsData)))
    }
  } else {
    #PDMI
    for (i in 1:M) {
      if (is.null(margins)==FALSE) {
        theta <- mix::dabipf.mix(s,margins=margins,design=design,start=thetahat,steps=steps)
      } else {
        theta <- mix::da.mix(s,thetahat,steps=steps)
      }
      imps[[i]] <- as.data.frame(mix::imp.mix(s,theta))
    }
  }

  if (M==1) {
    imps <- data.frame(imps[[1]])
  }
  attr(imps, "pd") <- pd

  imps
}
