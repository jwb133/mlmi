#' Categorical data model imputation using log linear models
#'
#' This function performs multiple imputation under a log-linear model
#' as described by Schafer (1997), using his \code{cat} package. By default
#' \code{catImp} will impute using a log-linear model allowing for all two-way
#' associations, but not higher order associations. This can be modified through
#' use of the \code{type} and \code{margins} arguments.
#'
#' @param obsData The data frame to be imputed. Variables must be coded taking positive
#' integer values only.
#' @param M Number of imputations to generate.
#' @param pd Specify whether to use posterior draws (\code{TRUE})
#' or not (\code{FALSE}).
#' @param type An integer specifying what type of log-linear model to impute using.
#' \code{type=1}, the default, allows for all two-way associations in the log-linear model.
#' \code{type=2} allows for all three-way associations (plus lower).
#' \code{type=3} fits a saturated model.
#' @param margins An optional argument that can be used instead of \code{type} to specify
#' the desired log-linear model.
#' @param steps If \code{pd} is \code{TRUE}, the \code{steps} argument specifies
#' how many MCMC iterations to perform.
#' @param rseed The value to set the \code{cat} package's random number seed to,
#' using the \code{rngseed} function of \code{cat}. This function must be called at least
#' once before imputing using \code{cat}. If the user wishes to set the seed using
#' \code{rngseed} before calling \code{catImp}, set \code{rseed=NULL}.
#' @return A list of imputed datasets, or if \code{M=1}, just the imputed data frame.
#' @export
catImp <- function(obsData, M=10, pd=FALSE, type=1, margins=NULL, steps=100, rseed) {

  if (is.data.frame(obsData)==FALSE) {
    stop("obsData argumment must be a data frame.")
  }

  if (is.null(rseed)==FALSE) {
    cat::rngseed(rseed)
  }

  imps <- vector("list", M)

  s <- cat::prelim.cat(as.matrix(obsData))

  if (is.null(margins)==FALSE) {
    print("Imputing using log-linear model specified by margins argument")
  } else if (type==1) {
    print("Imputing using log-linear model with 2-way associations")
    #create corresponding margins argument
    margins <- NULL
    for (i in 1:(ncol(obsData)-1)) {
      for (j in (i+1):ncol(obsData)) {
        margins <- c(margins, i,j,0)
      }
    }
    margins <- head(margins, -1)
  } else if (type==2) {
    if (ncol(obsData)<3) {
      stop("You cannot impute with 3-way associations with fewer than 3 variables.")
    } else {
      print("Imputing using log-linear model with 3-way (and lower order) associations")
      #create corresponding margins argument
      margins <- NULL
      for (i in 1:(ncol(obsData)-2)) {
        for (j in (i+1):(ncol(obsData)-1)) {
          for (k in (j+1):ncol(obsData)) {
            margins <- c(margins, i,j,k,0)
          }
        }
      }
      margins <- head(margins, -1)
    }
  }

  #find MLE
  if (is.null(margins)==FALSE) {
    thetahat <- cat::ecm.cat(s, margins=margins)
  } else if (type==3) {
    print("Imputing using saturated log-linear model")
    thetahat <- cat::em.cat(s)
  } else {
    stop("You must specify a valid type value or specify the margins argument.")
  }

  if (pd==FALSE) {
    #MLMI
    for (i in 1:M) {
      imps[[i]] <- as.data.frame(cat::imp.cat(s,thetahat))
    }
  } else {
    #PDMI
    for (i in 1:M) {
      if (is.null(margins)==FALSE) {
        theta <- cat::dabipf(s,margins=margins,start=thetahat,steps=steps)
      } else {
        theta <- cat::da.cat(s,thetahat,steps=steps)
      }
      imps[[i]] <- as.data.frame(cat::imp.cat(s,theta))
    }
  }

  if (M==1) {
    imps <- data.frame(imps[[1]])
  }
  attr(imps, "pd") <- pd

  imps
}
