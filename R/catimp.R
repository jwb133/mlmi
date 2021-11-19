#' Imputation for categorical variables using log linear models
#'
#' This function performs multiple imputation under a log-linear model
#' as described by Schafer (1997), using his \code{cat} package, either with
#' or without posterior draws.
#'
#' By default \code{catImp} will impute using a log-linear model allowing for all two-way
#' associations, but not higher order associations. This can be modified through
#' use of the \code{type} and \code{margins} arguments.
#'
#' With \code{pd=FALSE}, all imputed datasets are generated conditional on the MLE
#' of the model parameter, referred to as maximum likelihood multiple imputation
#' by von Hippel and Bartlett (2021).
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
#' @param obsData The data frame to be imputed. Variables must be coded such that
#' they take consecutive positive integer values, i.e. 1,2,3,...
#' @param M Number of imputations to generate.
#' @param pd Specify whether to use posterior draws (\code{TRUE})
#' or not (\code{FALSE}).
#' @param type An integer specifying what type of log-linear model to impute using.
#' \code{type=1}, the default, allows for all two-way associations in the log-linear model.
#' \code{type=2} allows for all three-way associations (plus lower).
#' \code{type=3} fits a saturated model.
#' @param margins An optional argument that can be used instead of \code{type} to specify
#' the desired log-linear model. See the documentation for the \code{margins} argument
#' in \code{\link[cat]{ecm.cat}} and Schafer (1997) on how to specify this.
#' @param steps If \code{pd} is \code{TRUE}, the \code{steps} argument specifies
#' how many MCMC iterations to perform in order to generate the model parameter value for
#' each imputation.
#' @param rseed The value to set the \code{cat} package's random number seed to,
#' using the \code{rngseed} function of \code{cat}. This function must be called at least
#' once before imputing using \code{cat}. If the user wishes to set the seed using
#' \code{rngseed} before calling \code{catImp}, set \code{rseed=NULL}.
#' @return A list of imputed datasets, or if \code{M=1}, just the imputed data frame.
#'
#' @references Schafer J.L. (1997). Analysis of incomplete multivariate data.
#' Chapman & Hall, Boca Raton, Florida, USA.
#'
#' @references von Hippel P.T. and Bartlett J.W. Maximum likelihood multiple imputation: faster,
#' more efficient imputation without posterior draws. Statistical Science 2021; 36(3) 400-420 \doi{10.1214/20-STS793}.
#'
#' @example data-raw/catExample.r
#'
#' @import stats
#'
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
    message("Imputing using log-linear model specified by margins argument")
  } else if (type==1) {
    message("Imputing using log-linear model with 2-way associations")
    #create corresponding margins argument
    margins <- NULL
    for (i in 1:(ncol(obsData)-1)) {
      for (j in (i+1):ncol(obsData)) {
        margins <- c(margins, i,j,0)
      }
    }
    margins <- utils::head(margins, -1)
  } else if (type==2) {
    if (ncol(obsData)<3) {
      stop("You cannot impute with 3-way associations with fewer than 3 variables.")
    } else {
      message("Imputing using log-linear model with 3-way (and lower order) associations")
      #create corresponding margins argument
      margins <- NULL
      for (i in 1:(ncol(obsData)-2)) {
        for (j in (i+1):(ncol(obsData)-1)) {
          for (k in (j+1):ncol(obsData)) {
            margins <- c(margins, i,j,k,0)
          }
        }
      }
      margins <- utils::head(margins, -1)
    }
  }

  #find MLE
  if (is.null(margins)==FALSE) {
    thetahat <- cat::ecm.cat(s, margins=margins)
  } else if (type==3) {
    message("Imputing using saturated log-linear model")
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
