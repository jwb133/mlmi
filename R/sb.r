#' Score based variance estimation
#'
#' This function implements the score based variance estimation approach described by von Hippel
#' and Bartlett, based on Wang and Robins (1998).
#'
#' @param imps A list of imputed datasets produced by one of the imputation functions
#' in \code{mlmi} or another package.
#' @param analysisFun A function to analyse the imputed datasets that when applied to
#' a dataset returns a list containing a vector \code{est}.
#' @param scoreFun A function whose first argument is a dataset and whose second argument is
#' a vector of parameter values. It should return a matrix of subject level scores
#' evaluated at the parameter value passed to it.
#' @param dfComplete The complete data degrees of freedom. If \code{analysisFun} returns a vector
#' of parameter estimates, \code{dfComplete} should be a vector of the same length. If not
#' specified, it is assumed that the complete data degrees of freedom is effectively infinite (1e+05).
#' @param ... Other parameters that are to be passed through to \code{analysisFun}.
#' @return A list containing the overall parameter estimates, its corresponding covariance matrix, and
#' degrees of freedom for each parameter.
#'
#' @references von Hippel PT, Bartlett JW. Maximum likelihood multiple imputation: can multiple imputation
#' work without posterior draws? under review
#'
#' @example data-raw/wbExample.r
#'
#'
#' @export
scoreBased <- function(imps, analysisFun, scoreFun, pd=NULL, dfComplete=NULL, ...) {
  M <- length(imps)
  if ("pd" %in% names(attributes(imps)) == TRUE) {
    pd <- as.logical(attributes(imps)['pd'])
  } else {
    if (is.null(pd)==TRUE) {
      stop("Your imputed datasets doesn't have a pd attribute stored, so you must use specify the pd argument when calling wb.")
    }
  }

  #run analysis on first imputation to find out length of parameter vector
  result <- analysisFun(imps[[1]],...)
  numParms <- length(result$est)
  N <- dim(imps[[1]])[1]

  #analyse each imputed datasets
  ests <- array(0, dim=c(M,numParms))
  for (m in 1:M) {
    #result <- analysisFun(imps[[m]],...)
    result <- analysisFun(imps[[m]])
    ests[m,] <- result$est
  }
  thetaHat <- colMeans(ests)
  Bhat <- var(ests)
  scores <- array(0, dim=c(M,N,numParms))
  Vcom_inv <- array(0, dim=c(numParms,numParms))
  for (m in 1:M) {
    scores[m,,] <- scoreFun(imps[[m]], thetaHat)
    Vcom_inv <- Vcom_inv + t(scores[m,,]) %*% scores[m,,]
  }
  Vcom_inv <- Vcom_inv / M
  scom_bar <- apply(scores, c(2,3), mean)
  Vmis_inv <- array(0, dim=c(numParms,numParms))
  for (m in 1:M) {
    Vmis_inv <- Vmis_inv + t(scores[m,,]-scom_bar) %*% (scores[m,,]-scom_bar)
  }
  Vmis_inv <- Vmis_inv/(M-1)
  VML_inv <- Vcom_inv - Vmis_inv
  Vcom <- solve(Vcom_inv)

  gammaHatMis <- Vmis_inv %*% Vcom
  gammaHatObs <- diag(numParms) - gammaHatMis
  gammaTildeMis <- H(gammaHatMis, (M-1)*N)
  gammaTildeObs <- diag(numParms) - gammaTildeMis
  print(gammaTildeMis)
  VTildeML <- Vcom %*% solve(gammaTildeObs)
  totalVar <- VTildeML + Bhat/M

  if (is.null(dfComplete)) {
    #assume effectively infinite degrees of freedom for each parameter
    dfComplete <- rep(100000,numParms)
  }

  nuObs <- dfComplete*diag(gammaTildeObs)*(dfComplete+3)/(dfComplete+1)
  MIdf <-  totalVar^2 / (totalVar^2/nuObs + (diag(Bhat)/M)^2/(M-1))

  list(est=thetaHat, var=totalVar, df=MIdf)
}

