#' Normal regression imputation of a single variable
#'
#' Performs multiple imputation of a single continuous variable using a normal
#' linear regression model. The covariates in the imputation model must be fully
#' observed. By default \code{normUniImp} imputes every dataset using the
#' maximum likelihood estimates of the imputation model parameters, which here
#' coincides with the OLS estimates. If \code{pd=TRUE} is specified, it instead
#' performs posterior draw Bayesian imputation.
#'
#' @param obsData The data frame to be imputed.
#' @param impFormula The linear model formula.
#' @param M Number of imputations to generate.
#' @param pd Specify whether to use posterior draws (\code{TRUE})
#' or not (\code{FALSE}).
#' @return A list of imputed datasets, or if \code{M=1}, just the imputed data frame.
#' @export
normUniImp <- function(obsData, impFormula, M=5, pd=FALSE) {
  #put some checks in to check there are no missing values in covariates of imp model

  #fit desired imputation model
  impMod <- lm(impFormula, data=obsData)

  outcomeVar <- all.vars(impFormula)[1]

  if (M>1) {
    imps <- vector("list", M)
  }

  if (pd==FALSE) {
    fittedVals <- predict(impMod, obsData)
    for (i in 1:M) {
      newimp <- obsData
      #impute
      newimp[is.na(obsData[,outcomeVar]),outcomeVar] <- fittedVals[is.na(obsData[,outcomeVar])] +
        rnorm(sum(is.na(obsData[,outcomeVar])), mean=0, sd=sigma(impMod))
      if (M>1) {
        imps[[i]] <- newimp
      }
    }
  } else {
    beta <- impMod$coef
    sigmasq <- summary(impMod)$sigma^2
    varcov <- vcov(impMod)
    for (i in 1:M) {
      newimp <- obsData
      outcomeModResVar <- (sigmasq*impMod$df) / rchisq(1,impMod$df)
      covariance <- (outcomeModResVar/sigmasq)*varcov
      outcomeModBeta <- beta + MASS::mvrnorm(1, mu=rep(0,ncol(covariance)), Sigma=covariance)
      impMod$coefficients <- outcomeModBeta
      fittedVals <- predict(impMod, obsData)
      newimp[is.na(obsData[,outcomeVar]),outcomeVar] <- fittedVals[is.na(obsData[,outcomeVar])] +
        rnorm(sum(is.na(obsData[,outcomeVar])), mean=0, sd=outcomeModResVar^0.5)
      if (M>1) {
        imps[[i]] <- newimp
      }
    }
  }
  if (M>1) {
    attr(imps, "pd") <- pd
    imps
  } else {
    attr(newimp, "pd") <- pd
    newimp
  }
}
