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
#' @return A list of imputed datasets.
#' @export
normUniImp <- function(obsData, impFormula, M=5, pd=FALSE) {
  #put some checks in to check there are no missing values in covariates of imp model

  #fit desired imputation model
  impMod <- lm(impFormula, data=obsData)

  outcomeVar <- all.vars(impFormula)[1]

  imps <- vector("list", M)

  if (pd==FALSE) {
    fittedVals <- predict(impMod, obsData)
    print(fittedVals)
    for (i in 1:M) {
      imps[[i]] <- obsData
      #impute
      imps[[i]][is.na(obsData[,outcomeVar]),outcomeVar] <- fittedVals[is.na(obsData[,outcomeVar])] +
        rnorm(sum(is.na(obsData[,outcomeVar])), mean=0, sd=sigma(impMod))
    }
  } else {
    beta <- impMod$coef
    sigmasq <- summary(impMod)$sigma^2
    varcov <- vcov(impMod)
    for (i in 1:M) {
      imps[[i]] <- obsData
      outcomeModResVar <- (sigmasq*impMod$df) / rchisq(1,impMod$df)
      covariance <- (outcomeModResVar/sigmasq)*varcov
      outcomeModBeta <- beta + MASS::mvrnorm(1, mu=rep(0,ncol(covariance)), Sigma=covariance)
      impMod$coefficients <- outcomeModBeta
      fittedVals <- predict(impMod, obsData)
      imps[[i]][is.na(obsData[,outcomeVar]),outcomeVar] <- fittedVals[is.na(obsData[,outcomeVar])] +
        rnorm(sum(is.na(obsData[,outcomeVar])), mean=0, sd=outcomeModResVar^0.5)
    }
  }
  attributes(imps) <- list(pd=pd)
  imps
}
