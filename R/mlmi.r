normUniImp <- function(inputData, impFormula, M=5, pd=FALSE) {
  #put some checks in to check there are no missing values in covariates of imp model

  #fit desired imputation model
  impMod <- lm(impFormula, data=inputData)

  outcomeVar <- all.vars(impFormula)[1]
  print(outcomeVar)

  imps <- vector("list", M)

  if (pd==FALSE) {
    fittedVals <- predict(impMod, inputData)
    print(fittedVals)
    for (i in 1:M) {
      imps[[i]] <- inputData
      #impute
      imps[[i]][is.na(inputData[,outcomeVar]),outcomeVar] <- fittedVals[is.na(inputData[,outcomeVar])] +
        rnorm(sum(is.na(inputData[,outcomeVar])), mean=0, sd=sigma(impMod))
    }
  } else {
    beta <- impMod$coef
    sigmasq <- summary(impMod)$sigma^2
    varcov <- vcov(impMod)
    for (i in 1:M) {
      imps[[i]] <- inputData
      outcomeModResVar <- (sigmasq*impMod$df) / rchisq(1,impMod$df)
      covariance <- (outcomeModResVar/sigmasq)*varcov
      outcomeModBeta <- beta + MASS::mvrnorm(1, mu=rep(0,ncol(covariance)), Sigma=covariance)
      impMod$coefficients <- outcomeModBeta
      fittedVals <- predict(impMod, inputData)
      imps[[i]][is.na(inputData[,outcomeVar]),outcomeVar] <- fittedVals[is.na(inputData[,outcomeVar])] +
        rnorm(sum(is.na(inputData[,outcomeVar])), mean=0, sd=outcomeModResVar^0.5)
    }
  }

  imps
}
