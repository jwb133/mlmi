normUniImp <- function(inputData, impFormula, M=5) {
  #put some checks in to check there are no missing values in covariates of imp model

  #fit desired imputation model
  impMod <- lm(impFormula, data=inputData)
  fittedVals <- predict(impMod, inputData)

  outcomeVar <- all.vars(impFormula)[1]
  print(outcomeVar)

  imps <- vector("list", M)
  for (i in 1:M) {
    imps[[i]] <- inputData
    #impute
    imps[[i]][is.na(inputData[,outcomeVar]),outcomeVar] <- fittedVals[is.na(inputData[,outcomeVar])] +
      rnorm(sum(is.na(inputData[,outcomeVar])), mean=0, sd=sigma(impMod))
  }

  imps
}
