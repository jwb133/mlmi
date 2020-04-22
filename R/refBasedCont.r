#' Reference based imputation of repeated measures continuous data
#'
#' Performs multiple imputation of a repeatedly measured continuous endpoint
#' using reference based imputation as proposed by Carpenter et al. This approach
#' can be used for imputation of missing data in randomised clinical trials.
#'
#' Due to uncongeniality between the imputation and analysis models, the resulting
#' imputed datasets should be analysed using the
#' \href{https://cran.r-project.org/package=bootImpute}{bootImpute} package.
#' Unlike most implementations of reference based imputation, this implementation
#' imputes conditional on the maximum likelihood estimates of the model parameters,
#' rather than a posterior draw. This is ok provided the aforementioned bootstrapping
#' approach is used for inference, rather than Rubin's rules.
#'
#' THINGS TO DO:
#' 1) J2R IMPLEMENTATION
#' 2) NON-MONOTONE MISSINGNESS
#' 3) option for INTERACTIONS BETWEEN BASELINE AND VISIT
#'
#' @param obsData The data frame to be imputed.
#' @param outcomeVarStem String for stem of outcome variable name, e.g. y if y1, y2, y3 are the outcome columns
#' @param nVisits The integer number of visits (not including baseline)
#' @param trtVar A string giving the name of the randomised treatment group variable
#' @param M Number of imputations to generate.
#' @return A list of imputed datasets, or if \code{M=1}, just the imputed data frame.
#'
#' @references von Hippel P.T. (2018) Maximum likelihood multiple imputation: faster,
#' more efficient imputation without posterior draws. \href{https://arxiv.org/abs/1210.0870v9}{arXiv:1210.0870v9}.
#'
#' @example data-raw/normUniExample.r
#'
#' @export
refBasedCts <- function(obsData, outcomeVarStem, nVisits, trtVar, baselineVars=NULL, baselineVisitInt=FALSE, M=5) {
  imps <- vector("list", M)

  trtCol <- which(colnames(obsData)==trtVar)
  outcomeCols <- which(colnames(obsData)==paste(outcomeVarStem, 1:nVisits, sep=""))
  controlWide <- obsData[obsData[,trtCol]==0,]
  controlN <- dim(controlWide)[1]
  #reshape to long form
  controlLong <- reshape(controlWide, varying=paste(outcomeVarStem, 1:nVisits, sep=""),
                         direction="long", sep="", idvar="mlmiId")
  controlLong <- controlLong[order(controlLong$mlmiId,controlLong$time),]

  #fit mixed model assuming MAR
  if (is.null(baselineVars)) {
    mixedFormula <- formula(paste(outcomeVarStem, "~ factor(time)", sep=""))
  } else {
    mixedFormula <- formula(paste(outcomeVarStem, "~ factor(time)+", paste(baselineVars, sep="+"), sep=""))
  }
  controlMod <- nlme::gls(mixedFormula,
                    na.action=na.omit, data=controlLong,
                    correlation=corSymm(form=~time | mlmiId),
                    weights=varIdent(form=~1|time))
  #print(summary(controlMod))

  controlCov <- mmrmCov(controlMod)
  #create matrix of covariate conditional means
  meanMat <- array(0, dim=c(controlN,nVisits))
  if (baselineVisitInt==FALSE) {
    if (length(baselineVars)==1) {
      meanMat <- controlWide[,baselineVars]*tail(coef(controlMod),1)
    } else {
      meanMat <- subset(controlWide, select=baselineVars) %*% array(tail(coef(controlMod),length(baselineVars)),
                                                                    dim=c(length(baselineVars), 1))
    }
    meanMat <- array(rep(meanMat, nVisits),dim=c(controlN,nVisits))
  } else {
    stop("Interactions between baseline covariates and visit is not yet coded up")
  }
  #now add in visit effects
  meanMat[,1] <- meanMat[,1] + coef(controlMod)[1]
  for (i in 2:nVisits) {
    meanMat[,i] <- meanMat[,i] + coef(controlMod)[1] + coef(controlMod)[i]
  }

  #create data frame of just the outcome columns
  yObs <- subset(controlWide, select=paste(outcomeVarStem, 1:nVisits, sep=""))
  for (imp in 1:M) {
    yImp <- yObs

    #impute for those who are missing all measurements
    yImp[is.na(yObs[,1]),] <- mvrnorm(n=sum(is.na(yObs[,1])), mu=rep(0,nVisits), Sigma=controlCov)
    yImp[is.na(yObs[,1]),] <- yImp[is.na(yObs[,1]),] + meanMat[is.na(yObs[,1]),]

    #impute for each drop-out pattern up to those who just missed final time point
    for (i in 2:nVisits) {
      #calculate conditional covariance matrix
      condVar <- controlCov[i:nVisits,i:nVisits] - array(controlCov[i:nVisits,1:(i-1)], dim=c(nVisits-i+1,i-1)) %*% solve(controlCov[1:(i-1),1:(i-1)]) %*% t(array(controlCov[i:nVisits,1:(i-1)], dim=c(nVisits-i+1,i-1)))
      yImp[is.na(yObs[,i]),i:nVisits] <- mvrnorm(n=sum(is.na(yObs[,i])), mu=rep(0,(nVisits-i+1)), Sigma=condVar)
      #calculate mean conditional on observed outcomes
      condMean <- meanMat[is.na(yObs[,i]),i:nVisits] + array(as.matrix(yImp[is.na(yObs[,i]),1:(i-1)]) - meanMat[is.na(yObs[,i]),1:(i-1)],dim=c(sum(is.na(yObs[,i])),i-1)) %*%
         solve(controlCov[1:(i-1),1:(i-1)]) %*% t(array(controlCov[i:nVisits,1:(i-1)], dim=c(nVisits-i+1,i-1)))
      yImp[is.na(yObs[,i]),i:nVisits] <- yImp[is.na(yObs[,i]),i:nVisits] + condMean
    }

    imps[[imp]] <- obsData
    imps[[imp]][obsData[,trtCol]==0,outcomeCols] <- yImp
  }

  #now impute active arm
  activeWide <- obsData[obsData[,trtCol]==1,]
  activeN <- dim(activeWide)[1]
  #reshape to long form
  activeLong <- reshape(activeWide, varying=paste(outcomeVarStem, 1:nVisits, sep=""),
                         direction="long", sep="", idvar="mlmiId")
  activeLong <- activeLong[order(activeLong$mlmiId,activeLong$time),]

  activeMod <- nlme::gls(mixedFormula,
                          na.action=na.omit, data=activeLong,
                          correlation=corSymm(form=~time | mlmiId),
                          weights=varIdent(form=~1|time))
  #print(summary(activeMod))

  activeCov <- mmrmCov(activeMod)
  #create matrix of covariate conditional means
  meanMat <- array(0, dim=c(activeN,nVisits))
  if (baselineVisitInt==FALSE) {
    if (length(baselineVars)==1) {
      meanMat <- activeWide[,baselineVars]*tail(coef(activeMod),1)
    } else {
      meanMat <- subset(activeWide, select=baselineVars) %*% array(tail(coef(activeMod),length(baselineVars)),
                                                                    dim=c(length(baselineVars), 1))
    }
    meanMat <- array(rep(meanMat, nVisits),dim=c(activeN,nVisits))
  } else {
    stop("Interactions between baseline covariates and visit is not yet coded up")
  }
  #now add in visit effects
  meanMat[,1] <- meanMat[,1] + coef(activeMod)[1]
  for (i in 2:nVisits) {
    meanMat[,i] <- meanMat[,i] + coef(activeMod)[1] + coef(activeMod)[i]
  }

  #create data frame of just the outcome columns
  yObs <- subset(activeWide, select=paste(outcomeVarStem, 1:nVisits, sep=""))
  for (imp in 1:M) {
    yImp <- yObs

    #impute for those who are missing all measurements
    yImp[is.na(yObs[,1]),] <- mvrnorm(n=sum(is.na(yObs[,1])), mu=rep(0,nVisits), Sigma=activeCov)
    yImp[is.na(yObs[,1]),] <- yImp[is.na(yObs[,1]),] + meanMat[is.na(yObs[,1]),]

    #impute for each drop-out pattern up to those who just missed final time point
    for (i in 2:nVisits) {
      #calculate conditional covariance matrix
      condVar <- activeCov[i:nVisits,i:nVisits] - array(activeCov[i:nVisits,1:(i-1)], dim=c(nVisits-i+1,i-1)) %*% solve(activeCov[1:(i-1),1:(i-1)]) %*% t(array(activeCov[i:nVisits,1:(i-1)], dim=c(nVisits-i+1,i-1)))
      yImp[is.na(yObs[,i]),i:nVisits] <- mvrnorm(n=sum(is.na(yObs[,i])), mu=rep(0,(nVisits-i+1)), Sigma=condVar)
      #calculate mean conditional on observed outcomes
      condMean <- meanMat[is.na(yObs[,i]),i:nVisits] + array(as.matrix(yImp[is.na(yObs[,i]),1:(i-1)]) - meanMat[is.na(yObs[,i]),1:(i-1)],dim=c(sum(is.na(yObs[,i])),i-1)) %*%
        solve(activeCov[1:(i-1),1:(i-1)]) %*% t(array(activeCov[i:nVisits,1:(i-1)], dim=c(nVisits-i+1,i-1)))
      yImp[is.na(yObs[,i]),i:nVisits] <- yImp[is.na(yObs[,i]),i:nVisits] + condMean
    }

    imps[[imp]][obsData[,trtCol]==1,outcomeCols] <- yImp
  }

  if (M>1) {
    imps
  } else {
    imps[[1]]
  }
}


  impMod <- lm(impFormula, data=obsData)

  outcomeVar <- all.vars(impFormula)[1]

  if (M>1) {
    imps <- vector("list", M)
  }

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
  if (M>1) {
    imps
  } else {
    newimp
  }
}

#this function takes an MMRM fitted using GLS and creates the marginal covariance matrix
mmrmCov <- function(glsMod) {
  #construct correlation matrix
  nVisits <- length(coef(glsMod$modelStruct$varStruct, uncons = FALSE, allCoef = TRUE))
  mmrmCor <- matrix(1, nrow=nVisits, nVisits)
  mmrmCor[lower.tri(mmrmCor, diag=FALSE)] <- coef(glsMod$modelStruct$corStruct, uncons = FALSE, allCoef = TRUE)
  mmrmCor <- t(mmrmCor)
  mmrmCor[lower.tri(mmrmCor, diag=FALSE)] <- coef(glsMod$modelStruct$corStruct, uncons = FALSE, allCoef = TRUE)
  mmrmCor

  #construct vector of visit specific SDs
  mmrmSDs <- (glsMod$sigma) * sqrt(coef(glsMod$modelStruct$varStruct, uncons = FALSE, allCoef = TRUE))

  #construct covariance matrix
  diag(mmrmSDs) %*% mmrmCor %*% diag(mmrmSDs)
}
