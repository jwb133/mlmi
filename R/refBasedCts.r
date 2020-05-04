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
#' Baseline covariates must be numeric variables. If you have factor variables you must
#' code these into suitable dummy indicators and pass these to the function.
#'
#' THINGS TO DO:
#' 1) option for INTERACTIONS BETWEEN BASELINE AND VISIT
#' 2) implement copy reference method
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
refBasedCts <- function(obsData, outcomeVarStem, nVisits, trtVar, baselineVars=NULL, type="MAR", baselineVisitInt=FALSE, M=5) {
  if (nVisits<2) {
    stop("You must have at least 2 post-baseline visits.")
  }
  if (!is.null(baselineVars)) {
    #check all baseline covariates are not factors
    for (i in 1:length(baselineVars)) {
      if (is.factor(obsData[[baselineVars[i]]])) {
        stop("Factor baseline variables are not allowed currently. Please code dummy variables and pass these.")
      }
    }
  }
  imps <- vector("list", M)

  trtCol <- which(colnames(obsData)==trtVar)
  outcomeCols <- which(colnames(obsData) %in% paste(outcomeVarStem, 1:nVisits, sep=""))
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
    mixedFormula <- formula(paste(outcomeVarStem, "~ factor(time)+", paste(baselineVars, collapse="+"), sep=""))
  }
  controlMod <- nlme::gls(mixedFormula,
                    na.action=na.omit, data=controlLong,
                    correlation=corSymm(form=~time | mlmiId),
                    weights=varIdent(form=~1|time))

  controlCov <- mmrmCov(controlMod)
  #create matrix of covariate conditional means
  if (baselineVisitInt==FALSE) {
    if (length(baselineVars)==1) {
      meanMat <- controlWide[,baselineVars]*tail(coef(controlMod),1)
    } else if (length(baselineVars)>1) {
      meanMat <- as.matrix(subset(controlWide, select=baselineVars)) %*% array(tail(coef(controlMod),length(baselineVars)),
                                                                    dim=c(length(baselineVars), 1))
    } else {
      meanMat <- rep(0, controlN)
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

  #identify missingness patterns
  controlMissPat <- missingnessPatterns(yObs)

  for (imp in 1:M) {
    yImp <- yObs

    if (is.null(controlMissPat$nonMonotonePatterns)==FALSE) {
      #impute intermediate missingness assuming MAR
      for (pat in 1:nrow(controlMissPat$nonMonotonePatterns)) {
        #determine which are the intermediate missing values in the pattern
        visitsToImpute <- as.numeric(intermediateDetect(controlMissPat$nonMonotonePatterns[pat,]))
        visitsObs <- as.numeric(which(controlMissPat$nonMonotonePatterns[pat,]==1))
        #identify individuals with this pattern
        rowsToImpute <- which(apply(1*(!is.na(yObs)), 1, function(x) identical(x, controlMissPat$nonMonotonePatterns[pat,])))
        #calculate conditional covariance matrix
        condVar <- controlCov[visitsToImpute,visitsToImpute] - array(controlCov[visitsToImpute,visitsObs],
          dim=c(length(visitsToImpute),length(visitsObs))) %*%
          solve(controlCov[visitsObs,visitsObs]) %*% t(array(controlCov[visitsToImpute,visitsObs], dim=c(length(visitsToImpute),length(visitsObs))))

        yImp[rowsToImpute,visitsToImpute] <- mvrnorm(n=length(rowsToImpute), mu=rep(0,length(visitsToImpute)), Sigma=condVar)
        #calculate mean conditional on observed outcomes
        condMean <- meanMat[rowsToImpute,visitsToImpute] + array(as.matrix(yImp[rowsToImpute,visitsObs]) - meanMat[rowsToImpute,visitsObs],dim=c(length(rowsToImpute),length(visitsObs))) %*%
          solve(controlCov[visitsObs,visitsObs]) %*% t(array(controlCov[visitsToImpute,visitsObs], dim=c(length(visitsToImpute),length(visitsObs))))
        yImp[rowsToImpute,visitsToImpute] <- yImp[rowsToImpute,visitsToImpute] + condMean
      }
    }

    #now impute any monotone missingness
    if (is.null(controlMissPat$monotonePatterns)==FALSE) {
      for (pat in 1:nrow(controlMissPat$monotonePatterns)) {
        if (sum(controlMissPat$monotonePatterns[pat,])==0) {
          #every visit missing
          yImp[is.na(yObs[,1]),] <- mvrnorm(n=sum(is.na(yObs[,1])), mu=rep(0,nVisits), Sigma=controlCov)
          yImp[is.na(yObs[,1]),] <- yImp[is.na(yObs[,1]),] + meanMat[is.na(yObs[,1]),]
        } else if (sum(controlMissPat$monotonePatterns[pat,]) < nVisits) {
          visitsToImpute <- as.numeric(which(controlMissPat$monotonePatterns[pat,]==0))
          visitsObs <- as.numeric(which(controlMissPat$monotonePatterns[pat,]==1))
          #identify individuals with this pattern
          rowsToImpute <- which(apply(1*(!is.na(yObs)), 1, function(x) identical(x, controlMissPat$monotonePatterns[pat,])))
          #calculate conditional covariance matrix
          condVar <- controlCov[visitsToImpute,visitsToImpute] - array(controlCov[visitsToImpute,visitsObs],
                                                                       dim=c(length(visitsToImpute),length(visitsObs))) %*%
            solve(controlCov[visitsObs,visitsObs]) %*% t(array(controlCov[visitsToImpute,visitsObs], dim=c(length(visitsToImpute),length(visitsObs))))

          yImp[rowsToImpute,visitsToImpute] <- mvrnorm(n=length(rowsToImpute), mu=rep(0,length(visitsToImpute)), Sigma=condVar)
          #calculate mean conditional on observed outcomes
          condMean <- meanMat[rowsToImpute,visitsToImpute] + array(as.matrix(yImp[rowsToImpute,visitsObs]) - meanMat[rowsToImpute,visitsObs],dim=c(length(rowsToImpute),length(visitsObs))) %*%
            solve(controlCov[visitsObs,visitsObs]) %*% t(array(controlCov[visitsToImpute,visitsObs], dim=c(length(visitsToImpute),length(visitsObs))))
          yImp[rowsToImpute,visitsToImpute] <- yImp[rowsToImpute,visitsToImpute] + condMean
        }
      }
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
  #create matrix of covariate conditional means under both active and control model fits
  activeMeanMat <- array(0, dim=c(activeN,nVisits))
  controlMeanMat <- array(0, dim=c(activeN,nVisits))
  if (baselineVisitInt==FALSE) {
    if (length(baselineVars)==1) {
      activeMeanMat <- activeWide[,baselineVars]*tail(coef(activeMod),1)
      controlMeanMat <- activeWide[,baselineVars]*tail(coef(controlMod),1)
    } else if (length(baselineVars)>1) {
      activeMeanMat <- as.matrix(subset(activeWide, select=baselineVars)) %*% array(tail(coef(activeMod),length(baselineVars)),
                                                                    dim=c(length(baselineVars), 1))
      controlMeanMat <- as.matrix(subset(activeWide, select=baselineVars)) %*% array(tail(coef(controlMod),length(baselineVars)),
                                                                         dim=c(length(baselineVars), 1))
    } else {
      #no baseline covariates
      activeMeanMat <- rep(0, activeN)
      controlMeanMat <- rep(0, activeN)
    }
    activeMeanMat <- array(rep(activeMeanMat, nVisits),dim=c(activeN,nVisits))
    controlMeanMat <- array(rep(controlMeanMat, nVisits),dim=c(activeN,nVisits))
  } else {
    stop("Interactions between baseline covariates and visit is not yet coded up")
  }
  #now add in visit effects
  activeMeanMat[,1] <- activeMeanMat[,1] + coef(activeMod)[1]
  controlMeanMat[,1] <- controlMeanMat[,1] + coef(controlMod)[1]
  for (i in 2:nVisits) {
    activeMeanMat[,i] <- activeMeanMat[,i] + coef(activeMod)[1] + coef(activeMod)[i]
    controlMeanMat[,i] <- controlMeanMat[,i] + coef(controlMod)[1] + coef(controlMod)[i]
  }

  #create data frame of just the outcome columns
  yObs <- subset(activeWide, select=paste(outcomeVarStem, 1:nVisits, sep=""))

  #identify missingness patterns
  activeMissPat <- missingnessPatterns(yObs)

  for (imp in 1:M) {
    yImp <- yObs

    if (is.null(activeMissPat$nonMonotonePatterns)==FALSE) {
      #impute intermediate missingness assuming MAR
      for (pat in 1:nrow(activeMissPat$nonMonotonePatterns)) {
        #determine which are the intermediate missing values in the pattern
        visitsToImpute <- as.numeric(intermediateDetect(activeMissPat$nonMonotonePatterns[pat,]))
        visitsObs <- as.numeric(which(activeMissPat$nonMonotonePatterns[pat,]==1))
        #identify individuals with this pattern
        rowsToImpute <- which(apply(1*(!is.na(yObs)), 1, function(x) identical(x, activeMissPat$nonMonotonePatterns[pat,])))
        #calculate conditional covariance matrix
        condVar <- activeCov[visitsToImpute,visitsToImpute] - array(activeCov[visitsToImpute,visitsObs],
                                                                     dim=c(length(visitsToImpute),length(visitsObs))) %*%
          solve(activeCov[visitsObs,visitsObs]) %*% t(array(activeCov[visitsToImpute,visitsObs], dim=c(length(visitsToImpute),length(visitsObs))))

        yImp[rowsToImpute,visitsToImpute] <- mvrnorm(n=length(rowsToImpute), mu=rep(0,length(visitsToImpute)), Sigma=condVar)
        #calculate mean conditional on observed outcomes
        condMean <- activeMeanMat[rowsToImpute,visitsToImpute] + array(as.matrix(yImp[rowsToImpute,visitsObs]) - activeMeanMat[rowsToImpute,visitsObs],dim=c(length(rowsToImpute),length(visitsObs))) %*%
          solve(activeCov[visitsObs,visitsObs]) %*% t(array(activeCov[visitsToImpute,visitsObs], dim=c(length(visitsToImpute),length(visitsObs))))
        yImp[rowsToImpute,visitsToImpute] <- yImp[rowsToImpute,visitsToImpute] + condMean
      }
    }

    #now impute any monotone missingness
    if (is.null(activeMissPat$monotonePatterns)==FALSE) {
      for (pat in 1:nrow(activeMissPat$monotonePatterns)) {
        if (sum(activeMissPat$monotonePatterns[pat,])==0) {
          #every visit missing
          if (type=="J2R") {
            yImp[is.na(yObs[,1]),] <- mvrnorm(n=sum(is.na(yObs[,1])), mu=rep(0,nVisits), Sigma=controlCov)
            yImp[is.na(yObs[,1]),] <- yImp[is.na(yObs[,1]),] + controlMeanMat[is.na(yObs[,1]),]
          } else {
            #MAR
            yImp[is.na(yObs[,1]),] <- mvrnorm(n=sum(is.na(yObs[,1])), mu=rep(0,nVisits), Sigma=activeCov)
            yImp[is.na(yObs[,1]),] <- yImp[is.na(yObs[,1]),] + activeMeanMat[is.na(yObs[,1]),]
          }
        } else if (sum(activeMissPat$monotonePatterns[pat,]) < nVisits) {
          visitsToImpute <- as.numeric(which(activeMissPat$monotonePatterns[pat,]==0))
          visitsObs <- as.numeric(which(activeMissPat$monotonePatterns[pat,]==1))
          #identify individuals with this pattern
          rowsToImpute <- which(apply(1*(!is.na(yObs)), 1, function(x) identical(x, activeMissPat$monotonePatterns[pat,])))
          if (type=="J2R") {
            #calculate conditional covariance matrix
            condVar <- controlCov[visitsToImpute,visitsToImpute] - array(controlCov[visitsToImpute,visitsObs],
                                                                        dim=c(length(visitsToImpute),length(visitsObs))) %*%
              solve(controlCov[visitsObs,visitsObs]) %*% t(array(controlCov[visitsToImpute,visitsObs], dim=c(length(visitsToImpute),length(visitsObs))))

            #calculate mean conditional on observed outcomes
            condMean <- controlMeanMat[rowsToImpute,visitsToImpute] + array(as.matrix(yImp[rowsToImpute,visitsObs]) - activeMeanMat[rowsToImpute,visitsObs],dim=c(length(rowsToImpute),length(visitsObs))) %*%
              solve(controlCov[visitsObs,visitsObs]) %*% t(array(controlCov[visitsToImpute,visitsObs], dim=c(length(visitsToImpute),length(visitsObs))))
          } else {
            #calculate conditional covariance matrix
            condVar <- activeCov[visitsToImpute,visitsToImpute] - array(activeCov[visitsToImpute,visitsObs],
                                                                         dim=c(length(visitsToImpute),length(visitsObs))) %*%
              solve(activeCov[visitsObs,visitsObs]) %*% t(array(activeCov[visitsToImpute,visitsObs], dim=c(length(visitsToImpute),length(visitsObs))))

            #calculate mean conditional on observed outcomes
            condMean <- activeMeanMat[rowsToImpute,visitsToImpute] + array(as.matrix(yImp[rowsToImpute,visitsObs]) - activeMeanMat[rowsToImpute,visitsObs],dim=c(length(rowsToImpute),length(visitsObs))) %*%
              solve(activeCov[visitsObs,visitsObs]) %*% t(array(activeCov[visitsToImpute,visitsObs], dim=c(length(visitsToImpute),length(visitsObs))))
          }
          yImp[rowsToImpute,visitsToImpute] <- mvrnorm(n=length(rowsToImpute), mu=rep(0,length(visitsToImpute)), Sigma=condVar)
          yImp[rowsToImpute,visitsToImpute] <- yImp[rowsToImpute,visitsToImpute] + condMean
        }
      }
    }

    imps[[imp]][obsData[,trtCol]==1,outcomeCols] <- yImp
  }

  if (M>1) {
    imps
  } else {
    imps[[1]]
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

#function which when passed a data matrix, finds the missing data patterns
#which occur
missingnessPatterns <- function(yObs) {
  rMat <- 1-1*(is.na(yObs))
  uniquePatterns <- unique(rMat)
  nonMonotoneInd <- rep(0, nrow(uniquePatterns))
  for (i in 1:nrow(uniquePatterns)) {
    if (identical(intermediateDetect(uniquePatterns[i,]),NA)) {
      nonMonotoneInd[i] <- FALSE
    } else {
      nonMonotoneInd[i] <- TRUE
    }
  }
  monotoneInd <- !nonMonotoneInd

  if (sum(nonMonotoneInd)==0) {
    #no intermediate missingness
    nonMonotonePatterns <- NULL
  } else {
    #some intermediate missingness
    if (sum(nonMonotoneInd)>1) {
      nonMonotonePatterns <- uniquePatterns[nonMonotoneInd==TRUE,]
    } else {
      nonMonotonePatterns <- matrix(uniquePatterns[nonMonotoneInd==TRUE,], nrow=1)
    }
  }
  if (sum(monotoneInd)==0) {
    #no monotone missingness
    monotonePatterns <- NULL
  } else {
    #some monotone missingness
    if (sum(monotoneInd)>1) {
      monotonePatterns <- uniquePatterns[monotoneInd==TRUE,]
    } else {
      monotonePatterns <- matrix(uniquePatterns[monotoneInd==TRUE,], nrow=1)
    }
  }
  list(monotonePatterns=monotonePatterns, nonMonotonePatterns=nonMonotonePatterns)
}

#function which when passed a vector of 0s and 1s detects intemediate missingness and returns
#positions on intermediate missing values
intermediateDetect <- function(obsVec) {
  missPos <- which(obsVec==0)
  lastObsPos <- tail(which(obsVec==1),1)
  if (length(missPos[missPos<lastObsPos])>0) {
    missPos[missPos<lastObsPos]
  } else {
    NA
  }
}
