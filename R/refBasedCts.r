#' Reference based imputation of repeated measures continuous data
#'
#' Performs multiple imputation of a repeatedly measured continuous endpoint in a
#' randomised clinical trial
#' using reference based imputation as proposed by \doi{10.1080/10543406.2013.834911}{Carpenter et al (2013)}. This approach
#' can be used for imputation of missing data in randomised clinical trials.
#'
#' Unlike most implementations of reference based imputation, this implementation
#' imputes conditional on the maximum likelihood estimates of the model parameters,
#' rather than a posterior draw. If one is interested in frequentist valid inferences,
#' this is ok provided the bootstrapping used, for example with using the
#' \href{https://cran.r-project.org/package=bootImpute}{bootImpute} package.
#'
#' Intermediate missing values are imputed assuming MAR, based on the mixed model fit
#' to that patient's treatment arm. Monotone missing values are imputed using the specified
#' imputation type.
#'
#' Baseline covariates must be numeric variables. If you have factor variables you must
#' code these into suitable dummy indicators and pass these to the function.
#'
#' @param obsData The data frame to be imputed.
#' @param outcomeVarStem String for stem of outcome variable name, e.g. y if y1, y2, y3 are the outcome columns
#' @param nVisits The integer number of visits (not including baseline)
#' @param trtVar The string variable name of the randomised treatment group variable. The reference arm is assumed
#' to correspond to \code{trtVar==0}.
#' @param baselineVars A string or vector of strings specfying the baseline variables. Often this will include
#' the baseline measurement of the outcome
#' @param baselineVisitInt TRUE/FALSE indicating whether to allow for interactions between each baseline variable
#' and visit. Default is TRUE.
#' @param type A string specifying imputation type to use. Valid options are "MAR", "J2R"
#'
#' @param M Number of imputations to generate.
#' @return A list of imputed datasets, or if \code{M=1}, just the imputed data frame.
#'
#' @references Carpenter JR, Roger JH, Kenward MG. Analysis of Longitudinal Trials with Protocol Deviation:
#' A Framework for Relevant, Accessible Assumptions, and Inference via Multiple Imputation. (2013) 23(6) 1352-1371
#' @references von Hippel PT & Bartlett JW (2019) Maximum likelihood multiple imputation: Faster imputations and
#' consistent standard errors without posterior draws \href{https://arxiv.org/abs/1210.0870v10}{arXiv:1210.0870v10}.
#'
#' @example data-raw/refBasedCtsExample.R
#'
#' @export
refBasedCts <- function(obsData, outcomeVarStem, nVisits, trtVar, baselineVars=NULL,
                        baselineVisitInt=TRUE, type="MAR", M=5) {
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
  } else if (baselineVisitInt==FALSE) {
    mixedFormula <- formula(paste(outcomeVarStem, "~ factor(time)+", paste(baselineVars, collapse="+"), sep=""))
  } else {
    #interactions between time and baseline covariates
    mixedFormula <- formula(paste(outcomeVarStem, "~ ",
                                  paste("factor(time)*", baselineVars, collapse="+", sep=""), sep=""))
  }
  controlMod <- nlme::gls(mixedFormula,
                    na.action=na.omit, data=controlLong,
                    correlation=nlme::corSymm(form=~time | mlmiId),
                    weights=nlme::varIdent(form=~1|time))

  controlCov <- mmrmCov(controlMod)
  #create matrix of covariate conditional means
  if (baselineVisitInt==FALSE) {
    if (length(baselineVars)==1) {
      meanMat <- controlWide[,baselineVars]*utils::tail(coef(controlMod),1)
    } else if (length(baselineVars)>1) {
      meanMat <- as.matrix(subset(controlWide, select=baselineVars)) %*% array(utils::tail(coef(controlMod),length(baselineVars)),
                                                                    dim=c(length(baselineVars), 1))
    } else {
      meanMat <- rep(0, controlN)
    }
    meanMat <- array(rep(meanMat, nVisits),dim=c(controlN,nVisits))

    #now add in visit effects
    meanMat[,1] <- meanMat[,1] + coef(controlMod)[1]
    for (i in 2:nVisits) {
      meanMat[,i] <- meanMat[,i] + coef(controlMod)[1] + coef(controlMod)[i]
    }
  } else {
    #interactions between time and baseline variables
    if (length(baselineVars)==1) {
      #intercept + covariate effect at visit 1
      meanMat <- coef(controlMod)[1] + controlWide[,baselineVars]*coef(controlMod)[nVisits+1]
      meanMat <- array(rep(meanMat, nVisits),dim=c(controlN,nVisits))

      for (i in 2:nVisits) {
        #visit effect plus extra effect of baseline variable at this visit
        meanMat[,i] <- meanMat[,i] + coef(controlMod)[i] + controlWide[,baselineVars]*coef(controlMod)[nVisits+i]
      }
    } else if (length(baselineVars)>1) {
      #intercept + covariate effects at visit 1
      meanMat <- rep(coef(controlMod)[1], controlN)
      for (i in 1:length(baselineVars)) {
        meanMat <- meanMat + controlWide[,baselineVars[i]]*coef(controlMod)[nVisits+i]
      }
      meanMat <- array(rep(meanMat, nVisits),dim=c(controlN,nVisits))

      for (i in 2:nVisits) {
        #add in visit effect first
        meanMat[,i] <- meanMat[,i] + coef(controlMod)[i]
        #add in baseline visit interaction effects
        for (j in 1:length(baselineVars)) {
          #figure out which coefficient we need for this visit and baseline covariate
          coefIndexNeeded <- nVisits + length(baselineVars) + (nVisits-1)*(j-1) + (i-1)
          meanMat[,i] <- meanMat[,i] + controlWide[,baselineVars[j]]*coef(controlMod)[coefIndexNeeded]
        }
      }
    } else {
      stop("You have specified interactions between baseline variables and time, but you have specified any baseline variables.")
    }
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

        yImp[rowsToImpute,visitsToImpute] <- MASS::mvrnorm(n=length(rowsToImpute), mu=rep(0,length(visitsToImpute)), Sigma=condVar)
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
          yImp[is.na(yObs[,1]),] <- MASS::mvrnorm(n=sum(is.na(yObs[,1])), mu=rep(0,nVisits), Sigma=controlCov)
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

          yImp[rowsToImpute,visitsToImpute] <- MASS::mvrnorm(n=length(rowsToImpute), mu=rep(0,length(visitsToImpute)), Sigma=condVar)
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
                          correlation=nlme::corSymm(form=~time | mlmiId),
                          weights=nlme::varIdent(form=~1|time))
  #print(summary(activeMod))

  activeCov <- mmrmCov(activeMod)
  #create matrix of covariate conditional means for active patients under both active and control model fits
  if (baselineVisitInt==FALSE) {
    if (length(baselineVars)==1) {
      activeMeanMat <- activeWide[,baselineVars]*utils::tail(coef(activeMod),1)
      controlMeanMat <- activeWide[,baselineVars]*utils::tail(coef(controlMod),1)
    } else if (length(baselineVars)>1) {
      activeMeanMat <- as.matrix(subset(activeWide, select=baselineVars)) %*% array(utils::tail(coef(activeMod),length(baselineVars)),
                                                                    dim=c(length(baselineVars), 1))
      controlMeanMat <- as.matrix(subset(activeWide, select=baselineVars)) %*% array(utils::tail(coef(controlMod),length(baselineVars)),
                                                                         dim=c(length(baselineVars), 1))
    } else {
      #no baseline covariates
      activeMeanMat <- rep(0, activeN)
      controlMeanMat <- rep(0, activeN)
    }
    activeMeanMat <- array(rep(activeMeanMat, nVisits),dim=c(activeN,nVisits))
    controlMeanMat <- array(rep(controlMeanMat, nVisits),dim=c(activeN,nVisits))

    #now add in visit effects
    activeMeanMat[,1] <- activeMeanMat[,1] + coef(activeMod)[1]
    controlMeanMat[,1] <- controlMeanMat[,1] + coef(controlMod)[1]
    for (i in 2:nVisits) {
      activeMeanMat[,i] <- activeMeanMat[,i] + coef(activeMod)[1] + coef(activeMod)[i]
      controlMeanMat[,i] <- controlMeanMat[,i] + coef(controlMod)[1] + coef(controlMod)[i]
    }
  } else {
    #interactions between time and baseline variables
    if (length(baselineVars)==1) {
      #intercept + covariate effect at visit 1
      controlMeanMat <- coef(controlMod)[1] + activeWide[,baselineVars]*coef(controlMod)[nVisits+1]
      controlMeanMat <- array(rep(controlMeanMat, nVisits),dim=c(activeN,nVisits))
      activeMeanMat <- coef(activeMod)[1] + activeWide[,baselineVars]*coef(activeMod)[nVisits+1]
      activeMeanMat <- array(rep(activeMeanMat, nVisits),dim=c(activeN,nVisits))

      for (i in 2:nVisits) {
        #visit effect plus extra effect of baseline variable at this visit
        controlMeanMat[,i] <- controlMeanMat[,i] + coef(controlMod)[i] + activeWide[,baselineVars]*coef(controlMod)[nVisits+i]
        activeMeanMat[,i] <- activeMeanMat[,i] + coef(activeMod)[i] + activeWide[,baselineVars]*coef(activeMod)[nVisits+i]
      }
    } else if (length(baselineVars)>1) {
      #intercept + covariate effects at visit 1
      controlMeanMat <- rep(coef(controlMod)[1], activeN)
      activeMeanMat <- rep(coef(activeMod)[1], activeN)

      for (i in 1:length(baselineVars)) {
        controlMeanMat <- controlMeanMat + activeWide[,baselineVars[i]]*coef(controlMod)[nVisits+i]
        activeMeanMat <- activeMeanMat + activeWide[,baselineVars[i]]*coef(activeMod)[nVisits+i]
      }
      controlMeanMat <- array(rep(controlMeanMat, nVisits),dim=c(activeN,nVisits))
      activeMeanMat <- array(rep(activeMeanMat, nVisits),dim=c(activeN,nVisits))

      for (i in 2:nVisits) {
        #add in visit effect first
        controlMeanMat[,i] <- controlMeanMat[,i] + coef(controlMod)[i]
        activeMeanMat[,i] <- activeMeanMat[,i] + coef(activeMod)[i]
        #add in baseline visit interaction effects
        for (j in 1:length(baselineVars)) {
          #figure out which coefficient we need for this visit and baseline covariate
          coefIndexNeeded <- nVisits + length(baselineVars) + (nVisits-1)*(j-1) + (i-1)
          controlMeanMat[,i] <- controlMeanMat[,i] + activeWide[,baselineVars[j]]*coef(controlMod)[coefIndexNeeded]
          activeMeanMat[,i] <- activeMeanMat[,i] + activeWide[,baselineVars[j]]*coef(activeMod)[coefIndexNeeded]
        }
      }
    }
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

        yImp[rowsToImpute,visitsToImpute] <- MASS::mvrnorm(n=length(rowsToImpute), mu=rep(0,length(visitsToImpute)), Sigma=condVar)
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
            yImp[is.na(yObs[,1]),] <- MASS::mvrnorm(n=sum(is.na(yObs[,1])), mu=rep(0,nVisits), Sigma=controlCov)
            yImp[is.na(yObs[,1]),] <- yImp[is.na(yObs[,1]),] + controlMeanMat[is.na(yObs[,1]),]
          } else {
            #MAR
            yImp[is.na(yObs[,1]),] <- MASS::mvrnorm(n=sum(is.na(yObs[,1])), mu=rep(0,nVisits), Sigma=activeCov)
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
          yImp[rowsToImpute,visitsToImpute] <- MASS::mvrnorm(n=length(rowsToImpute), mu=rep(0,length(visitsToImpute)), Sigma=condVar)
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
      nonMonotonePatterns <- t(as.matrix(uniquePatterns[nonMonotoneInd==TRUE,]))
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
      monotonePatterns <- t(as.matrix(uniquePatterns[monotoneInd==TRUE,]))
    }
  }
  list(monotonePatterns=monotonePatterns, nonMonotonePatterns=nonMonotonePatterns)
}

#function which when passed a vector of 0s and 1s detects intemediate missingness and returns
#positions on intermediate missing values
intermediateDetect <- function(obsVec) {
  missPos <- which(obsVec==0)
  lastObsPos <- utils::tail(which(obsVec==1),1)
  if (length(missPos[missPos<lastObsPos])>0) {
    missPos[missPos<lastObsPos]
  } else {
    NA
  }
}
