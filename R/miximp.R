#' Imputation for a mixture of continuous and categorical variables using
#' the general location model.
#'
#' This function performs multiple imputation under a general location model
#' as described by Schafer (1997), using the \code{mix} package.
#'
#' @param obsData The data frame to be imputed. The categorical variables must be
#' in the first \code{nCat} columns, and they must be coded using consecutive positive
#' integers.
#' @param nCat The number of categorical variables in \code{obsData}.
#' @param M Number of imputations to generate.
#' @param pd Specify whether to use posterior draws (\code{TRUE})
#' or not (\code{FALSE}).
#' @param marginsType An integer specifying what type of log-linear model to use for the
#' categorical variables. \code{marginsType=1}, the default, allows for all two-way associations
#' in the log-linear model. \code{marginsType=2} allows for all three-way associations (plus lower).
#' \code{marginsType=3} assumes a saturated log-linear model for the categorical variables.
#' @param margins If \code{marginsType} is not specified, \code{margins} must be
#' supplied to specify the margins of the log-linear model for the categorical variable.
#' See the help for \code{mix::ecm.mix} for details on specifying \code{margins}.
#' @param designType An integer specifying how the continuous variables' means should depend
#' on the categorical variables. \code{designType=1}, the default, assumes the mean of each continuous
#' variable is a linear function with main effects of the categorical variables.
#' \code{designType=2} assumes each continuous variables has a separate mean for each
#' combination of the categorical variables.
#' @param design If \code{designType} is not specified, \code{design} must be supplied
#' to specify how the mean of the continuous variables depends on the categorical variables.
#' See the help for \code{mix::ecm.mix} for details on specifying \code{design}.
#' @param steps If \code{pd} is \code{TRUE}, the \code{steps} argument specifies
#' how many MCMC iterations to perform.
#' @param rseed The value to set the \code{mix} package's random number seed to,
#' using the \code{rngseed} function of \code{mix}. This function must be called at least
#' once before imputing using \code{mix}. If the user wishes to set the seed using
#' \code{rngseed} before calling \code{mixImp}, set \code{rseed=NULL}.
#' @return A list of imputed datasets, or if \code{M=1}, just the imputed data frame.
#' @export
mixImp <- function(obsData, nCat, M=10, pd=FALSE,
                   marginsType=1, margins=NULL, designType=1, design=NULL,
                   steps=100, rseed) {

  if (is.data.frame(obsData)==FALSE) {
    stop("obsData argumment must be a data frame.")
  }

  if (is.null(rseed)==FALSE) {
    mix::rngseed(rseed)
  }

  imps <- vector("list", M)

  s <- mix::prelim.mix(as.matrix(obsData), p=nCat)

  if ((marginsType==3) & (designType==2)) {
    #a completely unrestricted general location model is to be used
    print("Imputing using unrestricted general location model")
    restrictedMod <- 0
  } else {
    #a restricted model is to be used
    restrictedMod <- 1
    #figure out margins specification
    if (is.null(margins)==FALSE) {
      print("Imputing using specified margins argument")
    } else if (marginsType==1) {
      print("Imputing categorical variables using all 2-way associations.")
      #create corresponding margins argument
      margins <- NULL
      for (i in 1:(nCat-1)) {
        for (j in (i+1):nCat) {
          margins <- c(margins, i,j,0)
        }
      }
      margins <- head(margins, -1)
    } else if (marginsType==2) {
      if (nCat<3) {
        stop("You cannot impute with 3-way associations with fewer than 3 categorical variables.")
      } else {
        print("Imputing categorical variables using all three-way interactions.")
        #create corresponding margins argument
        margins <- NULL
        for (i in 1:(nCat-2)) {
          for (j in (i+1):(nCat-1)) {
            for (k in (j+1):nCat) {
              margins <- c(margins, i,j,k,0)
            }
          }
        }
        margins <- head(margins, -1)

      }
    } else {
      print("Imputing categorical variables using saturated log-linear model")
      margins <- 1:nCat
    }

    #figure out design specification
    if (is.null(design)==FALSE) {
      print("Imputing using specified design argument")
    } else if (designType==1) {
      print("Imputing continuous variables assuming main effects of categorical variables")
      #create corresponding design argument
      #create dummy data with one row for each combination of categorical variables
      dummyData <- data.frame(expand.grid(mySeq(s$d)))
      for (i in 1:ncol(dummyData)) {
        dummyData[,i] <- as.factor(dummyData[,i])
      }
      #create formula of main effects of each variable
      dummyFormula <- "~Var1"
      for (i in 2:ncol(dummyData)) {
        dummyFormula <- paste(dummyFormula, "+Var", i, sep="")
      }
      #use model.matrix to create design matrix
      design <- model.matrix(as.formula(dummyFormula), data=dummyData)

    } else {
      print("Imputing continuous variables assuming a separate mean for each combination of categorical variable values")
      design <- diag(s$ncells)
    }
  }

  #find MLE
  if (restrictedMod==1) {
    thetahat <- mix::ecm.mix(s, margins=margins, design=design)
  } else {
    thetahat <- mix::em.mix(s)
  }

  if (pd==FALSE) {
    #MLMI
    for (i in 1:M) {
      imps[[i]] <- as.data.frame(mix::imp.mix(s,thetahat,as.matrix(obsData)))
    }
  } else {
    #PDMI
    for (i in 1:M) {
      if (restrictedMod==1) {
        theta <- mix::dabipf.mix(s,margins=margins,design=design,start=thetahat,steps=steps)
      } else {
        theta <- mix::da.mix(s,thetahat,steps=steps)
      }
      imps[[i]] <- as.data.frame(mix::imp.mix(s,theta))
    }
  }

  if (M==1) {
    imps <- data.frame(imps[[1]])
  }
  attr(imps, "pd") <- pd

  imps
}

#vectorised version of seq
mySeq <- function(arg) {
  myList <- vector("list", length(arg))
  for (i in 1:length(arg)) {
    myList[[i]] <- 1:arg[i]
  }
  myList
}
