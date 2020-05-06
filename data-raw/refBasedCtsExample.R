#take a look at ctsTrialWide data
head(ctsTrialWide)

#impute the missing outcome values twice assuming MAR
imps <- refBasedCts(ctsTrialWide, outcomeVarStem="y", nVisits=3, trtVar="trt",
                    baselineVars=c("v", "y0"), type="MAR", M=2)

#now impute using jump to reference method
imps <- refBasedCts(ctsTrialWide, outcomeVarStem="y", nVisits=3, trtVar="trt",
                    baselineVars=c("v", "y0"), type="J2R", M=2)

#for frequentist valid inferences we use bootstrapping from the bootImpute package
\dontrun{
  #bootstrap 10 times using 2 imputations per bootstrap. Note that to do this
  #we specify nImp=2 to bootImpute by M=1 to the refBasedCts function.
  #Also, 10 bootstraps is far too small to get reliable inferences. To do this
  #for real you would want to use a lot more (e.g. at least nBoot=1000).
  library(bootImpute)
  bootImps <- bootImpute(ctsTrialWide, refBasedCts, nBoot=10, nImp=2,
                         outcomeVarStem="y", nVisits=3, trtVar="trt",
                         baselineVars=c("v", "y0"), type="J2R", M=1)

  #write a small wrapper function to perform an ANCOVA at the final time point
  ancova <- function(inputData) {
    coef(lm(y3~v+y0+trt, data=inputData))
  }
  ests <- bootImputeAnalyse(bootImps, ancova)
  ests
}
