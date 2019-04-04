withinBetween <- function(imps, analysisfun, dfComplete=100000, ...) {
  M <- length(imps)
  pd <- as.logical(attributes(imps)['pd'])

  #run analysis on first imputation to find out length of parameter vector
  result <- analysisfun(imps[[1]],...)
  numParms <- length(result$est)

  #analyse each imputed datasets
  ests <- array(0, dim=c(M,numParms))
  vars <- array(0, dim=c(M,numParms,numParms))
  for (m in 1:M) {
    result <- analysisfun(imps[[m]],...)
    ests[m,] <- result$est
    vars[m,,] <- result$var
  }
  thetaHat <- colMeans(ests)
  What <- apply(vars, c(2,3), mean)
  Bhat <- var(ests)

  if (pd==FALSE) {
    #implement von Hippel's WB variance
    gammaHatMis <- solve(What) %*% Bhat
    gammaTildeMis <- H(gammaHatMis, M-1)
    gammaTildeObs <- diag(length(thetaHat))-gammaTildeMis
    VTildeML_MLMI <- What %*% solve(gammaTildeObs)
    vTildeMLMI <- VTildeML_MLMI + Bhat/M

    #degrees of freedom
    nuTildeML_WB <- (M-1)*(mean(diag(gammaTildeObs))/mean(diag(gammaTildeMis)))^2 - 4
    nuHatMLMI_WB <- VTildeMLMI^2 / (VTildeML_MLMI^2/nuTildeML_WB + (Bhat/M)^2 / (M-1))
    nuTildeObs <- dfComplete*gammaTildeObs*(dfComplete+3)/(dfComplete+1)
    nuTildeMLMI_WB <- max(3,((1/nuHatMLMI_WB) + (1/nuTildeObs))^(-1))

  } else {
    #Rubin's rules
    VhatPDMI_WB <- What + (1+1/M)*Bhat
    totalVar <- VhatPDMI_WB
    gammaTildeMis <- (1+1/M)*sum(diag(Bhat %*% solve(VhatPDMI_WB)))/numParms
    nuHatPDMI_WB <- (M-1)/gammaTildeMis^2
    nuTildeObs <- dfComplete*(1-gammaTildeMis)*(dfComplete+1)/(dfComplete+3)
    nuTildePDMI_WB <- max(3, (1/nuHatPDMI_WB + 1/nuTildeObs)^(-1))
    MIdf <- nuTildePDMI_WB
  }

  list(est=thetaHat, var=totalVar, df=MIdf)
}

#shrinkage function h(,)
h <- function(gamma,nu) {
  (nu/2)*gamma*gsl::gamma_inc((nu-2)/2, nu*gamma/2) / gsl::gamma_inc(nu/2, nu*gamma/2)
}

#matrix shrinkage function H(,)
H <- function(gamma,nu) {
  eigenValsVec <- base::eigen(gamma)$vectors
  Q <- eigenValsVec$vectors
  DeltaTilde <- diag(h(eigenValsVec$values, nu))
  Q %*% DeltaTilde %*% solve(Q)
}