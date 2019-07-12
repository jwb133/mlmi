#shrinkage functions

#shrinkage function h(,)
h <- function(gamma,nu) {
  (nu/2)*gamma*exp(gsl::lngamma((nu-2)/2) + log(gsl::gamma_inc_Q((nu-2)/2, nu*gamma/2)) -
    (gsl::lngamma(nu/2) + log(gsl::gamma_inc_Q(nu/2, nu*gamma/2))))
}

#matrix shrinkage function H(,)
H <- function(gamma,nu) {
  eigenValsVec <- base::eigen(gamma)
  #print(eigenValsVec)
  Q <- eigenValsVec$vectors
  DeltaTilde <- diag(h(eigenValsVec$values, nu))
  Q %*% DeltaTilde %*% solve(Q)
}
