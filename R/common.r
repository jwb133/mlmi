#shrinkage functions

#shrinkage function h(,)
h <- function(gamma,nu) {
  (nu/2)*gamma*gsl::gamma_inc((nu-2)/2, nu*gamma/2) / gsl::gamma_inc(nu/2, nu*gamma/2)
}

#matrix shrinkage function H(,)
H <- function(gamma,nu) {
  eigenValsVec <- base::eigen(gamma)
  Q <- eigenValsVec$vectors
  DeltaTilde <- diag(h(eigenValsVec$values, nu))
  Q %*% DeltaTilde %*% solve(Q)
}