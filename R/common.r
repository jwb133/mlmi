#shrinkage functions

#shrinkage function h(,)
#we calculate log(h,) first, to avoid numerical overflow problems. This avoids calculating
#the two incomplete gamma function values directly, which can be very large
h <- function(gamma,nu) {
  exp(log(nu/2) + log(gamma) + gsl::lngamma((nu-2)/2) + log(gsl::gamma_inc_Q((nu-2)/2, nu*gamma/2)) -
    (gsl::lngamma(nu/2) + log(gsl::gamma_inc_Q(nu/2, nu*gamma/2))))
}

#matrix shrinkage function H(,)
H <- function(gamma,nu) {
  eigenValsVec <- base::eigen(gamma)
  Q <- eigenValsVec$vectors
  DeltaTilde <- diag(h(eigenValsVec$values, nu))
  Q %*% DeltaTilde %*% solve(Q)
}
