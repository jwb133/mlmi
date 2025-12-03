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

#check if a vector takes consectutive integers starting from 1
is_pos_consecutive <- function(x) {
  #remove any NAs
  x <- x[!is.na(x)]
  #unique values taken
  ux <- sort(unique(x))
  if (length(ux)==max(x)) {
    #values should be consecutive integers
    all(ux == seq_len(max(ux)))
  } else {
    FALSE
  }
}
