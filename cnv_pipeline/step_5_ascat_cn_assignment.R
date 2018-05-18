
library(nloptr)

# play around with Ascat:
get_r <- function(rho, psi, n) {
  log2((2 * (1-rho) + rho * n) / psi)
}

compare <- function(observed, ideal) {
  sq <- (observed - ideal)^2
  sum(sq)
}

objective <- function(x) {
  rho <- x[1]
  psi <- x[2]
  n <- x[3]
  ideal <- get_r(rho, psi, n)
  compare(observed, ideal)
}
