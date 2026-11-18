vmf.kernel <- function(y, h = "rot", ncores = 1) {
  ## y is the data
  ## h is "rot" (Garcia-Portugues, 2013) or a number
  p <- dim(y)[2]  ## dimensionality of the data
  n <- dim(y)[1]  ## sample size of the data
  if ( is.character(h) )  h <- .vmf.bw(y, h)
  ## exponentially scaled, vmf_kernel_cpp() works with exp( (y'y - 1) / h^2 )
  cpk <- (1/h^2)^(p/2 - 1) / ( (2 * pi)^(p/2) * besselI(1/h^2, p/2 - 1, expon.scaled = TRUE) )
  com <- cpk / (n - 1)

  f <- vmf_kernel_cpp(t(y), h, com, as.integer(ncores))
  list( h = h, f = f )
}



## y: matrix of unit vectors, one per row. h: "rot" (Garcia-Portugues, 2013).
## The value returned is h, where the concentration of the kernel is 1 / h^2.
.vmf.bw <- function(y, h) {
  if ( h != "rot" )  stop( 'h must be "rot" or numeric' )
  n <- dim(y)[1]
  q <- dim(y)[2] - 1
  k <- Rfast::vmf.mle(y)$kappa
  up <- 4 * pi^0.5 * besselI(k, (q - 1)/2, expon.scaled = TRUE)^2
  down <- k^( (q + 1)/2 ) * n * ( 2 * q * besselI(2 * k, (q + 1)/2, expon.scaled = TRUE) +
                                    (2 + q) * k * besselI(2 * k, (q + 3)/2, expon.scaled = TRUE) )
  ( up / down )^( 1/(4 + q) )
}
