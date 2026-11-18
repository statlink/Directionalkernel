vm.kernel <- function(y, h = "tay", rads = TRUE, ncores = 1) {
  ## y is the data
  ## h is "tay" (Taylor, 2008), "rot" (Garcia-Portugues, 2013) or a number
  ## if the data are in degrees we transform them into radians
  if ( !rads )  y <- y / 180 * pi
  n <- length(y)  ## sample size
  if ( is.character(h) )  h <- .vm.bw(y, h)
  ## exponentially scaled, vm_kernel_cpp() works with exp( h * (cos - 1) )
  com <- (n - 1) * 2 * pi * besselI(h, 0, expon.scaled = TRUE)

  f <- vm_kernel_cpp( rbind( cos(y), sin(y) ), h, com, as.integer(ncores) )

  list( h = h, f = f )
}


## u: angles in radians. h: "tay" (Taylor, 2008) or "rot" (Garcia-Portugues, 2013).
## The value returned is the concentration of the von Mises kernel.
.vm.bw <- function(u, h) {
  if ( !h %in% c("tay", "rot") )  stop( 'h must be "tay", "rot" or numeric' )
  n <- length(u)
  est <- Rfast::rowmeans( rbind( cos(u), sin(u) ) )
  Rbar <- sqrt( sum(est^2) )  ## mean resultant length
  lik <- function(k)  k * Rbar - log( besselI(k, 0, expon.scaled = TRUE) ) - k
  k <- optimize(lik, c(0, 50000), maximum = TRUE, tol = 1e-6)$maximum
  i0 <- besselI(k, 0, expon.scaled = TRUE)^2
  i2 <- besselI(2 * k, 2, expon.scaled = TRUE)
  if ( h == "tay" ) {
    res <- 3 * n * k^2 * i2 / ( 4 * pi^0.5 * i0 )
  } else {
    i1 <- besselI(2 * k, 1, expon.scaled = TRUE)
    res <- k * n * ( 2 * i1 + 3 * k * i2 ) / ( 4 * pi^0.5 * i0 )
  }
  res^0.4
}
