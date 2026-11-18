cv.kernreg <- function(y, x, rads = TRUE, h = seq(0.1, 1, length = 10), type = "gauss",
                         nfolds = 10, folds = NULL, graph = FALSE, ncores = 1) {
  if ( is.null(dim(y)) ) {
    if ( !rads )  y <- y / 180 * pi
    y <- cbind( cos(y), sin(y) )
  } else y <- as.matrix(y)
  n <- dim(y)[1]
  runtime <- proc.time()
  ina <- 1:n
  if ( is.null(folds) )  folds <- .makefolds(ina, nfolds = nfolds, stratified = FALSE)
  nfolds <- length(folds)
  msp <- matrix( nrow = length(folds), ncol = length(h) )
  lh <- length(h)
  y <- as.matrix(y)   ;   x <- as.matrix(x)

  for ( i in 1:nfolds ) {
    nu <- folds[[ i ]]
    ytrain <- y[-nu, , drop = FALSE]
    ytest <- y[nu, , drop = FALSE]
    xtest <- x[nu, , drop = FALSE]
    xtrain <- x[-nu, , drop = FALSE]
    est <- Directionalkernel::kernreg(xtest, ytrain, xtrain, rads = TRUE, h, type = type, ncores = ncores)
    for ( j in 1:lh )  msp[i, j] <- 0.5 - 0.5 * mean( Rfast::rowsums(est[[ j ]] * ytest) )
  }

  runtime <- proc.time() - runtime
  mspe <- Rfast::colmeans(msp)
  names(mspe) <- paste("h=", h, sep = "")
  if (graph) {
    plot(h, mspe, xlab = "Bandwidth parameter (h)", ylab = "MSPE", cex.lab = 1.2, cex.axis = 1.2)
    abline(v = h, col = "lightgrey", lty = 2)
    abline(h = seq(min(mspe), max(mspe), length = 10), col = "lightgrey", lty = 2)
    points(h, mspe, pch = 19)
    lines(h, mspe, lwd = 2)
  }
  list(mspe = mspe, h = h[ which.min(mspe) ], performance = min(mspe), runtime = runtime)
}





.makefolds <- function(ina, nfolds = 10, stratified = TRUE) {
  n <- length(ina)
  if ( stratified ) {
    ## shuffle within each class, then deal the observations to the folds in turn
    id <- lapply( unique(ina), function(lab) {
      w <- which(ina == lab)
      w[ sample.int( length(w) ) ]
    } )
    id <- unlist(id)
  } else  id <- sample.int(n)
  runs <- split( id, rep_len(1:nfolds, n) )
  names(runs) <- paste("Fold", seq_along(runs))
  runs
}
