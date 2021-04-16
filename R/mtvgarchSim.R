mtvgarchSim <- function (n, m = 2, order.g = c(1,1), order.h = c(1,1,0, 1,1,0), order.x = NULL, 
                         intercept.g = c(1.2,1), size = c(0.8,1), speed = c(1.5,2), 
                         location = c(0.5,0.8), intercept.h = c(0.2,0.3), arch = c(0.10,0.05), 
                         garch = c(0.80,0.90), asym = NULL, xtv = NULL, xreg = NULL, par.xreg = NULL, 
                         R = c(1,0.6,0.6,1), dcc = FALSE, par.dcc = NULL, 
                         opt = 2, verbose = FALSE, innovations = NULL)
{
  if(m == 1) stop("For univariate processes use tvgarchSim.")
  names.ID <- paste("y", 1:m, sep = "")
  if (is.null(innovations)) innovations <- matrix(rnorm(n * m), n, m, byrow=T)
  if (!is.null(order.g)) {
    order.g <- as.matrix(order.g)
    max.s <- ncol(order.g)
    if (sum(order.g) != length(location)) stop("Mismatch between the number of location parameters and the order of the long-term component.")
    names.g <- c("intercept", paste("size", 1:max.s, sep = ""), paste("speed", 1:max.s, sep = ""))
    for (j in 1:max.s) {
      names.g <- c(names.g, paste("location", j, 1:max(order.g[,j]), sep = ""))
    }
    par.g <- matrix(NA, m, 1+2*max.s+max(rowSums(order.g)))
    colnames(par.g) <- names.g
    rownames(par.g) <- names.ID
    intercept.g <- matrix(intercept.g, m, 1)
    size <- matrix(size, m, max.s, byrow = TRUE)
    speed <- matrix(speed, m, max.s, byrow = TRUE)
    location <- matrix(location, m, max(rowSums(order.g)), byrow = TRUE)
  }
  if (is.null(order.g)) par.g <- NULL
  order.h <- matrix(order.h, m, 3, byrow = TRUE)
  max.p <- max(order.h[,1])
  max.q <- max(order.h[,2])
  max.r <- max(order.h[,3])
  par.h <- matrix(NA, m, 1+max.q)
  names.h <- c("intercept", paste("arch", 1:max.q, sep = ""))
  intercept.h <- matrix(intercept.h, m, 1)
  arch <- matrix(arch, m, max.q, byrow = TRUE)
  if (max.p != 0) {
    par.h <- cbind(par.h, matrix(NA, m, max.p))
    names.h <- c(names.h, paste("garch", 1:max.p, sep = ""))
    garch <- matrix(garch, m, max.p, byrow = TRUE)
  }
  if (max.r != 0) {
    par.h <- cbind(par.h, matrix(NA, m, max.r))
    names.h <- c(names.h, paste("asym", 1:max.r, sep = ""))  
    asym <- matrix(asym, m, max.r, byrow = TRUE)
  }
  if (!is.null(order.x)) {
    if (is.null(xreg)) {
      par.h <- cbind(par.h, matrix(NA, m, m))
      par.xreg <- matrix(par.xreg, m, m, byrow = TRUE)
      order.x <- matrix(order.x, m, m, byrow = TRUE)
      diag(order.x) <- 0
      names.h <- c(names.h, paste("xreg", 1:m, sep = ""))  
    } 
    if (!is.null(xreg)) {
      par.h <- cbind(par.h, matrix(NA, m, ncol(xreg)))
      if (length(order.x) != m*ncol(xreg)) warning("The number of par.xreg should be equal to m*ncol(xreg).")
      par.xreg <- matrix(par.xreg, m, ncol(xreg), byrow = TRUE)
      order.x <- matrix(order.x, m, ncol(xreg), byrow = TRUE)
      names.h <- c(names.h, paste("xreg", 1:ncol(xreg), sep = ""))  
    }
  }
  colnames(par.h) <- names.h
  rownames(par.h) <- names.ID
  '
    Constructing correlated innovations
  '
  R <- matrix(R, m, m, byrow = TRUE)
  if (sum( diag( R ) ) != m || any(R > 1) || any(R < (-1)))  stop("R is not a valid correlation matrix.")
  ID <- 1:m
  namesR <- numeric(0)
  for (i in 1:m) {          
    for (j in 1:m) {    
      namesR <- c(namesR, paste("r", paste(paste(ID[j]), ID[i], sep = ""), sep = ""))
    }
  }
  namesR <- matrix(namesR, m, m)
  namesR <- namesR[lower.tri(namesR)]
  '
    Constructing constant conditional correlations (cc)
  '
  if (dcc == FALSE) { 
    chol_P <- chol(R)
    innovations <- innovations %*% chol_P  
    ccc <- matrix(R[lower.tri(R)], n, m*(m-1)/2)
    colnames(ccc) <- namesR
  }
  '
    Constructing dynamic conditional correlations (cc)
  '
  if (dcc == TRUE) {   
    dcc <- matrix(0, n, m*(m-1)/2)
    Qbar <- (1 - par.dcc[1] - par.dcc[2]) * R
    Qt <- R  
    zt2 <- R
    for (t in 1:n) {
      Qt <- Qbar + par.dcc[1] * zt2 + par.dcc[2] * Qt
      invQt <- diag( 1 / sqrt( diag( Qt ) ) )
      Rt  <- invQt %*% Qt %*% invQt
      dcc[t,] <- Rt[lower.tri(Rt)]
      cholRt <- t(chol(Rt))
      zt <- drop(cholRt %*% innovations[t,])
      zt2 <- zt %o% zt
      innovations[t,] <- zt
    }
    colnames(dcc) <- namesR
  }
  if (is.null(colnames(innovations))) colnames(innovations) <- paste("innovations", 1:m, sep = "")
  '
    Constructing the g component and standardising residuals by sqrt(g)
  '
  if (!is.null(order.g)) { 
    if (is.null(xtv)) xtv <- (1:n) / n
    g <- matrix(1, n, m)
    for (i in 1:m) {
      if (order.g[i,1] != 0) {
        s.i <- length(which(order.g[i,] > 0))
        if (!is.null(intercept.g)) intercept.g.i <- intercept.g[i]
        if (!is.null(size)) size.i <- size[i,1:s.i]
        if (!is.null(speed)) speed.i <- speed[i,1:s.i]
        if (!is.null(location)) location.i <- location[i,1:sum(order.g[i,])]
        G <- matrix(1, n, s.i)
        for (j in 1:s.i) {
          G[,j] <- tv(speed = speed.i[j], location = location.i[(sum(order.g[i,1:j])-order.g[i,j]+1):sum(order.g[i,1:j])], xtv = xtv, opt = opt, order.g = order.g[i,j])
        }
        g[,i] <- intercept.g.i[1] + G %*% size.i[1:s.i]
        par.g[i,1] <- intercept.g.i
        par.g[i,(1+1):(1+length(order.g[i,]))] <- size.i
        par.g[i,(1+length(order.g[i,])+1):(1+2*length(order.g[i,]))] <- speed.i
        par.g[i,(1+2*length(order.g[i,])+1):(1+2*length(order.g[i,])+sum(order.g[i,]))] <- location.i
      }
    }
    colnames(g) <- paste("g", 1:m, sep = "")
  }
  if (is.null(order.g)) g <- matrix(1, n, m)
  '
    Constructing the h component
  '
  spillovers <- FALSE
  if (!is.null(order.x) && is.null(xreg)) {
    spillovers <- TRUE
    xreg <- rbind(colMeans(innovations^2), (innovations^2)[-n,])
    colnames(xreg) <- colnames(innovations)
  }
  h <- matrix(1, n, m)
  if (is.null(order.x)) xreg.i <- NULL
  for (i in 1:m){
    if (!is.null(intercept.h)) intercept.h.i <- intercept.h[i]
    if (!is.null(arch)) arch.i <- arch[i,1:order.h[i,2]]
    if (order.h[i,1] != 0) {
      if (!is.null(garch)) garch.i <- garch[i,1:order.h[i,1]]
    }
    if (order.h[i,1] == 0) garch.i <- NULL
    if (order.h[i,3] != 0) {
      if (!is.null(asym)) asym.i <- asym[i,1:order.h[i,3]]
    }
    if (order.h[i,3] == 0) asym.i <- NULL
    if (!is.null(order.x)) {
      if (any(order.x[i,] != 0)) xreg.i <- as.matrix(xreg[,which(order.x[i,] == 1)]) %*% par.xreg[i,which(order.x[i,] == 1)]
      if (sum(order.x[i,] == 0)) xreg.i <- NULL
    }
    h[,i] <- garchxSim(n = n, intercept = intercept.h.i, arch = arch.i, garch = garch.i, asym = asym.i, xreg = xreg.i, innovations = innovations[,i], as.zoo = FALSE, verbose = TRUE)[,"sigma2"]
    par.h[i,1] <- intercept.h.i
    par.h[i,(1+1):(1+order.h[i,2])] <- arch.i
    if (order.h[i,1] != 0) par.h[i,(1+max.q+1):(1+max.q+order.h[i,1])] <- garch.i
    if (order.h[i,3] != 0) par.h[i,(1+max.q+max.p+1):(1+max.q+max.p+order.h[i,3])] <- asym.i
    if (!is.null(order.x)) par.h[i,1+max.q+max.p+max.r+which(order.x[i,] == 1)] <- par.xreg[i,which(order.x[i,] == 1)]
  }
  colnames(h) <- paste("h", 1:m, sep = "")
  '
    Output
  '
  sigma2 <- h * g
  colnames(sigma2) <- paste("sigma2", 1:m, sep = "")
  y <- sqrt(sigma2) * innovations
  colnames(y) <- names.ID
  if (verbose == TRUE) {
    results <- list(y = zoo(y), innovations = zoo(innovations), sigma2 = zoo(sigma2), spillovers = spillovers,
                    g = zoo(g), order.g = order.g, par.g = par.g, 
                    h = zoo(h), order.h = order.h, par.h = par.h, xreg = xreg)
    if (m > 1) {
      results$ccc <- zoo(ccc)
      if(!is.null(par.dcc)) results$dcc <- zoo(dcc)
    }
    return(results)
  }
  if(verbose == FALSE) return(zoo(y))
}
