tvgarch <- function (y, order.g = 1, order.h = c(1,1,0), xtv = NULL, xreg = NULL, 
                     initial.values = list(), opt = 2, turbo = FALSE, trace = FALSE)
{
  n <- length(y)
  if (order.g == 0) order.g <- NULL
  npar.h <- 1 + sum(order.h)
  if (!is.null(xreg)) {
    xreg <- as.matrix(xreg)
    npar.h <- npar.h + ncol(xreg)
  } 
  robse.par.h <- numeric(npar.h)
  maxpqr <- max(order.h)
  if (!is.null(order.g)) {
    if (!is.null(xtv)){
      xtv <- matrix(xtv, n, 1)
      if (is.null(colnames(xtv))) colnames(xtv) <- "xtv"
    }
    if (is.null(xtv)) {
      xtv <- matrix((1:n)/n, n, 1)
      colnames(xtv) <- "time"
    }
    s <- length(order.g)
    npar.g <- 1 + 2 * s + sum(order.g)
  }
  '
    Initial values
  '
  if (!is.null(order.g)) {
    if (is.null(initial.values$intercept.g)) intercept.g <- 1
    else intercept.g <- initial.values$intercept.g
    if (is.null(initial.values$size)) size <- rep(0.1, s)
    else size <- initial.values$size
    if (is.null(initial.values$speed)) speed <- rep(1, s)
    else speed <- initial.values$speed
    if (is.null(initial.values$location)) {
      location <- NULL
      for (i in 1:s) {
        if (order.g[i] == 1) location <- c(location, mean(xtv))
        if (order.g[i] > 1) location <- c(location, seq(from = min(xtv)+0.5*sd(xtv), to = max(xtv)-0.5*sd(xtv), by = (max(xtv)-min(xtv)-sd(xtv))/(order.g[i]-1)))
      }
    }
    else location <- initial.values$location
    if (s != length(size)) stop("Mismatch between the number of transition functions, s, and the number of initial size values.")
    if (s != length(speed)) stop("Mismatch between the number of transition functions, s, and the number of initial speed values.")
    if (sum(order.g) != length(location)) stop("Mismatch between order.gand the initial location values.")
  }
  if (is.null(initial.values$intercept.h)) intercept.h <- 0.1
  if (!is.null(initial.values$intercept.h))  intercept.h <- initial.values$intercept.h
  if (is.null(initial.values$arch)) arch <- rep(0.1/order.h[2], order.h[2])
  if (!is.null(initial.values$arch))  arch <- initial.values$arch
  if (order.h[2] != length(arch)) stop("Mismatch between the number of ARCH-type parameters, q, and the arch initial values.")
  if (order.h[1] != 0) {
    if (is.null(initial.values$garch)) garch <- rep(0.7/order.h[1], order.h[1])
    else garch <- initial.values$garch
    if (order.h[1] != length(garch)) stop("Mismatch between the number of GARCH-type parameters, p, and the garch initial values.")
  }
  if (order.h[1] == 0) garch <- NULL
  if (order.h[3] != 0) {
    if (is.null(initial.values$asym)) asym <- rep(0.02/order.h[3], order.h[3])
    if (!is.null(initial.values$asym)) asym <- initial.values$asym
    if (order.h[3] != length(asym)) stop("Mismatch between the number of GJR-type parameters, r, and the asym initial values.")
  }
  if (order.h[3] == 0)  asym <- NULL
  if (!is.null(xreg)) {
    if (is.null(initial.values$par.xreg)) par.xreg <- rep(0.01, ncol(xreg))
    if (!is.null(initial.values$par.xreg)) par.xreg <- initial.values$par.xreg
    if (ncol(xreg) != length(par.xreg)) stop("Mismatch between the number of covariates, X, and the par.xreg initial values.")
  }
  if (is.null(xreg)) par.xreg <- NULL
  '
    Vector of initial parameters for h component
  '
  names.h <- c("intercept.h", paste("arch", paste(seq(1:order.h[2]), sep = ""), sep = ""))
  ini.par.h <- c(intercept.h, arch)
  if (order.h[1] != 0) {
    names.h <- c(names.h, paste("garch", paste(seq(1:order.h[1]), sep = ""), sep = "")) 
    ini.par.h <- c(ini.par.h, garch)
  }
  if (order.h[3] != 0) {
    names.h <- c(names.h, paste("asym", paste(seq(1:order.h[3]), sep = ""), sep = "")) 
    ini.par.h <- c(ini.par.h, asym)
  }
  if (!is.null(xreg)) { 
    if(!is.null(colnames(xreg))) names.h <- c(names.h, colnames(xreg))
    else names.h <- c(names.h, paste("xreg", paste(seq(1:ncol(xreg)), sep = ""), sep = ""))
    ini.par.h <- c(ini.par.h, par.xreg)
  }
  if (length(ini.par.h) != npar.h) stop("Check initial parameters in the h component.")
  initial.h <- matrix(ini.par.h, 1, npar.h)
  colnames(initial.h) <- names.h
  rownames(initial.h) <- "Value:"
  if (trace == TRUE) {
    cat("\n Initial parameters in the h component\n\n")
    print(round(initial.h, 4))
  }
  '
   g component and maximization by parts
  '
  g <- y
  g[1:n] <- 1
  if (!is.null(order.g)) {
    '
      Vector of initial parameters for g component
    ' 
    ini.par.g <- c(intercept.g, size, speed, location)
    names.g <- c("intercept.g", paste("size", paste(seq(1:s), sep = ""), sep = ""), paste("speed", paste(seq(1:s), sep = ""), sep = ""))
    for (i in 1:s) {
      if (s == 1) names.g <- c(names.g, paste("location", paste(seq(1:max(order.g[i])), sep = ""), sep = ""))
      else names.g <- c(names.g, paste("location", paste(i, sep = ""), paste(seq(1:max(order.g[i])), sep = ""), sep = ""))
    }
    if (length(ini.par.g) != npar.g) stop("Check initial parameters in the g component.")
    initial.g <- matrix(ini.par.g, 1, npar.g)
    colnames(initial.g) <- names.g
    rownames(initial.g) <- "Value:"
    if (trace == TRUE) {
      cat("\n Initial parameters in the g component:\n\n")
      print(round(initial.g, 4))
    }
    '
      Parameter constraints for g component and iter = 0
    '
    r <- diag(npar.g)
    r1 <- r[1,]
    ui.g <- r1
    ci.g <- 1e-5
    if (s == 1) {
      r2 <- r[2,]
      r2[1] <- 1
      ui.g <- rbind(ui.g, r2)
      ci.g <- c(ci.g, 1e-5)
    } 
    if (s > 1) {
      r.size <- combos(s)$binary
      r3 <- matrix(0, nrow(r.size), npar.g)
      r3[,1] <- 1
      r3[,2:(s+1)] <- r.size
      ui.g <- rbind(ui.g, r3)
      ci.g <- c(ci.g, rep(1e-5, nrow(r.size))) 
    }
    r4 <- r[(2+s):(2*s+1),]
    r5 <- r[(2*s+2):npar.g,]
    ui.g <- rbind(ui.g, r4, -r4, r5, -r5)
    if (opt != 2) ci.g <- c(ci.g, rep(1e-5, s), rep(-250, s), rep(min(xtv)+1e-3, sum(order.g)), rep(-max(xtv)+1e-3, sum(order.g))) 
    if (opt == 2) ci.g <- c(ci.g, rep(log(1e-5), s), rep(-log(250), s), rep(min(xtv)+1e-3, sum(order.g)), rep(-max(xtv)+1e-3, sum(order.g))) 
    if (any(order.g > 1)) {
      for (i in which(order.g > 1)) {
        r6 <- -r[(1+2*s+sum(order.g[1:i])-order.g[i]+1):(1+2*s+sum(order.g[1:i])-1),] + r[(1+2*s+sum(order.g[1:i])-order.g[i]+2):(1+2*s+sum(order.g[1:i])),]
        ui.g <- rbind(ui.g, r6)
        ci.g <- c(ci.g, rep(0, order.g[i]-1)) 
      }
    }
    '
      Maximization by parts: estimating initial TV parameters
    '
    iter <- 0
    h <- rep(1, n)
    iter0.fit.g <- constrOptim(theta = ini.par.g, f = tvObj, grad = NULL, ui = ui.g, ci = ci.g, y = y, order.g = order.g, xtv = xtv, opt = opt, fixed.h = h, iter0 = TRUE, flag = 1)
    par.hat0.g <- iter0.fit.g$par
    g[1:n] <- tvObj(par.g = iter0.fit.g$par, fixed.par.g = NULL, xtv = xtv, opt = opt, order.g = order.g, fixed.h = h, y = y, iter0 = TRUE, flag = 2)
    estimates <- matrix(par.hat0.g, 1, length(names.g))
    colnames(estimates) <- names.g
    rownames(estimates) <- "Value:"
    if (trace == TRUE) {
      cat("\n Initial estimation for the g component:\n\n")
      if (turbo == TRUE) print(round(estimates[1,], 4))
      else print(round(estimates, 4))
    }
    '
      Parameter constraints for g component and iter > 0
    '
    par.hat.g <- par.hat0.g[-c(1,(s+2):(2*s+1))]
    par.hat0.g <- par.hat0.g[c(1,(s+2):(2*s+1))]
    r <- diag(length(par.hat.g))
    ui.g <- rbind(r[(s+1):(npar.g-s-1),], -r[(s+1):(npar.g-s-1),])
    ci.g <- c(rep(min(xtv)+1e-3, sum(order.g)), rep(-max(xtv)+1e-3, sum(order.g))) 
    if (s == 1) {
      r1 <- r[1,]
      ui.g <- rbind(ui.g, r1)
      ci.g <- c(ci.g, -par.hat0.g[1]+1e-5)
    } 
    if (s > 1) {
      r.size <- combos(s)$binary
      r1 <- matrix(0, nrow(r.size), length(par.hat.g))
      r1[,1:s] <- r.size
      ui.g <- rbind(ui.g, r1)
      ci.g <- c(ci.g, rep(-par.hat0.g[1]+1e-5, nrow(r.size))) 
    }
    if (any(order.g > 1)) {
      for(i in which(order.g > 1)){
        r2 <- -r[(s+sum(order.g[1:i])-order.g[i]+1):(s+sum(order.g[1:i])-1),] + r[(s+sum(order.g[1:i])-order.g[i]+2):(s+sum(order.g[1:i])),]
        ui.g <- rbind(ui.g, r2)
        ci.g <- c(ci.g, rep(0, order.g[i]-1)) 
      }
    }
    '
      Iterative estimation
    '
    iter <- iter + 1
    maxiter <- 1000 
    phi <- y/sqrt(g)
    par.hat.h <- ini.par.h
    repeat{
      if (trace == TRUE) {
        cat(" \n\n")
        cat("Iteration round:",iter, "\n\n")
      }
      '
        Estimating GJR-GARCH-X parameters:
      '
      iter.fit.h <- garchx(y = phi, order = order.h, xreg = xreg, initial.values = par.hat.h)
      h <- garchxRecursion(pars = as.numeric(iter.fit.h$par), aux = iter.fit.h)
      if (iter > 1) conv.h <- max(abs(par.hat.h-iter.fit.h$par))
      par.hat.h <- iter.fit.h$par
      if (trace == TRUE) {
        cat("Estimates (h component):", par.hat.h, "\n")
        cat("Likelihood:", logLik.garchx(iter.fit.h), "\n")
      }
      '
        Estimating TV parameters:
      '
      iter.fit.g <- constrOptim(theta = par.hat.g, f = tvObj, grad = NULL, ui = ui.g, ci = ci.g, fixed.par.g = par.hat0.g, y = y, order.g = order.g, xtv = xtv, opt = opt, fixed.h = h, iter0 = FALSE, flag = 1)
      conv.g <- max(abs(par.hat.g - iter.fit.g$par))
      par.hat.g <- iter.fit.g$par
      g[1:n] <- tvObj(par.g = par.hat.g, fixed.par.g = par.hat0.g, xtv = xtv, opt = opt, order.g = order.g, fixed.h = h, y = y, iter0 = FALSE, flag = 2)
      phi <- y/sqrt(g)
      if (trace == TRUE) {
        cat("Estimates (g component):", par.hat.g, "\n")
        cat("Likelihood:", - iter.fit.g$value, "\n")
        if (iter > 1) cat("Max. convergence:", max(conv.h, conv.g), "\n")
      }
      '
        Next round
      '
      iter <- iter + 1
      if (iter > 2) {
        if (iter >= maxiter || max(conv.g, conv.h) < 1e-5){
          if (iter == maxiter) warning("Convergence has not been reached after ", paste(maxiter, sep = "")," iterations.")
          break
        }
      }
    }
    par.hat.g <- c(par.hat0.g[1], iter.fit.g$par[1:s], par.hat0.g[-1], iter.fit.g$par[-(1:s)])
    if (turbo == TRUE) {
      vcov.g <- NULL
      robse.par.g <- rep(NA, length(iter.fit.g$par))
    } 
    if (turbo == FALSE) {
      '
          TV variance-covariance matrix:
      '
      vcov.g <- matrix(NA, npar.g, npar.g)
      rownames(vcov.g) <- names.g
      colnames(vcov.g) <- names.g
      jac.g <- jacobian(func = tvObj, x = iter.fit.g$par, fixed.par.g = par.hat0.g, xtv = xtv, opt = opt, order.g = order.g, fixed.h = h, y = y, iter0 = FALSE, flag = 0)
      J.g <- crossprod(jac.g)  
      H.g <-  optimHess(par = iter.fit.g$par, fn = tvObj, fixed.par.g = par.hat0.g, xtv = xtv, opt = opt, order.g = order.g, fixed.h = h, y = y, iter0 = FALSE, flag = 1)
      vcov.iter <- solve(-H.g) %*% J.g %*% solve(-H.g)
      vcov.g[c(2:(s+1),(2*s+2):npar.g),c(2:(s+1),(2*s+2):npar.g)] <- vcov.iter
      robse.par.g <- sqrt(diag(vcov.g))
    }
    '
      TV final parameter estimates:
    '
    estimates.g <- as.matrix(rbind(par.hat.g, robse.par.g))
    rownames(estimates.g) <- c("Estimate:", "Std. Error:")
    colnames(estimates.g) <- names.g
  }
  '
    Final GJR-GARCH-X estimation
  '
  iter.fit.h <- garchx(y = y/sqrt(g), order = order.h, xreg = xreg, initial.values = ini.par.h)
  h <- garchxRecursion(pars = as.numeric(iter.fit.h$par), aux = iter.fit.h)
  if (turbo == TRUE) {
    vcov.h <- NULL
    robse.par.h <- rep(NA, length(iter.fit.h$par))
  }
  if (turbo == FALSE) {
    '
      GJR-GARCH-X variance-covariance matrix:
    '
    vcov.h <- vcov.garchx(iter.fit.h)
    robse.par.h <- sqrt(diag(vcov.h))
  }
  '
    GJR-GARCH-X final parameter estimates:
  '
  estimates.h <- as.matrix(rbind(iter.fit.h$par, robse.par.h))
  rownames(estimates.h) <- c("Estimate:", "Std. Error:")
  colnames(estimates.h) <- names.h
  logLik <- logLik.garchx(iter.fit.h)
  '
    TV-GJR-GARCH-X output:
  '
  if (!is.null(order.g)) {
    sigma2 <- h*g
    residuals <- y / sqrt(sigma2)
    results <- list(par.g = estimates.g[1,], se.g = estimates.g[2,], par.h = estimates.h[1,], se.h = estimates.h[2,], names.g = names.g, names.h = names.h, sigma2 = sigma2, 
                    residuals = residuals, h = h, g = g, logLik = logLik, vcov.g = vcov.g, vcov.h = vcov.h, message.g = iter.fit.g$message, message.h = iter.fit.h$message,
                    order.g = order.g, order.h = order.h, xtv = xtv, xreg = xreg, opt = opt, y = y, date = date(), iter = iter, turbo = turbo)
    class(results) <- "tvgarch"
  }
  if (is.null(order.g)){
    results <- iter.fit.h
    results$order.g <- order.g
  }
  return(results)
} 



