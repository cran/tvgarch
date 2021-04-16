mtvgarch <-  function (y, order.g = c(1, 1), order.h = c(1,1,0, 1,1,0), order.x = NULL,
                       initial.values = list(), xtv = NULL, xreg = NULL, opt = 2, dcc = FALSE, 
                       turbo = TRUE, trace = FALSE)
{
    y <- as.matrix(y)
    n <- nrow(y)
    m <- ncol(y)
    if (is.null(colnames(y)) ) colnames(y) <- paste("y", 1:m, sep = "")
    names.y <- colnames(y)
    if (!is.null(order.g)) {
      if (!is.null(xtv)) {
        xtv <- matrix(xtv, n, 1)
        if (is.null(colnames(xtv))) 
          colnames(xtv) <- "xtv1"
      }
      if (is.null(xtv)) {
        xtv <- as.matrix((1:n)/n)
        colnames(xtv) <- "time"
      }
      order.g <- as.matrix(order.g)
      max.s <- ncol(order.g)
      names.g <- c("intercept.g", paste("size", 1:max.s, sep = ""), paste("speed", 1:max.s, sep = ""))
      for(i in 1:max.s) {
        names.g <- c(names.g, paste("location", i, 1:max(order.g[,i]), sep = ""))
      }
      npar.g <- 1 + 2 * max.s + sum(matrixStats::colMaxs(order.g))
      par.g <- matrix(NA, m, npar.g)
      colnames(par.g) <- names.g
      rownames(par.g) <- names.y
      se.g <- par.g
    }
    if (is.null(order.g)) {
      par.g <- NULL
      se.g <- par.g
    } 
    order.h <- matrix(order.h, m, 3, byrow = TRUE)
    max.p <- max(order.h[,1])
    max.q <- max(order.h[,2])
    max.r <- max(order.h[,3])
    names.h <- c("intercept.h", paste("arch", 1:max.q, sep = ""))
    if (max.p != 0) names.h <- c(names.h, paste("garch", 1:max.p, sep = ""))
    if (max.r != 0) names.h <- c(names.h, paste("asym", 1:max.r, sep = ""))  
    npar.h <- 1 + max.p + max.q + max.r
    if (!is.null(order.x)) {
      diag(order.x) <- 1
      if (sum(order.x) == m) order.x <- NULL
      if (sum(order.x) > m) {
        if (is.null(xreg)) {
          order.x <- matrix(order.x, m, m, byrow = TRUE)
          diag(order.x) <- 0
          y2 <- as.matrix(rbind(colMeans(y^2), as.matrix(y[-n,]^2)))
          names.x <- paste("xreg", 1:m, sep = "")
          colnames(y2) <- names.x 
          colnames(order.x) <- names.y
          rownames(order.x) <- names.y
        } 
        if (!is.null(xreg)) {
          order.x <- matrix(order.x, m, ncol(xreg))
          rownames(order.x) <- names.y
          if (is.null(colnames(xreg))) colnames(xreg) <- paste("xreg", 1:m, sep = "")
          names.x <- colnames(xreg)
          colnames(order.x) <- names.x
        }
        max.x <- ncol(order.x)
        npar.h <- npar.h + max.x
        names.h <- c(names.h, names.x)  
      }
    }
    par.h <- matrix(NA, m, npar.h)
    colnames(par.h) <- names.h
    rownames(par.h) <- names.y
    se.h <- par.h
    if (!is.null(initial.values)) {
      if (!is.null(initial.values$intercept.g)) initial.values$intercept.g <- matrix(initial.values$intercept.g, m, 1)
      if (!is.null(initial.values$size)) initial.values$size <- matrix(initial.values$size, m, max.s, byrow = TRUE)
      if (!is.null(initial.values$speed)) initial.values$speed <- matrix(initial.values$speed, m, max.s, byrow = TRUE)
      if (!is.null(initial.values$location)) initial.values$location <- matrix(initial.values$location, m, max(rowSums(order.g)), byrow = TRUE)
      if (!is.null(initial.values$intercept.h)) initial.values$intercept.h <- matrix(initial.values$intercept.h, m, 1)
      if (!is.null(initial.values$arch)) initial.values$arch <- matrix(initial.values$arch, m, max.q, byrow = TRUE)
      if (!is.null(initial.values$garch)) initial.values$garch <- matrix(initial.values$garch, m, max.p, byrow = TRUE)
      if (!is.null(initial.values$asym)) initial.values$asym <- matrix(initial.values$asym, m, max.r, byrow = TRUE)
      if (!is.null(initial.values$par.xreg)) initial.values$par.xreg <- matrix(initial.values$par.xreg, m, max.x, byrow = TRUE)
    } 
    h <- matrix(1, n, m)
    g <- matrix(1, n, m)
    robse.par.h <- par.h
    robse.par.g <- par.g
    conv.i <- numeric(m)
    iter <- 1
    '
      Estimating volatilities equation by equation 
    '
    nx <- 0
    if (is.null(order.g)) {
      order.g.i <- NULL
      intercept.g.i <- NULL
      size.i <- NULL
      speed.i <- NULL
      location.i <- NULL
    }
    Objs <- list()
    for(i in 1:m) {
      if (!is.null(order.g)) {
        if (order.g[i,1] != 0) {
          s.i <- length(which(order.g[i,] > 0))
          if (is.null(initial.values$intercept.g)) intercept.g.i <- NULL
          else intercept.g.i <- initial.values$intercept.g[i]
          if (is.null(initial.values$size)) size.i <- NULL
          else size.i <- initial.values$size[i,1:s.i]
          if (is.null(initial.values$speed)) speed.i <- NULL
          else speed.i <- initial.values$speed[i,1:s.i]
          if (is.null(initial.values$location)) location.i <- NULL
          else location.i <- initial.values$location[i,1:sum(order.g[i,])]
          par.g.i <- c(intercept.g.i, size.i, speed.i, location.i)
          order.g.i <- order.g[i,1:s.i]
        }
        else{
          order.g.i <- NULL
          intercept.g.i <- NULL
          size.i <- NULL
          speed.i <- NULL
          location.i <- NULL
        } 
      }
      if (is.null(initial.values$intercept.h)) intercept.h.i <- NULL
      else intercept.h.i <- initial.values$intercept.h[i]
      if (is.null(initial.values$arch)) arch.i <- NULL
      else arch.i <- initial.values$arch[i,1:order.h[i,2]]
      par.h.i <- c(intercept.h.i, arch.i)
      if (order.h[i,1] != 0) {
        if (is.null(initial.values$garch)) garch.i <- NULL
        else garch.i <- initial.values$garch[i,1:order.h[i,1]]
        par.h.i <- c(par.h.i, garch.i)
      }
      if (order.h[i,1] == 0) garch.i <- NULL
      if (order.h[i,3] != 0) {
        if (is.null(initial.values$asym)) asym.i <- NULL
        else asym.i <- initial.values$asym[i,1:order.h[i,3]]
        par.h.i <- c(par.h.i, asym.i)
      }
      if (order.h[i,3] == 0) asym.i <- NULL
      if (is.null(order.x)) {
        par.xreg.i <- NULL
        xreg.i <- NULL
      }
      if (!is.null(order.x)) {
        if (sum(order.x[i,]!=0)) {
          if (is.null(initial.values$par.xreg)) par.xreg.i <- NULL
          else par.xreg.i <- initial.values$par.xreg[i,which(order.x[i,] == 1)]
          par.h.i <- c(par.h.i, par.xreg.i)
          if (!is.null(xreg)) {
            xreg.i <- as.matrix(xreg[,which(order.x[i,] == 1)])
            colnames(xreg.i) <- names.x[which(order.x[i,] == 1)]
          }
          if (is.null(xreg)) {
            xreg.i <- as.matrix(y2[,which(order.x[i,] == 1)])
            colnames(xreg.i) <- names.x[which(order.x[i,] == 1)]
          }
        }
        if (sum(order.x[i,]) == 0) {
          par.xreg.i <- NULL
          xreg.i <- NULL
        }
      }
      yi.i <- as.matrix(y[,i])
      colnames(yi.i) <- names.y[i]
      if (!is.null(order.g) && !is.null(order.x)) {
        if (sum(order.x[i,] != 0) && is.null(xreg) && order.g[which(order.x[i,] != 0),1] != 0) { 
          turbo2 <- TRUE
          nx <- nx + 1 
        }
        else turbo2 <- FALSE
      }
      else turbo2 <- turbo
      tvgarch.i <- tvgarch(y = yi.i, order.g = order.g.i, order.h = order.h[i,], xtv = xtv, xreg = xreg.i, 
                           opt = opt, initial.values = list(intercept.g = intercept.g.i, size = size.i, speed = speed.i, location = location.i,
                                                            intercept.h = intercept.h.i, arch = arch.i, garch = garch.i, asym = asym.i, par.xreg = par.xreg.i),  
                           turbo = turbo2, trace = trace)
      Objs[[paste("obj", i, sep = "")]] <- tvgarch.i
      '
        Saving estimated parameters
      '       
      if (!is.null(order.g)) {
        if (order.g[i,1] != 0) {
          par.g[i,1] <- tvgarch.i$par.g[1]
          par.g[i,(1+1):(1+s.i)] <- tvgarch.i$par.g[(1+1):(1+s.i)]
          par.g[i,(1+max.s+1):(1+max.s+s.i)] <- tvgarch.i$par.g[(1+s.i+1):(1+2*s.i)]
          par.g[i,(1+2*max.s+1):(1+2*max.s+sum(order.g[i,]))] <- tvgarch.i$par.g[(1+2*s.i+1):(1+2*s.i+sum(order.g[i,]))]
          h[,i] <- tvgarch.i$h
          g[,i] <- tvgarch.i$g
        }
      }
      if (is.null(order.g) || order.g[i,1] == 0) {
        tvgarch.i$par.h <- coef.garchx(object = tvgarch.i)
        h[,i] <- garchxRecursion(pars = as.numeric(tvgarch.i$par.h), aux = tvgarch.i)
        if (!is.null(tvgarch.i$fitted)) Objs[[paste("obj", i, sep = "")]]$fitted <- h[,i]
        if (!is.null(tvgarch.i$residuals)) Objs[[paste("obj", i, sep = "")]]$residuals <- y[,i]/sqrt(h[,i])
      }
      par.h[i,1] <- tvgarch.i$par.h[1]
      par.h[i,(1+1):(1+order.h[i,2])] <- tvgarch.i$par.h[(1+1):(1+order.h[i,2])]
      if (order.h[i,1] != 0) par.h[i,(1+max.q+1):(1+max.q+order.h[i,1])] <- tvgarch.i$par.h[(1+order.h[i,2]+1):(1+sum(order.h[i,1:2]))]
      if (order.h[i,3] != 0) par.h[i,(1+max.q+max.p+1):(1+max.q+max.p+order.h[i,3])] <- tvgarch.i$par.h[(1+sum(order.h[i,1:2])+1):(1+sum(order.h[i,]))]
      if (sum(order.x[i,]) != 0) par.h[i,(1+max.q+max.p+max.r+which(order.x[i,] == 1))] <- tvgarch.i$par.h[(1+sum(order.h[i,])+1):(1+sum(order.h[i,])+sum(order.x[i,]))]
      '
        Saving final standard errors
      ' 
      if (turbo == FALSE) {
        if (!is.null(order.g)) {
          if (order.g[i,1] != 0) {
            se.g[i,1] <- tvgarch.i$se.g[1]
            se.g[i,(1+1):(1+length(order.g[i,]))] <- tvgarch.i$se.g[(1+1):(1+length(order.g[i,]))]
            se.g[i,(1+length(order.g[i,])+1):(1+2*length(order.g[i,]))] <- tvgarch.i$se.g[(1+length(order.g[i,])+1):(1+2*length(order.g[i,]))]
            se.g[i,(1+2*length(order.g[i,])+1):(1+2*length(order.g[i,])+sum(order.g[i,]))] <- tvgarch.i$se.g[(1+2*length(order.g[i,])+1):(1+2*length(order.g[i,])+sum(order.g[i,]))]
          }
        }
        else tvgarch.i$se.h <- sqrt(diag(vcov.garchx(object = tvgarch.i)))
        se.h[i,1] <- tvgarch.i$se.h[1]
        se.h[i,(1+1):(1+order.h[i,2])] <- tvgarch.i$se.h[(1+1):(1+order.h[i,2])]
        if (order.h[i,1] != 0) se.h[i,(1+max.q+1):(1+max.q+order.h[i,1])] <- tvgarch.i$se.h[(1+order.h[i,2]+1):(1+sum(order.h[i,1:2]))]
        if (order.h[i,3] != 0) se.h[i,(1+max.q+max.p+1):(1+max.q+max.p+order.h[i,3])] <- tvgarch.i$se.h[(1+sum(order.h[i,1:2])+1):(1+sum(order.h[i,]))]
        if (sum(order.x[i,]) != 0) se.h[i,(1+max.q+max.p+max.r+which(order.x[i,] == 1))] <- tvgarch.i$se.h[(1+sum(order.h[i,])+1):(1+sum(order.h[i,])+sum(order.x[i,]))]
      }
    }
    if (!is.null(order.g) && !is.null(order.x) && is.null(xreg)) {
      phi2 <- rbind(colMeans(y^2/g), (y^2/g)[-n,])
      colnames(phi2) <- names.x
      if (trace == TRUE) {
        cat("\n")
        cat(" Iterative estimation for the volatility spillovers \n")
      }
      repeat{
        for(i in 1:m) {
          if (sum(order.x[i,]) != 0 && is.null(xreg) && order.g[which(order.x[i,] != 0),1] != 0) {
            if (order.g[i,1] != 0) {
              intercept.g.i <- par.g[i,1]
              size.i <- par.g[i,(1+1):(1+s.i)]
              speed.i <- par.g[i,(1+max.s+1):(1+max.s+s.i)]
              location.i <- par.g[i,(1+2*max.s+1):(1+2*max.s+sum(order.g[i,]))]
              par.g.i <- c(intercept.g.i,size.i,speed.i,location.i)
              order.g.i <- order.g[i,1:s.i]
            }
            else{
              order.g.i <- NULL
              intercept.g.i <- NULL
              size.i <- NULL
              speed.i <- NULL
              location.i <- NULL
            }
            intercept.h.i <- par.h[i,1]
            par.h.i <- intercept.h.i
            arch.i <- par.h[i,(1+1):(1+order.h[i,2])]
            par.h.i <- c(par.h.i, arch.i)
            if (order.h[i,1] != 0) {
              garch.i <- par.h[i,(1+max.q+1):(1+max.q+order.h[i,1])]
              par.h.i <- c(par.h.i, garch.i)
            }
            else garch.i <- NULL
            if (order.h[i,3] != 0) {
              asym.i <- par.h[i,(1+max.q+max.p+1):(1+max.q+max.p+order.h[i,3])]
              par.h.i <- c(par.h.i, asym.i)
            } 
            else asym.i <- NULL
            par.xreg.i <- par.h[i,(1+max.q+max.p+max.r+which(order.x[i,] == 1))]
            par.h.i <- c(par.h.i,par.xreg.i)
            xreg.i <- as.matrix(phi2[,which(order.x[i,] == 1)])
            colnames(xreg.i) <- colnames(phi2)[which(order.x[i,] == 1)]
            if (trace == TRUE) {
              cat("\n")
              cat("Estimating volatility for ", colnames(y)[i], "\n", sep = "")
            }
            yi.i <- as.matrix(y[,i])
            colnames(yi.i) <- names.y[i]
            if (nx == 1) turbo2 <- FALSE
            if (nx > 1) turbo2 <- TRUE
            tvgarch.i <- tvgarch(y = yi.i, order.g = order.g.i, order.h = order.h[i,], xtv = xtv, xreg = xreg.i, 
                                 opt = opt, initial.values = list(intercept.g = intercept.g.i, size = size.i, speed = speed.i, location = location.i,
                                                                  intercept.h = intercept.h.i, arch = arch.i, garch = garch.i, asym = asym.i, par.xreg = par.xreg.i), 
                                 trace = trace, turbo = turbo2)
            Objs[[paste("obj", i, sep = "")]] <- tvgarch.i
            if (order.g[i,1] != 0) conv.i[i] <- max(par.g.i - tvgarch.i$par.g, par.h.i - tvgarch.i$par.h) 
            else conv.i[i] <- max(par.h.i - tvgarch.i$par.h) 
            '
              Saving final estimated parameters
            '       
            if (order.g[i,1] != 0) {
              par.g[i,1] <- tvgarch.i$par.g[1]
              par.g[i,(1+1):(1+s.i)] <- tvgarch.i$par.g[(1+1):(1+s.i)]
              par.g[i,(1+max.s+1):(1+max.s+s.i)] <- tvgarch.i$par.g[(1+s.i+1):(1+2*s.i)]
              par.g[i,(1+2*max.s+1):(1+2*max.s+sum(order.g[i,]))] <- tvgarch.i$par.g[(1+2*s.i+1):(1+2*s.i+sum(order.g[i,]))]
              h[,i] <- tvgarch.i$h
              g[,i] <- tvgarch.i$g
            }
            if (order.g[i,1] == 0) {
              tvgarch.i$par.h <- coef.garchx(object = tvgarch.i)
              h[,i] <- garchxRecursion(pars = as.numeric(tvgarch.i$par.h), aux = tvgarch.i)
              if (!is.null(tvgarch.i$fitted)) Objs[[paste("obj", i, sep = "")]]$fitted <- h[,i]
              if (!is.null(tvgarch.i$residuals)) Objs[[paste("obj", i, sep = "")]]$residuals <- y[,i]/sqrt(h[,i])
            }
            par.h[i,1] <- tvgarch.i$par.h[1]
            par.h[i,(1+1):(1+order.h[i,2])] <- tvgarch.i$par.h[(1+1):(1+order.h[i,2])]
            if (order.h[i,1] != 0) par.h[i,(1+max.q+1):(1+max.q+order.h[i,1])] <- tvgarch.i$par.h[(1+order.h[i,2]+1):(1+sum(order.h[i,1:2]))]
            if (order.h[i,3] != 0) par.h[i,(1+max.q+max.p+1):(1+max.q+max.p+order.h[i,3])] <- tvgarch.i$par.h[(1+sum(order.h[i,1:2])+1):(1+sum(order.h[i,]))]
            par.h[i,(1+max.q+max.p+max.r+which(order.x[i,] == 1))] <- tvgarch.i$par.h[(1+sum(order.h[i,])+1):(1+sum(order.h[i,])+sum(order.x[i,]))]
            if (!is.null(tvgarch.i$message.h)) message.h[i] <- tvgarch.i$message.h
            '
              Saving final standard errors
            '        
            if (turbo == FALSE) {
              if (order.g[i,1] != 0) {
                se.g[i,1] <- tvgarch.i$se.g[1]
                se.g[i,(1+1):(1+length(order.g[i,]))] <- tvgarch.i$se.g[(1+1):(1+length(order.g[i,]))]
                se.g[i,(1+length(order.g[i,])+1):(1+2*length(order.g[i,]))] <- tvgarch.i$se.g[(1+length(order.g[i,])+1):(1+2*length(order.g[i,]))]
                se.g[i,(1+2*length(order.g[i,])+1):(1+2*length(order.g[i,])+sum(order.g[i,]))] <- tvgarch.i$se.g[(1+2*length(order.g[i,])+1):(1+2*length(order.g[i,])+sum(order.g[i,]))]
              }
              else tvgarch.i$se.h <- sqrt(diag(vcov.garchx(object = tvgarch.i)))
              se.h[i,1] <- tvgarch.i$se.h[1]
              se.h[i,(1+1):(1+order.h[i,2])] <- tvgarch.i$se.h[(1+1):(1+order.h[i,2])]
              if (order.h[i,1] != 0) se.h[i,(1+max.q+1):(1+max.q+order.h[i,1])] <- tvgarch.i$se.h[(1+order.h[i,2]+1):(1+sum(order.h[i,1:2]))]
              if (order.h[i,3] != 0) se.h[i,(1+max.q+max.p+1):(1+max.q+max.p+order.h[i,3])] <- tvgarch.i$se.h[(1+sum(order.h[i,1:2])+1):(1+sum(order.h[i,]))]
              se.h[i,(1+max.q+max.p+max.r+which(order.x[i,] == 1))] <- tvgarch.i$se.h[(1+sum(order.h[i,])+1):(1+sum(order.h[i,])+sum(order.x[i,]))]
            }
          }
          # third condition: xreg is mtv-garch and so estimates need to be updated
        }
        phi2 <- rbind(colMeans(y^2/g), (y^2/g)[-n,])
        colnames(phi2) <- names.x
        if (nx == 1) break
        if (nx != 1) { 
          iter <- iter + 1
          if (max(conv.i) < 1e-3) break
          if (trace == TRUE) {
            cat("\n")
            cat("Iteration round: ",iter,"\n",sep = "")
          }
        }
      }
    }
    '
      Computing Constant Conditional Correlations
    '
    sigma2 <- h * g
    colnames(sigma2) <- paste("sigma2.", 1:m, sep = "")
    colnames(h) <- paste("h", 1:m, sep = "")
    colnames(g) <- paste("g", 1:m, sep = "")
    z <- y/sqrt(sigma2)
    R <- cor(z)
    par.r <- R[lower.tri(R)]
    ccc <- matrix(par.r, n, m*(m-1)/2)
    ID <- 1:m
    namesR <- numeric(0)
    for (i in 1:m) {          
      for (j in 1:m) {    
        namesR <- c(namesR, paste(paste(colnames(y)[j]), colnames(y)[i], sep = "."))
      }
    }
    namesR <- matrix(namesR,m,m)
    namesR <- namesR[lower.tri(namesR)]
    colnames(ccc) <- namesR
    if (dcc == TRUE) {
      '
        Estimating Dynamic Conditional Correlations
      '
      if (is.null(initial.values$par.dcc)) initial.values$par.dcc = c(0.05, 0.80)   
      ui.dcc <- rbind(diag(2), -c(1,1))
      ci.dcc <- c(5e-5, 5e-5, -1+5e-5)
      fit.dcc <- constrOptim(theta = initial.values$par.dcc, f = dccObj, grad = NULL, ui = ui.dcc, ci = ci.dcc, z = z, sigma2 = sigma2, flag = 1)
      if (fit.dcc$convergence != 0 && is.null(fit.dcc$message) == FALSE) warning(paste(fit.dcc$convergence))
      if (turbo == FALSE) {
        jac.dcc <- jacobian(func = dccObj, x = fit.dcc$par, z = z, sigma2 = sigma2, flag = 0)
        J.dcc <- crossprod(jac.dcc)  
        H.dcc <- optimHess(par = fit.dcc$par, fn = dccObj, z = z, sigma2 = sigma2, flag = 1)
        vcov.dcc <- solve(-H.dcc) %*% J.dcc %*% solve(-H.dcc)
      }
      if (turbo == TRUE) vcov.dcc <- matrix(NA, 2, 2)
      names.dcc <- c("alpha", "beta")
      rownames(vcov.dcc) <- names.dcc
      colnames(vcov.dcc) <- names.dcc
      dcc <- dccObj(par.dcc = fit.dcc$par, z = z, sigma2 = sigma2, flag = 2)
      colnames(dcc) <- namesR
      par.dcc <- matrix(fit.dcc$par, 1, 2)
      se.dcc <- matrix(sqrt(diag(vcov.dcc)), 1, 2)
      colnames(par.dcc) <- names.dcc
      rownames(par.dcc) <- ""
      colnames(se.dcc) <- names.dcc
      rownames(se.dcc) <- ""
    }
    else par.dcc <- NULL
    if (!is.null(order.x) && is.null(xreg)) {
      spillovers <- TRUE
      if (nx > 0) xreg <- phi2
      else xreg <- y2
    }
    else spillovers <- FALSE
    results <- list(Objs = Objs, sigma2 = sigma2, residuals = z, par.h = par.h, se.h = se.h, h = h, g = g, R = R,
                    order.g = order.g, order.h = order.h, order.x = order.x, xtv = xtv, xreg = xreg, spillovers = spillovers, 
                    y = y, date = date(), turbo = turbo, trace = trace)
    if (is.null(par.dcc)) results$ccc = ccc
    if (!is.null(par.dcc)) {
      results$dcc <- dcc
      results$par.dcc <- par.dcc
      results$se.dcc <- se.dcc
      results$vcov.dcc <- vcov.dcc
      results$logLik.dcc <- - fit.dcc$value
    }
    if (!is.null(order.g)) {
      results$par.g <- par.g
      results$se.g <- se.g
      results$opt <- opt
    }
    class(results) <- "mtvgarch"
    return(results)
  }
