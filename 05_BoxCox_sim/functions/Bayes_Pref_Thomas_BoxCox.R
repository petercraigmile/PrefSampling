


##lambda.eta.shape <- 15
##lambda.eta.shape <- 100


calc.nu <- function (lambda, alpha, kappa, cell.area) {

    alpha * Box.Cox(lambda, kappa) 
}



init.Pref.Thomas <- function (z, locs, grid, parents, cell.area,
                              locs.pp = locs,
                              beta.mu=c(0,0),
                              beta.prec=diag(c(1, 1)),
                              delta.shape = 5,
                              delta.rate  = 1,
                              lambda.shape.W=30,
                              lambda.rate.W=100) {
    
    x <- MCMC()

    x$z <- z
        
    x$D      <- Edist(locs)    
    x$D.grid <- Edist(grid)
    
    x$parents   <- parents
    
    x$delta.shape <- delta.shape
    x$delta.rate  <- delta.rate

    x$lambda.shape.W <- lambda.shape.W
    x$lambda.rate.W  <- lambda.rate.W

    x$omega2.shape <- 16
    x$omega2.rate  <- 1
    
    x$prop.log.delta.sd  <- 0.07
    x$prop.log.omega2.sd <- 0.05
    x$prop.kappa.sd      <- 0.02
#    x$prop.alpha         <- 0.02
    
    x$locs.in.grid <- find.locs.in.grid(locs, grid)

    x$counts    <- create.grid.of.counts(locs.pp, grid)
    x$cell.area <- cell.area   
    x$grid      <- grid

    x$beta.mu   <- beta.mu
    x$beta.prec <- beta.prec
    
    create(x, "theta.W", c(var(z)/2, 1))
    create(x, "theta.W.jumps", FALSE)
    create(x, "sigma2", var(z)/2)
    create(x, "beta", c(mean(z)))
    create(x, "y", rnorm(length(z), z, sqrt(x$sigma2)))
    
    create(x, "alpha", 1)

    create(x, "kappa", 0)
    create(x, "kappa.jumps", FALSE)
    
    x$Sigma.W     <- GSP.exp(x$D, x$theta.W)
    x$Sigma.W.inv <- solve(x$Sigma.W)
    
    create(x, "delta", 1)
    create(x, "omega2", 1)
    create(x, "delta.omega2.jumps", FALSE)
    
    return(x)
}



fit.Pref.Thomas  <- function (x, max.iter, thin = 10, every = 500, burn.in=FALSE) {
    
    to.upd <- c("sigma2.y", "kappa", "theta.W.alpha",  "delta.omega2")

    run.MCMC(x, to.upd, max.iter=max.iter, every=every, thin=thin, burn.in=burn.in)

    return(x)
}



fit.Pref.Thomas.fix.kappa  <- function (x, max.iter, thin = 10, every = 500, burn.in=FALSE) {
    
    to.upd <- c("sigma2.y", "kappa.dummy", "theta.W.alpha",  "delta.omega2")

    run.MCMC(x, to.upd, max.iter=max.iter, every=every, thin=thin, burn.in=burn.in)

    return(x)
}



update.sigma2.y <- function (x) {

    new.sigma2 <- update.var.invgamma(x$z, x$y)
    
    update(x, "sigma2", new.sigma2)

    lambda.y <- Thomas.intensity(x$grid, x$parents, x$delta, x$omega2)[x$locs.in.grid]

    mu  <- rep(x$beta, nrow(x$D)) + calc.nu(lambda.y, x$alpha, x$kappa, x$cell.area)

    PP <- diag(rep(1/x$sigma2, nrow(x$D))) + x$Sigma.W.inv
    mm <- x$z/x$sigma2 + x$Sigma.W.inv %*% mu

    new.y <- drop(rmvnorm.cond.precision(mm, PP))

    update(x, "y", new.y)
}





update.theta.W.alpha <- function (x) {

    lambda.y <- Thomas.intensity(x$grid, x$parents, x$delta, x$omega2)[x$locs.in.grid]

    nu <- calc.nu(lambda.y, x$alpha, x$kappa, x$cell.area)

    X <- cbind(rep(1, nrow(x$D)), nu)
    
    mh <- update.GSP.pars(x$theta.W[1], x$theta.W[2], x$y, X, x$D,
                          GSP.exp, 0.07,
                          beta.mu=x$beta.mu, beta.prec=x$beta.prec,
                          tau2.shape=1, tau2.rate=1,
                          lambda.shape=x$lambda.shape.W,
                          lambda.rate=x$lambda.rate.W)

    update(x, "beta", mh$beta[1])
    update(x, "alpha", mh$beta[2])
    update(x, "theta.W", c(mh$tau2, mh$lambda))
    update(x, "theta.W.jumps", mh$jump)

    if (mh$jump) {

        x$Sigma.W     <- mh$tau2 * mh$R
        x$Sigma.W.inv <- solve(x$Sigma.W)
    }
}




update.delta.omega2 <- function (x) {

    cur.delta  <- x$delta
    cur.omega2 <- x$omega2

    res <- x$y - x$beta
    
    new.log.delta <- rnorm(1, log(cur.delta), x$prop.log.delta.sd)
    new.delta     <- exp(new.log.delta)   
    
    new.log.omega2 <- rnorm(1, log(cur.omega2), x$prop.log.omega2.sd)
    new.omega2     <- exp(new.log.omega2)
    
    new.lambda <- Thomas.intensity(x$grid, x$parents, new.delta, new.omega2)
    cur.lambda <- Thomas.intensity(x$grid, x$parents, cur.delta, cur.omega2)
    
    new.lambda.y <- new.lambda[x$locs.in.grid]
    cur.lambda.y <- cur.lambda[x$locs.in.grid]

    new.nu <- calc.nu(new.lambda.y, x$alpha, x$kappa, x$cell.area)
    cur.nu <- calc.nu(cur.lambda.y, x$alpha, x$kappa, x$cell.area)

    new.r <- res - new.nu
    cur.r <- res - cur.nu
    
    new.term <- -0.5 * sum(crossprod(x$Sigma.W.inv, new.r) * new.r)
    cur.term <- -0.5 * sum(crossprod(x$Sigma.W.inv, cur.r) * cur.r)

    lp.new <- sum(dpois(x$counts, x$cell.area * new.lambda, log=TRUE)) +
        new.term +
        dgamma(new.delta, x$delta.shape, x$delta.rate, log=TRUE) + new.log.delta +
        dgamma(new.omega2, x$omega2.shape, x$omega2.rate, log=TRUE) + new.log.omega2
    
    lp.cur <- sum(dpois(x$counts, x$cell.area * cur.lambda, log=TRUE)) +
        cur.term +
        dgamma(cur.delta, x$delta.shape, x$delta.rate, log=TRUE) + log(cur.delta) +
        dgamma(cur.omega2, x$omega2.shape, x$omega2.rate, log=TRUE) + log(cur.omega2)
    
    jump <- log(runif(1)) <= (lp.new - lp.cur)
    
    update(x, "delta.omega2.jumps", jump)
    
    if (jump) {
        
        update(x, "delta",  new.delta)
        update(x, "omega2", new.omega2)
    } else {
        
        update(x, "delta" , cur.delta)
        update(x, "omega2", cur.omega2)
    }
}




update.kappa <- function (x) {

    lambda.y <- Thomas.intensity(x$grid, x$parents, x$delta, x$omega2)[x$locs.in.grid]

    cur.kappa <- x$kappa
    new.kappa <- rnorm.truncated(1, cur.kappa, x$prop.kappa.sd, -1, 1)

    res <- x$y - x$beta

    new.nu <- calc.nu(lambda.y, x$alpha, new.kappa, x$cell.area)
    cur.nu <- calc.nu(lambda.y, x$alpha, cur.kappa, x$cell.area)

    new.r <- res - new.nu
    cur.r <- res - cur.nu
    
    new.term <- -0.5 * sum(crossprod(x$Sigma.W.inv, new.r) * new.r)
    cur.term <- -0.5 * sum(crossprod(x$Sigma.W.inv, cur.r) * cur.r)

    lp.new <- new.term + dnorm.truncated(new.kappa, 0, 0.3, -1, 1, log=TRUE) -
        dnorm.truncated(new.kappa, cur.kappa, x$prop.kappa.sd, -1, 1, log=TRUE)
    lp.cur <- cur.term + dnorm.truncated(cur.kappa, 0, 0.3, -1, 1, log=TRUE) -
        dnorm.truncated(cur.kappa, new.kappa, x$prop.kappa.sd, -1, 1, log=TRUE)

    jump <- MH(lp.new, lp.cur)

    update(x, "kappa.jumps", jump)

    if (jump) {

        update(x, "kappa", new.kappa)
    } else {
        
        update(x, "kappa", cur.kappa)
    }
}



update.kappa.dummy <- function (x) {

    update(x, "kappa", x$kappa)
    update(x, "kappa.jumps", FALSE)
}





Pref.Thomas.predict.GP <- function (x, plocs, every=50, CORES=5, K) {

    tt <- x["theta.W"]
    bb <- x["beta"]
    aa <- x["alpha"]
    yy <- x["y"]
    dd <- x["delta"]
    oo <- x["omega2"]
    kk <- x["kappa"]

#    X  <- cbind(rep(1, length(x$z)))
#    PX <- cbind(rep(1, length(plocs)))
    
    X <- x$X 

    PD <- Edist(plocs)
    
    B <- Edist(locs, plocs)

    ks <- 1:ncol(tt)

    if (!missing(K)) {

        ks <- sort(sample(ks, K))
    }
    
    ws <- mclapply(ks, function (k) {
        
        if (k%%every==0) cat(k, " ")
        
        Sig  <- GSP.exp(x$D, tt[,k])
        BSig <- GSP.exp(B,   tt[,k])
        PSig <- GSP.exp(PD,  tt[,k])
        
        CC <- solve(Sig, BSig)
        
        mu.P <- bb[k]
        mu   <- bb[k]
        
        cond.mean <- mu.P + crossprod(CC, yy[,k] - mu)
        cond.cov  <- PSig - crossprod(BSig, CC)
        
        ws <- drop(rmvnorm(1, cond.mean, cond.cov))
        
        ll <- Thomas.intensity(Plocs, x$parents, dd[k], oo[k])

        ws + calc.nu(ll, aa[k], kk[k], x$cell.area)
    
        ##    list(cond.mean=cond.mean, cond.cov=cond.cov)
    }, mc.cores=CORES)

    simplify2array(ws)
}

