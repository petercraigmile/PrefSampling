


est.intensity <- function (locs, W, eps, as.im=FALSE, use.diggle=FALSE) {

    z <- ppp(x=locs[,1], y=locs[,2], window=W)

    if (use.diggle) {
        dens.est <- density.ppp(z, at="pixels", sigma=bw.diggle, diggle=TRUE,
                                eps=eps)
    } else {
        dens.est <- density.ppp(z, at="pixels", eps=eps)
    }

    if (as.im) {

        lambda.hat <- pmax(t(dens.est$v), 1e-8)
    } else {
    
        lambda.hat <- as.numeric(pmax(t(dens.est$v), 1e-8))
    }

    return(lambda.hat)
}



create.grid.of.counts <- function (locs, grid) {
    
    tt <- table(find.locs.in.grid(locs, grid))
    
    counts <- rep(0, nrow(grid))
    counts[as.numeric(names(tt))] <- tt
    counts
}



LGCP.sim <- function (grid, window,
                      mu.eta, theta.eta, cov.fun.eta) {

    ## Simulate the log intensity Gaussian process
    eta <- GSP.sim(n=1, rep(mu.eta, grid$M),
                   cov.fun=cov.fun.eta,
                   dists=Edist(grid$grid), theta=theta.eta)

    ## Calculate the first order intensity
    lambda    <- exp(eta)
    lambda.im <- as.im(t(matrix(lambda, grid$nx)), W=window)

    ## Simulate the LGCP locations
    pp   <- rpoispp(lambda.im)
    locs <- cbind(pp$x, pp$y)

    list(grid=grid, window=window,
         mu.eta=mu.eta,
         theta.eta=mu.eta,
         cov.fun.eta=cov.fun.eta,
         eta=eta,
         lambda=lambda, lambda.im=lambda.im,
         pp=pp, locs=locs)    
}




Thomas.intensity <- function (locs, parents, delta, omega2) {

    D <- Edist(locs, parents)
    b <- 2 * omega2

    delta/(pi*b) * rowSums(exp(-D^2 / b))
}


Thomas.sim.parents <- function (lambda0, window) {

    parents0 <- rpoispp(lambda0, win=window)

    cbind(parents0$x, parents0$y)
}


Thomas.sim <- function (grid, window,
                        parents, delta, omega2) {

    ## Calculate the first order intensity
    lambda    <- Thomas.intensity(grid$grid, parents, delta, omega2)
    lambda.im <- as.im(t(matrix(lambda, grid$nx)), W=window)

    ## Simulate the LGCP locations
    pp   <- rpoispp(lambda.im)
    locs <- cbind(pp$x, pp$y)

    list(grid=grid, window=window,
         parents=parents, delta=delta, omega2=omega2,
         lambda=lambda, lambda.im=lambda.im,
         pp=pp, locs=locs)    
}

