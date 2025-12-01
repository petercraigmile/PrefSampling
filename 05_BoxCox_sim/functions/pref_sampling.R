

Box.Cox <- function (x, kappa) {

    if (abs(kappa) < 1e-8) {

        log(x)

    } else {
        
        (x^kappa-1.0) / kappa
    }
}



inv.Box.Cox <- function (x, kappa) {

    if (abs(kappa) < 1e-8) {

        exp(x)

    } else {
        
        (kappa * x + 1.0)^(1/kappa)
    }
}



GSP.pref.sim <- function (pp,
                          mu.W, theta.W, W.cov.fun,
                          alpha, sigma2.Z, kappa=0) {
    
    sites <- rbind(pp$locs, pp$grid$grid)

    ## Simulate W at the point process and grid locations
    W.all <- GSP.sim(1, rep(mu.W, nrow(sites)), cov.fun=W.cov.fun,
                     dists=Edist(sites), theta=theta.W)

    sel <- 1:nrow(pp$locs)

    in.grid <- find.locs.in.grid(pp$locs, pp$grid$grid)

    ## W at the point process locations    
    W      <- W.all[sel]

    ## W at the grid locations
    W.grid <- W.all[-sel]

    nu <- alpha * Box.Cox(exp(pp$eta), kappa)
        
    ## Y at the point process locations
    Y <- W + nu[in.grid]

    ## Y at the grid locations
    Y.grid <- W.grid + nu
    
    ## Simulate the Gaussian observation process
    Z <- rnorm(length(Y), Y, sqrt(sigma2.Z))
       
    list(pp=pp,
         mu.W=mu.W,
         theta.W=theta.W,
         W.cov.fun=W.cov.fun,
         sigma2.Z=sigma2.Z,
         alpha=alpha, kappa=kappa,
         W=W, W.grid=W.grid,
         Y=Y, Y.grid=Y.grid,
         Z=Z)
}




pred.Y <- function (Z, locs, pred.grid) {
    ## Assumes constant mean

    D <- Edist(locs)    
    X <- cbind(rep(1, nrow(D)))
    
    est <- GSP.mle(c(var(Z)/2,1,var(Z)/2), z=Z, X=X, cov.fun=GSP.exp, dists=D)
   
    that <- est$theta.hat
    mhat <- est$beta.hat
    
    B <- Edist(locs, pred.grid)
    
    CC <- solve(GSP.exp(D, that), GSP.exp(B, that[1:2]))

    Y.hat <- drop(mhat + crossprod(CC, Z - mhat))

    list(est=est, Y.hat=Y.hat)
}




pred.Y.pref <- function (Z, locs, pred.grid, W, eps, use.diggle=FALSE) {
    ## Assumes constant mean

    lambda.hat <- est.intensity(locs, W, eps, use.diggle=use.diggle)

    eta.hat <- log(lambda.hat)

    pxs <- seq(W$xrange[1]+eps/2, W$xrange[2]-eps/2, eps)
    pys <- seq(W$yrange[1]+eps/2, W$yrange[2]-eps/2, eps)

    pg <- make.grid(list(pxs, pys))

    D <- Edist(locs)

    sel1 <- find.locs.in.grid(locs,      pg)
    sel2 <- find.locs.in.grid(pred.grid, pg)
    
    X  <- cbind(rep(1, nrow(locs)),      eta.hat[sel1])
    PX <- cbind(rep(1, nrow(pred.grid)), eta.hat[sel2])
    
    est <<- GSP.mle(c(var(Z)/2,1,var(Z)/2), z=Z, X=X, cov.fun=GSP.exp, dists=D)
    
    that <- est$theta.hat
    bhat <- est$beta.hat

    mhat.obs  <- X  %*% bhat
    mhat.pred <- PX %*% bhat
    
    B <- Edist(locs, pred.grid)
    
    CC <- solve(GSP.exp(D, that), GSP.exp(B, that[1:2]))

    Y.hat <- drop(mhat.pred + crossprod(CC, Z - mhat.obs))

    list(est=est, Y.hat=Y.hat)
}


