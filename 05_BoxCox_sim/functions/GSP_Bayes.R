
update.GSP.pars <- function (tau2, lambda, z, X, dists,
                             cov.fun, log.lambda.prop.sd,
                             R, chol.R,
                             beta.mu, beta.prec,
                             tau2.shape, tau2.rate,
                             lambda.shape, lambda.rate,
                             other.pars) {
    
    if (missing(chol.R)) {
        
        if (missing(R)) {

            if (missing(other.pars)) {
                R <- cov.fun(dists, c(1, lambda))
            } else {
                R <- cov.fun(dists, c(1, lambda, other.pars))
            }
        }

        chol.R <- tryCatch(base::chol(R), error = function(e) e)        
    }

    new.beta <- update.reg.beta.normal(z, X, sqrt(tau2) * chol.R, chol=TRUE,
                                       beta.mu=beta.mu, beta.prec=beta.prec)

    new.mu <- drop(X %*% new.beta)
    
    new.tau2 <- update.var.invgamma(z, new.mu, chol.R, chol=TRUE,
                                    shape=tau2.shape, rate=tau2.rate)
    
    log.lambda <- log(lambda)
    
    new.log.lambda <- rnorm(1, log.lambda, log.lambda.prop.sd)
    new.lambda     <- exp(new.log.lambda)
    
    if (missing(other.pars)) {
        new.R      <- cov.fun(dists, c(1, new.lambda))
    } else {
        new.R      <- cov.fun(dists, c(1, new.lambda, other.pars))
    }
    new.chol.R <- tryCatch(base::chol(new.R), error = function(e) e)

    res <- z - new.mu
    
    ll.new <- log.dmvnorm(res, sigma = new.chol.R * sqrt(new.tau2), chol=TRUE) +
        dgamma(new.lambda, lambda.shape, lambda.rate, log=TRUE) +
        new.log.lambda ## log Jacobian
    
    ll.old <- log.dmvnorm(res, sigma = chol.R * sqrt(new.tau2), chol=TRUE) +
        dgamma(lambda, lambda.shape, lambda.rate, log=TRUE) +
        log.lambda ## log Jacobian
    
    if (MH(ll.new, ll.old)) {

        list(beta=new.beta,
             mu=new.mu,
             tau2=new.tau2,
             lambda=new.lambda,
             jump=TRUE,
             R=new.R,
             chol.R=new.chol.R)
    } else {

        list(beta=new.beta,
             mu=new.mu,
             tau2=new.tau2,
             lambda=lambda,
             jump=FALSE)
    }
}



update.GSP.pars.no.beta <- function (tau2, lambda, z, dists,
                                     cov.fun, log.lambda.prop.sd,
                                     R, chol.R,
                                     tau2.shape, tau2.rate,
                                     lambda.shape, lambda.rate,
                                     other.pars) {

    if (missing(chol.R)) {
        
        if (missing(R)) {

            if (missing(other.pars)) {
                R <- cov.fun(dists, c(1, lambda))
            } else {
                R <- cov.fun(dists, c(1, lambda, other.pars))
            }
        }

        chol.R <- tryCatch(base::chol(R), error = function(e) e)        
    }
    
    new.tau2 <- update.var.invgamma(z, corr=chol.R, chol=TRUE,
                                    shape=tau2.shape, rate=tau2.rate)
    
    log.lambda <- log(lambda)
    
    new.log.lambda <- rnorm(1, log.lambda, log.lambda.prop.sd)
    new.lambda     <- exp(new.log.lambda)

    if (missing(other.pars)) {
        new.R      <- cov.fun(dists, c(1, new.lambda))
    } else {
        new.R      <- cov.fun(dists, c(1, new.lambda, other.pars))
    }
    new.chol.R <- tryCatch(base::chol(new.R), error = function(e) e)
    
    ll.new <- log.dmvnorm(z, sigma = new.chol.R * sqrt(new.tau2), chol=TRUE) +
        dgamma(new.lambda, lambda.shape, lambda.rate, log=TRUE) +
        new.log.lambda ## log Jacobian
    
    ll.old <- log.dmvnorm(z, sigma = chol.R * sqrt(new.tau2), chol=TRUE) +
        dgamma(lambda, lambda.shape, lambda.rate, log=TRUE) +
        log.lambda ## log Jacobian
    
    if (MH(ll.new, ll.old)) {

        list(tau2=new.tau2,
             lambda=new.lambda,
             jump=TRUE,
             R=new.R,
             chol.R=new.chol.R)
    } else {

        list(tau2=new.tau2,
             lambda=lambda,
             jump=FALSE)
    }
}



##update.GSP.pars(1, 2, chs$beta[,1], X, D, cov.fun=GSP.exp,
##                beta.mu=0, beta.prec=1/10, log.lambda.prop.sd=0.3)




update.GSP.pars2 <- function (tau2, lambda, z, X, dists,
                              cov.fun, log.theta.prop.sd,
                              R, chol.R,
                              beta.mu, beta.prec,
                              tau2.shape, tau2.rate,
                              lambda.shape, lambda.rate) {
    
    if (missing(chol.R)) {
        
        if (missing(R)) {
            
            R <- cov.fun(dists, c(1, lambda))
        }

        chol.R <- tryCatch(base::chol(R), error = function(e) e)        
    }

    new.beta <- update.reg.beta.normal(z, X, sqrt(tau2) * chol.R, chol=TRUE,
                                       beta.mu=beta.mu, beta.prec=beta.prec)

    new.mu <- drop(X %*% new.beta)
    
    log.lambda <- log(lambda)
    log.tau2   <- log(tau2)

    log.theta <- c(log.tau2, log.lambda)
    
    new.log.theta <- rnorm(2, log.theta, log.theta.prop.sd)
    new.theta     <- exp(new.log.theta)
    
    new.R      <- cov.fun(dists, c(1, new.theta[2]))
    new.chol.R <- tryCatch(base::chol(new.R), error = function(e) e)

    res <- z - new.mu
    
    ll.new <- log.dmvnorm(res, sigma = new.chol.R * sqrt(new.theta[1]), chol=TRUE) +
        dgamma(new.theta[2], lambda.shape, lambda.rate, log=TRUE) +
        dinvgamma(new.theta[1], tau2.shape, tau2.rate, log=TRUE) +
        sum(new.log.theta) ## log Jacobian
    
    ll.old <- log.dmvnorm(res, sigma = chol.R * sqrt(tau2), chol=TRUE) +
        dgamma(lambda, lambda.shape, lambda.rate, log=TRUE) +
        dinvgamma(tau2, tau2.shape, tau2.rate, log=TRUE) +
        sum(log.theta) ## log Jacobian
    
    if (MH(ll.new, ll.old)) {

        list(beta=new.beta,
             tau2=new.theta[1],
             lambda=new.theta[2],
             jump=TRUE,
             R=new.R,
             chol.R=new.chol.R)
    } else {

        list(beta=new.beta,
             tau2=tau2,
             lambda=lambda,
             jump=FALSE)
    }
}
