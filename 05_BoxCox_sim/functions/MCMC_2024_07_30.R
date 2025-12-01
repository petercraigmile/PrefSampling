
MCMC <- function () {
    ## ======================================================================
    ## Create a new MCMC object
    ## ======================================================================

    x <- new.env()
    class(x) <- "MCMC"
    x
}


"[.MCMC" <- function (x, var.name) {
    ## ======================================================================
    ## Return the chain for the variable 'var.name' converted to an
    ## array from the MCMC object x.
    ## ======================================================================

    ch <- paste(var.name, ".chain", sep = "")

    simplify2array(get(ch, env=x))
}


create <- function(x, ...) UseMethod("create")

chain <- function(x, ...) UseMethod("chain")


create.MCMC <- function (x, var.name, value) {
    ## ======================================================================
    ## Initialize the chain called 'var.name' with value 'value' in
    ## the MCMC object 'x'.
    ## ======================================================================

    assign(var.name, value, env=x)
    assign(paste(var.name, ".chain", sep=""), list(), envir=x)
}



update.MCMC <- function (x, var.name, value, save) {

    if (missing(save)) {

        if (exists("iter", x)) {
            
            save <- (!x$burn.in) & ((x$iter %% x$thin) == 0)
        } else {

            save <- TRUE
        }
    }
    
    assign(var.name, value, env=x)

    ## if we want to append to the chain (when save is TRUE)
    if (save) {

        chain <- paste(var.name, ".chain", sep = "")
        m <- length(get(chain, env=x)) + 1
        eval(substitute(a[[m]] <- av, list(a = as.name(chain), 
                                           av = value, m = m)), env=x)
    }
}



plot.MCMC <- function (x, var.name, ...) {

    y <- x[var.name]

    if (any(class(y)=="numeric")) {

        plot(y, type="l", ylab=var.name, lty=1, ...)
    } else {

        matplot(t(y), type="l", ylab=var.name, lty=1, ...)
    }
}



run.MCMC <- function (x, the.vars, max.iter, thin = 1, every = 1000,
                       burn.in = FALSE, how.long=TRUE) {

    if (how.long)
        t1 <- proc.time()

    x$burn.in <- burn.in
    x$thin    <- thin
    
    for (k in 1:max.iter) {

        x$iter <- k-1
        
        if (k%%every == 0) {
            
            if (how.long) {
                t2 <- proc.time()
                delta <- round(as.numeric((t2-t1)[3]) * (max.iter-k) / k, 1)
                cat(round(k/max.iter * 100, 1), "% (", delta, "s remaining) ", sep="")
            } else {
                cat(round(k/max.iter * 100, 1), "% ", sep="")
            }
        }
        
        sapply(the.vars, function(the.var)
            eval(call(paste("update.", the.var, sep = ""), x = x)))
    }
}




MCMC.adapt.cov <- function (theta) {

    if (any(class(theta)=="numeric")) {

        2.38^2 * var(theta)
        
    } else {

        2.38^2 * drop(cov(t(theta))) / nrow(theta)
    }
}




MH <- function (ll.new, ll.old) {

    log(runif(1)) <= (ll.new - ll.old)
}


## ======================================================================
## Inverse gamma functions
## Defined so that it is reciprocal of the gamma function in R
## ======================================================================

dinvgamma <- function (x, shape, rate=1, scale=1/rate, log=FALSE) {

    if (log) {
        
        -rate/x + shape*log(rate) - lgamma(shape) - (shape+1.0)*log(x)
        
    } else {
        
        exp(-rate/x) * rate^shape / (gamma(shape) * x^(shape+1.0))
    }
}



rinvgamma <- function (n, shape, rate=1.0, scale=1.0/rate) {

  1.0 / rgamma(n=n, shape=shape, rate=rate)
}



log.dmvnorm <- function (x, mean = rep(0, p), sigma = diag(p), chol=FALSE) {
    ## ======================================================================
    ## Adapted from dmvnorm in the mvtnorm R library.
    ## If chol=TRUE, then sigma is the Cholesky factor.
    ##
    ## The checkSymmetry option has been removed.  This is assumed to be true
    ## if the covariance matrix is passed to the function.
    ## ======================================================================

    if (is.vector(x)) 
        x <- matrix(x, ncol = length(x))
    
    p <- ncol(x)
    
    if (!missing(mean)) {
        if (!is.null(dim(mean))) 
            dim(mean) <- NULL
        if (length(mean) != p) 
            stop("x and mean have non-conforming size")
    }
    
    if (!missing(sigma)) {
        if (p != ncol(sigma)) 
            stop("x and sigma have non-conforming size")
    }

    if (chol) {
        dec <- sigma
    } else {
        dec <- tryCatch(base::chol(sigma), error = function(e) e)
    }

    if (inherits(dec, "error")) {
        x.is.mu <- colSums(t(x) != mean) == 0
        logretval <- rep.int(-Inf, nrow(x))
        logretval[x.is.mu] <- Inf
    }
    else {
        tmp <- base::backsolve(dec, t(x) - mean, transpose = TRUE)
        rss <- colSums(tmp^2)
        logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 * pi) - 0.5 * rss
    }
    
    names(logretval) <- rownames(x)
    logretval
}



rmvnorm.cond.precision <- function (m, P, P.chol=chol(P)) {
    ## ======================================================================
    ## Purpose: Sample from N( P^{-1} m, P^{-1} )
    ## Assumes: m and P are real-valued.
    ##          The length of m is equal to the dimension of the square matrix P.
    ## ======================================================================

    a     <- forwardsolve(t(P.chol), m)
    draw  <- backsolve(P.chol, a + rnorm(nrow(P)))
  
    return(draw)
}




##methods(class="MCMC")


update.var.invgamma <- function (x, mu=0, corr, chol=FALSE,
                                 shape=0.01, rate=0.01, scale=1/rate) {
    ## ======================================================================
    ## Update tau^2 for X ~ N(mu, tau^2 corr),
    ##        where tau^2 ~ inverse gamma(shape, rate)
    ##
    ## If 'chol'=TRUE, then 'corr' is the Cholesky factor of the correlation
    ## matrix of 'x'; otherwise 'corr' is the correlation structure for 'x'.
    ## If the correlation matrix 'corr' is missing, assume an uncorrelated
    ## correlation structure.
    ##
    ## Assumes 'x' is a vector.
    ## ======================================================================

    if (missing(corr)) {

        RSS <- sum( (x - mu)^2 )
            
    } else {
    
        if (chol) {
            
            chol.corr <- corr
        } else {
            
            chol.corr <- tryCatch(base::chol(corr), error=function (e) e)
        }

        Lz  <- base::backsolve(chol.corr, matrix(x - mu, nrow=length(x)), transpose=TRUE)
        RSS <- sum(Lz^2)
    }
    
    1.0 / rgamma(1, shape + 0.5 * length(x), rate + 0.5 * RSS)
}



update.reg.beta.normal <- function (y, X, sigma, chol=FALSE,
                                    beta.mu, beta.prec) {
    ## ======================================================================
    ## Update beta for Y ~ N(X beta, sigma)
    ##        where beta ~ N(beta.mu, beta.prec^{-1})
    ##
    ## If 'chol'=TRUE, then 'sigma' is the Cholesky factor of the covariance
    ## matrix of 'y'; otherwise 'sigma' is the covariance structure for 'y'.
    ##
    ## Assumes 'y' is a vector.every
    ## ======================================================================

    if (chol) {
        
        chol.sigma <- sigma
    } else {
        
        chol.sigma <- tryCatch(base::chol(sigma), error=function (e) e)
    }

    LX <- base::backsolve(chol.sigma, X, transpose=TRUE)
    Ly <- base::backsolve(chol.sigma, matrix(y, nrow=length(y)), transpose=TRUE)

    P <- crossprod(LX)     + beta.prec
    m <- crossprod(LX, Ly) + beta.prec %*% beta.mu
    
    rmvnorm.cond.precision(m, P)
}






