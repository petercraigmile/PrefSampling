
## ======================================================================
## 'GSP' stands for Gaussian stationary process, defined by
## covariances through distances.
## ======================================================================


GSP.exp <- function (h, theta, transform.fun, eps=1e-12) {
    ## ======================================================================
    ## Exponential covariance at distance 'h' and parameter values 'theta'.
    ## (Matern with nu=1/2).
    ## If 'transform.fun' is given transform 'theta' using this function.
    ##
    ## On the original scale (untransformed),
    ## theta = (partial sill, range parameter) or
    ## theta = (partial sill, range parameter, nugget).
    ## ======================================================================
    
    p <- length(theta)

    if (p<2 || p>3) {

        stop("The length of 'theta' must be 2 or 3")
    }

    if (missing(transform.fun)) {

        drop(cppGSPexp(rbind(h), theta, eps))
        
    } else {
        
        drop(cppGSPexp(rbind(h), transform.fun(theta), eps))
    }
}



GSP.Diggle <- function (h, theta, transform.fun, eps=1e-12) {
    ## ======================================================================
    ## Diggle covariance at distance 'h' and parameter values 'theta'
    ## (Matern with nu=3/2).
    ## If 'transform.fun' is given transform 'theta' using this function.
    ##
    ## theta = (partial sill, range parameter) or
    ## theta = (partial sill, range parameter, nugget).
    ## ======================================================================

    p <- length(theta)

    if (p<2 || p>3) {

        stop("The length of 'theta' must be 2 or 3")
    }

    if (!missing(transform.fun)) {
        
         theta <- transform.fun(theta)
    }
    
    if (missing(transform.fun)) {
        
        drop(cppGSPDiggle(rbind(h), theta, eps))
        
    } else {
        
        drop(cppGSPDiggle(rbind(h), transform.fun(theta), eps))
    }
}



GSP.Gaussian <- function (h, theta, transform.fun, eps=1e-12) {
    ## ======================================================================
    ## Gaussian covariance at distance 'h' and parameter values 'theta'
    ## (Matern with nu=infinity).
    ## If 'transform.fun' is given transform 'theta' using this function.
    ##
    ## theta = (partial sill, range parameter) or
    ## theta = (partial sill, range parameter, nugget).
    ## ======================================================================

    p <- length(theta)

    if (p<2 || p>3) {

        stop("The length of 'theta' must be 2 or 3")
    }

    if (!missing(transform.fun)) {
        
         theta <- transform.fun(theta)
    }

    if (missing(transform.fun)) {

        drop(cppGSPGaussian(rbind(h), theta, eps))
        
    } else {
        
        drop(cppGSPGaussian(rbind(h), transform.fun(theta), eps))
    }
}




GSP.powerexp.fix.power <- function(h, theta, transformed.fun, esp=1e-12){
  if(length(theta) == 2){
    GSP.powerexp(h, c(theta[1], theta[2], 1.75), transformed.fun)
  }else{
    GSP.powerexp(h, c(theta[1], theta[2], 1.75, theta[3]), transformed.fun)
  }
}


GSP.powerexp <- function (h, theta, transform.fun, eps=1e-12) {
    ## ======================================================================
    ## Power exponential covariance at distance 'h' and parameter values 'theta'.
    ## If 'transform.fun' is given transform 'theta' using this function.
    ##
    ## On the original scale (untransformed),
    ## theta = (partial sill, range parameter, exponent) or
    ## theta = (partial sill, range parameter, exponent, nugget).
    ## ======================================================================

    p <- length(theta)

    if (p<3 || p>4) {

        stop("The length of 'theta' must be 3 or 4")
    }

    if (missing(transform.fun)) {
        
        drop(cppGSPpowerexp(rbind(h), theta, eps))
        
    } else {
        
        drop(cppGSPpowerexp(rbind(h), transform.fun(theta), eps))
    }
}



GSP.Matern <- function (h, theta, transform.fun, eps=1e-12) {
    ## ======================================================================
    ## Matern covariance at distance 'h' and parameter values 'theta'.
    ## If 'transform.fun' is given transform 'theta' using this function.
    ##
    ## On the original scale (untransformed),
    ## theta = (partial sill, range parameter, smoothness parameter) or
    ## theta = (partial sill, range parameter, smoothness parameter, nugget).
    ## ======================================================================

    p <- length(theta)

    if (p<3 || p>4) {

        stop("The length of 'theta' must be 3 or 4")
    }

    if (!missing(transform.fun)) {
        
         theta <- transform.fun(theta)
    }

    h.over.rp <- h / theta[2]
    nu        <- theta[3]   

    if (p==3) {
        
        Cov <- ifelse(h > eps,
                      theta[1] * 2.0^(1.0-nu) * h.over.rp^nu / gamma(nu) * besselK(h.over.rp, nu),
                      theta[1])
    } else {

        Cov <- ifelse(h > eps,
                      theta[1] * 2.0^(1.0-nu) * h.over.rp^nu / gamma(nu) * besselK(h.over.rp, nu),
                      theta[1] + theta[4])
    }

    return(Cov)
}



Edist <- function (x, y=x) {

    if (!any(class(x)=="matrix")) {
        stop("'x' must be a matrix")
    }
    if (!any(class(y)=="matrix")) {
        stop("'y' must be a matrix")
    }
    if (ncol(x) != ncol(y)) {
        stop("'x' and 'y' must have the same number of columns")
    }
    
    cppEdist(x,y)
}



make.grid <- function (list.of.coords) {

    cbind(apply(expand.grid(list.of.coords), 2, function (x) x))
}



GSP.sim <- function (n, mu, cov.fun, dists, theta, X, beta)
  ## ======================================================================
  ## By Peter F. Craigmile, pfc@stat.osu.edu
  ##
  ## Simulate 'n' realizations of a Gaussian stochastic process of
  ## length 'nrow(dists)' with mean 'mu' and a covariance function
  ## 'cov.fun'.  If 'X' is given calculate the mean using mu = X %*% beta.
  ## ======================================================================
{
    if (!missing(X)) {

        mu <- drop(X %*% beta)
    }

    ## calculate the covariance matrix
    Sigma <- cov.fun(dists, theta)

    if (missing(mu)) {

        drop(mvtnorm::rmvnorm(n=n, sigma=Sigma))
    } else {    
        
        drop(mvtnorm::rmvnorm(n=n, mean=mu, sigma=Sigma))
    }
}




GSP.m2l <- function (theta, z, cov.fun, dists, X=NULL, transform.fun) {
    ## if X is NULL then the mean is assumed to be zero.

    Sigma <- cov.fun(dists, theta, transform.fun)
    
    chol.Sigma <- tryCatch(chol(Sigma), error = function(e) e)

    if (inherits(chol.Sigma, "error")) {

        Inf

    } else {

        n   <- length(z)

        Lz <- backsolve(chol.Sigma, matrix(z, nrow=n), transpose = TRUE)
        
        if (!is.null(X)) {
            
            LX <- backsolve(chol.Sigma, X, transpose = TRUE)
            
            eps <- z - X %*% .lm.fit(LX, Lz)$coef
            
            Lz <- backsolve(chol.Sigma, matrix(eps, nrow=n), transpose = TRUE)
        }

        res <- 2*sum(log(diag(chol.Sigma))) + n * log(2 * pi) + colSums(Lz^2)

        return(res)
    }
}




    

GSP.beta.hat <- function (theta, z, cov.fun, dists, X=NULL) {
    ## if X is NULL then the mean is assumed to be zero.

    if (is.null(X)) {
        
        NA
        
    } else {

        Sigma <- cov.fun(dists, theta)

        chol.Sigma <- tryCatch(chol(Sigma), error = function(e) e)

        if (inherits(chol.Sigma, "error")) {
            
            Inf
            
        } else {

            LX <- backsolve(chol.Sigma, X, transpose = TRUE)
            Lz <- backsolve(chol.Sigma, matrix(z, nrow=length(z)), transpose = TRUE)

            .lm.fit(LX, Lz)$coef
        }
    }
}


GSP.beta.hat.cov <- function (theta, cov.fun, dists, X=NULL) {
    ## if X is NULL then the mean is assumed to be zero.

    if (is.null(X)) {
        
        NA
        
    } else {

        Sigma <- cov.fun(dists, theta)

        chol.Sigma <- tryCatch(chol(Sigma), error = function(e) e)

        if (inherits(chol.Sigma, "error")) {
            
            Inf
            
        } else {

            LX <- backsolve(chol.Sigma, X, transpose = TRUE)

            solve(crossprod(LX))
        }
    }
}



GSP.mle <- function (theta, z, X=NULL, cov.fun, dists, hessian=FALSE) {
    ## is cov.fun is missing using an independent normal model

    if (!missing(cov.fun)) {
        
        
        out <- optim(log(theta), GSP.m2l, method="BFGS",
                     z=z, cov.fun=cov.fun, dists=dists, X=X,
                     transform.fun=exp, hessian=hessian)

        theta.hat <- exp(out$par)

        if (!is.null(X)) {

            beta.hat <- GSP.beta.hat(theta.hat, z, cov.fun, dists, X)
        }
        
        m2l <- out$value

        neg.hessian <- 0.5 * out$hessian ## remember we calculation -2 * loglik!
        
    } else {

        ols <- lm.fit(X, z)

        theta.hat <- mean(ols$residuals^2)        

        beta.hat <- ols$coef

        m2l <- -2*sum(dnorm(z, drop(X %*% beta.hat), sqrt(theta.hat), log=TRUE))

        if (hessian) {
            
            neg.hessian <- length(z)/(2*theta.hat^2)
        } else {
            
            neg.hessian <- NULL
        }
    }

    if (is.null(X)) {

        list(theta.hat = theta.hat,
             AIC       = m2l+2*length(theta.hat),
             X         = X,
             Hessian   = neg.hessian)
        
    } else {

        list(theta.hat = theta.hat,
             beta.hat  = beta.hat,
             AIC       = m2l+2*(length(theta.hat)+ncol(X)),
             X         = X,
             Hessian   = neg.hessian)
    }
}







GSP.pred <- function (theta, z, X=NULL, cov.fun, dists, beta=NULL,
                      dists.between, dists.pred, X.pred=NULL,
                      calc.var=TRUE, calc.cov=FALSE, 
                      without.nugget=FALSE) {
    
    Sigma <- cov.fun(dists, theta)
    
    omega <- theta    
    if (without.nugget) {
        
        omega[length(omega)] <- 0
    }

    BSigma <- cov.fun(dists.between, omega)
        
    if (is.null(X)) {

        inner <- solve(Sigma, z)
        pred  <- drop(crossprod(BSigma, inner))
        
    } else {

        inner <- solve(Sigma, z - X %*% beta)
        pred  <- drop(X.pred %*% beta + crossprod(BSigma, inner))
    }

    if (calc.cov) {
        
        PSigma   <- cov.fun(dists.pred, omega)
        pred.cov <- PSigma - crossprod(solve(Sigma, BSigma), BSigma)

        if (calc.var) {

            pred.var <- diag(pred.cov)
            
        } else {

            pred.var <- NA
        }
        
    } else {
        
        pred.cov <- NA
        
        if (calc.var) {

            pred.var <- cov.fun(0, omega) - colSums(solve(Sigma, BSigma) * BSigma)
            
        } else {

            pred.var <- NA
        }
    }
    
    list(pred     = pred,
         pred.var = pred.var,
         pred.cov = pred.cov)
}






GSP.image <- function (z, xs, ys) {

    if (missing(xs) && missing(ys)) {

        image.plot(matrix(z, sqrt(length(z))))
        
    } else {

        image.plot(xs, ys, matrix(z, length(xs)))
    }
}

