##============================================================
## All the parameter update functions that utilizes Gibbs
## Sampling.
## 1. alpha, the preferential sampling parameter;
## 2. tau2, the measurement error;
## 3. sig2_W, the variance term for the GSP in Y;
##============================================================

dinvgamma <- function (x, shape, rate=1, scale=1/rate, log=FALSE)
{
  if (log)
    -rate/x + shape*log(rate) - lgamma(shape) - (shape+1)*log(x)
  else
    exp(-rate/x) * rate^shape / (gamma(shape) * x^(shape+1))
}

rinvgamma <- function (n, shape, rate=1, scale=1/rate)
{
  1 / rgamma(n=n, shape=shape, rate=rate)
}


rmvnorm.cond.precision <- function (m, P, eigen.P) {
## ======================================================================
## Purpose: Sample from N( P^{-1} m, P^{-1} )
## Assumes: m and P are real-valued.
##          The length of m is equal to the dimension of the square matrix P.
## ======================================================================

  if (missing(eigen.P)) {
    
    eigen.P <- eigen(P, symmetric=TRUE)
  }
  
  ev <- eigen.P$values

  if (!all(ev >= -sqrt(.Machine$double.eps) * abs(ev[1]))) {
    warning("P is numerically not positive definite")
  }
  
  A <- t(eigen.P$vectors) * sqrt(1/ev)
  
  drop(crossprod(A, A %*% m + rnorm(length(m))))
}


z_xy <- function(x,y, parents, omega2, delta){
  # Return the density at one location
  delta*sum(dnorm(sqrt((rep(x,parents$n)-parents$x)^2 + (rep(y,parents$n)-parents$y)^2), sd = sqrt(omega2)))
}

lambda <- function(S, parents, delta, omega2){
  # Return the density at all locations in S
  sapply(c(1:nrow(S)), function(i) z_xy(S[i,1], S[i,2], parents, omega2, delta))
}






# Total: n observations on 20x20 grid
# W: nx1 vector. the value of W at the location of each observation.
# Z: nx1 vector. observations
# Z_cell: nx1 vector, indicating which cell each observation belongs to


update_alpha <- function(log.int, Z, obs_cell, tau2, W, mu_W,
                         prior_mean = 0, prior_var = 50){
    p.inv <- solve(tau2*crossprod(log.int[obs_cell], log.int[obs_cell]) + 1/prior_var)
    q <- tau2*crossprod(log.int[obs_cell], (Z - W[obs_cell] - mu_W)) + prior_mean/prior_var
    new_alpha <- rnorm(1, mean = p.inv*q, sd = sqrt(p.inv))
    new_alpha
}


update_tau2 <- function(Z, log.int, alpha, W, mu_W, # parameters that are being updated
                        obs_cell, # fixed values 
                        prior_shape = 1, prior_rate = 1){ # priors
  Y_obs <- alpha*log.int[obs_cell] + W[obs_cell] + mu_W
  post_shape <- prior_shape + length(Z)/2
  post_rate <- prior_rate + sum((Z - Y_obs)^2)/2
  rinvgamma(1, shape = post_shape, rate = post_rate)
}


update_sig2_W <- function(W, phi_W, 
                          D_grid, cov.power,
                          prior_shape = 1, prior_rate = 1){ 
  t1 <- solve(exp(-(D_grid/phi_W)^cov.power) + diag(rep(1e-5, nrow(D_grid))), W)
  post_shape <- prior_shape + length(W)/2
  post_rate <- prior_rate + crossprod(t1, W)/2
  rinvgamma(1, shape = post_shape, rate = post_rate)
}


update_mu_w <- function(Z, log.int, alpha, W, tau2, obs_cell,
                        prior_mean = 0, prior_var = 25){
  temp <- Z - W[obs_cell] - alpha*log.int[obs_cell]
  p.inv <- solve(length(temp)/tau2 + 1/prior_var)
  q <- prior_mean/prior_var + sum(temp)/tau2
  rnorm(1, p.inv*q, sqrt(p.inv))
}



f_log_phi_w <- function(W, log_phi_W, sig2_W,
                    D, cov.power.w){
  Sigma <- sig2_W*exp(-(D/exp(log_phi_W))^cov.power.w)
  t1 <- dmvnorm(W, rep(0, length(W)), Sigma, log = T)
  t2 <- dgamma(exp(log_phi_W), 0.33*10, 10, log = T)
  t1 + t2
}


update_phi_w <- function(W, curr_phi_W, sig2_W,
                         D_grid, cov.power.w){
  D <- D_grid
  curr_log_phi_W <- log(curr_phi_W)
  new_log_phi_W <- rnorm(1, curr_log_phi_W, 0.1) # try a few sd values to get ~68% accept.prob
  
  p_curr <- f_log_phi_w(W, curr_log_phi_W, sig2_W,
                        D, cov.power.w)# + curr_log_phi_W
  p_new <- f_log_phi_w(W, new_log_phi_W, sig2_W,
                       D, cov.power.w)#  + new_log_phi_W
  if(log(runif(1)) < (p_new - p_curr)){
    exp(new_log_phi_W)
  }else{
    curr_phi_W
  }
}


library(mvtnorm)
update_W <- function(curr_W, alpha, Z, log.int,
                       mu_W, tau2,
                       D_grid, obs_cell, blocks, cov.power,
                       sum_pts, Sig,
                       prior_mean = 0, prior_var = 10){
  new_W <- curr_W
  sum_obs <- rep(0, nrow(D_grid))
  for(k in 1:nrow(D_grid)){
    id <- which(obs_cell == k)
    if(length(id) != 0){
      sum_obs[k] <- sum(Z[id] - alpha*log.int[k] - mu_W)
    }
  }
  for(j in c(1:ncol(blocks))){
    ks <- blocks[,j]
    
    A <- Sig[-ks,]
    B <- A[,ks]
    C <- solve(A[,-ks], B)

    cond.mean.k <- crossprod(C, new_W[-ks])
    cond.cov.k  <- Sig[ks,][,ks] - crossprod(B, C)

    cond.cov.k.inv <- solve(cond.cov.k)

    P <- cond.cov.k.inv +  diag(sum_pts[ks]/tau2, length(ks)) #+ diag(1/prior_var, length(ks))
    m <- cond.cov.k.inv %*% cond.mean.k + sum_obs[ks]/tau2

    new_W[ks] <- rmvnorm.cond.precision(m, P)
  }
  new_W
}


f_delta <- function(delta, omega2, W, mu_W, alpha, tau2, Z, S, parents,
                    N_pts, obs_cell, cell_size, 
                    prior_mean, prior_var){
  log.int <- log(lambda(S, parents, delta, omega2))
  # From the observations:
  t1 <- sum(dnorm(Z, alpha*log.int[obs_cell] + mu_W + W[obs_cell], sqrt(tau2), log = T))
  # From the locations (points):
  t2 <- sum(dpois(N_pts, cell_size*exp(log.int), log = T))
  # Prior
  t3 <- dnorm(log(delta), prior_mean, prior_var, log = T)
  t1 + t2 + t3
}
  
update_delta <- function(delta, omega2, W, mu_W, alpha, tau2, Z, S, parents,
                         N_pts, obs_cell, cell_size, prior_mean, prior_var){
  log_delta <- log(delta)
  new_log_delta <- rnorm(1, log_delta, 0.035)
  p_curr <- f_delta(delta, omega2, W, mu_W, alpha, tau2, Z, S, parents, N_pts, 
                    obs_cell, cell_size, prior_mean, prior_var) + log_delta
  p_new <- f_delta(exp(new_log_delta), omega2, W, mu_W, alpha, tau2, Z, S, parents, N_pts, 
                    obs_cell, cell_size, prior_mean, prior_var) + new_log_delta
  if(log(runif(1)) < (p_new - p_curr)){
    exp(new_log_delta)
  }else{
    delta
  }
}

f_omega2 <- function(omega2, delta, W, mu_W, alpha, tau2, Z, S, parents,
                    N_pts, obs_cell, cell_size, 
                    prior_mean, prior_var){
  log.int <- log(lambda(S, parents, delta, omega2))
  # From the observations:
  t1 <- sum(dnorm(Z, alpha*log.int[obs_cell] + mu_W + W[obs_cell], sqrt(tau2), log = T))
  # From the locations (points):
  t2 <- sum(dpois(N_pts, cell_size*exp(log.int), log = T))
  # Prior
  t3 <- dnorm(log(omega2), prior_mean, prior_var, log = T)
  t1 + t2 + t3
}


update_omega2 <- function(omega2, delta, W, mu_W, alpha, tau2, Z, S, parents,
                    N_pts, obs_cell, cell_size, 
                    prior_mean, prior_var){
  log_omega2 <- log(omega2)
  new_log_omega2 <- rnorm(1, log_omega2, 0.075)
  p_curr <- f_omega2(omega2, delta, W, mu_W, alpha, tau2, Z, S, parents, N_pts, 
                    obs_cell, cell_size, prior_mean, prior_var) + log_omega2
  p_new <- f_omega2(exp(new_log_omega2), delta, W, mu_W, alpha, tau2, Z, S, parents, N_pts, 
                    obs_cell, cell_size, prior_mean, prior_var) + new_log_omega2
  if(log(runif(1)) < (p_new - p_curr)){
    exp(new_log_omega2)
  }else{
    omega2
  }
}


