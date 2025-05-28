##============================================================
## All the parameter update functions that utilizes Gibbs
## Sampling.
## 1. alpha, the preferential sampling parameter;
## 2. tau2, the measurement error;
## 3. sig2_W, the variance term for the GSP in Y;
## 4. sig2_eta, the variance term for the GSP in PP.
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



# Total: n observations on 20x20 grid
# eta: 400x1 vector
# W: nx1 vector. the value of W at the location of each observation.
# Z: nx1 vector. observations
# Z_cell: nx1 vector, indicating which cell each observation belongs to

## Condition on W, (Z-W-mu_W) ~ N(alpha*(eta + mu_eta), tau2)
update_alpha <- function(eta, Z, obs_cell, tau2, W, mu_W, mu_eta,
                         prior_mean = 0, prior_var = 100){
    p.inv <- solve(tau2*crossprod((eta[obs_cell] + mu_eta), (eta[obs_cell] + mu_eta)) + 1/prior_var)
    #p.inv <- solve(tau2*crossprod(eta[obs_cell], eta[obs_cell]) + 1/prior_var)
    q <- tau2*crossprod((eta[obs_cell] + mu_eta), (Z - W[obs_cell] - mu_W)) + prior_mean/prior_var
    #q <- tau2*crossprod(eta[obs_cell], (Z - W[obs_cell] - mu_W)) + prior_mean/prior_var
    new_alpha <- rnorm(1, mean = p.inv*q, sd = sqrt(p.inv))
    new_alpha
}


update_tau2 <- function(Z, eta, alpha, W, mu_eta, mu_W, # parameters that are being updated
                        obs_cell, # fixed values 
                        prior_shape = 0.01, prior_rate = 0.01){ # priors
  #Y_obs <- alpha*(eta[obs_cell]) + W[obs_cell] + mu_W
  Y_obs <- alpha*(eta[obs_cell] + mu_eta) + W[obs_cell] + mu_W
  post_shape <- prior_shape + length(Z)/2
  post_rate <- prior_rate + sum((Z - Y_obs)^2)/2
  rinvgamma(1, shape = post_shape, rate = post_rate)
}


update_sig2_W <- function(W, phi_W, 
                          D_grid, cov.power.w,
                          prior_shape = 0.01, prior_rate = 0.01){ 
  t1 <- solve(exp(-(D_grid/phi_W)^cov.power.w) + diag(rep(1e-5, nrow(D_grid))), W)
  post_shape <- prior_shape + length(W)/2
  post_rate <- prior_rate + crossprod(t1, W)/2
  rinvgamma(1, shape = post_shape, rate = post_rate)
}

update_sig2_eta <- function(eta, phi_eta,
                            D_grid, cov.power.eta,
                            prior_shape = 0.01, prior_rate = 0.01){
  t1 <- solve(exp(-(D_grid/phi_eta)^cov.power.eta)+ diag(rep(1e-5, nrow(D_grid))), eta)
  post_shape <- prior_shape + length(eta)/2
  post_rate <- prior_rate + crossprod(t1, eta)/2
  rinvgamma(1, shape = post_shape, rate = post_rate)
}

update_mu_w <- function(Z, mu_eta, eta, alpha, W, tau2, obs_cell,
                        prior_mean = 0, prior_var = 100){
  #temp <- Z - W[obs_cell] - alpha*eta[obs_cell]
  temp <- Z - W[obs_cell] - alpha*(eta[obs_cell] + mu_eta)
  p.inv <- solve(length(temp)/tau2 + 1/prior_var)
  q <- prior_mean/prior_var + sum(temp)/tau2
  rnorm(1, p.inv*q, sqrt(p.inv))
}


f_mu_eta <- function(mu_eta, transformed_obs, alpha, eta, tau2, N_pts, cell_area,
                          obs_cell, prior_mean, prior_var){
  #t3 <- sum(dnorm(Z,))
  # likelihood from the points
  t1 <- sum(dpois(N_pts, exp(mu_eta + eta)*cell_area, log = T))
  # prior
  t2 <- dnorm(mu_eta, prior_mean, sqrt(prior_var), log = T)
  t3 <- sum(dnorm(transformed_obs, mu_eta, sqrt(tau2/abs(alpha)), log = T))
  #t3 <- sum(dnorm(transform.obs, alpha(mu_eta + eta) + mu_W + W, sd = sqrt(tau2), log = T))
  t1 + t2 + t3
}

update_mu_eta <- function(curr_mu_eta, Z, eta, alpha, mu_W, W, tau2, obs_cell, N_pts, cell_size,
                          prior_mean = 0, prior_var = 100){
  new_mu_eta <- rnorm(1, curr_mu_eta, 0.08)
  trans_obs <- (Z - W[obs_cell] - mu_W - alpha*(eta[obs_cell]))/alpha
  #trans_obs <- Z
  p_curr <- f_mu_eta(curr_mu_eta, trans_obs, alpha, eta, tau2, N_pts, cell_size, obs_cell, prior_mean, prior_var)
  p_new <- f_mu_eta(new_mu_eta, trans_obs, alpha, eta, tau2, N_pts, cell_size, obs_cell, prior_mean, prior_var)
  if(log(runif(1)) < (p_new - p_curr)){
    return(new_mu_eta)
  }else{
    return(curr_mu_eta)
  }
}



f_log_phi_w <- function(W, log_phi_W, sig2_W,
                    D, cov.power, prior_mean, prior_var){
  Sigma <- sig2_W*exp(-(D/exp(log_phi_W))^cov.power)
  t1 <- dmvnorm(W, rep(0, length(W)), Sigma, log = T)
  t2 <- dgamma(exp(log_phi_W), 0.5*5, 5, log = T)
  # t1 + t2 + t3
  t1 + t2
}

update_phi_w <- function(W, curr_phi_W, sig2_W,
                         D_grid, cov.power){
  D <- D_grid
  curr_log_phi_W <- log(curr_phi_W)
  new_log_phi_W <- rnorm(1, curr_log_phi_W, 0.1) # try a few sd values to get ~68% accept.prob
  
  p_curr <- f_log_phi_w(W, curr_log_phi_W, sig2_W,
                        D, cov.power)# + curr_log_phi_W
  p_new <- f_log_phi_w(W, new_log_phi_W, sig2_W,
                       D, cov.power)#  + new_log_phi_W
  if(log(runif(1)) < (p_new - p_curr)){
    exp(new_log_phi_W)
  }else{
    curr_phi_W
  }
}

f_log_phi_eta <- function(eta, log_phi_eta, sig2_eta,
                    D, cov.power){
  Sigma <- sig2_eta*exp(-(D/exp(log_phi_eta))^cov.power)
  t1 <- dmvnorm(eta, rep(0, length(eta)), Sigma, log = T)
  t2 <- dgamma(exp(log_phi_eta), 0.5*5, 5, log = T)  
  t1 + t2
}


update_phi_eta <- function(eta, curr_phi_eta, sig2_eta,
                         D_grid, cov.power){

  curr_log_phi_eta <- log(curr_phi_eta)
  new_log_phi_eta <- rnorm(1, curr_log_phi_eta, 0.025) # try a few sd values to get ~68% accept.prob
  p_curr <- f_log_phi_eta(eta, curr_log_phi_eta, sig2_eta,
                        D_grid, cov.power)# - (curr_log_phi_eta)
  p_new <- f_log_phi_eta(eta, new_log_phi_eta, sig2_eta,
                       D_grid, cov.power)#  - (new_log_phi_eta)
  if(log(runif(1)) < (p_new - p_curr)){
    exp(new_log_phi_eta)
  }else{
    curr_phi_eta
  }
}

f_eta_j <- function(j, eta_j, alpha, Z, W, D, tau2, mu_eta, mu_W,
                    N_pts, prior_mean, prior_var, cell_area, obs_cell, cond_mean, cond_var){
  # Likelihood from the points
  t1 <- dpois(N_pts[j], exp(mu_eta + eta_j)*cell_area, log = T)
  t2 <- 0
  id <- which(obs_cell == j) # check how many observations in this box
  # Likelihood from the observations
  if(length(id) > 0){
    t2 <- sum(dnorm(Z[id], (alpha*(eta_j + mu_eta) + mu_W + W[j]), sqrt(tau2), log = T))
    #t2 <- sum(dnorm(Z[id], (alpha*(eta_j) + mu_W + W[j]), sqrt(tau2), log = T))
  }
  
  t3 <- dnorm(eta_j, cond_mean, sqrt(cond_var), log = T)
  t4 <- dnorm(eta_j, prior_mean, sqrt(prior_var), log = T)
  t1 + t2 + t3 + t4
}

update_eta <- function(curr_eta, alpha, Z, W, mu_eta, mu_W,
                       Sig_eta, tau2,
                       D_grid, N_pts, obs_cell,
                       nbs, D_nbs, cell_size, cov.power.eta,
                       prior_mean = 0, prior_var = 100){
  new_eta <- curr_eta

  for(j in 1:length(new_eta)){
    curr_eta_j <- new_eta[j]
    new_eta_j <- rnorm(1, curr_eta_j, 0.5)
    Sigma_22 <- Sig_eta[nbs[[j]],][,nbs[[j]]]# + diag(rep(1e-5, length(nbs[[j]])))
    t1 <- Sig_eta[nbs[[j]],][,j]
    t2 <- solve(Sigma_22, t1) # Sigma22^-1%*%Sigma_12
    cond_mean <-  crossprod(t2, new_eta[nbs[[j]]])
    cond_var <- Sig_eta[j,j] - crossprod(t2, t1)
    p_curr <- f_eta_j(j, curr_eta_j, alpha, Z, W, D_grid, tau2, mu_eta, mu_W,
                     N_pts, prior_mean, prior_var, cell_size[j], obs_cell, cond_mean, cond_var)
    p_new <- f_eta_j(j, new_eta_j, alpha, Z, W, D_grid, tau2, mu_eta, mu_W,
                     N_pts, prior_mean, prior_var, cell_size[j], obs_cell, cond_mean, cond_var)
    if(log(runif(1)) < p_new - p_curr){
      new_eta[j] <- new_eta_j
    }
  }
  new_eta
} 


library(mvtnorm)
update_w_block_Gibbs_new <- function(curr_W, alpha, Z, eta,
                       mu_W, mu_eta, tau2,
                       D_grid, obs_cell, blocks, cov.power,
                       sum_pts, Sig,
                       prior_mean = 0, prior_var = 100){
  new_W <- curr_W
  sum_obs <- rep(0, nrow(D_grid))
  for(k in 1:nrow(D_grid)){
    id <- which(obs_cell == k)
    if(length(id) != 0){
      #sum_obs[k] <- sum(Z[id] - alpha*eta[k] - mu_W)
      sum_obs[k] <- sum(Z[id] - alpha*(eta[k]+mu_eta) - mu_W)
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

    P <- cond.cov.k.inv +  diag(sum_pts[ks]/tau2, length(ks)) # + diag(1/prior_var, length(ks))
    m <- cond.cov.k.inv %*% cond.mean.k + sum_obs[ks]/tau2

    new_W[ks] <- rmvnorm.cond.precision(m, P)
  }
  new_W
}



text_on_cell <- function(id, words = id, size = 0.5, col = 'white'){
  text(S[id,1], S[id,2], words, cex = size, col = col)  
}



