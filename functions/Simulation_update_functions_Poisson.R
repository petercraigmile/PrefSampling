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
# sampling effor parameter, to control the observation scale
# assumed to be known, e.g. 0.1

f_alpha_mu_w <- function(alpha_mu_w, eta, mu_eta, W, Z, obs_cell, kappa, prior_mean, prior_var){
  rate <- exp(W[obs_cell] + alpha_mu_w[2] + alpha_mu_w[1]*(eta[obs_cell] + mu_eta))
  - (sum(dpois(Z, kappa*rate, log = T)) + dmvnorm(alpha_mu_w, rep(prior_mean, 2), diag(rep(prior_var, 2)), log = T))
}

update_alpha_mu_w <- function(alpha, mu_W, eta, mu_eta, W, Z, obs_cell, kappa, prior_mean, prior_var){
  m <- optim(c(0, 0), f_alpha_mu_w, 
             eta = eta, mu_eta = mu_eta, W = W, 
             Z = Z, obs_cell = obs_cell, 
             kappa = kappa, prior_mean = prior_mean, 
             prior_var = prior_var)$par
  
  curr_alpha_mu_w <- c(alpha, mu_W)
  # alpha_mu_w <- new_alpha_mu_W
  rate_curr <- exp(W[obs_cell] + curr_alpha_mu_w[2] + curr_alpha_mu_w[1]*(eta[obs_cell] + mu_eta))
  var_curr <- solve(matrix(c(sum((mu_eta + eta[obs_cell])^2*rate_curr),
                             sum((mu_eta + eta[obs_cell])*rate_curr), 
                             sum((mu_eta + eta[obs_cell])*rate_curr), 
                             sum(rate_curr)), 2, byrow = T))
  
  new_alpha_mu_W <- rmvnorm(1, m, 2*var_curr)
  # new_alpha_mu_W
  rate_new <- exp(W[obs_cell] + new_alpha_mu_W[2] + new_alpha_mu_W[1]*(eta[obs_cell] + mu_eta))
  var_new <- solve(matrix(c(sum((mu_eta + eta[obs_cell])^2*rate_new), 
                            sum((mu_eta + eta[obs_cell])*rate_new), 
                            sum((mu_eta + eta[obs_cell])*rate_new), 
                            sum(rate_new)), 2, byrow = T))
  # new_alpha_mu_W <- rmvnorm(1, m, diag(c(0.001, 0.001)))
  
  p_curr <- f_alpha_mu_w(curr_alpha_mu_w, eta, mu_eta, W, Z, obs_cell, kappa, prior_mean, prior_var)*(-1)# + dmvnorm(new_alpha_mu_W, m, 1.25*var_curr, log = T) 
  p_new <- f_alpha_mu_w(new_alpha_mu_W, eta, mu_eta, W, Z, obs_cell, kappa, prior_mean, prior_var)*(-1)# + dmvnorm(curr_alpha_mu_w, m, 1.25*var_new, log = T) 
    if(log(runif(1)) < (p_new - p_curr)){
      new_alpha_mu_W
    }else{
      curr_alpha_mu_w
    }
}




f_three <- function(three_var, Z, eta, W,
                    obs_cell, cell_size, kappa, 
                    prior_mean, prior_var){
  # likelihood from the observations
  alpha <- three_var[1]
  mu_W <- three_var[2]
  mu_eta <- three_var[3]
  rate <- exp(W[obs_cell] + mu_W + alpha*(eta[obs_cell] + mu_eta))
  t1 <- sum(dpois(Z, kappa*rate, log = T))
  
  # likelihood from the points
  t2 <- sum(mu_eta + eta[obs_cell]) - sum(cell_size*exp(eta + mu_eta))
  # t2 <- dpois(N_pts, cell_size*exp(eta + mu_eta), log = T)
  # prior
  t3 <- dmvnorm(three_var, rep(prior_mean, 3), diag(rep(prior_var, 3)), log = T)
  -(t1 + t2 + t3)
}


update_three <- function(alpha, mu_W, mu_eta, eta, W, Z, obs_cell,cell_size, kappa, prior_mean, prior_var){
   m <- optim(c(0, 0, 0), f_three, 
             Z = Z, eta = eta, W = W, 
             obs_cell = obs_cell, 
             cell_size = cell_size,
             kappa = kappa, prior_mean = prior_mean, 
             prior_var = prior_var)$par
  curr_var <- c(alpha, mu_W, mu_eta)
  #curr_var <- new_var
  # alpha_mu_w <- new_alpha_mu_W
  rate_curr <- exp(W[obs_cell] + curr_var[2] + curr_var[1]*(eta[obs_cell] + curr_var[3]))
  var_curr <- solve(matrix(c(sum((mu_eta + eta[obs_cell])^2*rate_curr), # 1,1
                             sum((mu_eta + eta[obs_cell])*rate_curr), # 1,2 
                             sum(rate_curr) + sum(alpha*(mu_eta + eta[obs_cell])*rate_curr) - sum(Z), # 1,3
                             sum((mu_eta + eta[obs_cell])*rate_curr), # 2,1
                             sum(rate_curr), # 2,2
                             alpha*sum(rate_curr), # 2,3
                             sum(rate_curr) + sum(alpha*(mu_eta + eta[obs_cell])*rate_curr) - sum(Z), # 3,1
                             alpha*sum(rate_curr), # 3,2 
                             alpha^2*sum(rate_curr) + sum(cell_size*exp(eta + mu_eta))), 3, byrow = T)) # 3,3
  #var_curr <- diag(1e-4, 3)
  new_var <- rmvnorm(1, m, 2*var_curr)
  # new_var
  p_curr <- f_three(curr_var, Z, eta, W, obs_cell, cell_size, kappa, prior_mean, prior_var)*(-1)# + dmvnorm(new_alpha_mu_W, m, 1.25*var_curr, log = T) 
  p_new <- f_three(new_var, Z, eta, W, obs_cell, cell_size, kappa, prior_mean, prior_var)*(-1)# + dmvnorm(curr_alpha_mu_w, m, 1.25*var_new, log = T) 
    if(log(runif(1)) < (p_new - p_curr)){
      new_var
    }else{
      curr_var
    }
}


f_mu <- function(mu, alpha, Z, eta, W,
                    obs_cell, cell_size, kappa, 
                    prior_mean, prior_var){
  # likelihood from the observations
  mu_W <- mu[1]
  mu_eta <- mu[2]
  rate <- exp(W[obs_cell] + mu_W + alpha*(eta[obs_cell] + mu_eta))
  t1 <- sum(dpois(Z, kappa*rate, log = T))
  
  # likelihood from the points
  t2 <- sum(mu_eta + eta[obs_cell]) - sum(cell_size*exp(eta + mu_eta))
  # t2 <- dpois(N_pts, cell_size*exp(eta + mu_eta), log = T)
  # prior
  t3 <- dmvnorm(mu, rep(prior_mean, 2), diag(rep(prior_var, 2)), log = T)
  -(t1 + t2 + t3)
}


update_mu <- function(mu_W, mu_eta, alpha, eta, W, Z, obs_cell,cell_size, kappa, prior_mean, prior_var){
  m <- optim(c(0, 0), f_mu, 
             alpha = alpha, Z = Z, eta = eta, W = W, 
             obs_cell = obs_cell, 
             cell_size = cell_size,
             kappa = kappa, prior_mean = prior_mean, 
             prior_var = prior_var)$par
  curr_var <- c(mu_W, mu_eta)
  rate_curr <- exp(W[obs_cell] + curr_var[1] + alpha*(eta[obs_cell] + curr_var[2]))
  var_curr <- solve(matrix(c(sum(rate_curr), # 2,2
                             alpha*sum(rate_curr), # 2,3
                             alpha*sum(rate_curr), # 3,2 
                             alpha^2*sum(rate_curr) + sum(cell_size*exp(eta + mu_eta))), 2, byrow = T)) # 3,3
  new_var <- rmvnorm(1, m, 2*var_curr)
  p_curr <- f_mu(curr_var, alpha, Z, eta, W, obs_cell, cell_size, kappa, prior_mean, prior_var)*(-1)# + dmvnorm(new_alpha_mu_W, m, 1.25*var_curr, log = T) 
  p_new <- f_mu(new_var, alpha, Z, eta, W, obs_cell, cell_size, kappa, prior_mean, prior_var)*(-1)# + dmvnorm(curr_alpha_mu_w, m, 1.25*var_new, log = T) 
    if(log(runif(1)) < (p_new - p_curr)){
      new_var
    }else{
      curr_var
    }
}


f_alpha_mu_eta <- function(alpha_mu_eta, mu_W, Z, eta, W,
                    obs_cell, cell_size, kappa, 
                    prior_mean, prior_var){
  alpha <- alpha_mu_eta[1]
  mu_eta <- alpha_mu_eta[2]
  rate <- exp(W[obs_cell] + mu_W + alpha*(eta[obs_cell] + mu_eta))
  t1 <- sum(dpois(Z, kappa*rate, log = T))
  
  # likelihood from the points
  # t2 <- sum(mu_eta + eta[obs_cell]) - sum(cell_size*exp(eta + mu_eta))
  t2 <- sum(dpois(N_pts, cell_size*exp(eta + mu_eta), log = T))
  # prior
  t3 <- dmvnorm(alpha_mu_eta, rep(prior_mean, 2), diag(rep(prior_var, 2)), log = T)
  -(t1 + t2 + t3)
}

update_alpha_mu_eta <- function(alpha, mu_eta, mu_W, eta, W, Z, 
                                obs_cell, cell_size, kappa, prior_mean, prior_var){
  m <- optim(c(0, 0), f_alpha_mu_eta, 
             mu_W = mu_W, Z = Z, eta = eta, W = W, 
             obs_cell = obs_cell, 
             cell_size = cell_size,
             kappa = kappa, prior_mean = prior_mean, 
             prior_var = prior_var)$par
  curr_var <- c(alpha, mu_eta)
  rate_curr <- exp(W[obs_cell] + mu_W + curr_var[1]*(eta[obs_cell] + curr_var[2]))
  var_curr <- solve(matrix(c(sum((mu_eta + eta[obs_cell])^2*rate_curr), # 1,1
                             sum(rate_curr) + sum(alpha*(mu_eta + eta[obs_cell])*rate_curr) - sum(Z), # 1,3
                             sum(rate_curr) + sum(alpha*(mu_eta + eta[obs_cell])*rate_curr) - sum(Z), # 3,1
                             alpha^2*sum(rate_curr) + sum(cell_size*exp(eta + mu_eta))), 2, byrow = T)) # 3,3
  new_var <- rmvnorm(1, m, 2*var_curr)
  p_curr <- f_alpha_mu_eta(curr_var, mu_W, Z, eta, W, obs_cell, cell_size, kappa, prior_mean, prior_var)*(-1)# + dmvnorm(new_alpha_mu_W, m, 1.25*var_curr, log = T) 
  p_new <- f_alpha_mu_eta(new_var, mu_W, Z, eta, W, obs_cell, cell_size, kappa, prior_mean, prior_var)*(-1)# + dmvnorm(curr_alpha_mu_w, m, 1.25*var_new, log = T) 
    if(log(runif(1)) < (p_new - p_curr)){
      new_var
    }else{
      curr_var
    }
}




f_alpha <- function(alpha, eta, mu_eta, W, mu_W, Z, obs_cell, kappa, prior_mean, prior_var){
  rate <- exp(W[obs_cell] + mu_W + alpha*eta[obs_cell])
  sum(dpois(Z, kappa*rate, log = T)) + dnorm(alpha, prior_mean, sqrt(prior_var), log = T)
}

update_alpha <- function(curr_a, eta, Z, obs_cell, W, mu_W, mu_eta, kappa,
                         prior_mean = 0, prior_var = 50){
    new_a <- rnorm(1, curr_a, 0.00075)
    p_curr <- f_alpha(curr_a, eta, mu_eta, W, mu_W, Z, obs_cell, kappa, prior_mean, prior_var)
    p_new <- f_alpha(new_a, eta, mu_eta, W, mu_W, Z, obs_cell, kappa, prior_mean, prior_var)
    if(log(runif(1)) < (p_new - p_curr)){
      new_a
    }else{
      curr_a
    }
}

f_mu_w <- function(mu_W, W, alpha, eta, mu_eta, Z, obs_cell, kappa, prior_mean, prior_var){
  rate <- exp(W[obs_cell] + mu_W + alpha*eta[obs_cell])
  sum(dpois(Z, kappa*rate, log = T)) + dnorm(mu_W, prior_mean, sqrt(prior_var), log = T)
}

update_mu_w <- function(curr_mu_W, W, alpha, eta, mu_eta, Z, obs_cell, kappa,
                        prior_mean = 0, prior_var = 25){
  new_mu_W <- rnorm(1, curr_mu_W, 0.005)
  p_curr <- f_mu_w(curr_mu_W, W, alpha, eta, mu_eta, Z, obs_cell, kappa, prior_mean, prior_var)
  p_new <- f_mu_w(new_mu_W, W, alpha, eta, mu_eta, Z, obs_cell, kappa, prior_mean, prior_var)
  if(log(runif(1)) < (p_new - p_curr)){
    new_mu_W
  }else{
    curr_mu_W
  }
}


update_mu_w_Gibbs <- function(W, alpha, eta, mu_eta, Z, obs_cell, kappa, prior_a = 1, prior_b = 1){
  new_a <- prior_a + sum(Z)
  new_b <- prior_b + sum(exp(W[obs_cell] + alpha*eta[obs_cell]))
  log(rgamma(1, shape = new_a, rate = new_b))
}






f_mu_eta <- function(mu_eta, Z, alpha, eta, mu_W, W, N_pts, 
                          obs_cell, cell_size, kappa, prior_mean, prior_var){
  # likelihood from the observations
  #rate <- exp(W[obs_cell] + mu_W + alpha*(eta[obs_cell] + mu_eta))
  #t1 <- sum(dpois(Z, kappa*rate, log = T))
  # likelihood from the points
  # t2 <- sum(mu_eta + eta[obs_cell]) - mean(exp(eta + mu_eta))
  t2 <- sum(dpois(N_pts, cell_size*exp(eta + mu_eta), log = T))
  # prior
  t3 <- dnorm(mu_eta, prior_mean, sqrt(prior_var), log = T)
  t2 + t3
}

update_mu_eta <- function(curr_mu_eta, Z, eta, alpha, mu_W, W, obs_cell, N_pts, cell_size,
                          kappa, prior_mean = 0, prior_var = 100){
  new_mu_eta <- rnorm(1, curr_mu_eta, 0.005)
  p_curr <- f_mu_eta(curr_mu_eta, Z, alpha, eta, mu_W, W, N_pts, obs_cell, cell_size, kappa, prior_mean, prior_var)
  p_new <- f_mu_eta(new_mu_eta, Z, alpha, eta, mu_W, W, N_pts, obs_cell, cell_size, kappa, prior_mean, prior_var)
  if(log(runif(1)) < (p_new - p_curr)){
    new_mu_eta
  }else{
    curr_mu_eta
  }
}


update_sig2_W <- function(W, Sig_W, 
                          prior_shape = 0.01, prior_rate = 0.01){ 
  post_shape <- prior_shape + length(W)/2
  post_rate <- prior_rate + crossprod(solve(Sig_W, W), W)/2
  rinvgamma(1, shape = post_shape, rate = post_rate)
}

update_sig2_eta <- function(eta, Sig_eta,
                            prior_shape = 0.01, prior_rate = 0.01){
  post_shape <- prior_shape + length(eta)/2
  post_rate <- prior_rate + crossprod(solve(Sig_eta, eta), eta)/2
  rinvgamma(1, shape = post_shape, rate = post_rate)
}



f_log_phi_w <- function(W, log_phi_W, sig2_W,
                    D, cov.power, prior_mean, prior_var){
  Sigma <- sig2_W*exp(-(D/exp(log_phi_W))^cov.power)
  #t1 <- -0.5*det(Sigma)
  #t2 <- -0.5*crossprod(solve(Sigma, W), W)
  # t3 <- dnorm(log_phi_W, prior_mean, sqrt(prior_var), log = T)
  t1 <- dmvnorm(W, rep(0, length(W)), Sigma, log = T)
  t3 <- dgamma(exp(log_phi_W), 0.5*10, 10, log = T)
  # t1 + t2 + t3
  t1 + t3
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
  t2 <- dgamma(exp(log_phi_eta), 0.5*10, 10, log = T)  
  t1 + t2
}


update_phi_eta <- function(eta, curr_phi_eta, sig2_eta,
                         D_grid, cov.power){

  curr_log_phi_eta <- log(curr_phi_eta)
  new_log_phi_eta <- rnorm(1, curr_log_phi_eta, 0.1) # try a few sd values to get ~68% accept.prob
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



f_eta_j <- function(j, eta_j, alpha, Z, W, mu_eta, mu_W,
                    N_pts, prior_mean, prior_var, cell_area, obs_cell, 
                    kappa, cond_mean, cond_var){
  # likelihood from the locations
  t1 <- dpois(N_pts[j], exp(mu_eta + eta_j)*cell_area, log = T)
  # t1 <- N_pts[j]*(mu_eta + eta_j) - cell_area*exp(mu_eta + eta_j) 
    
  t2 <- 0
  id <- which(obs_cell == j) # check how many observations in this box
  
  # likelihood from the observation
  if(length(id) > 0){
    rate <- exp(alpha*eta_j + mu_W + W[j])
    t2 <- sum(dpois(Z[id], kappa*rate, log = T))
  }
  
  t3 <- dnorm(eta_j, cond_mean, sqrt(cond_var), log = T)
  t4 <- dnorm(eta_j, prior_mean, sqrt(prior_var), log = T)
  t1 + t2 + t3 + t4
}

update_eta <- function(curr_eta, alpha, Z, W, mu_eta, mu_W,
                       Sig_eta,
                       N_pts, obs_cell,
                       nbs, cell_size, kappa,
                       prior_mean = 0, prior_var = 50){
  new_eta <- curr_eta

  for(j in 1:length(new_eta)){
    curr_eta_j <- new_eta[j]
    new_eta_j <- rnorm(1, curr_eta_j, 0.25)
    Sigma_22 <- Sig_eta[nbs[[j]],][,nbs[[j]]]# + diag(rep(1e-5, length(nbs[[j]])))
    t1 <- Sig_eta[nbs[[j]],][,j]
    t2 <- solve(Sigma_22, t1) # Sigma22^-1%*%Sigma_12
    cond_mean <-  crossprod(t2, new_eta[nbs[[j]]])
    cond_var <- Sig_eta[j,j] - crossprod(t2, t1)
    p_curr <- f_eta_j(j, curr_eta_j, alpha, Z, W, mu_eta, mu_W,
                     N_pts, prior_mean, prior_var, cell_size[j], obs_cell, kappa, cond_mean, cond_var)
    p_new <- f_eta_j(j, new_eta_j, alpha, Z, W, mu_eta, mu_W,
                     N_pts, prior_mean, prior_var, cell_size[j], obs_cell, kappa, cond_mean, cond_var)
    if(log(runif(1)) < p_new - p_curr){
      new_eta[j] <- new_eta_j
    }
  }
  new_eta
} 

library(mvtnorm)


f_w_j <- function(j, W_j, alpha, Z, eta, mu_eta, mu_W,
                    prior_mean, prior_var, obs_cell, kappa, cond_mean, cond_var){

  t2 <- 0
  id <- which(obs_cell == j) # check how many observations in this box
  # if there is an observation, compute the log likelihood
  if(length(id) > 0){
    rate <- exp(alpha*eta[j] + mu_W + W_j)
    t2 <- sum(dpois(Z[id], kappa*rate, log = T))
  }
  
  t3 <- dnorm(W_j, cond_mean, sqrt(cond_var), log = T)
  t4 <- dnorm(W_j, prior_mean, sqrt(prior_var), log = T)
  t2 + t3 + t4
}

update_w <- function(curr_w, alpha, Z, eta, mu_eta, mu_W,
                       Sig_W, obs_cell, nbs, kappa,
                       prior_mean = 0, prior_var = 50){
  new_w <- curr_w

  for(j in 1:length(new_w)){
    curr_w_j <- new_w[j]
    new_w_j <- rnorm(1, curr_w_j, 0.1)
    Sigma_22 <- Sig_W[nbs[[j]],][,nbs[[j]]]# + diag(rep(1e-5, length(nbs[[j]])))
    t1 <- Sig_W[nbs[[j]],][,j]
    t2 <- solve(Sigma_22, t1) # Sigma22^-1%*%Sigma_12
    cond_mean <-  crossprod(t2, new_w[nbs[[j]]])
    cond_var <- Sig_W[j,j] - crossprod(t2, t1)
    p_curr <- f_w_j(j, curr_w_j, alpha, Z, eta, mu_eta, mu_W,
                     prior_mean, prior_var, obs_cell, kappa, cond_mean, cond_var)
    p_new <- f_w_j(j, new_w_j, alpha, Z, eta, mu_eta, mu_W,
                     prior_mean, prior_var, obs_cell, kappa, cond_mean, cond_var)
    if(log(runif(1)) < p_new - p_curr){
      new_w[j] <- new_w_j
    }
  }
  new_w
} 





text_on_cell <- function(id, words = id, size = 0.5, col = 'white'){
  text(S[id,1], S[id,2], words, cex = size, col = col)  
}



