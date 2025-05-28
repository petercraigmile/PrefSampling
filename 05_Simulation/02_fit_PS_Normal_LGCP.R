##============================================================
## Title: LGCP with parameter estimation under Exp Cov Fun
## Author: Rui Qiang
## Date: 5.2.2025
##============================================================

library(spatstat)
load('../Data/Normal_PS.RData')
source('../functions/Simulation_update_functions_LGCP.R')
Rcpp::sourceCpp("../functions/GSP_2024_04_11.cpp")
source("../functions/GSP_2024_04_11.R")

cov.power.eta <- theta.eta[3]
cov.power.w <- theta.w[3]

##========================================
## Chop the entire grid into blocks
## for faster update.
##========================================
blocks <- matrix(c(1:M^2), ncol = 5)

##========================================
## Set the number of iterations for the
## MCMC chain.
##========================================
N.it <- 30000

## Choose a simluated set to do MCMC on
which.sim <- sample(c(1:n.sim), 1)
plot(log.int.im[[which.sim]])
points(pp[[which.sim]])
plot(mu.im[[which.sim]])
points(pp[[which.sim]])


##========================================
## Set up a temporary distance matrix
## This determines the nearest neighbors
## Set diagonal values to 1000 to avoid
## having neighbors with itself.
## Cutoff determines the range of the 
## nearest neighbors.
##========================================
cutoff <- 0.25
D.temp <- D + diag(1000, nrow(D))
nbs <- lapply(c(1:nrow(D)), function(i) as.vector(which(D.temp[i,] < cutoff)))

D_nbs <- lapply(c(1:nrow(D)), function(i) D[nbs[[i]], nbs[[i]]])
D_grid <- D
D_obs <- D.obs[[which.sim]]
rm(D.temp)

Z <- dat[[which.sim]]$z

obs_cell <- dat[[which.sim]]$cell

N_pts <- grid.dat[[which.sim]]$n.pts

cell_size <- rep(1/(M^2), M^2)

##========================================
## Vectors for parameter updates.
## Randomize a starting value
##========================================
sig2_W <- rep(theta.w[1], N.it)
sig2_W[1] <- runif(1)

phi_W <- rep(theta.w[2], N.it) 
phi_W[1] <- runif(1)

alpha <- rep(alpha.true, N.it)
alpha[1] <- runif(1)

sig2_eta <- rep(theta.eta[1], N.it)
sig2_eta[1] <- runif(1)

phi_eta <- rep(theta.eta[2], N.it) 
phi_eta[1] <- 0.5

mu_W <- rep(mu.w, N.it)
mu_W[1] <- mean(Z) - log(sum(N_pts))

mu_eta <- rep(mu.eta, N.it)
mu_eta[1] <- log(sum(N_pts))

tau2 <- rep(tau2.true, N.it)
tau2[1] <- runif(1)

W <- matrix(rep(dat.pred[[which.sim]]$w - mu.w, N.it), nrow = N.it, byrow = T)

W[1,] <- rnorm(ncol(W))

eta <- matrix(rep(dat.pred[[which.sim]]$log.int - mu.eta, N.it), nrow = N.it, byrow = T)
eta[1,] <-  rnorm(ncol(eta))

##====================
## Run MCMC
##====================

time1 = Sys.time()
for(i in 2:N.it){
  # Gibbs
  alpha[i] <- update_alpha(eta[(i-1),], Z, obs_cell, tau2[i-1], W[(i-1),], mu_W[i-1], mu_eta[i-1])
  tau2[i] <- update_tau2(Z, eta[(i-1),], alpha[i], W[(i-1),], mu_eta[i-1], mu_W[i-1], obs_cell)
  
  sig2_W[i] <- update_sig2_W(W[(i-1),], phi_W[i-1], D_grid, cov.power.w)
  sig2_eta[i] <- update_sig2_eta(eta[(i-1),], phi_eta[i-1], D_grid, cov.power.eta) 
  mu_W[i] <- update_mu_w(Z, mu_eta[i-1], eta[(i-1),], alpha[i], W[(i-1),], tau2[i], obs_cell)
  
  # Metropolis Hastings
  phi_eta[i] <- update_phi_eta(eta[(i-1),], phi_eta[i-1], sig2_eta[i], D_grid, cov.power.eta)
  phi_W[i] <- update_phi_w(W[(i-1),], phi_W[i-1], sig2_W[i], D_grid, cov.power.w)
  mu_eta[i] <- update_mu_eta(mu_eta[i-1], Z, eta[(i-1),], alpha[i], mu_W[i], W[(i-1),], tau2[i], obs_cell,
                             N_pts, cell_size, prior_mean = log(pp[[which.sim]]$n), prior_var = 100)
  
  Sig <- GSP.powerexp(D_grid, c(sig2_W[i], phi_W[i], cov.power.w)) + diag(1e-5, nrow(D_grid))
  W[i,] <- update_w_block_Gibbs_new(W[(i-1),], alpha[i], Z, eta[(i-1),], 
                                   mu_W[i], mu_eta[i], tau2[i],
                                    D_grid, obs_cell, blocks, cov.power, 
                                    N_pts, Sig)
  
  Sig_eta <- GSP.powerexp(D_grid, c(sig2_eta[i], phi_eta[i], cov.power.eta)) + diag(1e-5, nrow(D_grid))
  eta[i,] <- update_eta(eta[(i-1),], alpha[i], Z, W[i,], mu_eta[i], mu_W[i],
                        Sig_eta, tau2[i],
                        D_grid, N_pts, obs_cell, nbs, D_nbs, cell_size, cov.power)
  
  if(i %% 50 == 0) cat(paste(i, ", "))
}
time2 = Sys.time()
dur = time2-time1

##====================================
## Check acceptance probability
##====================================
mean(sapply(c(1:ncol(eta)), function(i) mean(eta[(2:N.it),i] != eta[(1:(N.it-1)),i])))
mean(mu_eta[2:N.it] != mu_eta[1:(N.it-1)])
mean(phi_W[2:N.it] != phi_W[1:(N.it-1)])
mean(phi_eta[2:N.it] != phi_eta[1:(N.it-1)])


##============================================
## Check MCMC Chain and Posterior Estimates
##============================================
burn <- floor(N.it/6)
lag <- 5
est.log.int.im <- im(matrix(colMeans(eta[seq(burn, N.it, lag),]) + mean(mu_eta[seq(burn, N.it, lag)]), nrow = M, byrow = T), xrange = c(0,1), yrange = c(0,1))
est.eta.im <- im(matrix(colMeans(eta[seq(burn, N.it, lag),]), nrow = M, byrow = T), xrange = c(0,1), yrange = c(0,1))

est_mu <- sapply(c(seq(burn, N.it, 5)), function(i) alpha[i]*eta[i,] + W[i,] + mu_W[i])
est_mu_im <- im(matrix(rowMeans(est_mu), nrow = M, byrow = T), xrange = c(0,1), yrange = c(0,1))
est.y.im <- im(matrix(colMeans(W[seq(burn, N.it, lag),]) + mean(mu_W[seq(burn, N.it, lag)]), nrow = M, byrow = T), xrange = c(0,1), yrange = c(0,1))

##========================
## Metrics
##========================
IMSE <- mean((est_mu_im - mu.im[[which.sim]]$v)^2)
POP <- mean(est_mu_im > mu.im[[which.sim]]$v)
MAE <- mean(abs(est_mu_im - mu.im[[which.sim]]$v))

sd_est_mu <- matrix(sapply(c(1:nrow(est_mu)), function(i) sd(est_mu[i,])), nrow = M, byrow = T)
est_mu_low <- matrix(sapply(c(1:nrow(est_mu)), function(i) quantile(est_mu[i,], 0.025)), nrow = M, byrow = T)
est_mu_high <- matrix(sapply(c(1:nrow(est_mu)), function(i) quantile(est_mu[i,], 0.975)), nrow = M, byrow = T)

##========================================
## Coverage Probability
##========================================
Cov_Prob_95 <- mean((est_mu_high > mu.im[[which.sim]]$v) & (est_mu_low < mu.im[[which.sim]]$v))
