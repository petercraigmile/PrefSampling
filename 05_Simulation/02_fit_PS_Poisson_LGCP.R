##============================================================
## Title: Bayesian inference for LGCP based PS data with Poisson
## distributed observations.
## Author: Rui Qiang
## Date: 5.2.2025
##============================================================
library(spatstat)

##====================================================================
## Load the functions
##====================================================================
load('../Data/Poisson_PS.RData')
source('../functions/Simulation_update_functions_Poisson.R')
Rcpp::sourceCpp("../functions/GSP_2024_04_11.cpp")
source("../functions/GSP_2024_04_11.R")

##====================================================================
## Specify the spatial parameters
##====================================================================
cov.power.eta <- theta.eta[3]
cov.power.w <- theta.w[3]
kappa <- 1

##====================================================================
## Set the number of iterations for the MCMC chain
##====================================================================
N.it <- 30000

which.sim <- sample(c(1:n.sim), 1)

##====================================================================
## Set the distance matrix.
## Use the cutoff to determine the range of the nearest neighbors. 
##====================================================================
cutoff <- 0.25
D.temp <- D + diag(1000, nrow(D))
nbs <- lapply(c(1:nrow(D)), function(i) as.vector(which(D.temp[i,] < cutoff)))

D_nbs <- lapply(c(1:nrow(D)), function(i) D[nbs[[i]], nbs[[i]]])
D_grid <- D
D_obs <- D.obs[[which.sim]]
rm(D.temp)

##====================================================================
## Create a vector to store the observations
##====================================================================
Z <- dat[[which.sim]]$z

##====================================================================
## obs_cell: see which cells have observations
## N_pts: number of points in each cell
##====================================================================
obs_cell <- dat[[which.sim]]$cell
N_pts <- grid.dat[[which.sim]]$n.pts

cell_size <- rep(1/(M^2), M^2)

##====================================================================
## Initialize the parameters' vectors for update
##====================================================================
sig2_W <- rep(theta.w[1], N.it)
sig2_W[1] <- 3

phi_W <- rep(theta.w[2], N.it) 
phi_W[1] <- 0.1

alpha <- rep(alpha.true, N.it)
alpha[1] <- 1.5

sig2_eta <- rep(1, N.it)
sig2_eta[1] <- 3

phi_eta <- rep(theta.eta[2], N.it) 
phi_eta[1] <- 0.5

mu_W <- rep(mean(log(Z+1)), N.it)
mu_W[1] <- -3

mu_eta <- rep(log(sum(N_pts)), N.it)

W <- matrix(rep(dat.pred[[which.sim]]$w - mu.w, N.it), nrow = N.it, byrow = T)
W[1,] <- rnorm(ncol(W))

eta <- matrix(rep(dat.pred[[which.sim]]$log.int - mu.eta, N.it), nrow = N.it, byrow = T)
eta[1,] <- rnorm(ncol(eta))


##====================================================================
## Run the MCMC chain
##====================================================================
time1 = Sys.time()
for(i in 2:N.it){
  # Update individually
  ## can try updating hyper-parameters after the processes
  alpha[i] <- update_alpha(alpha[i-1], eta[i-1,], Z, obs_cell, 
                           W[i-1,], mu_W[i-1], mu_eta[i-1], kappa)
  
  mu_W[i] <- update_mu_w(mu_W[i-1], W[(i-1),], alpha[i], eta[(i-1),], 
                         mu_eta[i-1], Z, obs_cell, kappa)
  
  # Metropolis Hastings
  phi_eta[i] <- update_phi_eta(eta[(i-1),], phi_eta[i-1], sig2_eta[i-1], D_grid, cov.power.eta)
  phi_W[i] <- update_phi_w(W[(i-1),], phi_W[i-1], sig2_W[i-1], D_grid, cov.power.w)
  
  mu_eta[i] <- update_mu_eta(mu_eta[i-1], Z, eta[(i-1),], alpha[i], mu_W[i], 
                           W[(i-1),], obs_cell, N_pts, cell_size, kappa, 
                        prior_mean = log(pp[[which.sim]]$n), prior_var = 5)

  Sig_W <- GSP.powerexp(D_grid, c(sig2_W[i-1], phi_W[i], cov.power.w)) # + diag(1e-5, nrow(D_grid))
  W[i,] <- update_w(W[(i-1),], alpha[i], Z, eta[(i-1),], mu_eta[i], mu_W[i], Sig_W, obs_cell, nbs, kappa)
  
  Sig_eta <- GSP.powerexp(D_grid, c(sig2_eta[i-1], phi_eta[i], cov.power.eta)) # + diag(1e-5, nrow(D_grid))
  eta[i,] <- update_eta(eta[(i-1),], alpha[i], Z, W[i,], mu_eta[i], mu_W[i], Sig_eta, N_pts, obs_cell, nbs, cell_size, kappa)

  Sig_W <- GSP.powerexp(D_grid, c(1, phi_W[i], cov.power.w))
  sig2_W[i] <- update_sig2_W(W[i,], Sig_W)
  
  Sig_eta <- GSP.powerexp(D_grid, c(1, phi_eta[i], cov.power.eta))
  sig2_eta[i] <- update_sig2_eta(eta[i,], Sig_eta)
  
  
  if(i %% 50 == 0) cat(paste(i, ", "))
}
time2 = Sys.time()
dur = time2-time1

##====================================================================
## Check acceptance probability
##====================================================================
mean(sapply(c(1:ncol(eta)), function(i) mean(eta[(2:N.it),i] != eta[(1:(N.it-1)),i])))
mean(sapply(c(1:ncol(W)), function(i) mean(W[(2:N.it),i] != W[(1:(N.it-1)),i])))

mean(mu_eta[2:N.it] != mu_eta[1:(N.it-1)])
mean(phi_W[2:N.it] != phi_W[1:(N.it-1)])
mean(phi_eta[2:N.it] != phi_eta[1:(N.it-1)])

mean(alpha[2:N.it] != alpha[1:(N.it-1)])
mean(mu_W[2:N.it] != mu_W[1:(N.it-1)])

##====================================================================
## Trim the MCMC chain
##====================================================================
burn <- floor(N.it/6)
lag <- 5

##====================================================================
## Plot the estimated log intensity and the estimated random
## process \eta
##====================================================================
est.log.int.im <- im(matrix(colMeans(eta[seq(burn, N.it, lag),]) + mean(mu_eta[seq(burn, N.it, lag)]), nrow = M, byrow = T), xrange = c(0,1), yrange = c(0,1))
est.eta.im <- im(matrix(colMeans(eta[seq(burn, N.it, lag),]), nrow = M, byrow = T), xrange = c(0,1), yrange = c(0,1))
plot(est.eta.im)
plot(est.log.int.im)
points(pp[[which.sim]], cex = 0.25, pch = 19)
plot(log.int.im[[which.sim]])
points(pp[[which.sim]], cex = 0.25, pch = 19)

##====================================================================
## Calculate the posterior mean of 
##the mean of the observation process
##====================================================================
est_mu <- sapply(c(seq(burn, N.it, 5)), function(i) alpha[i]*eta[i,] + W[i,] + mu_W[i])
est_mu_im <- matrix(rowMeans(est_mu), nrow = M, byrow = T)
est.y.im <- im(matrix(colMeans(W[seq(burn, N.it, lag),]) + mean(mu_W[seq(burn, N.it, lag)]), nrow = M, byrow = T), xrange = c(0,1), yrange = c(0,1))
eta.im <- im(matrix(colMeans(eta[seq(burn, N.it, lag),]), nrow=M, byrow = T), xrange = c(0,1), yrange =c(0,1))
  
##====================================================================
## Metrics
##====================================================================
IMSE <- mean((est_mu_im - mu.im[[which.sim]]$v)^2)
inner <- (abs(S[,1] - 0.5) < 0.4) & (abs(S[,2] - 0.5) < 0.4)
POP <- mean(est_mu_im > mu.im[[which.sim]]$v)
MAE <- mean(abs(est_mu_im - mu.im[[which.sim]]$v))


##====================================================================
## Calculate the coverage probability
##====================================================================
sd_est_mu <- matrix(sapply(c(1:nrow(est_mu)), function(i) sd(est_mu[i,])), nrow = M, byrow = T)

est_mu_low <- matrix(sapply(c(1:nrow(est_mu)), function(i) quantile(est_mu[i,], 0.025)), nrow = M, byrow = T)
est_mu_high <- matrix(sapply(c(1:nrow(est_mu)), function(i) quantile(est_mu[i,], 0.975)), nrow = M, byrow = T)

Cov_Prob_95 <- mean((est_mu_high > mu.im[[which.sim]]$v) & (est_mu_low < mu.im[[which.sim]]$v))