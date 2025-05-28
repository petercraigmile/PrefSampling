##============================================================
## Title: PS under a Thomas Point Process with Gaussian observations
## Author: Rui Qiang
## Date: 4.19.2024
##============================================================
library(spatstat)

load("../Data/Thomas_PS.RData")
Rcpp::sourceCpp("../functions/GSP_2024_04_11.cpp")
source("../functions/GSP_2024_04_11.R")
source("../functions/Simulation_update_functions_Thomas.R")
delta.true <- delta
omega2.true <- omega2
library(spatstat)

cov.power.w <- 1 

blocks <- matrix(c(1:M^2), ncol = 5)

# Set the number of iterations
N.it <- 30000

which.sim <- sample(c(1:n.sim), 1)


##============================================
## Distance Matrices for this simulation
##============================================
D_grid <- D
D_obs <- D.obs[[which.sim]]

Z <- dat[[which.sim]]$z

obs_cell <- dat[[which.sim]]$cell

N_pts <- grid.dat[[which.sim]]$n.pts

cell_size <- rep(1/(M^2), M^2)

sig2_W <- rep(theta.w[1], N.it)
sig2_W[1] <- 1

phi_W <- rep(theta.w[2], N.it) 
phi_W[1] <- 0.1

alpha <- rep(alpha.true, N.it)
alpha[1] <- 0.5

delta <- rep(delta.true, N.it)
delta[1] <- 20

omega2 <- rep(omega2.true, N.it)
omega2[1] <- 0.05

mu_W <- rep(mu.w, N.it)
mu_W[1] <- 3

tau2 <- rep(tau2.true, N.it)
tau2[1] <- 2

W <- matrix(rep(dat.pred[[which.sim]]$w - mu.w, N.it), nrow = N.it, byrow = T)

W[1,] <- rnorm(ncol(W))

time1 = Sys.time()
for(i in 2:N.it){
  ## Gibbs steps
  delta[i] <- update_delta(delta[i-1], omega2[i-1], W[i-1,], mu_W[i-1], alpha[i-1], tau2[i-1], 
                       Z, S, parents[[which.sim]],
                        N_pts, obs_cell, cell_size, 0, 10)
  
  omega2[i] <- update_omega2(omega2[i-1], delta[i], W[i-1,], mu_W[i-1], alpha[i-1], tau2[i-1],
                        Z, S, parents[[which.sim]],
                       N_pts, obs_cell, cell_size, 0, 10)
  
  log.int <- log(lambda(S, parents[[which.sim]], delta[i], omega2[i]))
  alpha[i] <- update_alpha(log.int, Z, obs_cell, tau2[i-1], W[i-1,], mu_W[i-1], 0, 50)
  tau2[i] <- update_tau2(Z, log.int, alpha[i], W[i-1,], mu_W[i-1], obs_cell, 1, 1)
  
  sig2_W[i] <- update_sig2_W(W[(i-1),], phi_W[i-1], D_grid, cov.power.w)
  mu_W[i] <- update_mu_w(Z, log.int, alpha[i], W[i-1,], tau2[i], obs_cell, 0, 25)

  ## Metropolis Hastings
  phi_W[i] <- update_phi_w(W[(i-1),], phi_W[i-1], sig2_W[i], D_grid, cov.power)
  
  Sig <- GSP.powerexp(D_grid, c(sig2_W[i], phi_W[i], cov.power)) + diag(1e-5, nrow(D_grid))
  W[i,] <- update_W(W[(i-1), ], alpha[i], Z, log.int, mu_W[i], tau2[i],
                                    D_grid, obs_cell, blocks, cov.power, N_pts, Sig, 0, 10)

  if(i %% 50 == 0) cat(paste(i, ", "))
}
time2 = Sys.time()
dur = time2-time1

##============================================
## Check acceptance probability
##============================================
mean(delta[2:N.it] != delta[1:(N.it-1)])
mean(omega2[2:N.it] != omega2[1:(N.it-1)])
mean(phi_W[2:N.it] != phi_W[1:(N.it-1)])

##============================================
## Trim the MCMC chain and do some plots
##============================================
burn <- floor(N.it/6)
lag <- 5

##============================================
## Posterior mean of the log intensity
##============================================
est.log.int <- rowMeans(simplify2array(mclapply(seq(burn, N.it, lag), 
                                        function(i) log(lambda(S, parents[[which.sim]], delta[i], omega2[i])),
                                        mc.cores = 5)))
est.log.int.im <- im(matrix(est.log.int, M, byrow = T), xrange = c(0,1), yrange = c(0,1))

##============================================
## Posterior mean of the mean of the observation
## process (\mu), and the random process W
##============================================
est.mu <- rowMeans(simplify2array(mclapply(seq(burn, N.it, lag), 
                                        function(i) alpha[i]*log(lambda(S, parents[[which.sim]], delta[i], omega2[i])) + mu_W[i] + W[i,],
                                        mc.cores = 5)))
est_mu_im <- im(matrix(est.mu, M, byrow = T), xrange = c(0,1), yrange = c(0,1))

est.w.im <- im(matrix(colMeans(W[seq(burn, N.it, lag),]), nrow=M, byrow = T), 
               xrange = c(0,1), yrange =c(0,1))

##============================================
## Metrics
##============================================
IMSE <- mean((est_mu_im - mu.im[[which.sim]]$v)^2)
POP <- mean(est_mu_im > mu.im[[which.sim]]$v)
MAE <- mean(abs(est_mu_im - mu.im[[which.sim]]$v))
