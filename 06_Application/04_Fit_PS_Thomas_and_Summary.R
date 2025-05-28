##============================================================
## Title: LGCP with parameter estimation under Exp Cov Fun
## Author: Rui Qiang
## Date: 4.19.2024
##
##============================================================
# set.seed(123)
# N.it <- 4
library(spatstat)
library(sf)

# Entire four stat
load('../Data/prepped_data_40.RData')
Rcpp::sourceCpp("../functions/GSP_2024_04_11.cpp")
source("../functions/GSP_2024_04_11.R")
source("../functions/GHCN_update_functions_Thomas.R")
cov.power <- 1
delta.true <- 10
omega2.true <- 0.1
library(spatstat)

parents <- rbind(c(-111.64, 33.35), # AZ
                 c(-118.3, 33.68),  # CA
                 c(-122.42, 37.68), # CA
                 c(-119.94, 38.97), # NV
                 c(-111.68, 41.47), # UT
                 c(-111.68, 40.88)) # UT
cov.power.w <- 1 

library(spatstat)

# Set the number of iterations
N.it <- 500000 

##====================================
## Select individual states here
##====================================

full_dat <- dat
full_stations <- stations
sel <- unlist(st_intersects(NV, cells))
dat <- dat[sel,]
stations <- stations[stations$state == 'NV',]

blocks <- matrix(c(1:nrow(dat)), ncol = 3)
parents <- matrix(parents[4,], ncol = 2)



##====================================
## Initialize parameters
##====================================
D_grid <- Edist(cbind(dat$lon, dat$lat))
D_obs <- Edist(cbind(stations$lon, stations$lat))

#Z <- stations$tavg
Z <- stations$tavg

obs_cell <- match(stations$cell, dat$cell_id)

N_pts <- dat$npts

cell_size <- dat$area

sig2_W <- rep(1, N.it)
phi_W <- rep(1,5, N.it) 
alpha <- rep(1, N.it)
delta <- rep(20, N.it)
omega2 <- rep(0.1, N.it)
mu_W <- rep(15, N.it)
tau2 <- rep(1, N.it)
W <- matrix(rep(rnorm(nrow(dat)), N.it), nrow = N.it, byrow = T)
W[1,] <- rnorm(ncol(W))
S <- cbind(dat$lon, dat$lat)


time1 = Sys.time()
for(i in 2:N.it){
  # Gibbs
  ## can try updating hyper-parameters after the processes
  delta[i] <- update_delta(delta[i-1], omega2[i-1], W[i-1,], mu_W[i-1], alpha[i-1], tau2[i-1], 
                           Z, S, parents,
                           N_pts, obs_cell, cell_size, 0, 5) # prior is on the log of delta
  omega2[i] <- update_omega2(omega2[i-1], delta[i], W[i-1,], mu_W[i-1], alpha[i-1], tau2[i-1],
                             Z, S, parents,
                             N_pts, obs_cell, cell_size, 0, 5) # prior is on the log of omega2
  est.int <- (lambda(S, parents, delta[i], omega2[i]))
  alpha[i] <- update_alpha(est.int, Z, obs_cell, tau2[i-1], W[i-1,], mu_W[i-1], 0, 50)
  tau2[i] <- update_tau2(Z, est.int, alpha[i], W[i-1,], mu_W[i-1], obs_cell, 1, 1)
  sig2_W[i] <- update_sig2_W(W[(i-1),], phi_W[i-1], D_grid, cov.power.w)
  mu_W[i] <- update_mu_w(Z, est.int, alpha[i], W[i-1,], tau2[i], obs_cell, 0, 25)
  # Metropolis Hastings
  phi_W[i] <- update_phi_w(W[(i-1),], phi_W[i-1], sig2_W[i], D_grid, cov.power) # look at trace and accept.prob
  Sig <- GSP.powerexp(D_grid, c(sig2_W[i], phi_W[i], cov.power)) + diag(1e-5, nrow(D_grid))
  W[i,] <- update_W(W[(i-1), ], alpha[i], Z, est.int, mu_W[i], tau2[i],
                                    D_grid, obs_cell, blocks, cov.power, N_pts, Sig, 0, 10)
  if(i %% 100 == 0) cat(paste(i, ", "))
}
time2 = Sys.time()
dur = time2-time1


## Check acceptance probability
mean(delta[2:N.it] != delta[1:(N.it-1)])
mean(omega2[2:N.it] != omega2[1:(N.it-1)])
mean(phi_W[2:N.it] != phi_W[1:(N.it-1)])

library(spatstat)

burn <- floor(N.it/5)
lag <- 5

library(parallel)
est.int <- rowMeans(simplify2array(mclapply(seq(burn, N.it, lag), 
                                        function(i) (lambda(S, parents, delta[i], omega2[i])),
                                        mc.cores = 5)))

est.int.im <- rep(NA, x.len*y.len)
est.int.im[sel] <- est.int
est.int.im <- im(matrix(est.int.im, y.len, byrow = T), 
                 xcol = unique(full_dat$lon), 
                 yrow = unique(full_dat$lat))


est.mu <- rowMeans(simplify2array(mclapply(seq(burn, N.it, lag), 
                                          function(i) alpha[i]*(lambda(S, parents, delta[i], omega2[i])) + mu_W[i] + W[i,],
                                          mc.cores = 5)))
est.sd <- apply(simplify2array(mclapply(seq(burn, N.it, lag), 
                                          function(i) alpha[i]*(lambda(S, parents, delta[i], omega2[i])) + mu_W[i] + W[i,],
                                          mc.cores = 5)),1, sd)

est.mu.im <- rep(NA, x.len*y.len)
est.mu.im[sel] <- est.mu
est.mu.im <- im(matrix(est.mu.im, y.len, byrow = T), 
                xcol = unique(full_dat$lon), 
                yrow = unique(full_dat$lat))
est.sd.im <- rep(NA, x.len*y.len)
est.sd.im[sel] <- est.sd
est.sd.im <- im(matrix(est.sd.im, y.len, byrow = T), 
                xcol = unique(full_dat$lon), 
                yrow = unique(full_dat$lat))

#save.image("GHCN_PS_Thomas_NV_052725.RData")

keep <- seq(burn, N.it, lag)
summarize <- function(par){
  round(c(mean(par[keep]), quantile(par[keep], c(0.025, 0.975))),3)
}

round(rbind(summarize(tau2),
      summarize(alpha),
      summarize(mu_W),
      summarize(phi_W),
      summarize(sig2_W),
      summarize(omega2),
      summarize(delta)),2)
