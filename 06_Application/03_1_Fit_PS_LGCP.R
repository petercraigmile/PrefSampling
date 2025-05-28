##============================================================
## Title: LGCP with parameter estimation under Exp Cov Fun
## Author: Rui Qiang
## Date: 5.7.2025
##============================================================
library(spatstat)
library(sf)

##================================================================
## Load the packages
##================================================================
load('../Data/prepped_data_40.RData')
source('../functions/GHCN_update_functions_LGCP.R')
Rcpp::sourceCpp("../functions/GSP_2024_04_11.cpp")
source("../functions/GSP_2024_04_11.R")


N.it <- 50000
## Grid level ldata
full_dat <- dat
## State level data
full_stations <- stations

##====================================================
## If we want to look at each state individually
##====================================================
state <- 'UT'
sel <- st_intersects(get(state), full_dat$geometry)[[1]]
dat <- full_dat[sel,]
stations <- full_stations[full_stations$state == state,]

## Dividing the grid cells into blocks
blocks <- matrix(c(1:nrow(dat)), ncol = 4)  

## Set nearest neighbor cutoff distance
cutoff <- 2

## cell size
cell_size <- dat$area

## Calculate the distance matrix for the observations and the grid cells
D_grid <- Edist(cbind(dat$lon, dat$lat))
D_obs <- Edist(cbind(stations$lon, stations$lat))

## Compute the nearest neighbors  
D.temp <- D_grid + diag(1000, nrow(D_grid))
nbs <- lapply(c(1:nrow(D_grid)), function(i) as.vector(which(D.temp[i,] < cutoff)))
D_nbs <- lapply(c(1:nrow(D_grid)), function(i) D_grid[nbs[[i]], nbs[[i]]])
rm(D.temp)

## Observations
Z <- stations$tavg

## Remap the station to the current indices of the cell
obs_cell <- match(stations$cell, dat$cell_id)
## Number of observations in each cell
N_pts <- dat$npts
  
##===================================
## Initialize the MCMC chain
##===================================
sig2_W <- rep(2, N.it)

phi_W <- rep(2.5, N.it)

alpha <- rep(rnorm(1), N.it)

sig2_eta <- rep(runif(1), N.it)

phi_eta <- rep(2.5, N.it)

mu_W <- rep(rnorm(1), N.it)

mu_eta <- rep(rnorm(1), N.it)
              
tau2 <- rep(runif(1), N.it)

W <- matrix(rep(rnorm(nrow(dat)), N.it), nrow = N.it, byrow = T)

eta <- matrix(rep(rnorm(nrow(dat)), N.it), nrow = N.it, byrow = T)

cov.power.w <- 1
cov.power.eta <- 1.75

##===================================
## Run the MCMC
##===================================
time1 = Sys.time()
for(i in 2:N.it){
  # Gibbs
  ## can try updating hyper-parameters after the processes
  alpha[i] <- update_alpha(eta[(i-1),], Z, obs_cell, tau2[i-1], W[(i-1),], mu_W[i-1], mu_eta[i-1])
  tau2[i] <- update_tau2(Z, eta[(i-1),], alpha[i], W[(i-1),], mu_eta[i-1], mu_W[i-1], obs_cell)
  
  sig2_W[i] <- update_sig2_W(W[(i-1),], phi_W[i-1], D_grid, cov.power.w)
  sig2_eta[i] <- update_sig2_eta(eta[(i-1),], phi_eta[i-1], D_grid, cov.power.eta) 
  mu_W[i] <- update_mu_w(Z, mu_eta[i-1], eta[(i-1),], alpha[i], W[(i-1),], tau2[i], obs_cell)
  
  # Metropolis Hastings
  phi_eta[i] <- update_phi_eta(eta[(i-1),], phi_eta[i-1], sig2_eta[i], D_grid, cov.power.eta)
  phi_W[i] <- update_phi_w(W[(i-1),], phi_W[i-1], sig2_W[i], D_grid, cov.power.w)
  mu_eta[i] <- update_mu_eta(mu_eta[i-1], Z, eta[(i-1),], alpha[i], mu_W[i], W[(i-1),], tau2[i], obs_cell,
                             N_pts, cell_size, prior_mean = log(sum(N_pts)), prior_var = 100)
  
  Sig <- GSP.powerexp(D_grid, c(sig2_W[i], phi_W[i], cov.power.w)) + diag(1e-5, nrow(D_grid))
  W[i,] <- update_w_block_Gibbs_new(W[(i-1),], alpha[i], Z, eta[(i-1),], 
                                   mu_W[i], mu_eta[i], tau2[i],
                                    D_grid, obs_cell, blocks, cov.power.w, 
                                    N_pts, Sig)
  
  Sig_eta <- GSP.powerexp(D_grid, c(sig2_eta[i], phi_eta[i], cov.power.eta)) + diag(1e-5, nrow(D_grid))
  eta[i,] <- update_eta(eta[(i-1),], alpha[i], Z, W[i,], mu_eta[i], mu_W[i],
                        Sig_eta, tau2[i],
                        D_grid, N_pts, obs_cell, nbs, D_nbs, cell_size, cov.power.eta)
  
  if(i %% 50 == 0) cat(paste(i, ", "))
}
time2 = Sys.time()
dur = time2-time1



##===================================
## Check acceptance probability
##===================================

mean(sapply(c(1:ncol(eta)), function(i) mean(eta[(2:N.it),i] != eta[(1:(N.it-1)),i])))
mean(mu_eta[2:N.it] != mu_eta[1:(N.it-1)])
mean(phi_W[2:N.it] != phi_W[1:(N.it-1)])
mean(phi_eta[2:N.it] != phi_eta[1:(N.it-1)])

## Trim the MCMC chain
burn <- floor(N.it/5)
lag <- 5

keep <- seq(burn, N.it, lag)
## Calculate the posterior mean and the 95% posterior credible interval
summarize <- function(par){
  round(c(mean(par[keep]), quantile(par[keep], c(0.025, 0.975))),3)
}

rbind(summarize(tau2),
      summarize(alpha),
      summarize(mu_W),
      summarize(phi_W),
      summarize(sig2_W),
      summarize(mu_eta),
      summarize(phi_eta),
      summarize(sig2_eta))

mean(alpha[seq(burn, N.it, lag)])
sd(alpha[seq(burn, N.it, lag)])

est.log.int <- rep(NA, nrow(full_dat))
est.log.int[sel] <- colMeans(eta[seq(burn, N.it, lag),]) + mean(mu_eta[seq(burn, N.it, lag)])
est.log.int.im <- im(matrix(est.log.int, nrow = y.len, byrow = T), xcol = unique(full_dat$lon), yrow = unique(full_dat$lat))

plot(est.log.int.im)
points(stations$lon, stations$lat, pch = 19, cex = 0.25)

est_mu <- rep(NA, nrow(full_dat))
est_mu[sel] <- rowMeans(sapply(c(seq(burn, N.it, 5)), function(i) alpha[i]*(mu_eta[i] + eta[i,]) + W[i,] + mu_W[i]))
est_mu_im <- im(matrix(est_mu, nrow = y.len, byrow = T), xcol = unique(full_dat$lon), yrow = unique(full_dat$lat))
dat$mu <- rowMeans(sapply(c(seq(burn, N.it, 5)), function(i) alpha[i]*(eta[i,]) + W[i,] + mu_W[i]))

est_sd <- rep(NA, nrow(full_dat))
est_sd[sel] <- apply(sapply(c(seq(burn, N.it, 5)), function(i) alpha[i]*(mu_eta[i] + eta[i,]) + W[i,] + mu_W[i]), 1, sd)
est_sd_im <- im(matrix(est_sd, nrow = y.len, byrow = T), xcol = unique(full_dat$lon), yrow = unique(full_dat$lat))

save.image("Analyzed_PS_Normal_GHCN_Results.RData")
