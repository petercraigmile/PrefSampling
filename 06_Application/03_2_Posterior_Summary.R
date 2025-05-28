##============================================================
## Title: Summarize the posterior draws
## Author: Rui Qiang
## Date: 5.7.2025
##============================================================
library(spatstat)

## Load the saved data and functions
load('Analyzed_PS_Normal_GHCN_ALL_Results.RData')

Rcpp::sourceCpp("../functions/GSP_2024_04_11.cpp")
source(".../functions/GSP_2024_04_11.R")


burn <- floor(N.it/5)
lag <- 5
keep <- seq(burn, N.it, lag)
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


mean(alpha[keep])
sd(alpha[keep])

est.log.int <- rep(NA, nrow(full_dat))
est.log.int[full_dat$area > 0] <- colMeans(eta[keep,]) + mean(mu_eta[keep])
est.log.int.im <- im(matrix(est.log.int, nrow = y.len, byrow = T), 
                     xcol = unique(dat$lon), 
                     yrow = unique(dat$lat))

plot(est.log.int.im)
points(stations$lon, stations$lat, pch = 19, cex = 0.25)

est_mu <- rep(NA, nrow(full_dat))
est_mu[full_dat$area > 0] <- rowMeans(sapply(c(keep), function(i) alpha[i]*(eta[i,] + mu_eta[i]) + W[i,] + mu_W[i]))
est_mu_im <- im(matrix(est_mu, nrow = y.len, byrow = T), 
                     xcol = unique(dat$lon), 
                     yrow = unique(dat$lat))
plot(est_mu_im)
plot(states_win, add = T)

est_sd <- rep(NA, nrow(full_dat))
est_sd[full_dat$area > 0] <- apply((sapply(c(keep), function(i) alpha[i]*(eta[i,] + mu_eta[i]) + W[i,] + mu_W[i])), 1, sd)
est_sd_im <- im(matrix(est_sd, nrow = y.len, byrow = T), 
                     xcol = unique(dat$lon), 
                     yrow = unique(dat$lat))

plot(est_sd_im)

var(Z - est_mu[stations$cell_id])
plot(states_win)
plot(est_mu_im, add = T)
points(stations$lon, stations$lat, pch = 19, cex = 0.25)
save.image("Analyzed_PS_Normal_GHCN_Results_Unity.RData")


