##========================================##
## Title: Generate an LGCP based PS dataset
## with Poisson distributed observation 
## Author: Rui Qiang
## Date: 5.2.2025
##========================================##
library(spatstat)
library(parallel)
library(mvtnorm)

##========================================================
## Load the functions
##========================================================
Rcpp::sourceCpp("../functions/GSP_2024_04_11.cpp")
source("../functions/GSP_2024_04_11.R")
source("../functions/sim_functions.R")

##========================================================
## Specify number of simulations
##========================================================
cov.fun <- GSP.powerexp
n.sim <- 5

##========================================================
## Spatial parameters for the random processes eta and W
##========================================================
theta.eta <- c(0, 0.15, 1.75)  # partial sill, phi, power
theta.w <- c(1, 0.33, 1) # partial sill, range, power

##========================================================
## Mean of the random process W
##========================================================
mu.w <-  4 

##========================================================
## Preferential Sampling parameter \alpha
##========================================================
alpha <- 1 # PS parameter

##========================================================
## Additional scaling parameter \kappa
## Set to 1
##========================================================
kappa <- 1

##========================================================
## Trend component for x-y
## Set to 0 for no spatial trend in the mean
##========================================================
a <- 0
b <- 0

##========================================================
## Mean of the random process eta
##========================================================
mu.eta <- 5.25

##========================================================
## Resolution of the gridded spatial domain
##========================================================
M  <- 20
xs <- seq(.5/M, 1-.5/M, length.out = M)

##========================================================
## Coordinates and distance matrix for the grid
##========================================================
S  <- make.grid(list(xs, xs))
D <- Edist(S)

##========================================================
## Rename some parameters for checking later.
## This avoids accidental overwrite of the data generating parameters
##========================================================
phi.true <- theta.eta[2]
alpha.true <- alpha

##========================================================
## Generate the LGCP
##========================================================
trend.lgcp <- im(matrix(mu.eta + a*S[,1] + b*S[,2], nrow = M, byrow = T), xrange = c(0,1), yrange = c(0,1))
gp <- lapply(c(1:n.sim), function(i) GSP.sim(1, mu = rep(0, nrow(D)), cov.fun = cov.fun, dists=D, theta=theta.eta))
gp.im <- lapply(gp, function(x) im(matrix(x, nrow = M, byrow = T), xrange = c(0,1), yrange = c(0,1)))
log.int.im <- lapply(gp.im, function(x) trend.lgcp + x)

sapply(log.int.im, function(x) mean(exp(x)))

##========================================================
## Generate the points
##========================================================
pp <- lapply(c(1:n.sim), function(i) rpoispp(lambda = exp(log.int.im[[i]])))
sapply(pp, function(x) x$n)

##========================================================
## Generate the random effect W in the observation process
##========================================================
w <- lapply(c(1:n.sim), function(i) GSP.sim(1, mu = rep(mu.w, nrow(D)), 
                                            cov.fun = cov.fun, dists=D, theta=theta.w))
w.im <- lapply(c(1:n.sim), function(i) im(matrix(w[[i]], nrow = M, byrow = T), xrange = c(0,1), yrange = c(0,1)))

##========================================================
## Genrate the observations
##========================================================
mu.im <- lapply(c(1:n.sim), function(i) w.im[[i]] + alpha*(log.int.im[[i]]))

dat <- lapply(pp, function(pts) data.frame(x1 = pts$x, x2 = pts$y))
proxy <- lapply(dat, function(x) closest.coord(x, S))
for(i in c(1:n.sim)){
  dat[[i]]$w <- w[[i]][proxy[[i]]]
  dat[[i]]$log.int <- extract.im.2(cbind(dat[[i]]$x1, dat[[i]]$x2), log.int.im[[i]])
  dat[[i]]$mu <- dat[[i]]$w + alpha*(dat[[i]]$log.int)
  dat[[i]]$z <- rpois(length(dat[[i]]$mu), kappa*exp(dat[[i]]$mu))
  dat[[i]]$cell <- proxy[[i]]
}
D.obs <- lapply(pp, function(x) Edist(as.matrix(cbind(x$x, x$y))))

##============================================
## Map the underlying random process values
## to a set of prediction locations using
## the value from the nearest grid location
##============================================
dat.pred <- lapply(pp, function(pts) data.frame(x1 = S[,1], x2 = S[,2]))
for(i in c(1:n.sim)){
  dat.pred[[i]]$w <- w[[i]]
  dat.pred[[i]]$log.int <- extract.im.2(cbind(dat.pred[[i]]$x1, dat.pred[[i]]$x2), log.int.im[[i]])
  dat.pred[[i]]$mu <- dat.pred[[i]]$w + alpha*(dat.pred[[i]]$log.int) # underlying mu = y + alpha * log(lambda)
}

##============================================
## Create a data frame for Bayesian analysis
##============================================
grid.dat <- lapply(c(1:n.sim), function(i) data.frame(x1 = S[,1], x2 = S[,2]))
pp.grid <- lapply(c(1:n.sim), function(i) closest.coord(cbind(pp[[i]]$x,pp[[i]]$y), S, method = 'Manhattan'))

for(i in c(1:n.sim)){
  # The observed
  pts.cts <- tapply(pp.grid[[i]], pp.grid[[i]] , length)
  # grid.dat[[i]]$x1 <- S[,1]
  grid.dat[[i]]$n.pts <- 0
  grid.dat[[i]]$n.pts[as.numeric(names(pts.cts))] <- pts.cts
  
  # The unobserved
  grid.dat[[i]]$w <- w[[i]]
  grid.dat[[i]]$log.int <- extract.im.2(cbind(grid.dat[[i]]$x1, grid.dat[[i]]$x2), log.int.im[[i]])
  grid.dat[[i]]$mu <- grid.dat[[i]]$w + alpha*(grid.dat[[i]]$log.int) 

}


rm(GSP.beta.hat, GSP.beta.hat.cov, GSP.Diggle, GSP.exp, GSP.Gaussian,
   GSP.m2l, GSP.Matern, GSP.mle, GSP.powerexp, GSP.pred, GSP.sim,
   trend.lgcp, gp, gp.im)

save.image("../Data/Poisson_PS.Rdata")

