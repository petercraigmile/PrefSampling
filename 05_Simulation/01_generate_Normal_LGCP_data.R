##========================================##
## Title: Generate Data for the Simulation Experiment
## What: Generate 20 different point processes
## with observations 
## Author: Rui Qiang
## Date: 5.2.2025
##========================================##
library(spatstat)
library(parallel)
library(mvtnorm)
## all the necessary universal functions
Rcpp::sourceCpp("../functions/GSP_2024_04_11.cpp")
source("../functions/GSP_2024_04_11.R")

source("../functions/sim_functions.R")

cov.fun <- GSP.powerexp
n.sim <- 5

##================================================================
## Covariance power for power exponential correlation function
## Set to 1 for exponential. Set to 2 for Gaussian.
##================================================================
cov.power <- 1

##============================================================================
## theta.eta: Covariance parameters for random effect Eta (point process)
## theta.w: Covariance parameters for random effect W (observation process)
## In the order of "partial sill \sigma^2", 
## "range parameter \phi", 
## "power (for power exponential)"
##============================================================================
theta.eta <- c(2, 0.15, 1.75)  
theta.w <- c(2, 0.33, 1) 

##====================================
## Mean of the observation process
##====================================
mu.w <- 10

##============================================
## Preferential Sampling Parameter \alpha
##============================================
alpha <- 1 

##============================================
## Variance of observation error \tau^2
##============================================
tau2.true <- 1

##============================================
## Linear trend and mean of the random process \eta
##============================================
a <- 0
b <- 0
mu.eta <- 4.25 # 5.4 for homo

##============================================
## Resolution of the gridded spatial domain 
##============================================
M  <- 20

##============================================
## Make the grid
##============================================
xs <- seq(.5/M, 1-.5/M, length.out = M)
S  <- make.grid(list(xs, xs))

##============================================
## Distance matrix
##============================================
D <- Edist(S)

## For reference, so we don't accidently change
## anything down the road.
phi.true <- theta.eta[2]
alpha.true <- alpha

##============================================
## Generate the LGCP (as a list of images)
##============================================
trend.lgcp <- im(matrix(mu.eta + a*S[,1] + b*S[,2], nrow = M, byrow = T), xrange = c(0,1), yrange = c(0,1))
gp <- lapply(c(1:n.sim), function(i) GSP.sim(1, mu = rep(0, nrow(D)), cov.fun = cov.fun, dists=D, theta=theta.eta))
gp.im <- lapply(gp, function(x) im(matrix(x, nrow = M, byrow = T), xrange = c(0,1), yrange = c(0,1)))
log.int.im <- lapply(gp.im, function(x) trend.lgcp + x)

##============================================
## Check the mean of the intensity
##============================================
sapply(log.int.im, function(x) mean(exp(x)))
mean(sapply(log.int.im, function(x) mean(exp(x))))

##============================================
## Generate the points and check the numbers
##============================================
pp <- lapply(c(1:n.sim), function(i) rpoispp(lambda = exp(log.int.im[[i]])))
sapply(pp, function(x) x$n)
mean(sapply(pp, function(x) x$n))


##============================================
## Generate the random process W
##============================================
w <- lapply(c(1:n.sim), function(i) GSP.sim(1, mu = rep(mu.w, nrow(D)), 
                                            cov.fun = cov.fun, dists=D, theta=theta.w))
w.im <- lapply(c(1:n.sim), function(i) im(matrix(w[[i]], nrow = M, byrow = T), xrange = c(0,1), yrange = c(0,1)))
mu.im <- lapply(c(1:n.sim), function(i) w.im[[i]] + alpha*(log.int.im[[i]]))

##============================================
## Distance matrix for the points
## As a list of distance matrices
##============================================
D.obs <- lapply(pp, function(x) Edist(as.matrix(cbind(x$x, x$y))))

##============================================
## Generate the observed values
##============================================
dat <- lapply(pp, function(pts) data.frame(x1 = pts$x, x2 = pts$y))
proxy <- lapply(dat, function(x) closest.coord(x, S))
for(i in c(1:n.sim)){
  dat[[i]]$w <- w[[i]][proxy[[i]]]
  dat[[i]]$log.int <- extract.im.2(cbind(dat[[i]]$x1, dat[[i]]$x2), log.int.im[[i]])
  dat[[i]]$mu <- dat[[i]]$w + alpha*(dat[[i]]$log.int) 
  dat[[i]]$z <- rnorm(nrow(dat[[i]]), 0, sd = sqrt(tau2.true)) + dat[[i]]$mu # observations with error
  dat[[i]]$cell <- proxy[[i]]
}

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
  grid.dat[[i]]$n.pts <- 0
  grid.dat[[i]]$n.pts[as.numeric(names(pts.cts))] <- pts.cts
  
  # The unobserved
  grid.dat[[i]]$w <- w[[i]]
  grid.dat[[i]]$log.int <- extract.im.2(cbind(grid.dat[[i]]$x1, grid.dat[[i]]$x2), log.int.im[[i]])
  grid.dat[[i]]$mu <- grid.dat[[i]]$w + alpha*(grid.dat[[i]]$log.int)
}

##============================================
## Remove functions that's not used for the subsequent steps
##============================================
rm(GSP.beta.hat, GSP.beta.hat.cov, GSP.Diggle, GSP.exp, GSP.Gaussian,
   GSP.m2l, GSP.Matern, GSP.mle, GSP.powerexp, GSP.pred, GSP.sim,
   trend.lgcp, gp, gp.im, pp.grid)

save.image("../Data/Normal_PS.RData")

