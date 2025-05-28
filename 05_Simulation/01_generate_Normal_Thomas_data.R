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
Rcpp::sourceCpp("../functions/GSP_2024_04_11.cpp")
source("../functions/GSP_2024_04_11.R")
source("../functions/sim_functions.R")


##================================================================
## Parameter Specification
##================================================================
cov.fun <- GSP.powerexp
n.sim <- 20
cov.power <- 1
theta.w <- c(2, 0.33, 1) # partial sill, range, power
mu.w <- 10
alpha <- 1
tau2.true <- 1

##========================
## resolution
##========================
M  <- 20
xs <- seq(.5/M, 1-.5/M, length.out = M)

##========================
## coords
##========================
S  <- make.grid(list(xs, xs))

##========================
## distance matrix
##========================
D <- Edist(S)

##========================
## Generate the point process 
##========================
a <- 0
b <- 0
parents.lambda <- 30
delta <- 50 
omega2 <- 0.1^2
omega2.true <- omega2
delta.true <- delta
parents <- lapply(c(1:n.sim), function(i) rpoispp(parents.lambda))
Thomas.int <- lapply(parents, function(x) lambda(S, x, delta, omega2))
int.im <- lapply(Thomas.int, function(x) im(matrix(x, nrow = M, byrow = T), xrange = c(0,1), yrange = c(0,1)))
log.int.im <- lapply(int.im, function(x) log(x))
pp <- lapply(int.im, function(x) rpoispp(x))

##========================
## Check the number of points
##========================
sapply(pp, function(x) x$n)

alpha.true <- alpha

##========================================
## Generate the random process
## in the observation process
##========================================
w <- lapply(c(1:n.sim), function(i) GSP.sim(1, mu = rep(mu.w, nrow(D)), 
                                            cov.fun = cov.fun, dists=D, theta=theta.w))
w.im <- lapply(c(1:n.sim), function(i) im(matrix(w[[i]], nrow = M, byrow = T), xrange = c(0,1), yrange = c(0,1)))
mu.im <- lapply(c(1:n.sim), function(i) w.im[[i]] + alpha*(log.int.im[[i]]))

##========================================
## Distance matrix for the observations
##========================================
D.obs <- lapply(pp, function(x) Edist(as.matrix(cbind(x$x, x$y))))

##========================================
## Generate the observations
##========================================
dat <- lapply(pp, function(pts) data.frame(x1 = pts$x, x2 = pts$y))
proxy <- lapply(dat, function(x) closest.coord(x, S))
for(i in c(1:n.sim)){
  dat[[i]]$w <- w[[i]][proxy[[i]]]
  dat[[i]]$log.int <- extract.im.2(cbind(dat[[i]]$x1, dat[[i]]$x2), log.int.im[[i]])
  dat[[i]]$mu <- dat[[i]]$w + alpha*(dat[[i]]$log.int) 
  dat[[i]]$z <- rnorm(nrow(dat[[i]]), 0, sd = sqrt(tau2.true)) + dat[[i]]$mu 
  dat[[i]]$cell <- proxy[[i]]
}

##========================================
## Map the underlying random process values
## to a set of prediction locations using
## the value from the nearest grid location
##========================================
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
  grid.dat[[i]]$mu <- grid.dat[[i]]$w + alpha*(grid.dat[[i]]$log.int) # underlying mu = y + alpha * log(lambda)

}



rm(GSP.beta.hat, GSP.beta.hat.cov, GSP.Diggle, GSP.exp, GSP.Gaussian,
   GSP.m2l, GSP.Matern, GSP.mle, GSP.powerexp, GSP.pred, GSP.sim, pp.grid)

save.image("../Data/Thomas_PS.RData")

