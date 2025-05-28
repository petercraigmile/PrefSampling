##========================================##
## Title: sim functions
## all functions needed for simulating the points
## Author: Rui Qiang
## Date: 01.23.2024
##========================================##

closest.coord <- function(coords, grid, method = 'Euclidean'){
  if(method == 'Euclidean'){
    return(simplify2array(sapply(1:nrow(coords), function(i) (which.min((grid[,1]-coords[i,1])^2 + (grid[,2]-coords[i,2])^2))))) 
  }
  if(method == 'Manhattan'){
    return(simplify2array(sapply(1:nrow(coords), function(i) (which.min(abs(grid[,1]-coords[i,1]) + abs(grid[,2]-coords[i,2]))))))
  }
}

extract.im.2 <- function(coords, image){
  xs <- image$xcol
  ys <- image$yrow
  im.grid <- make.grid(list(xs, ys))
  v.vector <- as.vector(t(image$v))
  v.vector[closest.coord(coords, im.grid)]
}


z_xy <- function(x,y, parents, omega2, delta){
  # Return the density at one location
  delta*sum(dnorm(sqrt((rep(x,parents$n)-parents$x)^2 + (rep(y,parents$n)-parents$y)^2), sd = sqrt(omega2)))
}

lambda <- function(S, parents, delta, omega2){
  # Return the density at all locations in S
  sapply(c(1:nrow(S)), function(i) z_xy(S[i,1], S[i,2], parents, omega2, delta))
}

