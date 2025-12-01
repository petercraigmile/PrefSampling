
nearest.neighbor.list <- function (grid, delta, number.nbs) {

    D <- Edist(grid)

    if (missing(delta)) {
        
        lapply(1:nrow(D), function (k)
            sort(order(D[k,])[2:(number.nbs+1)]))
    } else {

        lapply(1:nrow(D), function (k) {

            sel <- which(D[k,] <= delta)
            sel[sel!=k]
            })
    }
}


POP <- function (pred, truth) {

    mean( pred > truth )
}



IMSPE <- function (pred, truth) {

    mean( (pred - truth)^2 )
}



find.locs.in.grid <- function (locs, grid) {

    sapply(1:nrow(locs), function (k)
        which.min( colSums((t(grid) - locs[k,])^2) ))
}



make.pred.grid <- function (nx, ny) {
    ## Make a prediction grid of size nx x ny in the domain [0,1]^2

    xs <- seq(.5/nx, 1-.5/nx, length.out = nx)
    ys <- seq(.5/ny, 1-.5/ny, length.out = ny)
    
    grid  <- make.grid(list(xs, ys))

    list(grid=grid,
         xs=xs, ys=ys,
         nx=nx, ny=ny, M=nrow(grid))
}



psumm <- function (x, name, k, factor=1, digits=3) {

    y <- x[name]

    if (!missing(k)) {
        y <- y[k,]
    }

    z <- y * factor

    xx <- c(mean(z), quantile(z, c(0.025, 0.975)))
    names(xx) <- NULL
    round(xx, 3)
}
