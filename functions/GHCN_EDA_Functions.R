

CSR.deviation.plot <- function (pp, K=100, begin = 0, end = 360) {
    ## ======================================================================
    ## Produces a plot of the deviation of the estimated K function for
    ## the point pattern 'pp' from CSR.  Also shows the deviations for K
    ## simulations from a CSR process with a value of the intensity
    ## estimated from 'pp'.
    ## Peter F. Craigmile, pfc@stat.osu.edu
    ## Modified by Rui Qiang:
    ## Now allows for directional K-function comparisons to assess potential
    ## non-stationarity
    ## ======================================================================

    require(spatstat)
    
    ## Estimate lambda from the data, assuming CSR 
    lam.hat <- (length(pp$x)-1)/area(pp)
    lam.hat
    
    ## Estimate the K function
    K.est <- Ksector(pp, begin, end, correction="Ripley")

    ## Calculate the K function deviation from CSR
    K.dev <- K.est$iso - K.est$theo

    ## Simulate a K Poisson processes in the same domain ('window') as the original data    
    sims <- lapply(1:K, function (k) rpoispp(lam.hat, win=pp))
    
    ## For each simulation estimate the K function
    Kfuns <- lapply(sims, function (z) Kest(z, begin, end, correction="Ripley"))

    ## For each simulation estimate the K function deviation
    Kdevs <- lapply(Kfuns, function (z) z$iso-z$theo)

    ## Calculate the range of the simulated deviations
    sim.range <- range(unlist(Kdevs))
    
    main.title <- ""
    if(end-begin != 360){
      main.title <- paste("Angles: ", begin, " to ", end)
    }
    ## Set up the plot
    #plot(K.est$r, K.dev,
    #     xlab="r", ylab=expression(hat(K)(r) - K[pois](r)),
    #     ylim=c(min(K.dev, sim.range[1]), max(K.dev, sim.range[2])), type="n",
    #     main = main.title)
    

    ## Add the simulated deviations
    sapply(1:K, function (k) lines(Kfuns[[k]]$r, Kdevs[[k]], col="gray"))

    ## Add a horizontal line 
    abline(h=0, lty=2)

    ## Show the K function estimated from the data    
    lines(K.est$r, K.dev, lwd=2)
}


## A generic function to fit any cluster process using kppm()
fit.cluster <- function(pp, cluster){
  if(cluster == 'Poisson'){
    ppm(pp ~ x + y + logpopden, data = list(logpopden = log_popden.image))
  }else{
   kppm(pp, trend = ~ x + y + logpopden, 
                   covariates = list(logpopden = log_popden.image),
                   clusters = cluster , method = 'palm',
                   statistic = 'K')
  }
}



fit.deviation.plot<- function (pp, fit, K=100, begin = 0, end = 360) {
    require(spatstat)
    
    ## Estimate the K function
    K.est <- Ksector(pp, begin, end, correction="Ripley")

    ## Calculate the K function deviation from CSR
    K.dev <- K.est$iso - K.est$theo

    ## Simulate a K point processes in the same domain ('window') as the original data    
    ## These simulated point processes have the same class as the fitted one.
    ## Both Poisson (inhomogeneous) and cluster point processes are allowed.
    sims <- sapply(1:K, function (k) simulate(fit))
    
    ## For each simulation estimate the K function
    Kfuns <- lapply(sims, function (z) Kest(z, begin, end, correction="Ripley"))

    ## For each simulation estimate the K function deviation
    Kdevs <- lapply(Kfuns, function (z) z$iso-z$theo)

    ## Calculate the range of the simulated deviations
    sim.range <- range(unlist(Kdevs))
    
    main.title <- ""
    if(end-begin != 360){
      main.title <- paste("Angles: ", begin, " to ", end)
    }
    ## Set up the plot
    plot(K.est$r, K.dev,
         xlab="r", ylab=expression(hat(K)(r) - K[pois](r)),
         ylim=c(min(K.dev, sim.range[1]), max(K.dev, sim.range[2])), type="n",
         main = main.title)

    ## Add the simulated deviations
    sapply(1:K, function (k) lines(Kfuns[[k]]$r, Kdevs[[k]], col="gray"))

    ## Add a horizontal line 
    abline(h=0, lty=2)

    ## Show the K function estimated from the data    
    lines(K.est$r, K.dev, lwd=2)
}


## Plot the estimated K function of the point process and 
## K functions of simulated point processes under the same fit.
plot.k <- function(nbreaks, pp, myfit, K = 100, process.name = myfit$clusters){
  if(length(process.name) == 0){
    process.name <- "Poisson"
  }
  pdf(file= paste0("figures/",process.name, "_Fit_Directional_K_Function.pdf"), width=6.5, height= 8)
  par(mfrow=c(2,2), cex=0.65, mar=c(5, 3, 2, 2), mgp=c(3, 1, 0), bty="L")
  range <- 360/nbreaks
  for(i in 0:(nbreaks-1)){
    fit.deviation.plot(pp, myfit, K, begin = i*range, end = (i+1)*range)
  }
  dev.off()
} 



## Simulate from a fitted point process and plot them
plot.sim <- function(pp, myfit, n.sim = 9, arrange = c(3,3), process.name = myfit$clusters, cutoff = 10){
  accept.sim <- 0
  if(length(process.name) == 0){
    process.name <- "Poisson"
  }
  pdf(file= paste0("figures/",process.name, "_Simulated.pdf"), width=6.5, height= 8)
  par(mfrow=arrange, cex=0.65, mar=c(0, 0.5, 1.5, 0.5), mgp=c(1, 1, 0), bty="L")
  plot(pp, pch = 19, cex = 0.6)
  while(accept.sim < n.sim){
    sim.myfit <- simulate(myfit)
    if(abs(sim.myfit[[1]]$n - pp$n) < cutoff){
      plot(sim.myfit[[1]], pch = 19, main = "", cex = 0.6)
      title(paste0("n = ", sim.myfit[[1]]$n), line = -4, cex.main = 1.5)
      accept.sim <- accept.sim + 1
      cat(accept.sim, "..")
      }  
  }
  dev.off()
}

sim.window <- function(pp, cls, fit.win, eval.win, K = 20, min.area, adj, ct_cls = 1, plot_map){
  require(spatstat)
    ## the point process
    pp.fit <- fit.cluster(as.ppp(ppp(pp$x, pp$y, fit.win)), cls)
    ## Estimate the K function on the windows
    K.est <- lapply(eval.win, function(w) Kest(as.ppp(ppp(pp$x, pp$y, window = w)), correction = "Ripley"))
    K.dev <- lapply(K.est, function(k) k$iso - k$theo)
    sims_full <- sapply(1:K, function (k) simulate(pp.fit))
    is <- which(sapply(eval.win, function(z) area(z)) > min.area)
    names <- is
    if(!is.null(adj)){
      if(adj == 'states' & length(eval.win) == 1){
        names <- 'AZ CA NV UT'
      }else{
       names <- c('AZ', 'CA', 'NV', 'UT') 
      }
    }
    map.ct <- 0
    for(i in is){
      map.ct <- map.ct + 1
      if(adj != 'states' & length(eval.win) == 1){
        j = which(names == adj)
        i = 1
      }else{
        j = i
      }
      if(ct_cls == 1 & plot_map){
        pp.win <- as.ppp(ppp(pp$x, pp$y, eval.win[[i]]))
        plot(pp.win, pch = 19, cex = 0.75, main = paste("Window", names[j], ", n = ", pp.win$n))
      }
      sims <- lapply(sims_full, function(z) as.ppp(ppp(z$x, z$y, window = eval.win[[i]])))
      avg.n.window <- mean(unlist(lapply(sims, function(z) z$n)))
      ## For each simulation estimate the K function in the same window
      Kfuns <- lapply(sims, function (z) Kest(ppp(z$x, z$y, window = eval.win[[i]]), correction="Ripley"))
      ## For each simulation estimate the K function deviation
      Kdevs <- lapply(Kfuns, function (z) z$iso-z$theo)
      ## Calculate the range of the simulated deviations
      sim.range <- range(unlist(Kdevs))
      ## Set up the plot
      plot(K.est[[i]]$r, K.dev[[i]], 
         xlab="r", ylab=expression(hat(K)(r) - K[pois](r)),
         ylim=c(-0.2, 4), type="n",
         main = paste("Average simulated points = ", round(avg.n.window)))
      legend('topright', legend = cls)
    ## Add the simulated deviations
    sapply(1:K, function (k) lines(Kfuns[[k]]$r, Kdevs[[k]], col="gray"))

    ## Add a horizontal line 
    abline(h=0, lty=2)

    ## Show the K function estimated from the data    
    lines(K.est[[i]]$r, K.dev[[i]], lwd=2)
    }
}

plot.K.window <- function(pp, cls, fit.win, eval.win, K = 20, min.area = 5, 
                          adj = NULL, more_adj = NULL, one.page = FALSE, plot_map = TRUE){
  process.name <- cls
  pdf(file= paste0("figures/",paste(process.name, collapse = '_'), '_',"_Windows_", adj, "_",more_adj,".pdf"), width=9 , height = 6.5)
  par(mfrow=c(1,1), cex=0.65, mar=c(5, 3, 2, 2), mgp=c(3, 1, 0), bty="L")
  plot(pp, pch = 19)
  for(i in 1:length(eval.win)){
    plot(eval.win[[i]], add = T, border = 'red')
    text(x = mean(eval.win[[i]]$xrange), y = mean(eval.win[[i]]$yrange), labels = i, cex = 2, col = 'red')
  }
  dev.off()
  if(one.page){
    pdf(file= paste0("figures/",paste(process.name, collapse = '_'), '_', length(eval.win) ,"_Windows_", adj, "_",more_adj,"_K_Functions_.pdf"), 
        width= 6.5 , height = 8)
    par(mfrow=c(length(eval.win),ifelse(plot_map, (length(cls)+1), length(cls))), cex=0.5, mar=c(5, 3, 2, 2), mgp=c(3, 1, 0), bty="L")
  }else{
    pdf(file= paste0("figures/",paste(process.name, collapse = '_'), '_', length(eval.win) ,"_Windows_", adj,"_",more_adj, "_K_Functions_.pdf"), 
        width= 7 + length(cls) , height = 4.5-length(cls)/2.5)
    par(mfrow=c(1,(length(cls)+1)), cex=0.65, mar=c(5, 3, 2, 2), mgp=c(3, 1, 0), bty="L")
  }
  for(i in 1:length(eval.win)){
    for(j in 1:length(cls)){
      if(length(fit.win) == 1){
          sim.window(pp, cls[j], fit.win[[1]], eval.win[i], K, min.area, adj, ct_cls = j, plot_map)
      }else{
          sim.window(pp, cls[j], eval.win[[i]], eval.win[i], K, min.area, adj, ct_cls = j, plot_map)
        }
    }
  }
  dev.off()
}






create.windows <- function(bigwindow, nrow, ncol){
  x <- seq(bigwindow$xrange[1]-0.1, bigwindow$xrange[2]+0.1,length.out = ncol + 1)
  y <- seq(bigwindow$yrange[1]-0.1, bigwindow$yrange[2]+0.1,length.out = nrow + 1)
  
  x1 <- rep(x[1:ncol], nrow)
  x2 <- rep(x[2:(ncol+1)], nrow)
  xs <- t(rbind(x1, x2, x2, x1, x1))
  y1 <- sort(rep(y[1:nrow], ncol))
  y2 <- sort(rep(y[2:(nrow+1)], ncol))
  ys <- t(rbind(y1, y1, y2, y2, y1))

  windows <- list()
  for(i in 1:(nrow*ncol)){
    mywin <- owin(poly = cbind(xs[i,], ys[i,]))
    mywin <- intersect.owin(mywin, bigwindow)
    windows[[i]] <- mywin
  }
  windows
}










