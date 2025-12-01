

source("init.R")

## Define the simulation numbers
ks <- 1:200


## Summarize the nu=0.25 sims

SIMS <- 10

POPs <- IMSPEs <- c()

for (k in 1:length(ks)) {

    fname <- paste0("output/nu0_sim", ks[k], ".RData")
    
    load(fname)
           
    pm <- rowMeans(chs["y"])
    
    POPs[k] <- POP(pm, sim$Y)
    IMSPEs[k] <- IMSPE(pm, sim$Y)
}


fix.POPs <- fix.IMSPEs <- c()

for (k in 1:length(ks)) {

    fname <- paste0("output/nu25_kappa0_sim", ks[k], ".RData")
    
    load(fname)
           
    pm <- rowMeans(chs["y"])
    
    fix.POPs[k] <- POP(pm, sim$Y)
    fix.IMSPEs[k] <- IMSPE(pm, sim$Y)
}

summ <- function (x, digits=3) {

    tt <- t.test(x)
    aa <- round(c(tt$est, tt$conf.int), digits)
    names(aa) <- NULL
    aa
}


cat(length(POPs), "sims", "\n")

print( rbind(summ(POPs),
             summ(fix.POPs),
             summ(IMSPEs),
             summ(fix.IMSPEs)) )
      

## Summarize the nu=0 sims

SIMS <- 10

POPs <- IMSPEs <- c()

for (k in 1:length(ks)) {

    fname <- paste0("output/nu0_sim", ks[k], ".RData")
    
    load(fname)
           
    pm <- rowMeans(chs["y"])
    
    POPs[k] <- POP(pm, sim$Y)
    IMSPEs[k] <- IMSPE(pm, sim$Y)
}


fix.POPs <- fix.IMSPEs <- c()

for (k in 1:length(ks)) {

    fname <- paste0("output/nu0_kappa0_sim", ks[k], ".RData")
    
    load(fname)
           
    pm <- rowMeans(chs["y"])
    
    fix.POPs[k] <- POP(pm, sim$Y)
    fix.IMSPEs[k] <- IMSPE(pm, sim$Y)
}

summ <- function (x, digits=3) {

    tt <- t.test(x)

    aa <- round(c(tt$est, tt$conf.int), digits)
    names(aa) <- NULL
    aa
}


cat(length(POPs), "sims", "\n")

print( rbind(summ(POPs),
             summ(fix.POPs),
             summ(IMSPEs),
             summ(fix.IMSPEs)) )
      
