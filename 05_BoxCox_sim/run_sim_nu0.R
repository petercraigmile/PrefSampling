
source("init.R")

source("functions/sim_settings_Thomas.R")
source("functions/grid_2020.R")

for (chain.num in 1:200) {

    print(chain.num)
    
    alpha <- 1
    kappa <- 0
    
    parents <- Thomas.sim.parents(lambda0, the.window)
    lambda  <- Thomas.intensity(G$grid, parents, delta, omega2)
    
    pp      <- Thomas.sim(G, the.window, parents, delta, omega2)
    pp$eta <- log(pp$lambda)
    
    sim     <- GSP.pref.sim(pp, mu.W, theta.W, cov.fun.W, alpha, sigma2.Z, kappa)
    
    cell.area <- 1/(G$nx*G$ny)        
    
    chs <- init.Pref.Thomas(sim$Z, pp$locs, G$grid, parents, cell.area)
    
    chs <- fit.Pref.Thomas(chs, 5000, every=2500, burn.in=TRUE)
    
    chs <- fit.Pref.Thomas(chs, 10000, every=2500, thin=10)
    
    pdf(file=paste0("figures/nu0_trace_plots", chain.num, ".pdf"), width=8, height=6.5)
    
    source("trace_plots_Pref_Thomas.R")
    
    pm <- rowMeans(chs["y"])
    
    print(IMSPE(pm, sim$Y))
    print(POP(pm, sim$Y))
    
    dev.off()
    
    save(list=c("chs", "pp", "sim", "cell.area"),
         file=paste0("output/nu0_sim", chain.num, ".RData"))
    
    chs <- init.Pref.Thomas(sim$Z, pp$locs, G$grid, parents, cell.area)
    chs$kappa <- 0
    
    chs <- fit.Pref.Thomas.fix.kappa(chs, 5000, every=2500, burn.in=TRUE)
    
    chs <- fit.Pref.Thomas.fix.kappa(chs, 10000, every=2500, thin=10)
    
    pdf(file=paste0("figures/nu0_kappa0_trace_plots", chain.num, ".pdf"), width=8, height=6.5)
    
    source("trace_plots_Pref_Thomas.R")
    
    pm <- rowMeans(chs["y"])
    
    print(IMSPE(pm, sim$Y))
    print(POP(pm, sim$Y))
    
    dev.off()
    
    save(list=c("chs", "pp", "sim", "cell.area"),
         file=paste0("output/nu0_kappa0_sim", chain.num, ".RData"))
    
    gc()
    
}
