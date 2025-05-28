

library(sf)

library(spatstat)
## The prepped data comes in various resolutions (40, 50, 60).
## The higher the number, the more densely distributed the grid cells.
## Computational complexity is higher when the resolution is higher.
load("Data/prepped_data_40.RData")
source("functions/GHCN_EDA_Functions.R")

sf_use_s2(F)
temp <- stations[,c(3,4)]
pdf(file = "figures/stations_tavg.pdf", width = 6.5, height = 4.5)
par(mfrow = c(1,1))
plot(states)
plot(AZ, add = T)
plot(CA, add = T)
plot(NV, add = T)
plot(temp, pch = 19, add = T)
Dev.off()


plot(stations$tavg)

## Assess the density
pdf(file = 'figures/density_estimates_four_states.pdf', width = 8, height = 6.5)
par(mfrow = c(1,1), mar = c(1,1,2,2))
plot(density(pp.states), asp = 1, box = F, main = "", ribsep = 0.05)
points(pp.states, pch =19, cex = 0.5, col = 'gray80')
title(main = 'Estimated Density of All Weather Stations', line = 0)
plot(AZ_win, add = T)
plot(CA_win, add = T)
plot(NV_win, add = T)
plot(UT_win, add = T)
dev.off()

pdf(file = 'figures/density_estimates_individual_states.pdf', width = 8, height = 6.5)
par(mfrow = c(2,2), mar = c(1,1,2,2))
plot(density(pp.AZ), asp = 1, box = F, main = "", ribsep = 0.05)
points(pp.AZ, pch =19, cex = 0.3, col = 'gray80')
points(density(pp.AZ)$xcol[71], density(pp.AZ)$yrow[46], pch = 17, cex = 2, col = 'black')
plot(AZ_win, add = T)
title(main = 'Estimated Density of Stations in AZ', line = 0)
plot(density(pp.CA), asp = 1, box = F, main = "", ribsep = 0.05)
points(pp.CA, pch =19, cex = 0.3, col = 'gray80')
points(density(pp.CA)$xcol[77], density(pp.CA)$yrow[16], pch = 17, cex = 1.75, col = 'black')
points(density(pp.CA)$xcol[26], density(pp.CA)$yrow[70], pch = 17, cex = 1.75, col = 'black')
plot(CA_win, add = T)
title(main = 'Estimated Density of Stations in CA', line = 0)
plot(density(pp.NV), asp = 1, box = F, main = "", ribsep = 0.05)
points(pp.NV, pch =19, cex = 0.3, col = 'gray80')
points(density(pp.NV)$xcol[2], density(pp.NV)$yrow[73], pch = 17, cex = 1.75, col = 'black')
plot(NV_win, add = T)
title(main = 'Estimated Density of Stations in NV', line = 0)
plot(density(pp.UT), asp = 1, box = F, main = "", ribsep = 0.05)
points(pp.UT, pch =19, cex = 0.3, col = 'gray80')
points(density(pp.UT)$xcol[61], density(pp.UT)$yrow[115], pch = 17, cex = 1.75, col = 'black')
points(density(pp.UT)$xcol[61], density(pp.UT)$yrow[100], pch = 17, cex = 1.75, col = 'black')
title(main = 'Estimated Density of Stations in UT', line = 0)
plot(UT_win, add = T)
dev.off()

pdf(file = 'figures/density_estimates_combined.pdf', width = 6.5, height = 8)
par(mfrow = c(3,2), mar = c(1,1,2,2))
plot(density(pp.states), asp = 1, box = F, main = "", ribsep = 0.05)
points(pp.states, pch =19, cex = 0.35, col = 'gray80')
title(main = 'Estimated Density of All Weather Stations', line = 0)
plot(AZ_win, add = T)
plot(CA_win, add = T)
plot(NV_win, add = T)
plot(UT_win, add = T)
plot(density(pp.AZ), asp = 1, box = F, main = "", ribsep = 0.05)
points(pp.AZ, pch =19, cex = 0.35, col = 'gray80')
points(density(pp.AZ)$xcol[71], density(pp.AZ)$yrow[46], pch = 17, cex = 1.75, col = 'black')
plot(AZ_win, add = T)
title(main = 'Estimated Density of Stations in AZ', line = 0)
plot(density(pp.CA), asp = 1, box = F, main = "", ribsep = 0.05)
points(pp.CA, pch =19, cex = 0.35, col = 'gray80')
points(density(pp.CA)$xcol[77], density(pp.CA)$yrow[16], pch = 17, cex = 1.75, col = 'black')
points(density(pp.CA)$xcol[26], density(pp.CA)$yrow[70], pch = 17, cex = 1.75, col = 'black')
plot(CA_win, add = T)
title(main = 'Estimated Density of Stations in CA', line = 0)
plot(density(pp.NV), asp = 1, box = F, main = "", ribsep = 0.05)
points(pp.NV, pch =19, cex = 0.35, col = 'gray80')
points(density(pp.NV)$xcol[2], density(pp.NV)$yrow[73], pch = 17, cex = 1.75, col = 'black')
plot(NV_win, add = T)
title(main = 'Estimated Density of Stations in NV', line = 0)
plot(density(pp.UT), asp = 1, box = F, main = "", ribsep = 0.05)
points(pp.UT, pch =19, cex = 0.35, col = 'gray80')
points(density(pp.UT)$xcol[61], density(pp.UT)$yrow[115], pch = 17, cex = 1.75, col = 'black')
points(density(pp.UT)$xcol[61], density(pp.UT)$yrow[100], pch = 17, cex = 1.75, col = 'black')
title(main = 'Estimated Density of Stations in UT', line = 0)
plot(UT_win, add = T)
dev.off()





##========================
##  Assess stationarity
##========================
library(geoR)
plot_vario <- function(which.state){
  vg.cloud <- variog(coords=cbind(stations$lon[which.state], 
                                  stations$lat[which.state]), 
                     data=stations$tavg[which.state])
  plot(vg.cloud, main = "", xlim = c(0, 15), ylim = c(0, 40), ylab = "Semivariogram",
       xlab = "Distance in degrees")
  #wls <- variofit(vg.cloud, ini.cov.pars = c(50, 5), cov.model="exp", fix.nugget=F, nugget=10,
   #               max.dist = vg.cloud$max.dist)
  # ML fit looks better.
  ml.fit <- likfit(coords = cbind(stations$lon[which.state], 
                                  stations$lat[which.state]),
                   data=stations$tavg[which.state],
                   cov.model="exp", fix.nugget=F,
                   nugget=10, ini.cov.pars=c(50,5))
  lines(ml.fit)
  #lines(wls, lty = 3)
  #ml.fit
}

stations.geodata <- as.geodata()

## look at the ML fit and think about how far the correlation goes.
## Think about the correlation within each state.
## How much correlation there is


pdf(file = 'figures/Variograms_Region.pdf', width = 6.5, height = 8)
par(mfrow = c(3,2), cex = 0.65)
#vg.cloud <- variog(coords=cbind(stations$lon, stations$lat), data=stations$tavg, op = "cloud")
#plot(vg.cloud)
fit.states <- plot_vario(1:nrow(stations))
title(main = "Estimated Semivariogram")
exp(-1/3.111)
exp(-2/3.111)
exp(-3/3.111)

plot_vario(which.AZ)
title(main = "Estimated Semivariogram for Arizona")
plot_vario(which.CA)
title(main = "Estimated Semivariogram for California")
plot_vario(which.NV)
title(main = "Estimated Semivariogram for Nevada")
plot_vario(which.UT)
title(main = "Estimated Semivariogram for Utah")
dev.off()



##========================
##  Plot K Functions
##========================

#pdf('figures/K_functions_States_Both.pdf', width = 12, height = 4.5)
#par(mfrow=c(1,2), cex=1, mar=c(3.1,3.5,1,0.5), mgp=c(1.8,0.5,0), bty="L")
#CSR.deviation.plot(pp.states)
#dev.off()

#pdf('../figures/K_functions_Individual_States.pdf', width = 12, height = 3.5)
#par(mfrow=c(1,4), cex=0.75, mar=c(3.1,3.5,1,0.5), mgp=c(1.8,0.5,0), bty="L")
#CSR.deviation.plot(AZ.pp)
#CSR.deviation.plot(CA.pp)
#CSR.deviation.plot(NV.pp)
#CSR.deviation.plot(UT.pp)
#dev.off()


Kest.AZ <- Kest(pp.AZ)
Kest.CA <- Kest(pp.CA)
Kest.NV <- Kest(pp.NV)
Kest.UT <- Kest(pp.UT)
Kest.All <- Kest(pp.states)

lgcp.AZ <- lgcp.estK(pp.AZ)
lgcp.CA <- lgcp.estK(pp.CA)
lgcp.NV <- lgcp.estK(pp.NV)
lgcp.UT <- lgcp.estK(pp.UT)
lgcp.All <- lgcp.estK(pp.states)


pdf('figures/K_functions_combined.pdf', width = 6.5, height = 8)
par(mfrow=c(3,2), cex=0.75, mar=c(3.1,3.5,1,0.5), mgp=c(1.8,0.5,0), bty="L")
ylim = c(-0.3, 2)
plot((Kest.AZ$iso - Kest.AZ$theo) ~ Kest.AZ$r, 
     xlim = c(0, 1.5), ylim = c(-0.2, 2), type = "l", lwd = 2,
     xlab = 'Distance in degrees', ylab = expression(hat(K)(r) - K[pois](r)), lty = 5, col = 'purple')
lines((Kest.CA$iso[1:330] - Kest.CA$theo[1:330]) ~ Kest.CA$r[1:330], col = 'gray40', lwd = 2, lty = 2)
lines((Kest.NV$iso - Kest.NV$theo) ~ Kest.NV$r, col = 'blue', lwd = 2, lty = 3)
lines((Kest.UT$iso - Kest.UT$theo) ~ Kest.UT$r, col = 'red', lwd = 2, lty = 4)
lines((Kest.All$iso - Kest.All$theo) ~ Kest.All$r, col = 'black', lwd = 2, lty = 1)
abline(h = 0, lty = 3)
legend('topleft', lty = c(1,2,3,4,5), col = c('black', 'gray40', 'blue', 'red','purple'), 
       legend = c('All', 'Arizona', 'California', 'Nevada', 'Utah'),box.lwd = 0)
plot((Kest.All$iso - Kest.All$theo) ~ Kest.All$r, 
     xlim = c(0, 1.5), ylim = c(-0.2, 2), type = "l", lwd = 2,
     xlab = 'Distance in degrees', ylab = expression(hat(K)(r) - K[pois](r)), main = "")
lines(lgcp.All$fit$fit - lgcp.All$fit$theo ~ lgcp.All$fit$r, lwd = 2, lty = 2, col = 'gray40')
CSR.deviation.plot(pp.states)
title("All Four States", line = 0)
plot((Kest.AZ$iso - Kest.AZ$theo) ~ Kest.AZ$r, 
     xlim = c(0, 1.5), ylim = c(-0.2, 2), type = "l", lwd = 2,
     xlab = 'Distance in degrees', ylab = expression(hat(K)(r) - K[pois](r)), main = "")
lines(lgcp.AZ$fit$fit - lgcp.AZ$fit$theo ~ lgcp.AZ$fit$r, lwd = 2, lty = 2, col = 'gray40')
title("Arizona", line = 0)
CSR.deviation.plot(pp.AZ)
plot((Kest.CA$iso - Kest.CA$theo) ~ Kest.CA$r, 
     xlim = c(0, 1.5), ylim = c(-0.2, 2), type = "l", lwd = 2,
     xlab = 'Distance in degrees', ylab = expression(hat(K)(r) - K[pois](r)), main = "")
lines(lgcp.CA$fit$fit - lgcp.CA$fit$theo ~ lgcp.CA$fit$r, lwd = 2, lty = 2, col = 'gray40')
title("California", line = 0)
CSR.deviation.plot(pp.CA)
plot((Kest.NV$iso - Kest.NV$theo) ~ Kest.NV$r, 
     xlim = c(0, 1.5), ylim = c(-0.2, 2), type = "l", lwd = 2,
     xlab = 'Distance in degrees', ylab = expression(hat(K)(r) - K[pois](r)), main = "")
lines(lgcp.NV$fit$fit - lgcp.NV$fit$theo ~ lgcp.NV$fit$r, lwd = 2, lty = 2, col = 'gray40')
title("Nevada", line = 0)
CSR.deviation.plot(pp.NV)
plot((Kest.UT$iso - Kest.UT$theo) ~ Kest.UT$r, 
     xlim = c(0, 1.5), ylim = c(-0.2, 2), type = "l", lwd = 2,
     xlab = 'Distance in degrees', ylab = expression(hat(K)(r) - K[pois](r)), main = "")
lines(lgcp.UT$fit$fit - lgcp.UT$fit$theo ~ lgcp.UT$fit$r, lwd = 2, lty = 2, col = 'gray40')
title("Utah", line = 0)
CSR.deviation.plot(pp.UT)
dev.off()





##========================
## Plot the K functions
##========================
Kest.pps <- lapply(pp.list, function(x) Kest(x))
pdf('../figures/K_functions_States.pdf', width = 6.3, height = 6.5)
par(mfrow=c(1,1), cex=0.75, mar=c(3.1,3.1,1,0.5), mgp=c(1.8,0.5,0), bty="L")
cols <- c('black', 'red', 'blue', 'green', 'purple')
plot(1, type = "n", xlab = 'r', 
     ylab = 'K(r)', xlim = c(0, 2.8),  
     ylim = c(0, 2.3)) 
for(i in 1:5){
  lines(Kest.pps[[i]]$iso-Kest.pps[[i]]$theo ~ Kest.pps[[i]]$r, lty = 1, col = cols[i])
  gw <- Gauss.weights(0, max(Kest.pps[[i]]$r), 256)
  absc <- c(0,gw$abscissa)
  wts <- c(0,gw$weights)  
  lines(K.theo.lgcp(Kest.pps[[i]]$r,theta = log(LGCP.pars[c(2:3),i])) - Kest.pps[[i]]$theo ~ Kest.pps[[i]]$r,
        lty = 3, col = cols[i], lwd = 2)
  lines(K.theo.Thomas(Kest.pps[[i]]$r,theta = log(Thomas.pars[c(4:3),i])) - Kest.pps[[i]]$theo ~ Kest.pps[[i]]$r,
        lty = 3, col = cols[i], lwd = 1)
}
legend(0, 2, legend = c('All', 'AZ', 'CA', 'NV', 'UT'), lty = 1, col = cols)
dev.off()

#write.csv(LGCP.pars, 'LGCP_Pars.csv')









