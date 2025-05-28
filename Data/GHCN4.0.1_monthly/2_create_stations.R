

stations <- read.table("Derived/US_stations.txt", header=FALSE)

library(maps)
map("usa")

US.id   <- stations[,1]
US.lat  <- stations[,2]
US.long <- stations[,3]
US.alt  <- stations[,4]
US.name <- stations[,5]

continental <- US.long < -50 & US.long > -130

cont.US.stations <- data.frame(id   = US.id,
                               lat  = US.lat,
                               long = US.long,
                               alt  = US.alt,
                               name = US.name)[continental,]

save(cont.US.stations, file="Derived/cont_US_stations.RData")


library(maps)

pdf(file="figures/US_cont_stations.pdf", width=4.5, height=3)
par(mfrow=c(1,1), cex=0.65, mar=c(3.1,3.1,1,0.5), mgp=c(1.8,0.5,0), bty="L")

plot(cont.US.stations$long, cont.US.stations$lat,
     pch=20, cex=0.15, xlab="", ylab="")
map("state", add=T, col="gray")

dev.off()
