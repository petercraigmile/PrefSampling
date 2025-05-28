##====================================
## Prepare the GHCN data for spatial analysis
## and model fitting.
## 
## Author: Rui Qiang
## Email: rui.qiang@outlook.com
## 
## Date: 05.05.2025
##====================================

##============================
## Load the packages
##============================
library(spatstat)
library(geoR)
library(rgeos)
library(fields)
library(parallel)
library(sf)
library(sp)

##============================
## Set planar geometry
##============================
sf_use_s2(FALSE)

##========================================================
## Load the GHCN data and the shapefiles for the states. 
## The RData files are a comprehensive collection of 
## all 8895 stations in the United States.
##========================================================
load('../Data/cont_US_stations.RData')
load('../Data/cont_US_monthly_tavg.RData')
Nevada <- st_transform(st_read("../Data/shapefiles/tl_rd22_32_cousub/tl_rd22_32_cousub.shp"), "+proj=longlat")
California <- st_transform(st_read("../Data/shapefiles/tl_rd22_06_cousub/tl_rd22_06_cousub.shp"), "+proj=longlat")
Arizona <- st_transform(st_read("../Data/shapefiles/tl_rd22_04_cousub/tl_rd22_04_cousub.shp"), "+proj=longlat")
Utah <- st_transform(st_read("../Data/shapefiles/tl_rd22_49_cousub/tl_rd22_49_cousub.shp"), "+proj=longlat")

## Combine the county shapes to form a state shape
CA <- st_union(California$geometry)
NV <- st_union(Nevada$geometry)
AZ <- st_union(Arizona$geometry)
UT <- st_union(Utah$geometry)

CA.sep <- st_cast(CA, "POLYGON")
CA <- CA.sep[[1]]
NV <- NV[[1]]
AZ <- AZ[[1]]
UT <- UT[[1]]
CANV <- st_union(CA, NV)
AZUT <- st_union(AZ,UT)

states <- st_union(CANV,AZUT)
rm(California, Nevada, Arizona, Utah, CANV, AZUT)

##========================================================
## Convert the stations data to sf geometry
##========================================================
stations <- st_as_sf(SpatialPoints(rbind(cbind(cont.US.stations$long, cont.US.stations$lat))))
stations <- stations[(!duplicated(stations$geometry)),]

##========================================================
## Create labels for which stations lie within 
## the 4 states.
## Extract the useful information accordingly.
##========================================================
inds <- st_intersects(stations, states)

sel <- which(sapply(inds, function (x) x[1]) == 1)
stations <- stations[sel,]
temp <- cont.US.monthly.tavg[sel]
coords <- st_coordinates(stations$geometry)

stations <- cbind(stations, coords)
stations.all <- stations
rm(cont.US.monthly.tavg, cont.US.stations)
rm(sel, inds)

##============================================
## Calculate decade-average of the temperature
##============================================
start <- 2005
n.years <- 10
end <- start + n.years

n.in <- rep(0, length(temp))
tavg <- rep(NA, length(temp))

for(i in 1:length(temp)){
  ## See which months it has some record
  sel <- (temp[[i]]$year >= start) & (temp[[i]]$year < end)
  
  ## See how many months of data is not NA in this period
  n.in[i] <- sum(!is.na(temp[[i]]$tavg[sel]))
  
  ## The months that has some record
  time.stamp <- temp[[i]]$year[sel]
  
  ## The indices (1-120) that these months correspond to
  fill <- round((time.stamp - round(start-1/12, 3))*12,1)
  t.v <- rep(NA, 120)
  t.v[fill] <- temp[[i]]$tavg[sel]
  t.m <- matrix(t.v, nrow = 12, ncol = n.years)
  tavg[i] <- mean(rowMeans(t.m, na.rm = T), na.rm = T)
}

## Only keep stations that have at least 50% of data complete
sel <- n.in/120 >= 0.5
sum(sel)

## Percentage of complete data over the decade (120 months)
n.in <- round(n.in[sel]/120,2)

## Select stations that has some data
stations <- stations[sel,]
coords.all <- coords
coords <- coords[sel,]

## The temperature averages
temperature <- tavg[sel]
rm(temp, t.m, fill, i, t.v, tavg, time.stamp, sel)


##========================================================
## Create the grid cells 
##========================================================

## Set the range of the spatial domain
y.range <- range(st_coordinates(states)[,2])
x.range <- range(st_coordinates(states)[,1])
## Set the number of grid cells on the horizontal axis
## y.len is calculated such that each grid cell 
## is roughly a square shape
x.len <-  60 
y.len <- round(x.len/(diff(x.range)/diff(y.range)))
width <-  diff(x.range)/x.len
height <- diff(y.range)/y.len

## Determine the centroids based on the domain and resolution
lon.centers <- seq(x.range[1]+width/2, x.range[2]-width/2, width)
lat.centers <- seq(y.range[1]+height/2, y.range[2]-height/2, height)
centroids <- cbind(rep(lon.centers,length(lat.centers)), sort(rep(lat.centers, length(lon.centers))))

## Set the top, right, left, bottom of the grid cells
yPlus <- centroids[,2] + height/2
xPlus <- centroids[,1] + width/2
yMinus <- centroids[,2] - height/2
xMinus <- centroids[,1] - width/2

## Determine the four corners of the grid cells
square=cbind(xMinus, yPlus,  # NW corner
             xPlus, yPlus,  # NE corner
             xPlus, yMinus,  # SE corner
             xMinus, yMinus, # SW corner
             xMinus, yPlus)  # NW corner again - close ploygon

## Create polygon shapes for the grid cells
cells <- NULL
for (i in 1:nrow(centroids)) {
  cells[[i]] <- st_polygon(list(matrix(square[i,],ncol = 2, byrow = T))) 
	 
}

## Transform accordingly
cells <- st_as_sfc(cells)


##========================================================
## Some plots
##========================================================

plot(cells, border = 'gray60', asp = 1)
plot(states, add = T)
plot(stations$geometry, add = T, pch = 19, col = 'red', cex = 0.2)

plot(cells, border = 'gray60', asp = 1)
plot(states, add = T)
plot(stations.all$geometry, add = T, pch = 19, col = 'red', cex = 0.2)

length(unlist(st_intersects(AZ, cells)))
length(unlist(st_intersects(CA, cells)))
length(unlist(st_intersects(NV, cells)))
length(unlist(st_intersects(UT, cells)))


##=======================
## map points to cells
##=======================
dat <- NULL

## Create a data frame to store our modeling data at the grid cell level
dat$lon <- centroids[,1]
dat$lat <- centroids[,2]
dat <- data.frame(dat)
dat$geometry <- cells

dat$area <- 0
dat$area_coords <- 0
sel <- which(sapply(st_intersects(cells, states), function (x) x[1]) == 1)
## Calculate the area in square kms
dat$area[sel] <- round(st_area(st_intersection(cells, states)),4)
## Calculate the area in degrees
dat$area_coords[sel] <- dat$area[sel]/(st_area(cells[sel]))*width^2
dat$npts <- rep(0, nrow(dat))
dat$npts.all <- rep(0, nrow(dat))
dat$cell_id <- c(1:nrow(dat))
which.cell <- unlist(st_intersects(stations, cells))
repeat.station <- which(sapply(st_intersects(stations, cells), function(x) length(x)) > 1)

repeat.station.all <- which(sapply(st_intersects(stations.all, cells), function(x) length(x)) > 1)
ind.all <- lapply(st_intersects(stations.all, cells), function(x) x)
ind <- lapply(st_intersects(stations, cells), function(x) x)
for(i in repeat.station){
  ind[[i]] <- min(ind[[i]])
}
for(i in repeat.station.all){
  ind.all[[i]] <- min(ind.all[[i]])
}

which.cell.all <- unlist(ind.all)
which.cell <- unlist(ind)
rm(ind, ind.all, repeat.station, repeat.station.all)

cts <- tapply(which.cell, which.cell, length)
cts.all <- tapply(which.cell.all, which.cell.all, length)
# inds <- match(names(cts), 1:nrow(centroids))
dat$npts[as.numeric(names(cts))] <- cts
dat$npts.all[as.numeric(names(cts.all))] <- cts.all


which.CA <- unlist(st_intersects(CA, stations))
which.AZ <- unlist(st_intersects(AZ, stations))
which.UT <- unlist(st_intersects(UT, stations))
which.NV <- unlist(st_intersects(NV, stations))
which.CA.all <- unlist(st_intersects(CA, stations.all))
which.AZ.all <- unlist(st_intersects(AZ, stations.all))
which.UT.all <- unlist(st_intersects(UT, stations.all))
which.NV.all <- unlist(st_intersects(NV, stations.all))

length(c(which.CA, which.AZ, which.UT, which.NV))
length(c(which.CA.all, which.AZ.all, which.UT.all, which.NV.all))

##=================================
## Store temperature as dataframe
##=================================

stations$tavg <- temperature
stations$cell <- which.cell
stations$complete <- n.in
stations$state <- rep(NA, nrow(stations))
stations$state[which.CA] <- 'CA'
stations$state[which.AZ] <- 'AZ'
stations$state[which.UT] <- 'UT'
stations$state[which.NV] <- 'NV'
stations$lon <- stations$X1
stations$lat <- stations$X2

dat$Z <- rep(NA, nrow(dat))
for(i in sort(unique(which.cell))){
  dat$Z[i] <- mean(stations$tavg[which(which.cell == i)])
}
rm(temperature, n.in, which.cell, which.cell.all, i)

##=======================
## Create point processes
##=======================

sf_use_s2(FALSE)

convert_to_owin <- function(myshape){
  sf_use_s2(FALSE)
  shape.bdry <- st_coordinates(myshape)[,c(1,2)]
  shape.bdry <- shape.bdry[nrow(shape.bdry):1,]
  owin(poly = list(x = shape.bdry[,1],
                   y = shape.bdry[,2]))
}

states_win <- convert_to_owin(states)
AZ_win <- convert_to_owin(AZ)
CA_win <- convert_to_owin(CA)
NV_win <- convert_to_owin(NV)
UT_win <- convert_to_owin(UT)


## Used to generate
pp.states <- ppp(coords[,1], coords[,2], window = states_win)
pp.CA <- ppp(coords[which.CA,1], coords[which.CA,2], window = CA_win)
pp.AZ <- ppp(coords[which.AZ,1], coords[which.AZ,2], window = AZ_win)
pp.NV <- ppp(coords[which.NV,1], coords[which.NV,2], window = NV_win)
pp.UT <- ppp(coords[which.UT,1], coords[which.UT,2], window = UT_win)


pp.states.all <- ppp(coords.all[,1], coords.all[,2], window = states_win)
pp.CA.all <- ppp(coords.all[which.CA.all,1], coords.all[which.CA.all,2], window = CA_win)
pp.AZ.all <- ppp(coords.all[which.AZ.all,1], coords.all[which.AZ.all,2], window = AZ_win)
pp.NV.all <- ppp(coords.all[which.NV.all,1], coords.all[which.NV.all,2], window = NV_win)
pp.UT.all <- ppp(coords.all[which.UT.all,1], coords.all[which.UT.all,2], window = UT_win)

CA_cell_id <- unlist(st_intersects(CA, cells))
AZ_cell_id <- unlist(st_intersects(AZ, cells))
NV_cell_id <- unlist(st_intersects(NV, cells))
UT_cell_id <- unlist(st_intersects(UT, cells))

dat.CA <- dat[unlist(st_intersects(CA, cells)),]
dat.AZ <- dat[unlist(st_intersects(AZ, cells)),]
dat.NV <- dat[unlist(st_intersects(NV, cells)),]
dat.UT <- dat[unlist(st_intersects(UT, cells)),]
stations$Population_Density <- dat$Population_Density[stations$cell]


parents <- rbind(c(-111.64, 33.35),
                 c(-118.3, 33.68),
                 c(-122.42, 37.68),
                 c(-119.94, 38.97),
                 c(-111.68, 41.47),
                 c(-111.68, 40.88))
      
colnames(stations) <- c('lon', 'lat', 'geometry', 'tavg', 'cell_id', 'complete', 'state', 'pop_den')
colnames(dat) <- c('lon', 'lat', 'geometry', 'area', 'area_coords_DONOTUSE', 'npts', 'npts.all', 'cell_id', 'pop_den', 'Z')

rm(CA.sep, cts.all, convert_to_owin)

save.image(paste0("../Data/prepped_data_", x.len ,".RData"))
