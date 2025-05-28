
load("Derived/cont_US_stations.RData")

US.temps <- readLines("Derived/US_monthly_temps.txt")


extract.monthly.values <- function (row.string, raw=FALSE) {

    elem <- substr(row.string, 16, 19)

    year <- as.numeric(substr(row.string, 12, 15))

    values <- DMFLAGs <- QCFLAGs <- DSFLAGs <- rep(NA, 12)

    pos <- 20
    
    for (k in 1:12) {
        
        values[k]  <- substr(row.string, pos, pos+4)
        DMFLAGs[k] <- substr(row.string, pos+5, pos+5)
        QCFLAGs[k] <- substr(row.string, pos+6, pos+6)
        DSFLAGs[k] <- substr(row.string, pos+7, pos+7)
        
        pos <- pos + 8
    }

    if (!raw) {
        
        values <- ifelse(values=="-9999", NA, as.numeric(values)/100)

        ## Values that fail any quality assurance check
        values <- ifelse(QCFLAGs!=" ", NA, values)

        ## E = Data has been estimated from neighbors within the PHA
        values <- ifelse(DMFLAGs=="E", NA, values)

        ## Check no "Global Summary of the Day" (should not be!)
        ## If so print "*" to output.
        if (any(DSFLAGs=="S")) {
            cat("*")
        }

        list(year=year+(0:11)/12,
             elem=rep(elem, 12),
             values=values)
        
    } else {

        list(year=year,
             elem=elem,
             values=values,
             DMFLAGs=DMFLAGs,
             QCFLAGs=QCFLAGs,
             DSFLAGs=QCFLAGs)
    }
}


##extract.monthly.values(sub[1], TRUE)
##extract.monthly.values(sub[1], FALSE)

monthly.tavg <- function (the.id, lat, long, alt, name) {

    ids <- substr(US.temps, 1, 11)
    
    sel <- which(ids==the.id)
    
    sub <- US.temps[sel]

    df <- lapply(1:length(sub), function (j) extract.monthly.values(sub[j]))
    
    list(id   = the.id,
         lat  = lat,
         long = long,
         alt  = alt,
         name = name,
         year = as.numeric(sapply(df, function (x) x$year)),
         tavg = as.numeric(sapply(df, function (x) x$values)))
}


## k <- 1
 
## monthly.tavg(cont.US.stations$id[k],
##              cont.US.stations$lat[k],
##              cont.US.stations$long[k],
##              cont.US.stations$alt[k],
##              cont.US.stations$name[k])

## rm(k)


library(parallel)

cont.US.monthly.tavg <- mclapply(1:nrow(cont.US.stations), function (k) {

    if (k%%1000==0) cat(round(k/nrow(cont.US.stations)*100), "% ")
    
    monthly.tavg(cont.US.stations$id[k],
                 cont.US.stations$lat[k],
                 cont.US.stations$long[k],
                 cont.US.stations$alt[k],
                 cont.US.stations$name[k])

}, mc.cores=6)


save(cont.US.monthly.tavg, file="Derived/cont_US_monthly_tavg.RData")

