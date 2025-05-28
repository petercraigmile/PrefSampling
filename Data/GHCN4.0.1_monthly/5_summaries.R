
##load(file="Derived/cont_US_stations.RData")
load(file="Derived/cont_US_monthly_tavg.RData")


pct.missing <- function (x) {

    mean(is.na(x$tavg))*100
}


first <- sapply(cont.US.monthly.tavg, function (x) floor(min(x$year)))
last  <- sapply(cont.US.monthly.tavg, function (x) floor(max(x$year)))


df <- data.frame(id          = sapply(cont.US.monthly.tavg, function (x) x$id),
                 pct.missing = round(sapply(cont.US.monthly.tavg, pct.missing), 2),
                 first.year  = first,
                 last.year   = last,
                 num.years   = last-first+1)[order(last),]

write.csv(df, file='Derived/summaries.csv')
