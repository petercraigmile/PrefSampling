

guide <- function (x) {

    abline(h=x, col="blue", lwd=2)
}

par(mfrow=c(4,2), cex=0.75, mar=c(3.1,3.1,1,0.5), mgp=c(1.8,0.5,0), bty="L")

ts.plot(chs["sigma2"])
guide(sigma2.Z)

ts.plot(chs["beta"])
guide(mu.W)

ts.plot(chs["alpha"])
guide(alpha)

ts.plot(chs["theta.W"][1,])
guide(theta.W[1])

ts.plot(chs["theta.W"][2,])
guide(theta.W[2])
mtext(paste(round(mean(chs["theta.W.jumps"]), 3), "acceptance prop"), side=3)

ts.plot(chs["delta"])
guide(delta)

ts.plot(chs["omega2"])
guide(omega2)
mtext(paste(round(mean(chs["delta.omega2.jumps"]), 3), "acceptance prop"), side=3)

if (length(chs["kappa"])>0) {
    ts.plot(chs["kappa"])
    guide(kappa)
    mtext(paste(round(mean(chs["kappa.jumps"]), 3), "acceptance prop"), side=3)
}
#mtext(paste(round(mean(chs["kappa.jumps"]), 3), "acceptance prop"), side=3)

par(mfrow=c(2,2), cex=0.75, mar=c(3.1,3.1,1,0.5), mgp=c(1.8,0.5,0), bty="L")

dd <- chs["delta"]
oo <- chs["omega2"]

lls <- rowMeans(sapply(1:length(dd), function (k)
    Thomas.intensity(G$grid, parents, dd[k], oo[k])))

RR <- range(c(lambda, lls))


GSP.image(lambda, zlim=RR)
mtext(side=3, "Truth")
points(pp$locs)
par(cex=0.75)

GSP.image(lls, zlim=RR)
mtext(side=3, "Post mean")
points(pp$locs)
par(cex=0.75)

#c.eta <- chs["eta"]

#pm <- rowMeans(c.eta)
#RR <- range(c(pp$eta, pm))

#GSP.image(pp$eta, zlim=RR)
#points(pp$locs, cex=0.5)

#GSP.image(pm, zlim=RR)
#points(pp$locs, cex=0.5)

matplot(t(chs["y"]), type="l", lty=1)

#plot(rowMeans(chs["eta.jumps"]))

# ml.est$beta.hat
#[1] 13.38649
# ml.est$theta.W.hat
#[1]  6.327569  3.139750 15.578875

bsumm <- function (x) {

    round(c(mean(x), quantile(x, c(0.025, 0.975))), 3)
}


print(bsumm(chs["beta"]))
print(bsumm(chs["theta.W"][1,]))
print(bsumm(chs["theta.W"][2,]))
print(bsumm(chs["sigma2"]))
