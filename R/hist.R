examples.deround.hist = function() {
  z.deround.hist(dat)
}


z.deround.hist = function(dat, mode = c("reported", "uniform","zda")[1], min.z=0, max.z=5, bin.width=0.05, z.pdf=NULL, repl=1,  ...) {
  restore.point("deround.hist")
  z = dat$z

  breaks = c(-Inf,seq(min.z,max.z, by=bin.width),Inf)
  hdat = hist(z, breaks=breaks, plot=FALSE)


}
