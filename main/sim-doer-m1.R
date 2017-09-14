## Synopsis: Get M0,M1,M2 goodness of fit test and, generating from M1
source("../main/sim-driver-multiplesteps.R")
la("../../binSegInf/binSegInf")

levs = c(0,1,3,5)
mns = lapply(levs, function(lev)c(rep(0,3), rep(lev,3)))
nsim = 300
sim.settings.list = list(list(mn = mns[[1]], sigma=1, nsim=nsim, type="binseg", mc.cores=6, maxsteps=1, ngen = 10000),
                         list(mn = mns[[2]], sigma=1, nsim=nsim, type="binseg", mc.cores=6, maxsteps=1, ngen = 10000),
                         list(mn = mns[[3]], sigma=1, nsim=nsim, type="binseg", mc.cores=6, maxsteps=1, ngen = 10000),
                         list(mn = mns[[4]], sigma=1, nsim=nsim, type="binseg", mc.cores=6, maxsteps=1, ngen = 10000))

resultlist = lapply(sim.settings.list,function(mysetting){ dosim2(mysetting)})
save(list=c("pvslist", "sim.settings.list"), file="../output/m1-sims.Rdata")

## ## Try one setting out
## sim.settings = list(mn = mns[[1]], sigma=1, nsim=nsim, type="binseg", mc.cores=6, maxsteps=1, ngen = 10000)
## result = dosim2(sim.settings)



## ## Analyze
## load(file="../output/results.Rdata")
## par(mfrow=c(3,1)){
##     qqunif(sapply(pvs, function(p)p$pvseq[1]))
##     qqunif(sapply(pvs, function(p)p$pvseq[2]))
##     qqunif(sapply(pvs, function(p)p$pvseq[3]))
## }

## sapply(pvs, function(p)p$changepoint.set.seq)
## sapply(all.results,function(myresult)myresult[["pvseq"]])
## sapply(pvs,function(myresult)myresult[["fdp"]])

## save(list=c("pvslist", "sim.settings.list"), file="../output/null-pvals.Rdata")
