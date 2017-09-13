## Synopsis: Run simulations for certain settings, save them in appropriate files
source("../main/sim-driver.R")

## Things to run
sim.settings.list = list(list(mn=rep(0,10), sigma=1, nsim=500, type="binseg", teststep = 1),
                         list(mn = rep(0,10), sigma=1, nsim=500, type="fusedlasso", teststep = 1),
                         list(mn=rep(0,10), sigma=1, nsim=500, type="binseg", teststep = 0),
                         list(mn=rep(0,10), sigma=1, nsim=500, type="fusedlasso", teststep = 0))

pvslist = mclapply(sim.settings.list,function(mysetting)dosim(mysetting), mc.cores=4)
save(list=c("pvslist", "sim.settings.list"), file="../output/null-pvals.Rdata")

## Optionally plot them
load("../output/null-pvals.Rdata")

pdf("../output/null-pvals.pdf",width=5,height=5)
par(mar=c(2,2,2,2))
par(mfrow=c(2,2))
titles = c("M1 Binary Segmentation",
           "M1 Fused Lasso",
           "M0 Binary Segmentation",
           "M0 Fused Lasso")

for(iplot in 1:length(pvslist)){
    qqunif(pvslist[[iplot]], xlab="Observed",ylab="Expected", main=titles[iplot])
}
graphics.off()
