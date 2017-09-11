## Synopsis: Run simulations for certain settings, save them in appropriate files
source("../main/sim-driver.R")

## Things to run
## pvs.1step.binseg = dosim(list(mn=rep(0,10), sigma=1, nsim=1000, type="binseg",
##                               teststep = 1))
## pvs.1step.fusedlasso = dosim(list(mn = rep(0,4), sigma=1, nsim=1000,
##                                   type="fusedlasso", teststep = 1))
## pvs.0step.binseg = dosim(list(mn=rep(0,4), sigma=1, nsim=10, type="binseg",
##                               teststep = 0))
## pvs.0step.fusedlasso = dosim(list( mn=rep(0,n), sigma=1, nsim=500,
##                                   type="fusedlasso", teststep = 0))

sim.settings.list = list(list(mn=rep(0,10), sigma=1, nsim=500, type="binseg", teststep = 1),
                         list(mn = rep(0,10), sigma=1, nsim=500, type="fusedlasso", teststep = 1),
                         list(mn=rep(0,10), sigma=1, nsim=500, type="binseg", teststep = 0),
                         list(mn=rep(0,10), sigma=1, nsim=500, type="fusedlasso", teststep = 0))

pvslist = mclapply(sim.settings.list,function(mysetting)dosim(mysetting), mc.cores=4)

