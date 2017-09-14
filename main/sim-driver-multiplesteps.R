## Synopsis: Try k-step binary segementation selected model tests.
## sim.settings = list(mn=c(rep(0,3), rep(5,3)), sigma=1, nsim=500, type="binseg", teststep = 1)
library(binSegInf)
library(genlassoinf)
source("../main/sim-helper.R")

#' Gets goodness of fit test p-values for models M_0 up to M_|maxsteps|.
dosim2 <- function(sim.settings){

    ## Extracting things
    sigma = sim.settings$sigma
    mn = sim.settings$mn
    n = length(mn)
    nsim = sim.settings$nsim
    type = sim.settings$type
    covariance <- diag(rep(sigma,n))
    mc.cores = sim.settings$mc.cores
    maxsteps = sim.settings$maxsteps

    ## Collect results
    all.results = lapply(1:nsim, function(isim){
        printprogress(isim, nsim)

        ## Generate data
        y0 <- mn + rnorm(n,0,sigma)

        pvseq = rep(NA,maxsteps+1)
        names(pvseq) = sapply(0:maxsteps,toString)
        changepoint.set.seq = list()

        models = get.all.models(y0, type=type, maxsteps+1)

        ngens = c(500, 500, 500000)
        for(istep in 0:maxsteps){

            ## Collect i'th model
            cp.curr <- models$cp.list[[istep+1]]
            cp.sign.curr <- models$cp.sign.list[[istep+1]]
            poly.curr = models$poly.list[[istep+1]]
            cp.next <- models$cp.list[[istep+2]]
            cp.sign.next = models$cp.sign.list[[istep+2]]


            ## Get conditional Gauss parameters
            myparam <- get.cond.gauss.param(cp.curr, y0, covariance=covariance)
            muorig <- myparam$muorig
            Sigmaorig <- myparam$Sigmaorig


            ## Generate new ys, then rejection sample
            ngen = ngens[istep+1]
            ys = (MASS::mvrnorm(n=ngen,mu=muorig,Sigma=Sigmaorig))
            in.polyhedron <- apply(ys, 1, function(myrow){
                return(all(poly.curr$gamma %*% myrow >= poly.curr$u))
            })
            accepted = which(in.polyhedron)
            if(length(accepted) > 1000){
                accepted = accepted[1:1000]
            }
            ys = ys[accepted,]

            ## Compute the quantile of the vTY|AY
            ## vtvec = apply(ys, 1, function(my.y){
            vtvec = mclapply(1:nrow(ys), function(irow){
                my.y = ys[irow,]
                newmodels = get.all.models(my.y, type=type, maxsteps=istep+1, TRUE)
                new.cp <- newmodels$cp.list[[istep+2]]
                new.cp.sign <- newmodels$cp.sign.list[[istep+2]]
                all.vs = make_all_segment_contrasts_from_cp(new.cp, new.cp.sign, n)
                my.v <- all.vs[[toString(tail(new.cp,1)*tail(new.cp.sign,1))]]
                return(sum(my.y*my.v))
            },mc.cores=mc.cores)

            ## Compute original v^TY
            all.vs = make_all_segment_contrasts_from_cp(cp.next,cp.sign.next,n)
            v <- all.vs[[toString(tail(cp.next,1)*tail(cp.sign.next,1))]]
            observed.vt = as.numeric(v%*%y0)

            ## pv = sum(vtvec > observed.vt | vtvec < -observed.vt)/length(vtvec)
            pv = sum(vtvec > observed.vt)/length(vtvec)
            pvseq[toString(istep)] = pv
            changepoint.set.seq[[toString(istep)]] = cp.curr
        }

        ## Compare observed verdicts to truth
        true.cp = n/2
        true.verdicts= sapply(changepoint.set.seq, function(mycpset){all(true.cp %in% mycpset)})
        true.stop.time = which.min(true.verdicts)

        ## See verdicts
        fdrtarget=0.1
        pvseq = pvseq[which(!is.na(pvseq))]
        stop.time = selectiveInference::forwardStop(pvseq, alpha=fdrtarget)


        ## False discovery proportion
        fdp = (stop.time-true.stop.time)/true.stop.time

        return(list(stop.time= stop.time, true.stop.time = true.stop.time, fdp = fdp,
                    changepoint.set.seq=changepoint.set.seq,
                    pvseq=pvseq))

    })

    return(all.results)
}


## Analyze results
## (sapply(all.results,function(myresult)myresult[["pvseq"]]))
## (sapply(all.results,function(myresult)myresult[["fdp"]]))

