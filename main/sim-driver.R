# Synopsis: Try one-step binary segementation selected model tests.
library(binSegInf)
library(genlassoinf)
source("../main/sim-helper.R")

sim.settings = list(mn=rep(0,10), sigma=1, nsim=500, type="binseg", teststep = 1)
## sim.settings = list(mn=rep(0,n), sigma=1, nsim=500, type="binseg", teststep = 0)
## sim.settings = list(mn=rep(0,4), sigma=1, nsim=500, type="fusedlasso", teststep = 0)

dosim <- function(sim.settings){

    ## Extracting things
    sigma = sim.settings$sigma
    mn = sim.settings$mn
    n=length(mn)
    nsim = sim.settings$nsim
    type = sim.settings$type
    teststep = sim.settings$teststep
    covariance <- diag(rep(sigma,n))

    pvs =  mclapply(1:nsim, function(isim){
    ## for(isim in 1:nsim){
        printprogress(isim,nsim)

        ## Generate data
        y0 <- mn + rnorm(n,0,sigma)

        ## Do things
        orig <- getstuff(y0, type=type, teststep)


        ## Form the sufficient statistic linear operator
        plateaus = get_plateaus(orig$cp1, n)
        suff.rows = do.call(rbind,lapply(plateaus, function(plt){
                                      v= rep(0,n); v[plt] = 1/length(plt); return(v)}))
        A = suff.rows

        ## Get null span
        w = A %*% y0
        S = svd(t(A),nu=ncol(A))
        nr = nrow(A)
        A_rest = t(S$u[,(nr+1):ncol(A)])

        Ag = rbind(A,A_rest)

        ## Partition matrix into four blocks
        Sigma = (Ag) %*% covariance %*% t(Ag)
        Si = nrow(A)
        S11 = Sigma[1:Si, 1:Si, drop=FALSE]
        S12 = Sigma[1:Si, (Si+1):n, drop=FALSE]
        S21 = Sigma[(Si+1):n, 1:Si, drop=FALSE]
        S22 = Sigma[(Si+1):n, (Si+1):n, drop=FALSE]

        ## Calculate distr of (A_rest y| Ay) using Schur's complement
        Sigmarest = S22 - S21%*%solve(S11, S12)
        murest = A_rest%*%mn + S21 %*% solve(S11, cbind(w - A %*% mn))

        ## Linear this back to to the original y scale, via left multiplication by Ag
        Aginv = solve(Ag)
        Sigmanew = matrix(0, nrow=n, ncol=n)
        Sigmanew[(Si+1):n, (Si+1):n] = Sigmarest
        Sigmaorig = Aginv %*% Sigmanew %*% t(Aginv)
        muorig = Aginv %*%  cbind(c(w, murest))


        ngen = 10000
        ys = (MASS::mvrnorm(n=ngen,mu=muorig,Sigma=Sigmaorig))

        ## Rejection sample
        in.polyhedron <- apply(ys, 1, function(myrow){
            all(orig$mypoly$gamma %*% myrow > orig$mypoly$u)
           return(all(orig$mypoly$gamma %*% myrow >= orig$mypoly$u))
        })

        for(irow in 1:nrow(ys)){
            myrow = ys[irow,]
            all(orig$mypoly$gamma %*% myrow > orig$mypoly$u)
        }

        ## Compute the quantile of the vTY|AY
        ys= ys[which(in.polyhedron),]

        vtvec = apply(ys, 1, function(my.y){
            nu = getstuff(my.y, type=type, teststep=teststep)
            my.v <- make_all_segment_contrasts(nu$g2)[[toString(nu$cp2*nu$cp2.sign)]] ## segment contrast
            return(sum(my.y*my.v))
        })

        ## Compute original v^TY
        v <- make_all_segment_contrasts(orig$g2)[[toString(orig$cp2*orig$cp2.sign)]]
        observed.vt = as.numeric(v%*%y0)

        ## pv = sum(vtvec > observed.vt | vtvec < -observed.vt)/length(vtvec)
        pv = sum(vtvec > observed.vt)/length(vtvec)

        return(pv)
    }, mc.cores=2)
    ## }
    return(unlist(pvs))
}



## Synopsis: Try k-step binary segementation selected model tests.

sim.settings = list(mn=rep(0,10), sigma=1, nsim=500, type="binseg", teststep = 1)

##' Gets first and second p-valuex.
dosim2 <- function(sim.settings){

    ## Extracting things
    sigma = sim.settings$sigma
    mn = sim.settings$mn
    n = length(mn)
    nsim = sim.settings$nsim
    type = sim.settings$type
    teststep = sim.settings$teststep
    covariance <- diag(rep(sigma,n))
    lev=sim.settings$lev

    ## Collect results
    all.results = mclapply(1:nsim, function(isim){

        printprogress(isim,nsim)

        ## Generate data
        y0 <- mn + rnorm(n,0,sigma)

        maxsteps = 2
        pvseq = rep(NA,maxsteps+1)
        names(pvseq) = sapply(0:maxsteps,toString)
        changepoint.set.seq = list()

        get.all.models(y0, type=type, maxsteps)
        get.all.models <- function(y, type, maxsteps){

            n = length(y)

            maxsteps = 5
            obj.final = binSeg_fixedSteps(y0, numSteps=maxsteps) ## Change to y
            poly.list = list()
            model.list = list()
            if(type=="binseg"){
                for(istep in 0:maxsteps){
                    if(istep==0){
                        poly.list[[istep+1]] = polyhedra(obj=rbind(rep(0,n)),u=0)
                        cp.list[[istep+1]] = cp.sign.list[[istep+1]] = NA
                    }
                    poly.list[[istep+1]] = polyhedra(obj.final, numsteps=istep)
                    cp.list[[istep+1]] = obj.final$cp[1:istep]
                    cp.sign.list[[istep+1]] = obj.final$cp.sign[1:istep]
                }
            }

            obj.final$cp
            things

            ## it would be great if I could collect smaller submodels from binSeg_fixedSteps.


             ## Fit binseg model, collect pvals BSFS.
             if(type=="binseg"){

                 if(teststep==0){
                     obj.curr=NULL
                     poly.curr = polyhedra(obj=rbind(rep(0,n)),
                                           u=0)
                     cp.curr = cp.curr.sign = c()
                 } else {
                     obj.curr = binSeg_fixedSteps(y, numSteps=teststep)
                     poly.curr = polyhedra(obj.curr)
                     cp.curr = obj.curr$cp
                     cp.curr.sign = obj.curr$cp.sign
                 }
                 obj.next = binSeg_fixedSteps(y, numSteps=teststep+1)
             } else {
                 stop("type not written yet!")
             }

             cp.next = obj.next$cp[which(!(obj.next$cp%in%cp.curr))]
             cp.next.sign = obj.next$cp.sign[which(!(obj.next$cp%in%cp.curr))]

             return(list(obj.curr=obj.curr,
                         obj.next=obj.next,
                         cp.curr=cp.curr,
                         cp.next=cp.next,
                         cp.curr.sign=cp.curr.sign,
                         cp.next.sign=cp.next.sign,
                       poly.curr=poly.curr))
        }

        ngens = c(1000,10000,500000)
        for(istep in 0:maxsteps){

            ## tryCatch({
            ## Collect i'th model
            source("../main/sim-helper.R")
            modelinfo <- get.two.models(y0, type=type, istep)
            myparam <- get.cond.gauss.param(modelinfo$obj.curr, y0)
            muorig <- myparam$muorig
            Sigmaorig <- myparam$Sigmaorig

            ## Generate new ys, then rejection sample
            ngen = ngens[istep+1]
            ys = (MASS::mvrnorm(n=ngen,mu=muorig,Sigma=Sigmaorig))
            in.polyhedron <- apply(ys, 1, function(myrow){
                return(all(modelinfo$poly.curr$gamma %*% myrow >= modelinfo$poly.curr$u))
            })
            ys= ys[which(in.polyhedron),]

            ## Compute the quantile of the vTY|AY
            vtvec = apply(ys, 1, function(my.y){
                new.modelinfo = get.ith.model(my.y, type, istep)
                modelinfo <- get.two.models(my.y, type, istep)
                all.vs = make_all_segment_contrasts(new.modelinfo$obj.next)
                my.v <- all.vs[[toString(new.modelinfo$cp.next*new.modelinfo$cp.next.sign)]]
                return(sum(my.y*my.v))
            })

            ## Compute original v^TY
            all.vs = make_all_segment_contrasts(modelinfo$obj.next)
            v <- all.vs[[toString(modelinfo$cp.next*modelinfo$cp.next.sign)]]
            observed.vt = as.numeric(v%*%y0)

            ## pv = sum(vtvec > observed.vt | vtvec < -observed.vt)/length(vtvec)
            pv = sum(vtvec > observed.vt)/length(vtvec)
            pvseq[toString(istep)] = pv
            changepoint.set.seq[[toString(istep)]] = modelinfo$cp.curr
            ## }, error = function(err) {
            ##     err$message = paste(err$message,"\n(rejection sampling rejected all samples..!)",sep="")
            ##     warning(err)})
        }

        ## ## Compare observed verdicts to truth
        true.cp = n/2
        true.verdicts= sapply(changepoint.set.seq, function(mycpset){all(true.cp %in% mycpset)})
        true.stop.time = which.min(true.verdicts)

        ## ## See verdicts
        verdicts = (pvseq<0.05)
        pvseq = pvseq[which(!is.na(pvseq))]
        stop.time = selectiveInference::forwardStop(pvseq, alpha=fdrtarget)


        ## False discovery proportion
        fdp = (stop.time-true.stop.time)/true.stop.time

        return(list(stop.time=stop.time, true.stop.time = true.stop.time, fdp = fdp,
                    changepoint.set.seq=changepoint.set.seq,
                    pvseq=pvseq))
    },mc.cores=4)




    ## Analyze results
    ## (sapply(all.results,function(myresult)myresult[["pvseq"]]))
    ## (sapply(all.results,function(myresult)myresult[["fdp"]]))

    return()


}
