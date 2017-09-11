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

