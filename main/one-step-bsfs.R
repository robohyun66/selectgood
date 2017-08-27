## Synopsis: Try one-step binary segementation selected model tests.
library(binSegInf)

## Generate some data
n=10; sigma=1; mn = rep(0,n)
covariance <- diag(rep(sigma,n))

nsim = 1000
pvs = rep(NA,nsim)
## for(isim in 1:nsim){
pvs = mclapply(1:nsim, function(isim){
    printprogress(isim,nsim)
    ## Generate data
    set.seed(isim)
    y <- mn + rnorm(n,0,sigma)

    ## Fit a 1-step BSFS.
    g1 = binSeg_fixedSteps(y, numSteps=1)
    g2 = binSeg_fixedSteps(y, numSteps=2)
    mypoly = polyhedra(g1)
    cp1 = g1$cp[1]
    cp1.sign = g1$cp.sign[1]
    cp2 = g2$cp[which(!(g2$cp==g1$cp))]
    cp2.sign = g2$cp.sign[which(!(g2$cp==g1$cp))]

    ## Form contrast
    ## v <- make_all_segment_contrasts(g1)[[1]] ## segment contrast
    v <- make_contrast(cp2, c(cp1,cp2), cp2.sign, n)

    ## Form a fixed contrast
    ## v = c(-2,+1,2.5,-1.3)

    ## Form the sufficient statistic linear operator
    plateaus = get_plateaus(cp1, n)
    ## plateaus = get_plateaus(c(),n) ## Oddly, this gives uniform p-values
    suff.rows = do.call(rbind,lapply(plateaus, function(plt){v=rep(0,n); v[plt] = 1/length(plt); return(v)}))
    A = suff.rows

    ## Get null span
    w = A %*% y
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
        return(all(mypoly$gamma %*% myrow > mypoly$u))
    })


    ## Compute the quantile of the vTY|AY
    vtvec = (ys[which(in.polyhedron),]%*%v)
    observed.vt = as.numeric(v%*%y)
    ## pv = sum(vtvec > observed.vt | vtvec < -observed.vt)/length(vtvec)
    pv = sum(vtvec > observed.vt)/length(vtvec)
    ## pvs[isim] <- pv
    return(pv)

}, mc.cores=3)
