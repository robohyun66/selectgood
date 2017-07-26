## Synopsis: Generate data from a one-jump mean, test the first changepoint
## model's goodness of fit.

## Generate data
n = 12
sigma = 1
delta = 0 ## signal size
mn = c(rep(0,n/2), rep(delta,n/2))

onesim <- function(ii=NULL, ngen=300){
    ## Generate data
    if(!is.null(ii))  set.seed(ii)
    y0 <- mn + rnorm(n, 0, sigma)

    ## Run two-step fused lasso on it
    maxsteps = 2
    D = genlassoinf::makeDmat(n, type='tf',ord=0)
    f0 = genlassoinf::dualpathSvd2(y0, D=D, maxsteps, approx=T)
    cp1 = f0$action[1]
    cp2 = f0$action[2]

    ## Get naive 1-step poyhedron
    Gobj.naive = genlassoinf::getGammat.naive(obj = f0, y = y0,
                                              condition.step = 1)
    G = Gobj.naive$G
    u = Gobj.naive$u
    cps = cp1
    inds = do.call(c,lapply(cps, function(cp)return(cp+c(-1,0,1))))
    inds = inds[inds!=0]
    other.inds = (1:n)[-inds]

    ## Get A = the map to the sufficient statistics
    plateaus = get_other_segment_indices(other.inds)
    suff.rows = do.call(rbind,lapply(plateaus, function(plt){v=rep(0,n); v[plt] = 1/length(plt); return(v)}))
    sat.rows = do.call(rbind, lapply(inds, function(ind){v=rep(0,n); v[ind] = 1;return(v)}))
    A = rbind(suff.rows, sat.rows)


    ## Get conditional distribution of y|Ay
    dst = get_condit_distr_y(A, mn, diag(rep(sigma^2,n)), y0)

    ## Generate new y's and accept-reject
    ys = t(MASS::mvrnorm(n=ngen,mu=dst$mu,Sigma=dst$cov))
    which.y.in.polyhedron =  apply(G%*%ys, 2, function(mycol){
        all(as.numeric(mycol) > u)
    })
    ys = ys[,which.y.in.polyhedron, drop=FALSE]

    ## Option 1: Make next segment test statistic (fixed, upward direction contrast)
    d = +1## (f0$pathobj$s)[2]

    if(cp2 > cp1){
        test.stat.row = make.contrast(cp1,cp2,n,n,d)
    } else {
        test.stat.row = make.contrast(1,cp2,cp1,n,d)
    }

    ## ## Option 2: Make /fixed/ test statistic
    ## test.stat.row = rep(0,n)
    ## test.stat.row[c(1,9)] = 1
    ## test.stat.row = rbind(test.stat.row)

    ## Calculate the p-value
    mystat = as.numeric(rbind(test.stat.row)%*%cbind(y0))
    newstat = apply(ys, 2, function(mycol){
        return(rbind(test.stat.row)%*%cbind(mycol))
    })
    pv = sum(newstat>abs(mystat) | newstat < -abs(mystat) )/length(newstat) ## two sided
    ## pv = sum(newstat>mystat)/length(newstat) ## one sided

    ## Also compute naive z-test pv
    munaive = rbind(test.stat.row)%*%cbind(mn)
    sigmanaivesquare = (rbind(test.stat.row)%*%t(test.stat.row)) * (sigma^2)
    pvnaive = 1-pnorm(abs(mystat), mean=munaive, sd=sqrt(sigmanaivesquare)) +
        pnorm(-abs(mystat), mean=munaive, sd=sqrt(sigmanaivesquare)) ## two sided
    ## pvnaive = 1-pnorm(mystat, mean=munaive, sd=sqrt(sigmanaivesquare)) ## one sided

    return(list(pv=pv,pvnaive=pvnaive, cps = c(cp1,cp2)))
}

## Main
nsim = 300
pvs = mclapply((1:nsim), function(isim){
        cat('\r', isim, 'out of', nsim)
        tryCatch({
            return(onesim(ngen=300))
        }, error=function(e) {print("error occurred")})
}, mc.cores=3)

## Extract things
pvsadj = sapply(pvs,function(a)a[["pv"]])
pvsnaive = sapply(pvs,function(a)a[["pvnaive"]])

## Make plot
qqunif(pvsnaive)
qqunif(pvsadj,add=TRUE,col='red')
legend("topleft", col=c('black','red'),pch=c(16,16), legend = c("naive", "adjusted"))
