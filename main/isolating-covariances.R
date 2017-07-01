## Generate data
library(genlassoinf)
n = 12
sigma = 1
mn = rep(1,n)
ngen=10000
check.poly=TRUE
nsim = 50

results = mclapply(1:nsim, function(isim){
    print(isim)
    set.seed(isim)
    y0 <- mn + rnorm(n, 0, sigma)

    ## Run two-step fused lasso on it
    maxsteps = 2
    D = makeDmat(n, type='tf',ord=0)
    f0 = dualpathSvd2(y0, D=D, maxsteps, approx=T)
    cp1 = f0$action[1]
    cp2 = f0$action[2]

    ## Get naive 1-step poyhedron
    Gobj.naive = getGammat.naive(obj = f0, y = y0, condition.step = 1)
    G = Gobj.naive$G
    u = Gobj.naive$u
    cps = cp1 ##f0$action[1]   ## cps = c(4,9)
    inds = do.call(c,lapply(cps, function(cp)return(cp+c(-1,0,1))))
    inds = inds[inds!=0]
    other.inds = (1:n)[-inds]

    ## Get A = the map to the sufficient statistics
    plateaus = get_other_segment_indices(other.inds)
    suff.rows = do.call(rbind,lapply(plateaus, function(plt){v=rep(0,n); v[plt] = 1/length(plt); return(v)}))
    sat.rows = do.call(rbind, lapply(inds, function(ind){v=rep(0,n); v[ind] = 1;return(v)}))
    ## A = rbind(suff.rows, sat.rows)
    A = suff.rows

    myrow = rep(1,n)
    test.stat.row = rbind(myrow) ## Temporarily added, nonrandom statistic


    #### Build A_rest
opts = lapply(2, function(opt){
    if(opt==1){
        ## Method 1: manual
        rest.inds = unlist(sapply(plateaus, function(plt)plt[-1]))
        A_rest= do.call(rbind, lapply(rest.inds, function(ind){v=rep(0,n); v[ind] = 1;return(v)}))

    } else if (opt==2){
        ## Method 2: svd (most generic!!!)
        S = svd(t(A),nu=ncol(A))
        nr = nrow(suff.rows)## + nrow(sat.rows) ##+ nrow(test.stat.row)
        A_rest = t(S$u[,(nr+1):ncol(A)])

    } else if(opt==3){

        ## Method 3: smarter manual (null space row bases)
        B = matrix(0, nrow = n-nrow(A), ncol = n)
        plateaus2 = plateaus[which(sapply(plateaus, length)!=1)]
        Bblocks <- lapply(plateaus2, function(plt){
            b = matrix(0, nrow = length(plt)-1, ncol = n)
            if(nrow(b)!=0){
                b[,plt] = makeDmat(m=length(plt), type="tf", ord=0)
            }
            return(b)
        })
        binds = c(0,cumsum(sapply(Bblocks, nrow)))
        Bblocks.rownums = lapply(1:(length(binds)-1), function(ii){
            as.numeric((binds[ii]+1):binds[ii+1])
        })
        for(ii in 1:length(plateaus2)){B[Bblocks.rownums[[ii]], ] = Bblocks[[ii]]}
        A_rest = B

    } else ( stop("not written!"))

    ## Build A_aug.
    Ag = rbind(A,A_rest)

    ## Make mean + covariance for W_rest | W
    w = A %*% y0

    ## Partition matrix into four blocks
    Sigma = (Ag) %*% diag(rep(sigma^2,n)) %*% t(Ag)
    Si = nrow(A)
    S11 = Sigma[1:Si, 1:Si, drop=FALSE]
    S12 = Sigma[1:Si, (Si+1):n, drop=FALSE]
    S21 = Sigma[(Si+1):n, 1:Si, drop=FALSE]
    S22 = Sigma[(Si+1):n, (Si+1):n, drop=FALSE]

    ## Calculate using Schur's complement
    Sigmarest = S22 - S21%*%solve(S11, S12)
    murest = S21 %*% solve(S11, cbind(w))

    ## Using the original way
    wrests = MASS::mvrnorm(n=ngen, mu=murest, Sigma=Sigmarest)
    ys = apply(wrests, 1, function(wrest){
        ysample = solve(Ag) %*% cbind(c(w,wrest))
        return(ysample)
    })
    if(check.poly){
        which.y.in.polyhedron =  apply(G%*%ys, 2, function(mycol){
            all(as.numeric(mycol) > u)
        })
        ys = ys[,which.y.in.polyhedron]
    }
    mystat = as.numeric(rbind(test.stat.row)%*%cbind(y0))

    ## Take ys, run the second step, calculate the linear contrast
    ## Calculate p-value
    newstat = apply(ys, 2, function(mycol){
        return(rbind(test.stat.row)%*%cbind(mycol))
    })
    pv = sum(newstat>mystat)/length(newstat)


    ## Transform to the conditional distribution
    Aginv = solve(Ag)
    Sigmanew = matrix(0, nrow=n, ncol=n)
    Sigmanew[(Si+1):n, (Si+1):n] = Sigmarest
    Sigmaorig = Aginv %*% Sigmanew %*% t(Aginv)
    muorig = solve(Ag, c(w, murest))

    ## Further transform to get distribution of test statistic
    sigmatest = rbind(test.stat.row) %*% Sigmaorig %*% t(test.stat.row)
    mutest = rbind(test.stat.row)%*%cbind(muorig)

    ys = t(MASS::mvrnorm(n=ngen,mu=muorig,Sigma=Sigmaorig))
    if(check.poly){
        which.y.in.polyhedron =  apply(G%*%ys, 2, function(mycol){
            all(as.numeric(mycol) > u)
        })
        ys = ys[,which.y.in.polyhedron, drop=FALSE]
    }
    mystat = as.numeric(rbind(test.stat.row)%*%cbind(y0))

    ## Take ys, run the second step, calculate the linear contrast
    ## Calculate p-value
    D = makeDmat(n,type='tf',ord=0)
    newstat2 = apply(ys, 2, function(mycol){
        return(rbind(test.stat.row)%*%cbind(mycol))
    })
    pv2 = sum(newstat2>mystat)/length(newstat2)

    ## is there /any/ meaningful variance?
    meaningful = (var(newstat2)>1E-10 | var(newstat)>1E-10)

    return(list(Sigmaorig=Sigmaorig, muorig=muorig, pv=pv, pv2=pv2, newstat,
                newstat2, meaningful=meaningful,sigmatest=sigmatest,mutest=mutest))
})

    ## lapply(opts, function(myopt)(myopt[c(3,4,7)]))
    ## return(any(sapply(opts, function(myopt)(myopt[c(7)]))))
    ## return(sapply(opts, function(myopt)(myopt[["sigmatest"]])))
    return(opts[[1]][c("sigmatest", "mutest")])
}, mc.cores=3)

### Investigate things
results
