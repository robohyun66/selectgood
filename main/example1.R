## Generate data
library(genlassoinf)
## set.seed(1)
n = 12
sigma = 1
ngen = 300
mn = rep(0,n)
## mn = c(2,rep(0,n-1))


onesim <- function(ii=NULL, opt=3){

    ## Generate data
    if(!is.null(ii))  set.seed(ii)
    y0 <- mn + rnorm(n,0,sigma)

    ## Run two-step fused lasso on it
    maxsteps = 2
    D = makeDmat(n,type='tf',ord=0)
    f0 = dualpathSvd2(y0, D=D, maxsteps, approx=T)
    cp1 = f0$action[1]
    cp2 = f0$action[2]

    ## Get naive 1-step poyhedron
    Gobj.naive = getGammat.naive(obj = f0, y = y0, condition.step = 1)
    G = Gobj.naive$G
    u = Gobj.naive$u
    cps = cp1##f0$action[1]   ## cps = c(4,9)
    inds = do.call(c,lapply(cps, function(cp)return(cp+c(-1,0,1))))
    inds = inds[inds!=0]
    other.inds = (1:n)[-inds]

    ## Helper function
    get_other_segment_indices <- function(other.inds){
        brks = which(diff(other.inds)!=1)
        other.inds[brks]
        ilist = Map(function(a,b)(a+1):b, c(0,brks),c(brks,length(other.inds)))
        return( lapply(ilist, function(i)other.inds[i]))
    }
    make.contrast <- function(l, b, r, n, direction=+1){
        stopifnot(l <= b & b < r)
        ind1 = l:b
        ind2 = (b+1):r
        v = rep(0,n)
        v[ind1] = -1/length(ind1)
        v[ind2] = 1/length(ind2)
        return(rbind(v*direction))
    }


    ## Get A = the map to the sufficient statistics
    plateaus = get_other_segment_indices(other.inds)
    suff.rows = do.call(rbind,lapply(plateaus, function(plt){v=rep(0,n); v[plt] = 1/length(plt); return(v)}))
    sat.rows = do.call(rbind, lapply(inds, function(ind){v=rep(0,n); v[ind] = 1;return(v)}))
    A = rbind(suff.rows, sat.rows)

    #### Build A_rest
    if(opt==1){
        ## Method 1: manual
        rest.inds = unlist(sapply(plateaus, function(plt)plt[-1]))
        A_rest= do.call(rbind, lapply(rest.inds, function(ind){v=rep(0,n); v[ind] = 1;return(v)}))
    } else if (opt==2){
        ## Method 2: svd (most generic!!!)
        S = svd(t(A),nu=ncol(A))
        nr = nrow(suff.rows) + nrow(sat.rows) ##+ nrow(test.stat.row)
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

    if(opt==3){
        Sigmarest = B%*%t(B) * sigma^2
        murest = rep(0,nrow(B))
    } else {
        ## Partition matrix into four blocks
        Sigma = (Ag) %*% diag(rep(sigma^2,n)) %*% t(Ag)
        Si = nrow(A)
        S11 = Sigma[1:Si, 1:Si]
        S12 = Sigma[1:Si, (Si+1):n]
        S21 = Sigma[(Si+1):n, 1:Si]
        S22 = Sigma[(Si+1):n, (Si+1):n]

        ## Calculate using Schur's complement
        Sigmarest = S22-S21%*%solve(S11)%*%S12
        murest = S21 %*% solve(S11) %*% cbind(w)
    }

    ## Draw ys
    set.seed(0)
    wrests  = MASS::mvrnorm(n=ngen, mu=murest, Sigma=Sigmarest)
    ys = apply(wrests, 1, function(wrest){
        ysample = solve(Ag) %*% cbind(c(w,wrest))
        return(ysample)
    })

    ## Check polyhedron
    check.poly=TRUE
    if(check.poly){
        which.y.in.polyhedron =  apply(G%*%ys, 2, function(mycol){
            all(as.numeric(mycol) > u)
        })
        ys = ys[,which.y.in.polyhedron]
    }

    ## See that the next variable matches
    ## cp2.matches = apply(ys, 2, function(mycol){
    ##     fnew = dualpathSvd2(mycol, D=D, 2, approx=T)
    ##     cp2new = fnew$action[2]
    ##     return(cp2 == cp2new)
    ## })
    ## ys = ys[,cp2.matches]


    ## Take these samples, run the second step, collect the regularization parameter
    maxsteps = 2
    D = makeDmat(n,type='tf',ord=0)
    newlams = apply(ys, 2, function(mycol){
        ynew = as.numeric(mycol)
        fnew = dualpathSvd2(ynew, D=D, maxsteps, approx=T)
        newlambda = fnew$lambda[2]
        return(newlambda)
    })

    ## Get the current test statistic value
    d = sign(cp2)
    ## d = 1
    if(cp2 > cp1){
        test.stat.row = make.contrast(cp1,cp2,n,n,d)
    } else {
        test.stat.row = make.contrast(1,cp2,cp1,n,d)
    }

    test.stat.row = rbind(c(1,rep(0,n-1))) ## Temporarily added, nonrandom statistic
    mystat = as.numeric(rbind(test.stat.row)%*%cbind(y0))

    ## Take ys, run the second step, calculate the linear contrast
    D = makeDmat(n,type='tf',ord=0)
    newstat = apply(ys, 2, function(mycol){
        fnew = dualpathSvd2(mycol, D=D, 2, approx=T)
        cp1 = fnew$action[1]
        cp2 = fnew$action[2]
        d = sign(cp2)
        ## d = 1
        if(cp2 > cp1){
            test.stat.row = make.contrast(cp1,cp2,n,n,d)
        } else {
            test.stat.row = make.contrast(1,cp2,cp1,n,d)
        }
        return(rbind(test.stat.row)%*%cbind(mycol))
    })



    ## Calculate p-value
    pv = sum(newstat>mystat)/ncol(ys)
    return(pv)
}


## Actually run simulations
nsim = 100
pvs = rep(NA,nsim)
## for(ii in 100+(1:nsim)){
pvs = mclapply((1:nsim), function(isim){
    tryCatch({
        ## cat('\r', isim, 'out of', nsim)
        return(onesim(isim,opt=2))
    }, error=function(e) {print("error occurred")})
}, mc.cores=3)

## sapply(pvs, is.na)
## which(sapply(pvs, function(pv)pv==1 | pv==0))
pvs2 = unlist(pvs)


## Plot
hist(pvs2)
qqplot(y=pvs2,x=seq(from=0,to=1,length=length(pvs2)),xlab = "expected", ylab="observed")
abline(0,1)

## Why do I have boundary cases?
