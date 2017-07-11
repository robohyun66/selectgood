## Generate data
library(genlassoinf)

onesim <- function(ii=NULL, ngen=300){

    ## Other simulation settings
    n = 16
    sigma = 1
    delta = 1
    ## mn = rep(0,n)
    mn = c(rep(0,n/2),rep(delta,n/2))

    ## Generate data
    if(!is.null(ii))  set.seed(ii)
    y0 <- mn + rnorm(n, 0, sigma)

    ## Run two-step fused lasso on it
    maxsteps = 3
    D = makeDmat(n, type='tf',ord=0)
    f0 = dualpathSvd2(y0, D=D, maxsteps, approx=T)


    nstep = 3
    for(istep in 1:nstep){
        ## TODO: run (variant of) innersim in here.
    }

    innersim = function(istep,f0,y0){## Get naive 1-step poyhedron
        istep=1
        Gobj.naive = getGammat.naive(obj = f0, y = y0, condition.step = istep)
        G = Gobj.naive$G
        u = Gobj.naive$u
        cps = f0$action[1:istep]
        inds = do.call(c,lapply(cps, function(cp)return(cp+c(-1,0,1))))
        inds = inds[inds!=0]
        other.inds = (1:n)[-inds]

        ## Get A = the map to the sufficient statistics
        plateaus = get_other_segment_indices(other.inds)
        suff.rows = do.call(rbind,lapply(plateaus, function(plt){v=rep(0,n); v[plt] = 1/length(plt); return(v)}))
        sat.rows = do.call(rbind, lapply(inds, function(ind){v=rep(0,n); v[ind] = 1;return(v)}))
        A = rbind(suff.rows, sat.rows)

        ## SVD to get the complement null space
        S = svd(t(A),nu=ncol(A))
        nr = nrow(suff.rows) + nrow(sat.rows)
        A_rest = t(S$u[,(nr+1):ncol(A)])

        ## Build A_aug.
        Ag = rbind(A,A_rest)

        ## Make mean + covariance for W_rest | W
        w = A %*% y0
        a = get_condit_distr_y(A,Ag,sigma,n)
        cond.cov <- a$cov
        cond.mn <- a$mn

        ## Generate new y's and accept-reject
        ys = t(MASS::mvrnorm(n=ngen,mu=cond.mn, Sigma=cond.cov))
        check.poly=TRUE
        if(check.poly){
            which.y.in.polyhedron =  apply(G%*%ys, 2, function(mycol){
                all(as.numeric(mycol) > u)
            })
            ys = ys[,which.y.in.polyhedron, drop=FALSE]
        }

        ## Make next segment test statistic
        d = (f0$pathobj$s)[istep+1]#sign(cp2)
        locs = f0$action[1:(istep+1)]
        newloc = f0$action[istep+1]
        my.v.lrt = make.v.tf.fp(test.knot = newloc,
                                adj.knot  = locs,
                                test.knot.sign = d,
                                D=D)

        ## Calculate the p-value
        mystat = as.numeric(rbind(test.stat.row)%*%cbind(y0))
        newstat = apply(ys, 2, function(mycol){
            return(rbind(test.stat.row)%*%cbind(mycol))
        })
        pvs[istep] = sum(newstat>mystat)/length(newstat)

        ## Naive pv
        munaive = rbind(test.stat.row)%*%cbind(mn)
        sigmanaivesquare = (rbind(test.stat.row)%*%t(test.stat.row)) * (sigma^2)
        pvsnaive[istep] = 1-pnorm(mystat, mean=munaive, sd=sqrt(sigmanaivesquare))
    }

    return(list(pv=pv,pvnaive=pvnaive, cps = c(cp1,cp2)))
}


## Actually run simulations
nsim = 1000
pvs = mclapply((1:nsim), function(isim){
        cat('\r', isim, 'out of', nsim)
    tryCatch({
        return(onesim(ngen=1000))
    }, error=function(e) {print("error occurred")})
}, mc.cores=3)

## Extract things
pvsadj = sapply(pvs,function(a)a[["pv"]])
pvsnaive = sapply(pvs,function(a)a[["pvnaive"]])
cps1 = sapply(pvs,function(a)a[["cps"]])[1,]
cps2 = sapply(pvs,function(a)a[["cps"]])[2,]

## Make plot
qqplot(x=pvsnaive,y=seq(from=0,to=1,length=length(pvsnaive)),xlab = "expected", ylab="observed",ylim=c(0,1),xlim=c(0,1))
par(new=TRUE)
qqplot(x=pvsadj,y=seq(from=0,to=1,length=length(pvsadj)),xlab = "expected", ylab="observed",ylim=c(0,1),xlim=c(0,1),col='red')
abline(0,1,col='black')
legend("topleft", col=c('black','red'),pch=c(16,16), legend = c("naive", "adjusted"))



