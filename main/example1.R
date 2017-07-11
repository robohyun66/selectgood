## Synopsis: Generate data from a particular mean, test the first changepoint
## model's goodness of fit. When

## Generate data
library(genlassoinf)
n = 12
sigma = 1
## mn = rep(0,n)
delta = 0
mn = c(rep(0,n/2),rep(delta,n/2))
## mn = c(rep(0,n/3),rep(delta,n/3),rep(0,n/3))

onesim <- function(ii=NULL, ngen=300){

    ## Generate data
    if(!is.null(ii))  set.seed(ii)
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
    cps = cp1
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

    ## Partition matrix into four blocks
    Sigma = (Ag) %*% diag(rep(sigma^2,n)) %*% t(Ag)
    Si = nrow(A)
    S11 = Sigma[1:Si, 1:Si]
    S12 = Sigma[1:Si, (Si+1):n]
    S21 = Sigma[(Si+1):n, 1:Si]
    S22 = Sigma[(Si+1):n, (Si+1):n]

    ## Calculate using Schur's complement
    Sigmarest = S22 - S21%*%solve(S11, S12)
    murest = S21 %*% solve(S11, cbind(w))

    ## Transform to the conditional distribution
    Aginv = solve(Ag)
    Sigmanew = matrix(0, nrow=n, ncol=n)
    Sigmanew[(Si+1):n, (Si+1):n] = Sigmarest
    Sigmaorig = Aginv %*% Sigmanew %*% t(Aginv)
    muorig = solve(Ag, c(w, murest))

    ## Generate new y's and accept-reject
    ys = t(MASS::mvrnorm(n=ngen,mu=muorig,Sigma=Sigmaorig))
    ## which.y.in.polyhedron =  apply(G%*%ys, 2, function(mycol){
    ##     all(as.numeric(mycol) > u)
    ## })
    ## ys = ys[,which.y.in.polyhedron, drop=FALSE]

    ## Do new accept-reject that checks for changepoint identity only
    which.y.has.same.changepoint = unlist(apply(ys, 2, function(my.y){
        a = genlasso::fusedlasso1d(my.y,maxsteps=2)
        gaps = D%*%a$beta[,2]
        yscp = which(abs(gaps)>1E-10)
        return(yscp==cp1)
    }))
    ## which.y.in.polyhedron =  apply(G%*%ys, 2, function(mycol){
    ##     all(as.numeric(mycol) > u)
    ## })
    ys = ys[,which.y.has.same.changepoint, drop=FALSE]



    ## Make next segment test statistic
    d = (f0$pathobj$s)[2]#sign(cp2)
    if(cp2 > cp1){
        test.stat.row = make.contrast(cp1,cp2,n,n,d)
    } else {
        test.stat.row = make.contrast(1,cp2,cp1,n,d)
    }

    ## ## Make /fixed/ test statistic
    ## test.stat.row = rep(0,n)
    ## test.stat.row[c(1,9)] = 1
    ## test.stat.row = rbind(test.stat.row)

    ## Calculate the p-value
    mystat = as.numeric(rbind(test.stat.row)%*%cbind(y0))
    newstat = apply(ys, 2, function(mycol){
        return(rbind(test.stat.row)%*%cbind(mycol))
    })
    pv = sum(newstat>mystat)/length(newstat)

    ## Naive pv
    munaive = rbind(test.stat.row)%*%cbind(mn)
    sigmanaivesquare = (rbind(test.stat.row)%*%t(test.stat.row)) * (sigma^2)
    pvnaive = 1-pnorm(mystat, mean=munaive, sd=sqrt(sigmanaivesquare))

    return(list(pv=pv,pvnaive=pvnaive, cps = c(cp1,cp2)))
}


## Actually run simulations
nsim = 200
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
