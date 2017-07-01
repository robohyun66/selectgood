## Generate data
library(genlassoinf)
## set.seed(1)
n = 12
sigma = 2
ngen = 300
mn = rep(0,n)


onesim <- function(ii){

    ## Generate data
    set.seed(ii)
    y0 <- mn + rnorm(n,0,sigma)

    ## Run two-step fused lasso on it
    maxsteps = 2
    D = makeDmat(n,type='tf',ord=0)
    f0 = dualpathSvd2(y0, D=D, maxsteps, approx=T)

    ## Get naive 1-step poyhedron
    Gobj.naive = getGammat.naive(obj = f0, y = y0, condition.step = 1)
    G = Gobj.naive$G
    u = Gobj.naive$u
    cps = f0$action[1]   ## cps = c(4,9)
    inds = do.call(c,lapply(cps, function(cp)return(cp+c(-1,0,1))))
    other.inds = (1:n)[-inds]

    ## Helper functions
    get_other_segment_indices <- function(other.inds){
        brks = which(diff(other.inds)!=1)
        other.inds[brks]
        ilist = Map(function(a,b)(a+1):b, c(0,brks),c(brks,length(other.inds)))
        return( lapply(ilist, function(i)other.inds[i]))
    }
    make.contrast <- function(l, b, r, n){
        stopifnot(l <= b & b < r)
        ind1 = l:b
        ind2 = (b+1):r
        v = rep(0,n)
        v[ind1] = 1/length(ind1)
        v[ind2] = 1/length(ind2)
        return(rbind(v))
    }

    ## Get A = the map to the sufficient statistics + test statistic.
    plateaus = get_other_segment_indices(other.inds)
    suff.rows = do.call(rbind,lapply(plateaus, function(plt){v=rep(0,n); v[plt] = 1/length(plt); return(v)}))
    sat.rows = do.call(rbind, lapply(inds, function(ind){v=rep(0,n); v[ind] = 1;return(v)}))

    cp1 = f0$action[1]
    cp2 = f0$action[2]
    if(cp2 < cp1){
        test.stat.row = make.contrast(cp1,cp2,n,n)
    } else {
        test.stat.row = make.contrast(1,cp1,cp2,n)
    }
    A = rbind(suff.rows, sat.rows, test.stat.row)
    Asub = rbind(suff.rows, sat.rows)

    ## partition
    Sigma = (A) %*% diag(rep(sigma^2,n)) %*% t(A)
    Si = nrow(A)-1
    S11 = Sigma[1:Si, 1:Si]
    S12 = Sigma[1:Si, (Si+1):nrow(A)]
    S21 = Sigma[(Si+1):nrow(A), 1:Si]
    S22 = Sigma[(Si+1):nrow(A), (Si+1):nrow(A)]

    ## Make covariance for W_rest | W
    w = Asub %*% y0
    Sigmarest = S22-S21%*%solve(S11)%*%S12
    murest = S21 %*% solve(S11) %*% cbind(w)

    ## See if it is in the polyhedron.
    S = svd(t(A),nu=ncol(A))
    nr = nrow(suff.rows) + nrow(sat.rows) + nrow(test.stat.row)
    Arest = t(S$u[,(nr+1):ncol(A)])
    Ag = rbind(A,Arest)

    ## Calculate p-value via the closed form conditional null distribution.
    return(pv)
}


nsim = 500
pvs = rep(NA,nsim)
for(ii in 1:nsim){
    tryCatch({
        print(ii)
        pvs[ii] <- onesim(ii)
    }, error=function(e) {print("error occurred")})
}

ii=453
        pvs[ii] <- onesim(ii)

hist(pvs)

##

