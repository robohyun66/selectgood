## Generate data
library(genlassoinf)
## set.seed(1)
n = 12
sigma = 1
ngen = 300
## mn = c(+5, +5, rep(0,n-2))
mn = rep(0,n)



## Fixed polyhedron
G = rbind(c(1,-1,rep(0,n-2)))
u = 2


## Fix the test statistic row.
test.stat.row = c(1,-1,rep(0,n-2))


onesim <- function(ii=NULL,check.poly=FALSE,G=NULL,u=NULL){
    ## Generate data
    if(!is.null(ii))  set.seed(ii)
    ## set.seed(0)
    y0 <- mn + rnorm(n,0,sigma)
    ## seed = 0
    if(check.poly){
        while(G%*%y0 < u){
            ## seed = seed + 1
            ## print(seed)
            ## set.seed(seed);
            y0 <- mn + rnorm(n,0,sigma)
        }
    }

    ## ## Run two-step fused lasso on it
    ## maxsteps = 2
    ## D = makeDmat(n,type='tf',ord=0)
    ## f0 = dualpathSvd2(y0, D=D, maxsteps, approx=T)

    ## ## Get naive 1-step poyhedron
    ## cps = f0$action[1]   ## cps = c(4,9)

    cps = 4
    inds = do.call(c,lapply(cps, function(cp)return(cp+c(-1,0,1))))
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

    ## Build the rest of the rows to be the span of the null space.
    B = matrix(0, nrow = 7, ncol = 12)
    B[1,1:2] = c(1,-1)
    for(jj in 6:11){
        B[jj-4,(jj):(jj+1)] = c(1,-1)
    }

    ## ## Build A_rest.
    ## rest.inds = unlist(sapply(plateaus, function(plt)plt[-1]))
    ## A_rest= do.call(rbind, lapply(rest.inds, function(ind){v=rep(0,n); v[ind] = 1;return(v)}))

    ## Build A_aug.
    Ag = rbind(A,B)#A_rest)
    ## pivot = apply(Ag,1,function(myrow){min(which(myrow!=0))}) ## Reordering is not needed
    ## Ag[order(pivot),]

    ## partition
    Sigma = (Ag) %*% diag(rep(sigma^2,n)) %*% t(Ag)
    Si = nrow(A)
    S11 = Sigma[1:Si, 1:Si]
    S12 = Sigma[1:Si, (Si+1):n]
    S21 = Sigma[(Si+1):n, 1:Si]
    S22 = Sigma[(Si+1):n, (Si+1):n]

    ## ## Make covariance for W_rest | W
    w = A %*% y0
    ## Sigmarest = S22-S21%*%solve(S11)%*%S12
    ## murest = S21 %*% solve(S11) %*% cbind(w)
    Sigmarest = B%*%t(B)
    murest = rep(0,nrow(B))

    ## Draw from Y|W=w via sampling and inverting
    wrests  = MASS::mvrnorm(n=ngen, mu=murest, Sigma=Sigmarest)
    ys = apply(wrests, 1, function(wrest){
        ysample = solve(Ag) %*% cbind(c(w,wrest))
        return(ysample)
    })

    ## Rejection sampling to only retain ys that are in polyhedron
    if(check.poly){
        which.y.in.polyhedron =  apply(G%*%ys, 2, function(mycol){
            all(as.numeric(mycol) > u)
        })
        ys = ys[,which.y.in.polyhedron]
    }

    ## Calculate T(y) for y in ys
    newstat = apply(ys, 2, function(mycol){
        return(as.numeric(rbind(test.stat.row)%*%cbind(mycol)))
    })

    ## Calculate T(y0)
    mystat = as.numeric(rbind(test.stat.row)%*%cbind(y0))

    ## Calculate p-value
    pv = sum(newstat>mystat)/ncol(ys)
    return(pv)
}


nsim = 1000
pvs = rep(NA,nsim)
for(ii in (1:nsim)){
    tryCatch({
        print(ii)
        pvs[ii] <- onesim(ii,FALSE,G,u)
    }, error=function(e) {print("error occurred")})
}

hist(pvs)
qqplot(x=pvs,y=seq(from=0,to=1,length=nsim))
abline(0,1)


## The polyhedron messes with this.
-
