## Generate data
library(genlassoinf)
n = 12
sigma = 1
## mn = rep(0,n)
delta = 1
mn = c(rep(0,n/2),rep(delta,n/2))

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

    ## Make next segment test statistic
    d = (f0$pathobj$s)[2]#sign(cp2)
    if(cp2 > cp1){
        test.stat.row = make.contrast(cp1,cp2,n,n,d)
    } else {
        test.stat.row = make.contrast(1,cp2,cp1,n,d)
    }

    ## Calculate the naive p-value
    mystat = as.numeric(rbind(test.stat.row)%*%cbind(y0))
    munaive = rbind(test.stat.row)%*%cbind(mn)
    sigmanaivesquare = (rbind(test.stat.row)%*%t(test.stat.row)) * (sigma^2)
    pvnaive = 1-pnorm(mystat, mean=munaive, sd=sqrt(sigmanaivesquare))

    return(list(pvnaive=pvnaive, cps = c(cp1,cp2)))
}


nsim = 1000
pvs = mclapply((1:nsim), function(isim){
        cat('\r', isim, 'out of', nsim)
    tryCatch({
        return(onesim(ngen=1000))
    }, error=function(e) {print("error occurred")})
}, mc.cores=3)
naive.pvs = unlist(sapply(pvs,function(a)a[1]))
qqplot(x=naive.pvs,y=seq(from=0,to=1,length=length(naive.pvs)),xlab = "expected", ylab="observed",ylim=c(0,1),xlim=c(0,1))
