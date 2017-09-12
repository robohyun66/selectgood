## Helper
getstuff <- function(y, type, teststep){

    n = length(y)

    ## Fit model, collect pvals BSFS.
    if(type=="binseg"){
        if(teststep!=0){
            g1 = binSeg_fixedSteps(y, numSteps=teststep)
            mypoly = polyhedra(g1)
        } else {
            g1 = NULL
            mypoly = NULL
        }
        g2 = binSeg_fixedSteps(y, numSteps=teststep+1)
    }
    if(type=="fusedlasso"){
        if(teststep!=0){
            g1 = dualpathSvd2(y, D = makeDmat(n, ord=0, type='tf'), maxsteps=teststep)
            mypoly = polyhedra(obj=g1$Gobj.naive$G,
                           u =  g1$Gobj.naive$u)
        } else {
            g1 = NULL
            mypoly = NULL
        }
        g2 = dualpathSvd2(y, D = makeDmat(n, ord=0, type='tf'), maxsteps=teststep+1)
    }

    ## Collect various things
    if(teststep==0){
        mypoly = polyhedra(obj=rbind(rep(0,n)),
                           u=0)
        cp1 = cp1.sign = c()
    } else {
        cp1 = g1$cp
        cp1.sign = g1$cp.sign
    }

    cp2 = g2$cp[which(!(g2$cp%in%cp1))]
    cp2.sign = g2$cp.sign[which(!(g2$cp%in%cp1))]

    ## Return everything
    return(list(g1=g1,g2= g2,cp1=cp1,
                cp2= cp2,cp1.sign=cp1.sign,
                cp2.sign=cp2.sign,
                mypoly=mypoly))
}


## Helper function
get.two.models <- function(y, type, teststep){

    n = length(y)

    ## Fit binseg model, collect pvals BSFS.
    if(type=="binseg"){

        if(teststep==0){
            obj.curr=NULL
            poly.curr = polyhedra(obj=rbind(rep(0,n)),
                                  u=0)
        cp.curr = cp.curr.sign = c()
        } else {
            obj.curr = binSeg_fixedSteps(y, numSteps=teststep)
            poly.curr = polyhedra(obj.curr)
            cp.curr = obj.curr$cp
            cp.curr.sign = obj.curr$cp.sign
        }
        obj.next = binSeg_fixedSteps(y, numSteps=teststep+1)
    } else {
        stop("type not written yet!")
    }

    cp.next = obj.next$cp[which(!(obj.next$cp%in%cp.curr))]
    cp.next.sign = obj.next$cp.sign[which(!(obj.next$cp%in%cp.curr))]

    return(list(obj.curr=obj.curr,
                obj.next=obj.next,
                cp.curr=cp.curr,
                cp.next=cp.next,
                cp.curr.sign=cp.curr.sign,
                cp.next.sign=cp.next.sign,
                poly.curr=poly.curr))
}


##' Get conditional Gaussian parameters, given current model object.
get.cond.gauss.param = function(obj.curr, y0, mn = rep(0,length(y0))){

    ## Form the sufficient statistic linear operator
    plateaus = get_plateaus(obj.curr$cp.curr, n)
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
    return(list(muorig=muorig, Sigmaorig=Sigmaorig))
}



## Apply forwardstop.
myforwardstop <- function(pvseq, alpha){
    fdp = sapply(1:length(pvseq), function(k){
        ((-1/k) * sum(log(1-pvseq)[1:k]))
    })
    return(which.max(fdp<=alpha))
}
