## Helper
getstuff <- function(y, type, teststep){
    ## Fit model, collect pvals BSFS.
    if(type=="binseg"){
        if(teststep!=0){
            g1 = binSeg_fixedSteps(y0, numSteps=teststep)
            mypoly = polyhedra(g1)
        }
        g2 = binSeg_fixedSteps(y, numSteps=teststep+1)
    }
    if(type=="fusedlasso"){
        if(teststep!=0){
            g1 = dualpathSvd2(y, D = makeDmat(n, ord=0, type='tf'), maxsteps=teststep)
            mypoly = polyhedra(obj=g1$Gobj.naive$G,
                           u =  g1$Gobj.naive$u)
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
