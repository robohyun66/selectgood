## Helper function
get_other_segment_indices <- function(other.inds){
    brks = which(diff(other.inds)!=1)
    ilist = Map(function(a,b)(a+1):b, c(0,brks),c(brks,length(other.inds)))
    return( lapply(ilist, function(i)other.inds[i]))
}


## Helper function
make.contrast <- function(l, b, r, n, direction=+1){
    stopifnot(l <= b & b < r)
    ind1 = l:b
    ind2 = (b+1):r
    v = rep(0,n)
    v[ind1] = -1/length(ind1)
    v[ind2] = 1/length(ind2)
    return(rbind(v*direction))
}


## Say vector y comes from N(mu,sm). This function obtains the distribution
## (mean and covariance matrix) of y|Ay=Ay0. This is done by obtaining a
## row-augmented version of A, Ag, that makes Ag full row rank (call the new
## rows A'); then obtaining samples from A'y|Ay the multivariate normal; then
## finally inverting the solutions back via left multiplication with A'^{-1}.
## @param A linear transformation matrix of interest
## @param mn mean vector of y
## @param cova covariance matrix of y
## @param y0 original drawn instance of y
get_condit_distr_y <- function(A, mn, cova, y0){

    ## Basic checks n = ncol(A)
    n = length(y0)
    stopifnot(n == length(mn))
    stopifnot(n == nrow(cova))
    stopifnot(n == ncol(cova))
    w = A %*% y0

    ## SVD to get the complement null space
    S = svd(t(A),nu=ncol(A))
    nr = nrow(A)
    A_rest = t(S$u[,(nr+1):ncol(A)])

    ## Build A_aug.
    Ag = rbind(A,A_rest)

    ## Partition matrix into four blocks
    Sigma = (Ag) %*% cova %*% t(Ag)
    Si = nrow(A)
    S11 = Sigma[1:Si, 1:Si]
    S12 = Sigma[1:Si, (Si+1):n]
    S21 = Sigma[(Si+1):n, 1:Si]
    S22 = Sigma[(Si+1):n, (Si+1):n]

    ## Calculate distr of (A_rest y| Ay) using Schur's complement
    Sigmarest = S22 - S21%*%solve(S11, S12)
    murest = A_rest%*%mn + S21 %*% solve(S11, cbind(w - A %*% mn))

    ## Transform this back to to the original y scale.
    Aginv = solve(Ag)
    Sigmanew = matrix(0, nrow=n, ncol=n)
    Sigmanew[(Si+1):n, (Si+1):n] = Sigmarest
    Sigmaorig = Aginv %*% Sigmanew %*% t(Aginv)
    muorig = Aginv %*%  cbind(c(w, murest))

    return(list(mu=muorig, cov=Sigmaorig))
}


## Function to plot QQ plot of p-values, against uniform(0,1) distribution. Use
## extra parameters for plotting.
## @param pp numeric vector of p-values.
## @param main label to plot as main title.
## @param add TRUE if you want to overlay on existing plots, by activating
##     \code{par(new=TRUE)}. Defaults to FALSE.
## @param ... other graphical parameters for \code{plot()}.
qqunif <- function(pp, main=NULL,add=FALSE,returnxy=FALSE,...){
    xy <- stats::qqplot(x=pp,
                 y=seq(from=0,to=1,length=length(pp)), plot.it=FALSE)
    if(add){ graphics::par(new=TRUE) }
    graphics::plot(xy, axes=FALSE, ylim=c(0,1),xlim=c(0,1),...)
    if(!add){
        graphics::axis(2); graphics::axis(1)
        graphics::abline(0,1)
    }
    if(!is.null(main)) graphics::title(main=main)
    if(returnxy) return(xy)
}
