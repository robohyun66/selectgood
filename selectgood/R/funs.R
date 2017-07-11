##' Helper function
get_other_segment_indices <- function(other.inds){
    brks = which(diff(other.inds)!=1)
    ilist = Map(function(a,b)(a+1):b, c(0,brks),c(brks,length(other.inds)))
    return( lapply(ilist, function(i)other.inds[i]))
}

##' Helper function
make.contrast <- function(l, b, r, n, direction=+1){
    stopifnot(l <= b & b < r)
    ind1 = l:b
    ind2 = (b+1):r
    v = rep(0,n)
    v[ind1] = -1/length(ind1)
    v[ind2] = 1/length(ind2)
    return(rbind(v*direction))
}


##' Gets the distribution of y|Ay, where Ag is a row-augmented version of A so
##' that it has full row rank.
get_condit_distr_y <- function(A, Ag, sigma, n){
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

    return(list(mn=mn, cov=cov))
}
