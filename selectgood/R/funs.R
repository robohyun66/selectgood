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
