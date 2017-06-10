#----------------------------------------------------------------
##' Computes the equilibrium of three types of games, given a matrix of objectives (or a set of matrices) and the structure of the strategy space.
##' @title Equilibrium computation of a discrete game for a given matrix with objectives values
##' @param Z is a matrix of size [\code{npts x nsim*nobj}] of objective values, see details,
##' @param equilibrium considered type, one of \code{"NE"}, \code{"NKSE"}, \code{"KSE"}
##' @param nobj nb of objectives (or players)
##' @param n.s scalar of vector. If scalar, total number of strategies (to be divided equally among players),
##'  otherwise number of strategies per player.
##' @param expanded.indices is a matrix containing the indices of the \code{integ.pts} on the grid, see \code{\link[GPGame]{generate_integ_pts}}
##' @param return.design Boolean; if \code{TRUE}, the index of the optimal strategy is returned (otherwise only the pay-off is returned)
##' @param sorted Boolean; if \code{TRUE}, the last column of expanded.indices is assumed to be sorted in increasing order. This provides a substantial efficiency gain.
##' @param cross Should the simulation be crossed? (May be dropped in future versions)
##' @details If \code{nsim=1}, each line of \code{Z} contains the pay-offs of the different players for a given strategy s: [obj1(s), obj2(s), ...].
##' The position of the strategy \code{s} in the grid is given by the corresponding line of \code{expanded.indices}. If \code{nsim>1}, (vectorized call) \code{Z} contains
##' different trajectories for each pay-off: each line is [obj1_1(x), obj1_2(x), ... obj2_1(x), obj2_2(x), ...].
##' @export
##' @importFrom Rcpp evalCpp
##' @examples
##' \dontrun{
##' ## Setup
##' fun <- function (x)
##' {
##'   if (is.null(dim(x)))    x <- matrix(x, nrow = 1)
##'   b1 <- 15 * x[, 1] - 5
##'   b2 <- 15 * x[, 2]
##'   return(cbind((b2 - 5.1*(b1/(2*pi))^2 + 5/pi*b1 - 6)^2 + 10*((1 - 1/(8*pi)) * cos(b1) + 1),
##'                -sqrt((10.5 - b1)*(b1 + 5.5)*(b2 + 0.5)) - 1/30*(b2 - 5.1*(b1/(2*pi))^2 - 6)^2-
##'                 1/3 * ((1 - 1/(8 * pi)) * cos(b1) + 1)))
##' }
##'
##' d <- nobj <- 2
##'
##' # Generate grid of strategies for Nash and Nash-Kalai-Smorodinsky
##' n.s <- c(11,11) # number of strategies per player
##' x.to.obj <- 1:2 # allocate objectives to players
##' integcontrol <- generate_integ_pts(n.s=n.s,d=d,nobj=nobj,x.to.obj=x.to.obj,gridtype="cartesian")
##' integ.pts <- integcontrol$integ.pts
##' expanded.indices <- integcontrol$expanded.indices
##'
##' # Compute the pay-off on the grid
##' response.grid <- t(apply(integ.pts, 1, fun))
##'
##' # Compute the Nash equilibrium (NE)
##' trueEq <- getEquilibrium(Z = response.grid, equilibrium = "NE", nobj = nobj, n.s = n.s,
##'                          return.design = TRUE, expanded.indices = expanded.indices,
##'                          sorted = !is.unsorted(expanded.indices[,2]))
##'
##' # Pay-off at equilibrium
##' print(trueEq$NEPoff)
##'
##' # Optimal strategy
##' print(integ.pts[trueEq$NE,])
##'
##' # Index of the optimal strategy in the grid
##' print(expanded.indices[trueEq$NE,])
##'
##' # Plots
##' par(mfrow = c(1,2))
##' plotGameGrid(fun = fun, n.grid = n.s, x.to.obj = x.to.obj, integcontrol=integcontrol,
##'              equilibrium = "NE")
##'
##' # Compute KS equilibrium (KSE)
##' trueKSEq <- getEquilibrium(Z = response.grid, equilibrium = "KSE", nobj = nobj,
##'                          return.design = TRUE, sorted = !is.unsorted(expanded.indices[,2]))
##'
##' # Pay-off at equilibrium
##' print(trueKSEq$NEPoff)
##'
##' # Optimal strategy
##' print(integ.pts[trueKSEq$NE,])
##'
##' plotGameGrid(fun = fun, n.grid = n.s, integcontrol=integcontrol,
##'              equilibrium = "KSE", fun.grid = response.grid)
##'
##' # Compute the Nash equilibrium (NE)
##' trueNKSEq <- getEquilibrium(Z = response.grid, equilibrium = "NKSE", nobj = nobj, n.s = n.s,
##'                          return.design = TRUE, expanded.indices = expanded.indices,
##'                          sorted = !is.unsorted(expanded.indices[,2]))
##'
##' # Pay-off at equilibrium
##' print(trueNKSEq$NEPoff)
##'
##' # Optimal strategy
##' print(integ.pts[trueNKSEq$NE,])
##'
##' # Index of the optimal strategy in the grid
##' print(expanded.indices[trueNKSEq$NE,])
##'
##' # Plots
##' plotGameGrid(fun = fun, n.grid = n.s, x.to.obj = x.to.obj, integcontrol=integcontrol,
##'              equilibrium = "NKSE")
##' }
getEquilibrium <- function(Z, equilibrium = c("NE", "NKSE", "KSE"), nobj=2, n.s, expanded.indices=NULL, return.design=FALSE,
                           sorted=FALSE, cross=FALSE){
  #### Choose appropriate function ###################
  if (equilibrium=="NE"){
    return(getNashEquilibrium(Z = Z, nobj = nobj, n.s = n.s, expanded.indices = expanded.indices, return.design = return.design, sorted = sorted, cross = cross))
  } else if (equilibrium=="KSE") {
    return(getKSequilibrium(Z = Z, nobj = nobj, n.s = n.s, expanded.indices = expanded.indices, return.design = return.design, sorted = sorted, cross = cross))
  } else if (equilibrium=="NKSE") {
    return(getNKSequilibrium(Z = Z, nobj = nobj, n.s = n.s, expanded.indices = expanded.indices, return.design = return.design, sorted = sorted, cross = cross))
  } else {
    cat("wrong crit \n")
    break;
  }
}

#----------------------------------------------------------------
## ' Computes the equilibrium of finite Kalai-Smorodinski games given a matrix of objectives (or a set of matrices) and the structure of the strategy space.
## ' @title Nash equilibrium computation
## ' @param Z is a matrix of size [npts x nsim*nobj] of objective values, see details,
## ' @param nobj nb of objectives (or players)
## ' @param n.s scalar of vector. If scalar, total nb of strategies (to be divided equally among players), otherwise nb of strategies per player.
## ' @param expanded.indices is a matrix containing the indices of the integ.pts on the grid
## ' @param return.design Boolean; if TRUE, the index of the optimal strategy is returned (otherwise only the pay-off is returned)
## ' @param sorted Boolean; if TRUE, the last column of expanded.indices is assumed to be sorted in increasing order. This provides a substantial efficiency gain.
## ' @param cross if TRUE, all the combinations of trajectories are used
## ' @param ... not used, for compatibility
## ' @details If \code{nsim}=1, each line of Z contains the pay-offs of the different players for a given strategy s: [obj1(s), obj2(s), ...].
## ' The position of the strategy s in the grid is given by the corresponding line of \code{expanded.indices}. If \code{nsim}>1, (vectorized call) Z contains
## ' different trajectories for each pay-off: each line is [obj1_1(x), obj1_2(x), ... obj2_1(x), obj2_2(x), ...].
## '
## ' @export

getNashEquilibrium <- function(Z, nobj=2, n.s, expanded.indices=NULL, return.design=FALSE, sorted=FALSE, cross=FALSE,...){

  if(cross){
    nsim <- ncol(Z)/nobj

    if(!sorted)
      print("Non sorted case not implemented with cross")

    combisim <- NULL
    for(i in 1:nobj)
      combisim <- c(combisim, list(0:(nsim - 1))) ## starts at 0 for Rcpp indices compatibility
    combisim <- as.matrix(expand.grid(combisim))
    NE <- PSNE_sparseMat_cross(n.s, Z, expanded.indices - as.integer(1), combisim = combisim, ncross = nsim^(nobj-1))

    if (return.design == FALSE) return(getPoffsCross(isNash = NE, Poffs = Z, combisim = combisim, nsim = nsim))
    else                     return(list(NEPoff = getPoffsCross(isNash = NE, Poffs = Z, combisim = combisim, nsim = nsim), NE = unlist(apply(NE, 2, which))))

  }


  if (sorted) {
    NE <- PSNE_sparseMat_sorted(n.s, Z, expanded.indices - as.integer(1)) ## -1 for compatibility with Rcpp
  } else {
    NE <- PSNE_sparseMat(n.s, Z, expanded.indices - as.integer(1)) ## -1 for compatibility with Rcpp
  }

  NEPoff <- matrix(NA, 1, nobj)

  if (!is.null(NE) && length(which(NE)) > 0) {
    if (!return.design) return(getPoffs(NE, Z, nsim = ncol(Z)/nobj, nobj))
    else return(list(NEPoff = getPoffs(NE, Z, nsim = ncol(Z)/nobj, nobj), NE = unlist(apply(NE, 2, which))))
  }else{
    if(return.design)
      return(list(NEPoff = NEPoff, NE = NA))
    return(NEPoff)
  }
}
#----------------------------------------------------------------
## ' Computes the equilibrium of finite Nash/Kalai-Smorodinski games given a matrix of objectives (or a set of matrices) and the structure of the strategy space.
## ' @title Nash/Kalai-Smorodinski equilibrium
## ' @param Z is a matrix of size [npts x nsim*nobj] of objective values, see details,
## ' @param nobj nb of objectives (or players)
## ' @param n.s scalar of vector. If scalar, total nb of strategies (to be divided equally among players), otherwise nb of strategies per player.
## ' @param expanded.indices is a matrix containing the indices of the integ.pts on the grid
## ' @param return.design Boolean; if TRUE, the index of the optimal strategy is returned (otherwise only the pay-off is returned)
## ' @param cross if TRUE, all the combinations of trajectories are used
## ' @param ... not used, for compatibility
## ' @details If \code{nsim}=1, each line of Z contains the pay-offs of the different players for a given strategy s: [obj1(s), obj2(s), ...].
## ' The position of the strategy s in the grid is given by the corresponding line of \code{expanded.indices}. If \code{nsim}>1, (vectorized call) Z contains
## ' different trajectories for each pay-off: each line is [obj1_1(x), obj1_2(x), ... obj2_1(x), obj2_2(x), ...].
## '
## ' @export
getNKSequilibrium <- function(Z, nobj=2, n.s, return.design=FALSE, expanded.indices=NULL, cross=FALSE, ...){

  # allShadow <- getNashEquilibrium(Z=Z, nobj=nobj, n.s=n.s, expanded.indices=expanded.indices, cross=cross)
  # if (any(is.na(allShadow))) {
  #   NEPoff <- rep(NA, nobj)
  #   NE <- NA
  # } else {

  nsim <- ncol(Z) / nobj

  if (cross) {
    indices <- list()
    for (u in 1:nobj) indices[[u]] <- (1:nsim)+(u-1)*nsim
    Jmat <- as.matrix(expand.grid(indices))
  } else {
    Jmat <- matrix(NA, nsim, nobj)
    for (u in 1:nsim) {
      Jmat[u,] <- seq(u, ncol(Z), nsim)
    }
  }

  # NEPoff <- matrix(NA, nrow(Jmat), nobj)
  # NE     <- rep(NA, nrow(Jmat))
  NE <- NEPoff <- NULL

  for (u in 1:nrow(Jmat)) {

    # Only look at first NE if there are several
    Shadow <- getNashEquilibrium(Z=Z, nobj=nobj, n.s=n.s, expanded.indices=expanded.indices, cross=cross)[1,]

    J <- Jmat[u,]
    I <- which(!is_dominated(t(Z[,J])))
    Zred <- Z[I, J, drop=FALSE]
    Nadir  <- apply(Zred, 2, min)

    # Shadow <- allShadow[u,]

    if (nrow(Zred)==1) {
      i <- 1
    } else {
      alldist2 <- rowSums((Zred - matrix(rep(Nadir, nrow(Zred)), ncol=nobj, byrow=T))^2) -
        as.numeric(((Zred - matrix(rep(Nadir, nrow(Zred)), ncol=nobj, byrow=T))%*%(Nadir - Shadow))^2) /
        drop(crossprod(Nadir - Shadow, Nadir - Shadow))
      i <- which.min(alldist2)
    }
    NEPoff <- rbind(NEPoff, Zred[i,])
    NE     <- c(NE, I[i])
  }

  if(is.null(NEPoff)){
    NEPoff <- matrix(NA, nrow(Jmat), nobj)
    NE     <- rep(NA, nrow(Jmat))
  }

  # }
  if (return.design==FALSE) return(NEPoff)
  else                     return(list(NEPoff=NEPoff, NE=I[i]))
}
#----------------------------------------------------------------
## ' Computes the equilibrium of finite Kalai-Smorodinski games given a matrix of objectives (or a set of matrices) and the structure of the strategy space.
## ' @title Kalai-Smorodinski equilibrium computation
## ' @param Z is a matrix of size [npts x nsim*nobj] of objective values, see details,
## ' @param nobj nb of objectives (or players)
## ' @param return.design Boolean; if TRUE, the index of the optimal strategy is returned (otherwise only the pay-off is returned)
## ' @param cross if TRUE, all the combinations of trajectories are used
## ' @param ... not used, for compatibility
## ' @details If \code{nsim}=1, each line of Z contains the pay-offs of the different players for a given strategy s: [obj1(s), obj2(s), ...].
## ' The position of the strategy s in the grid is given by the corresponding line of \code{expanded.indices}. If \code{nsim}>1, (vectorized call) Z contains
## ' different trajectories for each pay-off: each line is [obj1_1(x), obj1_2(x), ... obj2_1(x), obj2_2(x), ...].
## '
## ' @export
getKSequilibrium <- function(Z, nobj=2, return.design=FALSE, cross=FALSE, ...){

  nsim <- ncol(Z) / nobj

  if (cross) {
    indices <- list()
    for (u in 1:nobj) indices[[u]] <- (1:nsim)+(u-1)*nsim
    Jmat <- as.matrix(expand.grid(indices))
  } else {
    Jmat <- matrix(NA, nsim, nobj)
    for (u in 1:nsim) {
      Jmat[u,] <- seq(u, ncol(Z), nsim)
    }
  }

  NEPoff <- matrix(NA, nrow(Jmat), nobj)
  NE     <- rep(NA, nrow(Jmat))

  for (u in 1:nrow(Jmat)) {
    J <- Jmat[u,]
    I <- which(!is_dominated(t(Z[,J])))
    Zred <- Z[I, J, drop=FALSE]
    Nadir  <- apply(Zred, 2, min)
    Shadow <- apply(Zred, 2, max)

    if (nrow(Zred)==1) {
      i <- 1
    } else {
      alldist2 <- rowSums((Zred - matrix(rep(Nadir, nrow(Zred)), ncol=nobj, byrow=T))^2) -
        as.numeric(((Zred - matrix(rep(Nadir, nrow(Zred)), ncol=nobj, byrow=T))%*%(Nadir - Shadow))^2) /
        drop(crossprod(Nadir - Shadow, Nadir - Shadow))
      i <- which.min(alldist2)
    }
    NEPoff[u,] <- Zred[i,,drop = FALSE]
    NE[u]     <- I[i]
  }
  if (return.design==FALSE) return(NEPoff)
  else                     return(list(NEPoff=NEPoff, NE=I[i]))
}
