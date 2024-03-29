#' Computes exactlty the probability of non-domination for a set of points up to 3 objectives.
#' Then it is approximated by sampling
#' @noRd
#' @title Probability of non-domination for a set of points (exact up to 3 objectives)
#' @param paretoFront (optional) matrix corresponding to the Pareto Front (one output per column).
#' @param model list of objects of class \code{\link[DiceKriging]{km}}, one for each objective functions,
#' @param integration.points Points to compute the probability
#' @param predictions An optional list of predictions (using \code{\link[DiceKriging]{predict.km}}) at
#'        integration points for each model.
#' @param nsamp number of samples to estimate the probability of non-domination, based on a set of points sampled from the predictive distribution
#' @details TODO: improve the sampling method
#' @return A vector of probabilities
## ' @export
## ' @useDynLib DiceMOO
## ' @references
## ' V. Picheny (2014), Multiobjective optimization using Gaussian process emulators via stepwise uncertainty reduction,
## ' \emph{Statistics and Computing}
## ' @examples
## ' #---------------------------------------------------------------------------
## ' # Probability of non-domination on a simple 1D problem.
## ' #---------------------------------------------------------------------------
## ' \donttest{
## ' set.seed(12)
## ' fun1D <- function(x){return(c(fundet(x), -fundet(12+2*x)))}
## '
## ' # Generate initial data and models
## ' design.init <- data.frame(runif(7))
## ' response.init    <- t(apply(design.init, 1, fun1D))
## ' mf1 <- km(~., design = design.init, response = response.init[,1])
## ' mf2 <- km(~., design = design.init, response = response.init[,2])
## '
## ' # Compute predictions on a grid
## ' test.grid   <- matrix(seq(0, 1, length.out =100),ncol=1)
## ' p1  <- predict(mf1, newdata=test.grid, type="UK", checkNames=FALSE)
## ' p2  <- predict(mf2, newdata=test.grid, type="UK", checkNames=FALSE)
## '
## ' # Compute probability of non domination
## ' p.nd <- prob.of.non.domination(model=list(mf1,mf2), integration.points=test.grid, predictions=list(p1,p2))
## '
## ' # Some plots
## ' par(mfrow=c(3,1))
## ' plot(design.init[[1]], response.init[,1], lwd=3, xlim=c(0,1), ylim=c(.2,1.2), main="model 1")
## ' polygon(x=c(test.grid,rev(test.grid)), y=c(p1$mean + 2*p1$sd,
## ' rev(p1$mean - 2*p1$sd)), border=NA,col="lightgrey")
## ' lines(test.grid, p1$mean)
## ' points(design.init[[1]], response.init[,1], lwd=3)
## ' plot(design.init[[1]], response.init[,2], lwd=3, xlim=c(0,1), ylim=c(-5000,2000), main="model 2")
## ' polygon(x=c(test.grid,rev(test.grid)), y=c(p2$mean + 2*p2$sd,
## ' rev(p2$mean - 2*p2$sd)), border=NA,col="lightgrey")
## ' lines(test.grid, p2$mean)
## ' points(design.init[[1]], response.init[,2], lwd=3)
## ' plot(test.grid, p.nd, type="l", , main="Prob of ND")
## ' }
## ' @export

prob.of.non.domination <- function(paretoFront=NULL, model=NULL, integration.points=NULL, predictions=NULL, nsamp = 100, target=NULL){
  ###########################################################################################
  # Computes the probability that integration points are non dominated by a given paretoFront
  # 2 or 3 objectives only
  # Providing the kriging predictions may accelerate considerably the evaluation
  ###########################################################################################
  if (is.null(nsamp)) nsamp <- 100
  
  if (is.null(model) && is.null(predictions)){
    stop("Error: either a list of models or a list of km predictions must be provided")
  }
  if ((is.null(model) || is.null(integration.points)) && is.null(predictions)){
    stop("Error: either models + integration points or km predictions must be provided")
  }
  
  n.integration.points <- max(nrow(integration.points), length(predictions[[1]]$mean))
  n.obj <- max(c(ncol(paretoFront), length(model), length(predictions)))
  
  # Compute current Pareto front if missing
  if (is.null(paretoFront)){
    if (length(model[[1]]@noise.var)>0) {
      # If noise.var exists: get "denoised" observations
      pred <- lapply(model, FUN=predict, newdata = model[[1]]@X, checkNames = FALSE, type = "UK", light.return = TRUE)
      observations <- Reduce(cbind, lapply(pred, function(alist) alist$mean))
    } else {
      observations <- Reduce(cbind, lapply(model, slot, "y"))
    }
    if (!is.null(target)) observations <- (observations - matrix(rep(target, nrow(observations)), byrow=TRUE, nrow=nrow(observations)))^2
    
    # paretoFront <- matrix(t(nondominated_points(t(observations))),ncol=n.obj)
    paretoFront <- nonDom(observations)
  }
  n.pareto <- nrow(paretoFront)
  
  # Generate kriging predictions if missing
  if ( is.null(predictions) ){
    predictions <- vector("list",n.obj)
    for (i in 1:n.obj) predictions[[i]] <- predict(object=model[[i]], newdata=integration.points, type="UK", checkNames=FALSE, cov.compute = FALSE, light.return = TRUE)
  }
  if (!is.null(target)) {
    # Calibration mode
    samps <- NULL
    for(j in 1:n.obj){
      samps <- cbind(samps, (target[j] - rnorm(nsamp*n.integration.points, mean = predictions[[j]]$mean, predictions[[j]]$sd))^2)
    }
    idx <- nonDom(samps, paretoFront, return.idx=TRUE)
    pn <- rowSums(matrix(tabulate(idx, nbins=nsamp*n.integration.points), n.integration.points))
    return(pn/nsamp)
  } else {
    # Precompute important quantities
    if(n.obj <= 3){
      
      phi.x.tilde <- pnorm( (matrix(rep(paretoFront[,1],n.integration.points), ncol = n.integration.points) - t(matrix(rep(predictions[[1]]$mean, n.pareto), ncol=n.pareto)) ) /
                              t(matrix(rep(predictions[[1]]$sd, n.pareto), ncol=n.pareto)) )
      if (max(predictions[[2]]$sd) != 0) {
        # Regular case
        phi.y.tilde <- pnorm( (matrix(rep(paretoFront[,2],n.integration.points), ncol = n.integration.points) - t(matrix(rep(predictions[[2]]$mean, n.pareto), ncol=n.pareto)) ) /
                                t(matrix(rep(predictions[[2]]$sd, n.pareto), ncol=n.pareto)) )
      } else {
        # y is a fast function
        phi.y.tilde <- 1*( (matrix(rep(paretoFront[,2],n.integration.points), ncol = n.integration.points) > t(matrix(rep(predictions[[2]]$mean, n.pareto), ncol=n.pareto)) ) )
      }
    }
    
    if (n.obj == 3){
      if (max(predictions[[2]]$sd) != 0) {
        # Regular case
        phi.z.tilde <- pnorm( (matrix(rep(paretoFront[,3],n.integration.points), ncol = n.integration.points) - t(matrix(rep(predictions[[3]]$mean, n.pareto), ncol=n.pareto)) ) /
                                t(matrix(rep(predictions[[3]]$sd, n.pareto), ncol=n.pareto)) )
      } else {
        # z is a fast function
        phi.z.tilde <- 1*( (matrix(rep(paretoFront[,3],n.integration.points), ncol = n.integration.points) > t(matrix(rep(predictions[[3]]$mean, n.pareto), ncol=n.pareto)) ) )
      }
    }
    
    #### 2D CASE ##############################################################################
    pn <- rep(0, n.integration.points)
    if (n.obj == 2){
      # Check if Pareto front is sorted increasingly by its first component
      if (is.unsorted(paretoFront[,1])){
        nondominated.sorted <- sort(paretoFront[,1], index.return=TRUE)[[2]]
        phi.x.tilde <- phi.x.tilde[nondominated.sorted,,drop = FALSE]
        phi.y.tilde <- phi.y.tilde[nondominated.sorted,,drop = FALSE]
      }
      
      pn <- phi.x.tilde[1,]
      if (n.pareto > 1){
        for (j in 2:(n.pareto)){
          pn <- pn + (phi.x.tilde[j,] - phi.x.tilde[j-1,])*phi.y.tilde[j-1,]
        }
      }
      pn <- pn + (1 - phi.x.tilde[n.pareto,])*phi.y.tilde[n.pareto,]
      return(pn)
      
      ##### 3D CASE ###########################################################################
    } else if (n.obj == 3){
      # Check if Pareto front is sorted decreasingly by its third component
      if (is.unsorted(-paretoFront[,3])){
        nondominated.sorted <- sort(paretoFront[,3], decreasing=TRUE, index.return=TRUE)[[2]]
        paretoFront <- paretoFront[nondominated.sorted,,drop=FALSE]
        phi.x.tilde <- phi.x.tilde[nondominated.sorted,]
        phi.y.tilde <- phi.y.tilde[nondominated.sorted,]
        phi.z.tilde <- phi.z.tilde[nondominated.sorted,]
      }
      
      for (i in 1:(n.pareto)){
        if (i ==1){ p.obj.zi <- 1 - phi.z.tilde[i,]
        } else    { p.obj.zi <- phi.z.tilde[i-1,] - phi.z.tilde[i,] }
        
        # nondominated.sub <- which(!is_dominated(t(paretoFront[(i:n.pareto),1:2, drop=FALSE])))
        nondominated.sub <- nonDom(paretoFront[(i:n.pareto),1:2, drop=FALSE], return.idx=TRUE)
        
        #       nondominated.sub <- which(!is_dominated((nondominated_points(t(paretoFront[(i:n.pareto),1:2, drop=FALSE])))))
        # nondominated.sub <- which(!is_dominated(t(nondominated_points(t(paretoFront[(i:n.pareto),1:2, drop=FALSE])))))
        #       nondominated.sub <- paretocheck(paretoFront[(i:n.pareto),1:2, drop=FALSE])                         # Subset of the Pareto front
        nondominated.sub <- i-1+nondominated.sub[sort(paretoFront[(i:n.pareto)[nondominated.sub],1],index.return=TRUE)[[2]]] # Reorder by increasing first objective
        n.pareto.sub     <- length(nondominated.sub)
        
        phi.x.tilde.sub <- phi.x.tilde[nondominated.sub,,drop=FALSE]
        phi.y.tilde.sub <- phi.y.tilde[nondominated.sub,,drop=FALSE]
        
        pi <- phi.x.tilde.sub[1,]
        if (n.pareto.sub > 1){
          for (j in 2:(n.pareto.sub)){
            pi <- pi + (phi.x.tilde.sub[j,] - phi.x.tilde.sub[j-1,])*phi.y.tilde.sub[j-1,]
          }
        }
        pi <- pi + (1 - phi.x.tilde.sub[n.pareto.sub,])*phi.y.tilde.sub[n.pareto.sub,]
        
        pn <- pn + pi*p.obj.zi
      }
      p.obj.zi <- phi.z.tilde[n.pareto,]
      pn <- pn + p.obj.zi
      
      return(pn)
      ###########################################################################################
    } else {
      samps <- NULL
      for (j in 1:n.obj){
        samps <- cbind(samps, rnorm(nsamp*n.integration.points, mean = predictions[[j]]$mean, sd = predictions[[j]]$sd))
      }
      # for(j in 1:nsamp){
      #   pn[i] = pn[i] + !is_dominated(t(rbind(samps[j,], paretoFront)))[1]
      # }
      # pn <- rowSums(matrix(nonDomSet(samps, paretoFront), n.integration.points))
      
      # pn <- rowSums(matrix(nonDom(samps, paretoFront), n.integration.points))
      idx <- nonDom(samps, paretoFront, return.idx=TRUE)
      
      if (length(idx) == 0)  pn <- rep(1e-12, n.integration.points)
      else                   pn <- rowSums(matrix(tabulate(idx, nbins=nsamp*n.integration.points), n.integration.points))
      
      return(pn/nsamp)
    }
  }
}



