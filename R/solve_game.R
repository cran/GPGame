##' Main function to solve games.
##' @title Main solver
##' @details If \code{noise.var="given_by_fn"}, \code{fn} returns a list of two vectors, the first being the objective functions and the second
##' the corresponding noise variances.
##'
##' \code{integcontrol} controls the way the design space is discretized. One can directly provide a set of points \code{integ.pts} with
##' corresponding indices \code{expanded.indices} (for \code{NE}). Otherwise, the points are generated according to the number of strategies \code{n.s}.
##' If \code{n.s} is a scalar, it corresponds to the total number of strategies (to be divided equally among players),
##' otherwise it corresponds to the nb of strategies per player. In addition, one may choose the type of discretization with \code{gridtype}.
##' Options are '\code{lhs}' or '\code{cartesian}'. Finally, \code{lb} and \code{ub} are vectors specifying the bounds for the design variables.
##' By default the design space is \code{[0,1]^d}.
##'
##' \code{simucontrol} controls options on conditional GP simulations. Options are \code{IS}: if \code{TRUE}, importance sampling is used for \code{ynew};
##' \code{n.ynew} number of samples of \code{Y(x_{n+1})} and \code{n.sim} number of sample path generated.
##'
##' \code{plotcontrol} can be used to generate plots during the search. Options are \code{plots} (Boolean, \code{FALSE} by default), \code{compute.actual}
##' (Boolean, \code{FALSE} by default, to draw the actual problem, only for inexpensive \code{fun}), and \code{pbname} (string, for figure title and pdf export).
##'
##' \code{filtercontrol} controls filtering options. \code{filter} sets how to select a subset of simulation and candidate points,
##' either either a single value or a vector of two to use different filters for simulation and candidate points.
##' Possible values are '\code{window}', '\code{Pnash}' (for \code{NE}), '\code{PND}' (probability of non domination), '\code{none}'.
##' \code{nsimPoints} and \code{ncandPoints} set the maximum number of simulation/candidate points wanted
##' (use with filter '\code{Pnash}' for now). Default values are \code{800} and \code{200}, resp.
##' \code{randomFilter} (\code{TRUE} by default) sets whereas the filter acts randomly or deterministically.
##'
##' \code{kmcontrol} Options for handling \code{nobj} \code{\link[DiceKriging]{km}} models.
##' \code{cov.reestim} (Boolean, \code{TRUE} by default) specifies if the kriging hyperparameters
##' should be re-estimated at each iteration,
##'
##' \code{returncontrol} sets options for the last iterations and what is returned by the algorithm.
##' \code{return.Eq} (Boolean, \code{TRUE} by default) specifies
##' if a final search for the equilibrium is performed at the end. \code{finalcrit} sets a different criterion for the last iteration.
##' \code{track.Eq} allows to estimate the equilibrium at each iteration; options are '\code{none}' to do nothing,
##' "\code{mean}" (default) to compute the equilibrium of the prediction mean (all candidates),
##'  "\code{empirical}" (for \code{KSE}) and "\code{pex}"/"\code{psim}" (\code{NE} only)
##' for using \code{Pnash} estimate (along with mean estimate, on integ.pts only, NOT reestimated if \code{filter.simu} or \code{crit} is \code{Pnash}).
##'
##' @param fun fonction with vectorial output
##' @param equilibrium either '\code{NE}', '\code{KSE}' or '\code{NKSE}' for Nash/Kalai-Smoridinsky/Nash-Kalai-Smoridinsky equilibria
##' @param crit '\code{sur}' (default) is available for all equilibria, '\code{psim}' and '\code{pex}' are available for Nash
##' @param model list of \code{\link[DiceKriging]{km}} models
##' @param n.init number of points of the initial design of experiments if no model is given
##' @param n.ite number of iterations of sequential optimization
##' @param d variable dimension
##' @param nobj number of objectives (players)
##' @param x.to.obj for \code{NE} and \code{NKSE}, which variables for which objective
##' @param noise.var noise variance. Either a scalar (same noise for all objectives), a vector (constant noise, different for each objective),
##' a function (type closure) with vectorial output (variable noise, different for each objective) or \code{"given_by_fn"}, see Details.
##' If not provided, \code{noise.var} is taken as the average of \code{model@noise.var}.
##' @param integcontrol optional list for handling integration points. See Details.
##' @param simucontrol optional list for handling conditional simulations. See Details.
##' @param plotcontrol optional list for handling during-optimization plots. See Details.
##' @param filtercontrol optional list for handling filters. See Details.
##' @param kmcontrol optional list for handling \code{\link[DiceKriging]{km}} models. See Details.
##' @param returncontrol optional list for choosing return options. See Details.
##' @param ncores number of CPU available (> 1 makes mean parallel \code{TRUE})
##' @param trace controls the level of printing: \code{0} (no printing), \code{1} (minimal printing), \code{3} (detailed printing)
##' @param seed to fix the random variable generator
##' @param ... additional parameter to be passed to \code{fun}
##' @return
##' A list with components:
##' \itemize{
##' \item{\code{model}}{: a list of objects of class \code{\link[DiceKriging]{km}} corresponding to the last kriging models fitted.}
##' \item{\code{Jplus}}{: recorded values of the acquisition function maximizer}
##' \item{\code{integ.pts} and  \code{expanded.indices}}{: the discrete space used,}
##' \item{\code{predEq}}{: a list containing the recorded values of the estimated best solution,}
##' \item{\code{Eq.design, Eq.poff}}{: estimated equilibrium and corresponding pay-off (if \code{return.Eq==TRUE})}
##' }
##'
##' @export
##' @importFrom grDevices dev.off pdf rainbow
##' @importFrom graphics axis pairs par points title
##' @importFrom stats pnorm qnorm rnorm dnorm
##' @importFrom methods slot
##' @importFrom graphics filled.contour
##' @import DiceKriging DiceDesign parallel
##' @importFrom MASS mvrnorm
##' @importFrom grDevices terrain.colors
##' @importFrom graphics legend
##' @useDynLib GPGame, .registration = TRUE
##' @references
##' V. Picheny, M. Binois, A. Habbal (2016+), A Bayesian optimization approach to find Nash equilibria,
##' \emph{https://arxiv.org/abs/1611.02440}.
##' @examples
##' \dontrun{
##'
##' ##############################################
##' # Example 1: 2 variables, 2 players, no filter
##' ##############################################
##' # Define objective function (R^2 -> R^2)
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
##' # To use parallel computation (turn off on Windows)
##' library(parallel)
##' parallel <- FALSE #TRUE # 
##' if(parallel) ncores <- detectCores() else ncores <- 1
##'
##' # Simple configuration: no filter, discretization is a 21x21 grid
##'
##' # Grid definition
##' n.s <- rep(21, 2)
##' x.to.obj   <- c(1,2)
##' gridtype <- 'cartesian'
##'
##' # Run solver with 6 initial points, 4 iterations
##' # Increase n.ite to at least 10 for better results
##' res <- solve_game(fun, equilibrium = "NE", crit = "sur", n.init=6, n.ite=4,
##'                   d = 2, nobj=2, x.to.obj = x.to.obj,
##'                   integcontrol=list(n.s=n.s, gridtype=gridtype),
##'                   ncores = ncores, trace=1, seed=1)
##'
##' # Get estimated equilibrium and corresponding pay-off
##' NE <- res$Eq.design
##' Poff <- res$Eq.poff
##'
##' # Draw results
##' plotGame(res)
##'
##' ##############################################
##' # Example 2: 4 variables, 2 players, filtering
##' ##############################################
##' fun <- function(x, nobj = 2){
##'   if (is.null(dim(x)))     x <- matrix(x, 1)
##'   y <- matrix(x[, 1:(nobj - 1)], nrow(x))
##'   z <- matrix(x[, nobj:ncol(x)], nrow(x))
##'   g <- rowSums((z - 0.5)^2)
##'   tmp <- t(apply(cos(y * pi/2), 1, cumprod))
##'   tmp <- cbind(t(apply(tmp, 1, rev)), 1)
##'   tmp2 <- cbind(1, t(apply(sin(y * pi/2), 1, rev)))
##'   return(tmp * tmp2 * (1 + g))
##' }
##'
##' # Grid definition: player 1 plays x1 and x2, player 2 x3 and x4
##' # The grid is a lattice made of two LHS designs of different sizes
##' n.s <- c(44, 43)
##' x.to.obj   <- c(1,1,2,2)
##' gridtype <- 'lhs'
##'
##' # Set filtercontrol: window filter applied for integration and candidate points
##' # 500 simulation and 200 candidate points are retained.
##' filtercontrol <- list(nsimPoints=500, ncandPoints=200,
##'                    filter=c("window", "window"))
##'
##' # Set km control: lower bound is specified for the covariance range
##' # Covariance type and model trend are specified
##' kmcontrol <- list(lb=rep(.2,4), model.trend=~1, covtype="matern3_2")
##'
##' # Run solver with 20 initial points, 4 iterations
##' # Increase n.ite to at least 20 for better results
##' res <- solve_game(fun, equilibrium = "NE", crit = "psim", n.init=20, n.ite=2,
##'                   d = 4, nobj=2, x.to.obj = x.to.obj,
##'                   integcontrol=list(n.s=n.s, gridtype=gridtype),
##'                   filtercontrol=filtercontrol,
##'                   kmcontrol=kmcontrol,
##'                   ncores = 1, trace=1, seed=1)
##'
##' # Get estimated equilibrium and corresponding pay-off
##' NE <- res$Eq.design
##' Poff <- res$Eq.poff
##'
##' # Draw results
##' plotGame(res)
##'
##' }

solve_game <- function(
  fun, ..., equilibrium="NE", crit="sur", model=NULL, n.init=NULL, n.ite, d, nobj, x.to.obj=NULL, noise.var = NULL,
  integcontrol=NULL, simucontrol=NULL, plotcontrol=NULL, filtercontrol=NULL, kmcontrol=NULL, returncontrol=NULL,
  ncores=1, trace=1, seed=NULL) {

  t1 <- Sys.time()

  ### Initialize variables
  set.seed(seed)
  if(ncores > 1) parallel <- TRUE

  # Check critcontrols
  if(crit == 'pex') PnashMethod <- 'exact'
  else PnashMethod <- 'simu'

  # Check integcontrol
  if (!is.null(integcontrol$integ.pts)) integ.pts <- integcontrol$integ.pts else integ.pts <- NULL
  if (!is.null(integcontrol$expanded.indices)) expanded.indices <- integcontrol$expanded.indices else expanded.indices <- NULL
  if (!is.null(integcontrol$n.s)) n.s <- integcontrol$n.s else n.s <- NULL
  if (!is.null(integcontrol$gridtype)) gridtype <- integcontrol$gridtype else gridtype <- "lhs"
  if (!is.null(integcontrol$lb)) lb <- integcontrol$lb else lb <- rep(0, d)
  if (!is.null(integcontrol$ub)) ub <- integcontrol$ub else ub <- rep(1, d)

  # Check simucontrol
  if (!is.null(simucontrol$n.ynew)) n.ynew <- simucontrol$n.ynew else n.ynew <- 10
  if (!is.null(simucontrol$n.sim)) n.sim <- simucontrol$n.sim else n.sim <- 10
  if (!is.null(simucontrol$IS)) IS <- simucontrol$IS else IS <- TRUE
  cross <- FALSE

  # Check plotcontrol
  if (!is.null(plotcontrol$compute.actual)) compute.actual <- plotcontrol$compute.actual else compute.actual <- FALSE
  if (!is.null(plotcontrol$plots)) plots <- plotcontrol$plots else plots <- FALSE
  if (!is.null(plotcontrol$pbname)) pbname <- plotcontrol$pbname else pbname <- NULL
  if (!is.null(pbname)) exportPDF <- TRUE else exportPDF <- FALSE

  # Check filtercontrol
  if (!is.null(filtercontrol$filter)) filter <- filtercontrol$filter else  filter <- 'none'
  if (!is.null(filtercontrol$nsimPoints)) nsimPoints <- filtercontrol$nsimPoints else nsimPoints <- 800
  if (!is.null(filtercontrol$ncandPoints)) ncandPoints <- filtercontrol$ncandPoints else ncandPoints <- 200
  if (!is.null(filtercontrol$randomFilter)) randomFilter <- filtercontrol$randomFilter else randomFilter <- TRUE

  if (length(filter)==1) {
    simu.filter <- cand.filter <- filter
  } else {
    simu.filter <- filter[1]
    cand.filter <- filter[2]
  }

  if (equilibrium=="KSE") {
    if ("Pnash" %in% c(simu.filter, cand.filter)) cat("Pnash filter only available for NE; switching to PND \n")
    if (simu.filter == "Pnash") simu.filter <- 'PND'
    if (cand.filter == "Pnash") cand.filter <- 'PND'
  }

  if (length(randomFilter)==1) randomFilter <- rep(randomFilter, 2)

  # Check kmcontrol
  if (!is.null(kmcontrol$cov.reestim)) plots <- kmcontrol$cov.reestim else cov.reestim <- TRUE
  if (!is.null(kmcontrol$model.trend)) model.trend <- kmcontrol$model.trend else model.trend <- ~1
  if (!is.null(kmcontrol$lb)) kmlb <- kmcontrol$lb else kmlb <- rep(.1,d)
  if (!is.null(kmcontrol$ub)) kmub <- kmcontrol$ub else kmub <- rep(1,d)
  if (!is.null(kmcontrol$nugget)) kmnugget <- kmcontrol$nugget else kmnugget <- 1e-8
  if (!is.null(kmcontrol$control)) control <- kmcontrol$control else control <- list(trace=FALSE)
  if (!is.null(kmcontrol$covtype)) covtype <- kmcontrol$covtype else covtype <- "matern5_2"

  # Check returncontrol
  if (!is.null(returncontrol$return.Eq)) return.Eq <- returncontrol$return.Eq else return.Eq <- TRUE
  if (!is.null(returncontrol$finalcrit)) finalcrit <- returncontrol$finalcrit else finalcrit <- crit
  if (!is.null(returncontrol$track.Eq)) track.Eq <- returncontrol$track.Eq else track.Eq <- 'mean'
  if (equilibrium=="KSE") {
    if (!is.null(track.Eq)) {
      if (track.Eq %in% c("psim", "pex")) {
        cat("track.Eq must be set to none or mean for KSE - switching to mean \n")
        track.Eq <- "mean"
      }
    }
    finalcrit <- "sur"
  } else {
    if (track.Eq == "empirical") {
      cat("track.Eq= empirical option only available for KSE; switching to psim \n")
      track.Eq <- "psim"
    }
  }

  # Check Noise
  if (!is.null(noise.var)) {
    if (typeof(noise.var) == "closure") noise.var <- match.fun(noise.var)
    else if (typeof(noise.var) == "double" && length(noise.var==1)) noise.var <- rep(noise.var, nobj)
    kmnugget <- NULL
  }

  ### Initialize variables
  window <- NULL
  all.Jplus <- c()
  include.obs <- FALSE
  Eq.design.estimate <- Eq.poff.estimate <- trackPnash <- NULL
  predEq <- list()

  if (crit!="sur" && equilibrium!="NE"){
    cat("pex and psim available only for Nash equilibria; crit switched to sur \n")
    crit <- "sur"
  }

  if (trace>0) cat("--------------------------\n Starting", equilibrium, "search with:", "\n",
                   crit, "strategy,","\n",
                   simu.filter, "simulation point filter,", "\n",
                   cand.filter, "candidate point filter,", "\n",
                   "among (", n.s, ") strategies \n --------------------------\n " )

  #### Initial design and models ####################
  if (is.null(model)){
    if (is.null(n.init))     n.init <- 5*d
    design <- lhsDesign(n.init, d, seed = 42)$design
    response <- t(apply(design, 1, fun, ... = ...))

    if (!is.null(noise.var)) {
      if (typeof(noise.var) == "closure") {
        newnoise.var <- apply(design, 1, noise.var, ...)
      } else if (typeof(noise.var) == "double") {
        newnoise.var <- matrix(noise.var, nrow=n.init, ncol=ncol(response), byrow=TRUE)
      } else {#noise.var ="given_by_fn"
        # newnoise.var <- response[[2]]
        # ynew <- response[[1]]
        tmp <- newnoise.var <- NULL # initialization
        for (i in 1:length(response)) {
          tmp <- rbind(tmp, response[[i]][[1]])
          newnoise.var <- rbind(newnoise.var, response[[i]][[2]])
        }
        response <- tmp
      }
    } else {
      newnoise.var <- NULL
    }
    nobj <- ncol(response)

    my.km <- function(i) {
      km(model.trend, design = design, response = response[, i], covtype=covtype,
         control=control, lower=kmlb, upper=kmub, nugget=kmnugget, noise.var=newnoise.var[,i])
    }
    model <- mclapply(1:nobj, my.km, mc.cores=ncores)
  }

  #### Integration points ###########################
  if (is.null(integcontrol$integ.pts) || is.null(integcontrol$expanded.indices)) {
    res <- generate_integ_pts(n.s=n.s, d=d, nobj=nobj, x.to.obj=x.to.obj, gridtype=gridtype, lb = lb, ub = ub)
    integcontrol$integ.pts <- integ.pts <- res$integ.pts
    integcontrol$expanded.indices <- expanded.indices <- res$expanded.indices
  }
  n.integ.pts <- nrow(integ.pts)
  cand.pts <- integ.pts

  if (is.null(expanded.indices)) {
    sorted <- FALSE
  } else {
    sorted <- !is.unsorted(expanded.indices[,nobj])
  }

  if (compute.actual) {
    if (plots) {
      if (nobj==2) {
        if(exportPDF) pdf(paste0(pbname, "_obj.pdf"), width=5, height=5)
        par(mfrow=c(1,1), mar=c(3,3,2,1), mgp=c(2,1,0))
        groundTruth <- plotGameGrid(fun = fun, domain = matrix(rep(c(0,1), each = d), d), graphs = "objective",
                                    n.grid = 31, x.to.obj = x.to.obj, integcontrol=integcontrol, equilibrium = equilibrium, ... = ...)
        if(exportPDF) dev.off()
      }

      if(exportPDF) pdf(paste0(pbname, "_design.pdf"), width=5, height=5)
      par(mfrow=c(1,1), mar=c(3,3,2,1), mgp=c(2,1,0))
      plotGameGrid(fun = fun, domain = matrix(rep(c(0,1), each = d), d), graphs = "design", n.grid = 31, x.to.obj = x.to.obj,
                   integcontrol=my.integcontrol, type = equilibrium, ... = ...)
      if(exportPDF) dev.off()
    }
  }

  if (trace>2) cat("Time for initialization: \n")
  if (trace>2) print(Sys.time() - t1)
  t1 <- Sys.time()

  #############################################################
  #### MAIN LOOP STARTS HERE ##################################
  #############################################################
  if (plots) par(mfrow=c(min(n.ite+1, 3),2), mar=c(2,2,2,2))

  Jnres <- rep(NA, n.ite + 1)

  for (ii in 1:(n.ite+1)){
    if (ii < (n.ite+1) && trace>0){cat("--- Iteration #", ii, " ---\n")}
    t0 <- t1 <- Sys.time()

    ##################################################################
    # MODELS UPDATE
    ##################################################################
    if (ii > 1) {

      if (ii == n.ite){
        crit <- finalcrit
        cand.filter <- 'none'
      }

      newmodel <- model
      X.new <- matrix(xnew,nrow=1)

      my.update <- function(u) {
        try(update(object = model[[u]], newX = X.new, newy=ynew[u], newX.alreadyExist=FALSE, newnoise.var = newnoise.var[u],
                   cov.reestim = cov.reestim, kmcontrol = list(control = list(trace = FALSE))), silent = TRUE)
      }
      newmodel <- mclapply(1:nobj, my.update, mc.cores=ncores)

      for (u in 1:nobj){
        if (typeof(newmodel[[u]]) == "character" && cov.reestim) {
          cat("Error in hyperparameter estimation - old hyperparameter values used instead for model ", u, "\n")
          newmodel[[u]] <- try(update(object = model[[u]], newX = X.new, newy=ynew[u], newnoise.var = newnoise.var[u],
                                      newX.alreadyExist=FALSE, cov.reestim = FALSE), silent = TRUE)
        }
        if (typeof(newmodel[[u]]) == "character") {
          cat("Unable to udpate kriging model ", u, " at iteration", ii-1, "- optimization stopped \n")
          cat("lastmodel ", u, " is the model at iteration", ii-1, "\n")
          cat("par and values contain the ",ii , "th observation \n \n")
          # if (ii > 1) allX.new <- rbind(model[[u]]@X[(model[[u]]@n + 1):(model[[u]]@n+ii-1),, drop=FALSE], X.new)
          return(list(
            par    = X.new,
            values = ynew,
            nsteps = ii,
            model = model,
            Jplus = Jplus, integcontrol=integcontrol))
        } else {
          model[[u]] <- newmodel[[u]]
        }
      }
    }

    if (trace>2) cat("Time for models updates: \n")
    if (trace>2) print(Sys.time() - t1)
    t1 <- Sys.time()
    observations <- Reduce(cbind, lapply(model, slot, "y"))

    #--------------------------------------------------------#
    # Regular loop : find infill point
    #--------------------------------------------------------#
    if (ii < (n.ite+1) || return.Eq){
      if (ii==(n.ite+1)) include.obs <- TRUE

      pred <- mclapply(model, FUN=predict, newdata = integ.pts, checkNames = FALSE, type = "UK", light.return = TRUE, mc.cores=ncores)

      if (trace>2) cat("Time for predictions: \n")
      if (trace>2) print(Sys.time() - t1)
      t1 <- Sys.time()

      ##################################################################
      # FILTERING INTEGRATION POINTS
      ##################################################################
      if (simu.filter == "window") {

        if (is.null(window) || crit != "sur") {
          # Special case for initialization - otherwise do not change the old "window" value
          predmean <- Reduce(cbind, lapply(pred, function(alist) alist$mean))
          Eq_simu <- getEquilibrium(Z = predmean,  equilibrium = equilibrium, nobj = nobj, expanded.indices=expanded.indices,
                                    n.s=n.s, sorted = sorted, cross = cross)
          if(nrow(Eq_simu) < 1 || all(is.na(Eq_simu))){
            ## if no window has been defined in previous iterations, use the mean, otherwise keep the old one
            if(is.null(window)){
              window <- colMeans(observations)
            }
          } else if(nrow(Eq_simu) == 1) {
            window <- as.vector(Eq_simu)
          } else {
            window <- colMeans(Eq_simu, na.rm = TRUE)
          }
        }
        options <- list(window = window)
      } else if (simu.filter == "Pnash") {
        options <- list(method=PnashMethod, nsim=100)
      }
      if (trace > 1 && simu.filter == "window") cat("window for integ.pts", window, "\n")

      # Reduce the number of integration points first
      if (simu.filter != 'none') {
        filt <- filter_for_Game(n.s.target = pmin(n.s, round(nsimPoints^(1/nobj))), integcontrol=integcontrol,
                                model = model, predictions=pred, type=simu.filter, options = options, random=randomFilter[1],
                                include.obs=include.obs, ncores = ncores, equilibrium=equilibrium)
        I <- filt$I
      } else {
        I <- 1:nrow(integ.pts)
      }
      my.integ.pts <- integ.pts[I,]
      my.expanded.indices <- expanded.indices[I,]
      my.pred <- lapply(pred, function(alist) list(mean=alist$mean[I], sd=alist$sd[I]))

      my.n.s <- apply(my.expanded.indices, 2, function(x) length(unique(x)))
      if (trace>1) cat("Number of strategies after simu filtering: (", my.n.s, ")\n")

      my.expanded.indices2 <- my.expanded.indices
      for (i in 1:ncol(my.expanded.indices)) {
        unik_i <- unique(my.expanded.indices[,i])
        for (j in 1:length(unik_i)) {
          my.expanded.indices2[which(my.expanded.indices[,i]==unik_i[j]),i] <- j
        }
      }
      my.expanded.indices <- my.expanded.indices2

      my.integcontrol <- list(expanded.indices=my.expanded.indices, integ.pts=my.integ.pts, n.s=my.n.s)

      if (trace>2) cat("Time for filter 1 (simu): \n")
      if (trace>2) print(Sys.time() - t1)
      t1 <- Sys.time()

      ##################################################################
      # Precalculations and simulations (if not performed already)
      ##################################################################

      precalc.data <- mclapply(model, FUN=precomputeUpdateData, integration.points=my.integ.pts, mc.cores=ncores)

      if(crit == 'sur' || cand.filter=="window"){
        my.Simu <- t(Reduce(rbind, mclapply(model, simulate, nsim=n.sim, newdata=my.integ.pts, cond=TRUE,
                                            checkNames=FALSE, nugget.sim = 10^-8, mc.cores=ncores)))
      }

      if (trace>2) cat("Time for precalc + simu: \n")
      if (trace>2) print(Sys.time() - t1)
      t1 <- Sys.time()

      ##################################################################
      # FILTERING CANDIDATE POINTS
      ##################################################################
      cand.pts <- my.integ.pts

      if ((simu.filter == "window" && crit == 'sur') || cand.filter == "window" || plots) {
        Eq_simu <- getEquilibrium(my.Simu, equilibrium = equilibrium, nobj = nobj, expanded.indices=my.expanded.indices,
                                  n.s=my.n.s, sorted = sorted, cross = cross)
        Eq_simu <- Eq_simu[which(!is.na(Eq_simu[,1])),, drop = FALSE]
        Jnres[ii] <- det(cov(Eq_simu))
        temp <- apply(Eq_simu, 2, range, na.rm = TRUE)

        if (!any(is.na(temp))) window <- temp
        options <- list(window = window)
      } else if (cand.filter == "Pnash") {
        options <- list(method=PnashMethod, nsim = 100)
      }

      # Plots
      if (nobj==2 && plots) {
        if(exportPDF) pdf(paste0(pbname, "_objspace_ite_", ii, ".pdf"), width=5, height=5)
        par(mfrow=c(1,1), mar=c(3,3,2,1), mgp=c(2,1,0))

        if(compute.actual){
          xrange <- range(groundTruth$response.grid[,1])
          yrange <- range(groundTruth$response.grid[,2])
          xrange[1] <- xrange[1] - .1*(xrange[2] - xrange[1])
          yrange[1] <- yrange[1] - .1*(yrange[2] - yrange[1])
          plot(NA, xlim=xrange, ylim=yrange,
               main=paste0("Iteration ", ii), xlab=expression(f[1]), ylab=expression(f[2]))
          points(groundTruth$trueEqPoff[1], groundTruth$trueEqPoff[2], pch = 24, bg = "green", cex=2)
          # plotParetoEmp(trueParetoFront , col="green")
          points(observations[,1], observations[,2], col="blue", pch=19, cex=1.2)
          # plotParetoEmp(t(nondominated_points(t(observations))), col="blue")
          points(Eq_simu[,1], Eq_simu[,2], col="red", pch=3, lwd=2)
          legend("topright", legend = c("Actual Eq", "Simulated Eq", "Current obs", "New obs"), lwd = rep(NA,4),
                 col = c("black", "red", "blue", "red"), pch=c(24, 3, 19, 19), pt.cex= c(2, 1.2, 1, 2),
                 pt.lwd=c(1,2,1,1), pt.bg=c("green", NA, NA, NA))
        }else{
          UQ_eq <- TRUE
          if(crit == "sur") UQ_eq <- FALSE
          plotGame(res = list(model = model, Jplus = all.Jplus, integcontrol=integcontrol,
                              predEq = predEq, trackPnash = trackPnash),
                   equilibrium = equilibrium, add = FALSE, UQ_eq = UQ_eq, simucontrol = simucontrol)
          title(paste0("Iteration ", ii))
          if(crit == "sur") points(Eq_simu[,1], Eq_simu[,2], col="red", pch=3, lwd=2)
        }

        if (ii > n.ite) if(exportPDF) dev.off()
      }


      if (trace>1 && cand.filter == "window") cat("window for cand.pts", window, "\n")

      if (cand.filter != 'none') {
        if (cand.filter == "Pnash" && simu.filter == "Pnash") {
          # Special case if Pnash used twice: do not redo calculations
          if (randomFilter[2]) {
            J <- sample.int(length(I), size=min(length(I), ncandPoints), replace=FALSE, prob=filt$crit[I])
          } else {
            J <- I[1:min(length(I), ncandPoints)]
          }
        } else {
          # Regular case
          filt2 <- filter_for_Game(n.s.target = pmin(my.n.s, round(ncandPoints^(1/nobj))), integcontrol=my.integcontrol,
                                   model = model, predictions=my.pred, type=cand.filter, options = options,
                                   random=randomFilter[2], include.obs=include.obs, ncores = ncores, equilibrium=equilibrium)
          J <- filt2$I
        }
      } else {
        # No filter case
        J <- 1:nrow(cand.pts)
      }



      ## Remove already evaluated designs from candidates
      if((ii > 1) && (ii<(n.ite+1))){
        my.pred.s <- my.pred[[1]]$sd[J]
        J <- J[which(my.pred.s >= sqrt(model[[1]]@covariance@sd2)/1e6)]
      }

      my.cand.pts <- cand.pts[J,, drop=FALSE]
      cand.s <- my.expanded.indices[J,, drop=FALSE]
      # my.pred.s <- my.pred[[1]]$sd[J]

      if (trace>2) cat("Time for filter 2 (cand): \n")
      if (trace>2) print(Sys.time() - t1)
      t1 <- Sys.time()

      if (trace>0) cat("Nb of integration / candidate pts:", length(I), "/", length(J), "\n")

      ##################################################################
      # CRITERION MINIMIZATION
      ##################################################################
      my.integ.pts <- data.frame(my.integ.pts)
      if (crit=="sur") {
        Jplus <- unlist(mclapply(X=J, FUN=crit_SUR_Eq, model=model, integcontrol=my.integcontrol,
                                 Simu=my.Simu, equilibrium = equilibrium, precalc.data=precalc.data,
                                 n.ynew=n.ynew, IS=IS, cross=cross, mc.cores = ncores))
      } else if (crit == 'pex') {
        Jplus <- -crit_PNash(idx=J, integcontrol=my.integcontrol,
                            type = 'exact', model = model, ncores = ncores)
      } else if (crit=="psim") {
        Jplus <- -crit_PNash(idx=J, integcontrol=my.integcontrol,
                            type = 'simu', model = model, ncores = ncores, control = list(nsim = 100))
      } else {
        cat("wrong crit \n")
        break;
      }
      if(all(is.na(Jplus))){
        cat("No equilibrium for this problem \n")
        return(list(model = model, Jplus = Jplus, integcontrol=integcontrol))
      }
      Jplus[is.na(Jplus)] <- max(Jplus, na.rm = TRUE)

      if (trace>2) cat("Time for optim: \n")
      if (trace>2) print(Sys.time() - t1)

      if (trace > 0 && crit=="sur") cat("Jplus:", min(Jplus), '\n')
      if (trace > 0 && crit!="sur") cat("Pmax:", -min(Jplus), '\n')

      if (ii < n.ite+1) {
        # Get best solution
        i <- which.min(Jplus)
        xnew <- my.cand.pts[i,]
        ynew <- fun(xnew, ...)

        if (!is.null(noise.var)) {
          if (typeof(noise.var) == "closure") {
            newnoise.var <- noise.var(xnew)
          } else if (typeof(noise.var) == "double") {
            newnoise.var <- noise.var
          } else {#noise.var ="given_by_fn"
            newnoise.var <- ynew[[2]]
            ynew <- ynew[[1]]
          }
        } else {
          newnoise.var <- NULL
        }

        all.Jplus <- c(Jplus, min(Jplus))

        if (trace>1) cat("New design evaluated: \n")
        if (trace>1) cat(xnew, "\n")
        if (trace>1) cat("Corresponding new observation: \n")
        if (trace>1) cat(ynew, "\n")

        if (trace>0) cat("Total iteration time: \n")
        if (trace>0) print(Sys.time() - t0)
        ##################################################################
        # PLOTS
        ##################################################################
        if (nobj==2 && plots) {
          # Representation du nouveau point
          points(ynew[1], ynew[2], col="red", pch=19, cex=2)
          if(exportPDF) dev.off()
        }
        if (d==2 && plots) {

          # Representation du critere
          if (gridtype =="cartesian") {

            Jgrid <- rep(max(Jplus), prod(n.s))
            Jgrid[J] <- Jplus

            if(exportPDF) pdf(paste0(pbname, "_crit_ite_", ii, ".pdf"), width=5, height=5)
            par(mfrow=c(1,1), mar=c(3,3,2,1), mgp=c(2,1,0))

            # TO CHECK: that integ.pts corresponds to the structure needed in filled.contour
            filled.contour(seq(0, 1, length.out = n.s[1]), seq(0, 1, length.out = n.s[2]), nlevels = 20,
                           matrix(Jgrid, n.s[1], n.s[2]), main = paste0("Iteration ", ii),
                           xlab = expression(x[1]), ylab = expression(x[2]),
                           plot.axes = {axis(1); axis(2);
                             points(xnew[1], xnew[2], pch = 21, bg = "blue", cex=2);
                             points(model[[1]]@X[,1], model[[1]]@X[,2], pch = 21, bg = "white");
                           }
            )
            if(compute.actual) points(groundTruth$trueEqdesign[1], groundTruth$trueEqdesign[2], pch = 24, bg = "green", cex=1.5)
            if(exportPDF) dev.off()
          }
        }
      }
      # else {
      #   #--------------------------------------------------------#
      #   # If optimization is over: find best Eq estimate
      #   #--------------------------------------------------------#
      #   Eq.design.estimate <- cand.pts
      #   Eq.poff.estimate <- -Jplus
      # }

      ##################################################################
      # Tracking of estimated best solution
      ##################################################################

      if (is.null(track.Eq)) track.Eq <- "psim"

      if (track.Eq != 'none' || (ii == n.ite+1 && return.Eq)){

        if (track.Eq == "empirical") {
          observations <- Reduce(cbind, lapply(model, slot, "y"))
          currentEq <- getEquilibrium(Z=observations, equilibrium = equilibrium, nobj=nobj, return.design=TRUE, cross=cross)
          # if (ii == 1) {
          #   predEq <- list(Eq_emp)
          # } else {
          #   predEq <- c(predEq, list(Eq_emp))
          # }

        } else if (track.Eq == "psim" || track.Eq == "pex") {
          if (ii == 1) trackPnash <- list()

          if (simu.filter == "Pnash") {
            ImaxPnash <- which.max(filt$crit)
            maxPnash <- max(filt$crit)
            # trackPnash <- c(trackPnash, list(list(c(predEqPoff = unlist(lapply(my.pred, function(x) x$mean[ImaxPnash])),
            # Eq = my.integ.pts[ImaxPnash,, drop = F]), Pnash = maxPnash)))
            currentEq <- list(c(predEqPoff = unlist(lapply(my.pred, function(x) x$mean[ImaxPnash])),
                                Eq = my.integ.pts[ImaxPnash,, drop = F]), Pnash = maxPnash)
          } else {
            if (track.Eq == 'pex') {
              estPnash <- crit_PNash(idx = 1:nrow(my.integcontrol$expanded.indices), integcontrol=my.integcontrol,
                                     type = 'exact', model = model, ncores = ncores)
            } else {
              estPnash <- crit_PNash(idx = 1:nrow(my.integcontrol$expanded.indices), integcontrol=my.integcontrol,
                                     type = 'simu', model = model, ncores = ncores,
                                    control = list(nsim = 100))
            }
            ImaxPnash <- which.max(estPnash)
            maxPnash <- max(estPnash)
            currentEq <- list(predEqPoff = unlist(lapply(my.pred, function(x) x$mean[ImaxPnash])),
                              Eq = my.integ.pts[ImaxPnash,, drop = F], Pnash = maxPnash)
            # trackPnash <- c(trackPnash, list(list(predEqPoff = unlist(lapply(my.pred, function(x) x$mean[ImaxPnash])),
            #                                       Eq = my.integ.pts[ImaxPnash,, drop = F], Pnash = maxPnash)))
          }

          if (trace>2) cat("Maximum estimated Pnash:", maxPnash," \n")
          if (trace>2) print(Sys.time() - t1)
        } else {
          ## Eq of the GP predictive means
          predmean <- Reduce(cbind, lapply(pred, function(alist) alist$mean))
          currentEq <- getEquilibrium(Z = predmean,  equilibrium = equilibrium, nobj = nobj,
                                      expanded.indices=expanded.indices, n.s=n.s,
                                      sorted = sorted, cross = cross, return.design = T)
        }
        predEq <- c(predEq, list(currentEq))
      }

      if (trace>2) cat("Time for estimating best candidate: \n")
      if (trace>2) print(Sys.time() - t1)
      t1 <- Sys.time()
    }
  }
  if (return.Eq) {
    #--------------------------------------------------------#
    # When optimization is over: find best Eq estimate
    #--------------------------------------------------------#
    Eq.design.estimate <- integ.pts[predEq[[length(predEq)]]$NE,]
    Eq.poff.estimate <- predEq[[length(predEq)]]$NEPoff
    return(list(model = model, Eq.design=Eq.design.estimate, Eq.poff=Eq.poff.estimate,
                Jplus = all.Jplus, integcontrol=integcontrol,
                predEq = predEq))
  } else {
    return(list(model = model, Jplus = all.Jplus, integcontrol=integcontrol,
                predEq = predEq))
  }
}
