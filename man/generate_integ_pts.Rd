% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_integ_pts.R
\name{generate_integ_pts}
\alias{generate_integ_pts}
\title{Strategy generation}
\usage{
generate_integ_pts(
  n.s,
  d,
  nobj,
  x.to.obj = NULL,
  gridtype = "cartesian",
  equilibrium = "NE",
  lb = rep(0, d),
  ub = rep(1, d),
  include.obs = FALSE,
  model = NULL,
  init_set = NULL,
  include_set = NULL,
  seed = 42
)
}
\arguments{
\item{n.s}{scalar or vector. If scalar, total number of strategies (to be divided equally among players),
otherwise number of strategies per player.}

\item{d}{number of variables}

\item{nobj}{number of objectives (or players)}

\item{x.to.obj}{vector allocating variables to objectives. If not provided, default is \code{1:nobj}, assuming that \code{d=nobj}}

\item{gridtype}{either "\code{cartesian}" or "\code{lhs}", or a vector to define a different type for each player.}

\item{equilibrium}{either "\code{NE}", "\code{KSE}", "\code{CKSE}" or "\code{NKSE}"}

\item{lb, ub}{vectors specifying the bounds of the design space, by default \code{[0,1]^d}}

\item{include.obs}{Boolean, if TRUE observations given in \code{model@X} are added to the integration points (only for \code{KSE} and \code{CKSE})}

\item{model}{optional list of \code{km} models (used if \code{include.obs=TRUE})}

\item{init_set}{large grid to subsample from}

\item{include_set}{grid to be included in the larger one generated}

\item{seed}{random seed used by \code{\link[DiceDesign]{lhsDesign}}}
}
\value{
A list containing two matrices, \code{integ.pts} the design of experiments and \code{expanded.indices}
the corresponding indices (for \code{NE}), and the vector \code{n.s}
}
\description{
Preprocessing to link strategies and designs.
}
\examples{
\donttest{
##############################################
### 4 variables, 2 players, no filter
##############################################

# Create a 11x8 lattice based on 2 LHS designs
n.s <- c(11,8)
gridtype = "lhs"
# 4D space is split in 2
x.to.obj <- c(1,1,2,2)
integcontrol <- generate_integ_pts(n.s=n.s, d=4, nobj=2, x.to.obj = x.to.obj, gridtype=gridtype)
pairs(integcontrol$integ.pts)

# Create a simple 11x11 grid
integcontrol <- generate_integ_pts(n.s=11^2, d=2, nobj=2, gridtype="cartesian")
pairs(integcontrol$integ.pts)
}
}
