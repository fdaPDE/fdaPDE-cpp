#' FELSPLINE test function
#' 
#' @param x,y Points at which to evaluate the test function.
#' @param r0 The test domain is a sort of bent sausage. This is the radius of the inner bend.
#' @param r The radius of the curve at the centre of the sausage.
#' @param l The length of an arm of the sausage.
#' @param b The rate at which the function increases per unit increase in distance along the centre line of the sausage.
#' @description Implements a finite area test function based on one proposed by Tim Ramsay (2002) proposed by 
#' Simon Wood (2008).
#' @return Returns function evaluations, or NAs for points outside the horseshoe domain. 
#' @usage fs.test(x, y, r0 = 0.1, r = 0.5, l = 3, b = 1)  
#' @references 
#' \itemize{
#'    \item{Ramsay, T. 2002. Spline smoothing over difficult regions. J.R.Statist. Soc. B 64(2):307-319}
#'    \item{Wood, S. N., Bravington, M. V., & Hedley, S. L. (2008). Soap film smoothing. Journal of the Royal 
#' Statistical Society: Series B (Statistical Methodology), 70(5), 931-955.}
#' }
#' @examples 
#' library(fdaPDE)
#' 
#' ## Upload the horseshoe2D data
#' data(horseshoe2D)
#' boundary_nodes = horseshoe2D$boundary_nodes
#' boundary_segments = horseshoe2D$boundary_segments
#' locations = horseshoe2D$locations
#' 
#' ## Create the 2D mesh
#' mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)
#' ## Create the FEM basis
#' FEMbasis = create.FEM.basis(mesh)
#' ## Compute the coeff vector evaluating the desired function at the mesh nodes
#' ## In this case we consider the fs.test() function introduced by Wood et al. 2008
#' coeff = fs.test(mesh$nodes[,1], mesh$nodes[,2])
#' ## Create the FEM object
#' FEMfunction = FEM(coeff, FEMbasis)
#' ## Plot it
#' plot(FEMfunction)
#' @export

fs.test <- function(x,y,r0=.1,r=.5,l=3,b=1)
  ## test function based on Tim Ramsay (2002) J.R.Statist. Soc. B
  ## 64(2):307-319 "Spline smoothing over difficult regions"
  { q <- pi*r/2 ## 1/2 length of semi-circle part of centre curve
  a <- d <- x*0 ## along and distance to arrays
  
  ## convert x,y to along curve and distance to curve (a,d) 
  ## co-ordinates. 0 distance along is at (x=-r,y=0)  
  
  ind <- x>=0 & y>0
  a[ind] <- q + x[ind]
  d[ind] <- y[ind]-r
  
  ind <- x>=0 & y<=0
  a[ind] <- -q - x[ind]
  d[ind] <- -r - y[ind]
  
  
  ind <- x < 0 
  a[ind] <- -atan(y[ind]/x[ind])*r
  d[ind] <- sqrt(x[ind]^2+y[ind]^2) - r
  
  ## create exclusion index
  
  ind <- abs(d)>r-r0 | (x>l & (x-l)^2+d^2 > (r-r0)^2)
  
  # f <- a*b # the original
  f <- a*b+d^2
  
}


#' Covariate test function for the horseshoe domain
#' 
#' @param x,y Points at which to evaluate the test function.
#' @description Implements a finite area test function the horseshoe domain.
#' @return Returns function evaluations. 
#' @usage covs.test(x, y)  
#' @export

covs.test = function(x,y){
  eval = rep(NA, length(x))
  ind <- y > 0
  eval[ind] = 2*exp(1 - 0.75*((x[ind]-1.5)^2 + 12*((y[ind])-0.5)^2))
  ind <- y <= 0
  eval[ind] = 2*exp(1 - 0.75*((x[ind]-2.5)^2 + 3*((-y[ind])-0.5)^2))
  ifelse(eval > 10^{-16}, eval, 0)
}