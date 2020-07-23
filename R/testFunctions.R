#' FELSPLINE test function
#' 
#' @param x,y Points at which to evaluate the test function.
#' @param r0 The test domain is a sort of bent sausage. This is the radius of the inner bend.
#' @param r The radius of the curve at the centre of the sausage.
#' @param l The length of an arm of the sausage.
#' @param b The rate at which the function increases per unit increase in distance along the centre line of the sausage.
#' @param exclude Should exterior points be set to NA?
#' @description Implements a finite area test function based on one proposed by Tim Ramsay (2002) proposed by 
#' Simon Wood (2008).
#' @return Returns function evaluations, or NAs for points outside the horseshoe domain. 
#' @usage fs.test(x, y, r0 = 0.1, r = 0.5, l = 3, b = 1, exclude = FALSE)  
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
#' coeff = fs.test(mesh$nodes[,1], mesh$nodes[,2], exclude = FALSE)
#' ## Create the FEM object
#' FEMfunction = FEM(coeff, FEMbasis)
#' ## Plot it
#' plot(FEMfunction)
#' @export

fs.test <- function (x, y, r0 = 0.1, r = 0.5, l = 3, b = 1, exclude = FALSE) 
{
  
  q <- pi * r/2
  a <- d <- x * 0
  
  ind <- x >= 0 & y > 0
  a[ind] <- q + x[ind]
  d[ind] <- y[ind] - r
  
  ind <- x >= 0 & y <= 0
  a[ind] <- (-q - x[ind])
  d[ind] <- -r - y[ind]
  
  ind <- x < 0
  a[ind] <- -atan(y[ind]/x[ind]) * r 
  d[ind] <-( sqrt(x[ind]^2 + y[ind]^2) - r )* (y[ind]/r0*(as.numeric(abs(y[ind])<=r0 & x[ind]>-0.5))+(as.numeric(abs(y[ind])>r0 | x[ind]<(-0.5))))
  
  ind <- abs(d) > r - r0 | (x > l & (x - l)^2 + d^2 > (r - r0)^2)
  
  f <- a * b + d^2
  
  if (exclude) 
    f[ind] <- NA
  
  attr(f, "exclude") <- ind
  f
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




#' FELSPLINE 3D test function
#' 
#' @param x,y,z Points at which to evaluate the test function.
#' @param r0 The test domain is a sort of bent sausage. This is the radius of the inner bend.
#' @param r The radius of the curve at the centre of the sausage.
#' @param l The length of an arm of the sausage.
#' @param b The rate at which the function increases per unit increase in distance along the centre line of the sausage.
#' @param exclude Should exterior points be set to NA?
#' @description Implements a finite area test function based on one proposed by Tim Ramsay (2002) and by 
#' Simon Wood (2008) in 3D.
#' @return Returns function evaluations, or NAs for points outside the horseshoe domain. 
#' @usage fs.test.3D(x, y, z, r0 = 0.25, r = 1.25, l = 5, b = 1, exclude = FALSE)  
#' @examples 
#' library(fdaPDE)
#' 
#' data(horseshoe2.5D)
#' mesh = horseshoe2.5D
#' FEMbasis=create.FEM.basis(mesh)
#' 
#' # Evaluation at nodes
#' sol_exact=fs.test.3D(mesh$nodes[,1],mesh$nodes[,3],mesh$nodes[,2])
#' plot(FEM(sol_exact, FEMbasis))
#' @export
fs.test.3D <- function (x, y, z, r0 = 0.25, r = 1.25, l = 5, b = 1, exclude = FALSE) 
{
  q <- pi * r/2
  a <- d <- x * 0
  f <- rep(NA, length(x))
  
  ind <- x >= 0 & z > 0
  a[ind] <- q + x[ind]
  d[ind] <- z[ind] - r
  f[ind] <- a[ind] * b + d[ind]^2 + y[ind]^2
  
  ind <- x >= 0 & z <= 0
  a[ind] <- (-q - x[ind])
  d[ind] <- -r - z[ind]
  f[ind] <- a[ind] * b + d[ind]^2 + y[ind]^2
  
  ind <- x < 0
  a[ind] <- -atan(z[ind]/x[ind]) * r 
  d[ind] <-( sqrt(x[ind]^2 + z[ind]^2) - r )
  
  f[ind] <- a[ind] * b + d[ind]^2 + y[ind]^2
  
  ind <- abs(d) > r - r0 | (x > l & (x - l)^2 + d^2 > (r - r0)^2)
  if (exclude) 
    f[ind] <- NA
  
  attr(f, "exclude") <- ind
  f
}
