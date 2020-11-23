#' Spatial regression with differential regularization
#'
#' @param observations A vector of length #observations with the observed data values over the domain.
#' If the \code{locations} argument is left NULL the vector of the observations have to be of length #nodes of the
#' mesh in the FEMbasis. In this case, each observation is associated to the corresponding node in the mesh.
#' If the observations are observed only on a subset of the mesh nodes, fill with \code{NA} the values of the vector
#' \code{observations} in correspondence of unobserved data.
#' @param locations A #observations-by-2 matrix in the 2D case and #observations-by-3 matrix in the 2.5D and 3D case, where
#' each row specifies the spatial coordinates \code{x} and \code{y} (and \code{z} in 2.5D and 3D) of the corresponding
#' observation in the vector \code{observations}.
#' If the locations of the observations coincide with (or are a subset of) the nodes of the mesh in the \code{FEMbasis},
#' leave the parameter \code{locations = NULL} for a faster implementation.
#' @param FEMbasis A \code{FEMbasis} object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param covariates A #observations-by-#covariates matrix where each row represents the covariates associated with
#' the corresponding observed data value in \code{observations} and each column is a different covariate.
#' @param PDE_parameters A list specifying the parameters of the PDE in the regularizing term. Default is NULL, i.e.
#' regularization is by means of the Laplacian (stationary, isotropic case).
#' If the coefficients of the PDE are constant over the domain \code{PDE_parameters} must contain:
#' \itemize{
#'    \item{\code{K}, a 2-by-2 matrix of diffusion coefficients. This induces an anisotropic
#' smoothing with a preferential direction that corresponds to the first eigenvector of the diffusion matrix K;}
#'    \item{\code{b}, a vector of length 2 of advection coefficients. This induces a
#' smoothing only in the direction specified by the vector \code{b};}
#'    \item{\code{c}, a scalar reaction coefficient. \code{c} induces a shrinkage of the surface to zero.}
#' }
#' If the coefficients of the PDE are space-varying \code{PDE_parameters} must contain:
#' \itemize{
#' \item{\code{K}, a function that for each spatial location in the spatial domain (indicated by the vector of the 2
#' spatial coordinates) returns a 2-by-2 matrix of diffusion coefficients. The function must support recycling for
#' efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
#' an array with dimensions 2-by-2-by-#points.}
#' \item{\code{b}, a function that for each spatial location in the spatial domain returns
#' a vector of length 2 of transport coefficients. The function must support recycling for efficiency reasons, thus
#' if the input parameter is a #point-by-2 matrix, the output should be
#' a matrix with dimensions 2-by-#points;}
#' \item{\code{c}, a function that for each spatial location in the spatial domain  returns a scalar reaction coefficient.
#' The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
#' a vector with length #points;}
#' \item{\code{u}, a function that for each spatial location in the spatial domain  returns a scalar reaction coefficient.
#' \code{u} induces a reaction effect. The function must support recycling for efficiency reasons, thus if the input
#' parameter is a #point-by-2 matrix, the output should be
#' a vector with length #points.}
#' }
#' For 2.5D and 3D, only the Laplacian is available (\code{PDE_parameters=NULL}).
#' @param BC A list with two vectors:
#'  \code{BC_indices}, a vector with the indices in \code{nodes} of boundary nodes where a Dirichlet Boundary Condition should be applied;
#'  \code{BC_values}, a vector with the values that the spatial field must take at the nodes indicated in \code{BC_indices}.
#' @param incidence_matrix A #regions-by-#triangles/tetrahedrons matrix where the element (i,j) equals 1 if the j-th
#' triangle/tetrahedron is in the i-th region and 0 otherwise.
#' This is needed only for areal data. In case of pointwise data, this parameter is set to \code{NULL}.
#' @param areal.data.avg Boolean. It involves the computation of Areal Data. If \code{TRUE} the areal data are averaged, otherwise not.
#' @param search a flag to decide the search algorithm type (tree or naive or walking search algorithm).
#' @param bary.locations A list with three vectors:
#'  \code{locations}, location points which are same as the given locations options. (checks whether both locations are the same);
#'  \code{element ids}, a vector of element id of the points from the mesh where they are located;
#'  \code{barycenters}, a vector of barycenter of points from the located element.
#' @param family This parameter specify the distibution within exponential family used for GLM model.
#' The following distribution are implemented: "binomial", "exponential", "gamma", "poisson", "gaussian", "invgaussian".
#' The default link function for binomial is \code{logit} if you want either \code{probit} or \code{clogloc} set \code{family = "probit"}, \code{family = "cloglog"}.
#' @param mu0 This parameter is a vector that set the starting point for FPIRLS algorithm. It represent an initial guess of the location parameter.
#' Default is set to observation for non binary distribution while equal to \code{0.5(observations + 0.5)} for binary data.
#' @param scale.param Dispersion parameter of the chosen distribution. This is only required for "gamma", "gaussian", "invgaussian".
#' User may specify the parameter as a positive real number. If the parameter is not supplied, it is estimated from data according to Wilhelm Sangalli 2016.
#' @param threshold.FPIRLS This parameter is used for arresting algorithm iterations. Algorithm stops when two successive iterations lead to improvement in penalized log-likelihood smaller than threshold.FPIRLS.
#' Default value \code{threshold.FPIRLS = 0.0002020}.
#' @param max.steps.FPIRLS This parameter is used to limit the maximum number of iteration.
#' Default value \code{max.steps.FPIRLS=15}.
#' @param lambda.selection.criterion This parameter is used to select the optimization method related to the smoothing parameter \code{lambda}.
#' The following methods are implemented: 'grid', 'newton', 'newton_fd'.
#' The former is a pure evaluation method, therefore a vector of \code{lambda} testing penalizations must be provided.
#' The remaining two are optimization methods that automatically select the best penalization according to \code{lambda.selection.lossfunction} criterion.
#' They implement respectively a pure Newton method and a finite differences Newton method.
#' Default value \code{lambda.selection.criterion='grid'}
#' @param DOF.evaluation This parameter is used to identify if and how degrees of freedom computation has to be performed.
#' The following possibilities are allowed: NULL, 'exact' and 'stochastic'
#' In the former case no degree of freedom is computed, while the other two methods enable computation.
#' Stochastic computation of DOFs may be slightly less accurate than its deterministic counterpart, but is highly suggested for meshes of more than 5000 nodes, being fairly less time consuming.
#' Default value \code{DOF.evaluation=NULL}
#' @param lambda.selection.lossfunction This parameter is used to understand if some loss function has to be evaluated.
#' The following possibilities are allowed: NULL and 'GCV' (generalized cross validation)
#' The former case is that of \code{lambda.selection.criterion='grid'} pure evaluation, while the second can be employed for optimization methods.
#' Default value \code{lambda.selection.lossfunction=NULL}
#' @param lambda a vector of spatial smoothing parameters to be provided for evaluation if \code{lambda.selection.criterion='grid'}, an optional initialization otherwise
#' @param DOF.stochastic.realizations This parameter is considered only when \code{DOF.evaluation = 'stochastic'}.
#' It is a positive integer that represents the number of uniform random variables used in stochastic GCV computation.
#' Default value \code{DOF.stochastic.realizations=100}.
#' @param DOF.stochastic.seed This parameter is considered only when \code{DOF.evaluation = 'stochastic'}.
#' It is a positive integer that represents user defined seed employed in stochastic GCV computation.
#' Default value \code{DOF.stochastic.seed = 0} means random.
#' @param DOF.matrix Matrix of degrees of freedom. This parameter can be used if the DOF.matrix corresponding to \code{lambda} is available from precedent computation. This allows to save time
#' since the computation of the DOFs is the most expensive part of GCV.
#' @param GCV.inflation.factor Tuning parameter used for the estimation of GCV. Default value \code{GCV.inflation.factor = 1.0} or \code{1.8} in GAM.
#' It is advised to set it grather than 1 to avoid overfitting.
#' @param lambda.optimization.tolerance Tolerance parameter, a double between 0 and 1 that fixes how much precision is required by the optimization method: the smaller the parameter, the higher the accuracy.
#' Used only if \code{lambda.selection.criterion='newton'} or \code{lambda.selection.criterion='newton_fd'}.
#' Default value \code{lambda.optimization.tolerance=0.05}.
#' @return A list with the following variables in \code{family="gaussian"} case:
#' \itemize{
#'    \item{\code{fit.FEM}}{A \code{FEM} object that represents the fitted spatial field.}
#'    \item{\code{PDEmisfit.FEM}}{A \code{FEM} object that represents the Laplacian of the estimated spatial field.}
#'    \item{\code{solution}}{A list, note that all terms are matrices or row vectors: the \code{j}th column represents the vector of related to \code{lambda[j]} if \code{lambda.selection.criterion="grid"} and \code{lambda.selection.lossfunction=NULL}.
#'          In all the other cases is returned just the column related to the best smoothing parameter
#'          \item{\code{f}}{Matrix, estimate of function f, first half of solution vector}
#'          \item{\code{g}}{Matrix, second half of solution vector}
#'          \item{\code{z_hat}}{Matrix, prediction of the output in the locations}
#'          \item{\code{beta}}{If \code{covariates} is not \code{NULL}, a matrix with number of rows equal to the number of covariates and number of columns equal to length of lambda. It is the regression coefficients estimate}
#'          \item{\code{rmse}}{Estimate of the root mean square error in the locations}
#'          \item{\code{estimated_sd}}{Estiimate of the standard deviation of the error}
#'          }
#'    \item{\code{optimization}}{A detailed list of optimization related data:
#'          \item{\code{lambda_solution}}{numerical value of best lambda acording to \code{lambda.selection.lossfunction}, -1 if \code{lambda.selection.lossfunction=NULL}}
#'          \item{\code{lambda_position}}{integer, postion in \code{lambda_vector} of best lambda acording to \code{lambda.selection.lossfunction}, -1 if \code{lambda.selection.lossfunction=NULL}}
#'          \item{\code{GCV}}{numeric value of GCV in correspondence of the optimum}
#'          \item{\code{optimization_details}}{list containing further information about the optimization method used and the nature of its termination, eventual number of iterations}
#'          \item{\code{dof}}{numeric vector, value of DOFs for all the penalizations it has been computed, empty if not computed}
#'          \item{\code{lambda_vector}}{numeric value of the penalization factors passed by the user or found in the iterations of the optimization method}
#'          \item{\code{GCV_vector}}{numeric vector, value of GCV for all the penalizations it has been computed}
#'          }
#'    \item{\code{time}}{Duration of the entire optimization computation}
#'    \item{\code{bary.locations}}{A barycenter information of the given locations if the locations are not mesh nodes.}
#' }
#' A list with the following variables in others GAM case:
#' \itemize{
#'    \item{\code{fit.FEM}}{A \code{FEM} object that represents the fitted spatial field.}
#'    \item{\code{PDEmisfit.FEM}}{A \code{FEM} object that represents the Laplacian of the estimated spatial field.}
#'    \item{\code{beta}}{If covariates is not \code{NULL}, a matrix with number of rows equal to the number of covariates and numer of columns equal to length of lambda.  The \code{j}th column represents the vector of regression coefficients when
#' the smoothing parameter is equal to \code{lambda[j]}.}
#'    \item{\code{edf}}{If GCV is \code{TRUE}, a scalar or vector with the trace of the smoothing matrix for each value of the smoothing parameter specified in \code{lambda}.}
#'    \item{\code{stderr}}{If GCV is \code{TRUE}, a scalar or vector with the estimate of the standard deviation of the error for each value of the smoothing parameter specified in \code{lambda}.}
#'    \item{\code{GCV}}{If GCV is \code{TRUE}, a  scalar or vector with the value of the GCV criterion for each value of the smoothing parameter specified in \code{lambda}.}
#'    \item{\code{bary.locations}}{A barycenter information of the given locations if the locations are not mesh nodes.}
#'    \item{\code{fn_hat}}{ A matrix with number of rows equal to number of locations and number of columns equal to length of lambda. Each column contain the evaluaton of the spatial field in the location points.}
#'    \item{\code{J_minima}}{A vector of the same length of lambda, containing the reached minima for each value of the smoothing parameter.}
#'    \item {\code{variance.est}}{ A vector which return the variance estimates for the Generative Additive Models}
#' }
#' @description This function implements a spatial regression model with differential regularization.
#'  The regularizing term involves a Partial Differential Equation (PDE). In the simplest case the PDE involves only the
#'  Laplacian of the spatial field, that induces an isotropic smoothing. When prior information about the anisotropy or
#'  non-stationarity is available the PDE involves a general second order linear differential operator with possibly
#'  space-varying coefficients.
#'  The technique accurately handle data distributed over irregularly shaped domains. Moreover, various conditions
#'  can be imposed at the domain boundaries.
#' @usage smooth.FEM(locations = NULL, observations, FEMbasis,
#'  covariates = NULL, PDE_parameters = NULL, BC = NULL,
#'  incidence_matrix = NULL, areal.data.avg = TRUE,
#'  search = "tree", bary.locations = NULL,
#'  family = "gaussian", mu0 = NULL, scale.param = NULL, threshold.FPIRLS = 0.0002020, 
#'  max.steps.FPIRLS = 15, lambda.selection.criterion = "grid", DOF.evaluation = NULL, 
#'  lambda.selection.lossfunction = NULL, lambda = NULL, DOF.stochastic.realizations = 100,
#'  DOF.stochastic.seed = 0, DOF.matrix = NULL, GCV.inflation.factor = 1, 
#'  lambda.optimization.tolerance = 0.05)
#' @export

#' @references
#' \itemize{
#'    \item{Sangalli, L. M., Ramsay, J. O., Ramsay, T. O. (2013). Spatial spline regression models.
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 75(4), 681-703.}
#'    \item{Azzimonti, L., Sangalli, L. M., Secchi, P., Domanin, M., Nobile, F. (2015). Blood flow velocity field estimation
#' via spatial regression with PDE penalization. Journal of the American Statistical Association, 110(511), 1057-1071.}
#'    \item{Matthieu Wilhelm & Laura M. Sangalli (2016). Generalized spatial regression with differential regularization.
#'  Journal of Statistical Computation and Simulation, 86:13, 2497-2518.}
#' }
#' @examples
#' library(fdaPDE)
#'
#' #### No prior information about anysotropy/non-stationarity (laplacian smoothing) ####
#' data(horseshoe2D)
#' boundary_nodes = horseshoe2D$boundary_nodes
#' boundary_segments = horseshoe2D$boundary_segments
#' locations = horseshoe2D$locations
#'
#' mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)
#' FEMbasis = create.FEM.basis(mesh)
#' lambda = 10^-1
#' # no covariate
#' data = fs.test(mesh$nodes[,1], mesh$nodes[,2]) + rnorm(nrow(mesh$nodes), sd = 0.5)
#'
#' solution = smooth.FEM(observations = data, FEMbasis = FEMbasis, lambda = lambda)
#' plot(solution$fit.FEM)
#'
#' # with covariates
#' covariate = covs.test(mesh$nodes[,1], mesh$nodes[,2])
#' data = fs.test(mesh$nodes[,1], mesh$nodes[,2]) + 2*covariate + rnorm(nrow(mesh$nodes), sd = 0.5)
#'
#' solution = smooth.FEM(observations = data, covariates = covariate, 
#'                       FEMbasis = FEMbasis, lambda = lambda)
#'
#' # beta estimate:
#' solution$solution$beta
#' # non-parametric estimate:
#' plot(solution$fit.FEM)
#'
#' # Choose lambda with GCV - stochastic grid evaluation:
#' lambda = 10^(-2:0)
#' solution = smooth.FEM(observations = data,
#'                             covariates = covariate,
#'                             FEMbasis = FEMbasis,
#'                             lambda = lambda, DOF.evaluation = 'stochastic', 
#'                             lambda.selection.lossfunction = 'GCV')
#' bestLambda = solution$optimization$lambda_solution
#' # Choose lambda with GCV - Newton finite differences stochastic evaluation -:
#' solution = smooth.FEM(observations = data,
#'                             covariates = covariate,
#'                             FEMbasis = FEMbasis,
#'                             DOF.evaluation = 'stochastic', lambda.selection.lossfunction = 'GCV')
#' bestLambda = solution$optimization$lambda_solution
#'
#'
#' #### Smoothing with prior information about anysotropy/non-stationarity and boundary conditions ####
#' # See Azzimonti et al. for reference to the current exemple
#' data(quasicircle2D)
#' boundary_nodes = quasicircle2D$boundary_nodes
#' boundary_segments = quasicircle2D$boundary_segments
#' locations = quasicircle2D$locations
#' data = quasicircle2D$data
#'
#' mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)
#' FEMbasis = create.FEM.basis(mesh)
#' lambda = 10^-2
#'
#' # Set the PDE parameters
#' R = 2.8
#' K1 = 0.1
#' K2 = 0.2
#' beta = 0.5
#' K_func<-function(points)
#' {
#'   output = array(0, c(2, 2, nrow(points)))
#'   for (i in 1:nrow(points))
#'     output[,,i]=10*rbind(c(points[i,2]^2+K1*points[i,1]^2+K2*(R^2-points[i,1]^2-points[i,2]^2),
#'                            (K1-1)*points[i,1]*points[i,2]),
#'                          c((K1-1)*points[i,1]*points[i,2],
#'                            points[i,1]^2+K1*points[i,2]^2+K2*(R^2-points[i,1]^2-points[i,2]^2)))
#'   output
#' }
#'
#' b_func<-function(points)
#' {
#'   output = array(0, c(2, nrow(points)))
#'   for (i in 1:nrow(points))
#'     output[,i] = 10*beta*c(points[i,1],points[i,2])
#'   output
#' }
#'
#' c_func<-function(points)
#' {
#'   rep(c(0), nrow(points))
#' }
#'
#' u_func<-function(points)
#' {
#'   rep(c(0), nrow(points))
#' }
#' PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)
#'
#' # Set the boundary conditions
#' BC = NULL
#' BC$BC_indices = which(mesh$nodesmarkers == 1) # b.c. on the complete boundary
#' BC$BC_values = rep(0,length(BC$BC_indices)) # homogeneus b.c.
#'
#' # Since the data locations are a subset of the mesh nodes for a faster solution use:
#' dataNA = rep(NA, FEMbasis$nbasis)
#' dataNA[mesh$nodesmarkers == 0] = data
#' #grid evaluation
#' solution = smooth.FEM(observations = dataNA,
#'                             FEMbasis = FEMbasis,
#'                             lambda = lambda,
#'                             PDE_parameters = PDE_parameters,
#'                             BC = BC)
#' plot(solution$fit.FEM)
#' image(solution$fit.FEM)
#' # Newton's method
#' solution = smooth.FEM(observations = dataNA,
#'                             FEMbasis = FEMbasis,
#'                             PDE_parameters = PDE_parameters,
#'                             BC = BC)
#' plot(solution$fit.FEM)
#' image(solution$fit.FEM)
#'
#' #### Smoothing with areal data ####
#' # See Azzimonti et al. for reference to the current exemple
#' data(quasicircle2Dareal)
#' incidence_matrix = quasicircle2Dareal$incidence_matrix
#' data = quasicircle2Dareal$data
#' mesh = quasicircle2Dareal$mesh
#'
#' FEMbasis = create.FEM.basis(mesh)
#' lambda = 10^-4
#'
#' # Set the PDE parameters
#' R = 2.8
#' K1 = 0.1
#' K2 = 0.2
#' beta = 0.5
#' K_func<-function(points)
#' {
#'   output = array(0, c(2, 2, nrow(points)))
#'   for (i in 1:nrow(points))
#'     output[,,i]=10*rbind(c(points[i,2]^2+K1*points[i,1]^2+K2*(R^2-points[i,1]^2-points[i,2]^2),
#'                            (K1-1)*points[i,1]*points[i,2]),
#'                          c((K1-1)*points[i,1]*points[i,2],
#'                            points[i,1]^2+K1*points[i,2]^2+K2*(R^2-points[i,1]^2-points[i,2]^2)))
#'   output
#' }
#'
#' b_func<-function(points)
#' {
#'   output = array(0, c(2, nrow(points)))
#'   for (i in 1:nrow(points))
#'     output[,i] = 10*beta*c(points[i,1],points[i,2])
#'   output
#' }
#'
#' c_func<-function(points)
#' {
#'   rep(c(0), nrow(points))
#' }
#'
#' u_func<-function(points)
#' {
#'   rep(c(0), nrow(points))
#' }
#' PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)
#'
#' # Set the boundary conditions
#' BC = NULL
#' BC$BC_indices = which(mesh$nodesmarkers == 1) # b.c. on the complete boundary
#' BC$BC_values = rep(0,length(BC$BC_indices)) # homogeneus b.c.
#' #grid evaluation
#' solution = smooth.FEM(observations = data,
#'                             incidence_matrix = incidence_matrix,
#'                             FEMbasis = FEMbasis,
#'                             lambda = lambda,
#'                             PDE_parameters = PDE_parameters,
#'                             BC = BC)
#' plot(solution$fit.FEM)
#' image(solution$fit.FEM)
#' #Newton's method
#' solution = smooth.FEM(observations = data,
#'                             incidence_matrix = incidence_matrix,
#'                             FEMbasis = FEMbasis,
#'                             PDE_parameters = PDE_parameters,
#'                             BC = BC)
#' plot(solution$fit.FEM)
#' image(solution$fit.FEM)
#'

smooth.FEM<-function(locations = NULL, observations, FEMbasis,
                     covariates = NULL, PDE_parameters = NULL, BC = NULL,
                     incidence_matrix = NULL, areal.data.avg = TRUE,
                     search = "tree", bary.locations = NULL,
                     family = "gaussian", mu0 = NULL, scale.param = NULL, threshold.FPIRLS = 0.0002020, max.steps.FPIRLS = 15,
                     lambda.selection.criterion = "grid", DOF.evaluation = NULL, lambda.selection.lossfunction = NULL,
                     lambda = NULL, DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0, DOF.matrix = NULL, GCV.inflation.factor = 1, lambda.optimization.tolerance = 0.05)
{
  # Mesh identification
  if(class(FEMbasis$mesh) == "mesh.2D")
  {
    ndim = 2
    mydim = 2
  }else if(class(FEMbasis$mesh) == "mesh.2.5D")
  {
    ndim = 3
    mydim = 2
  }else if(class(FEMbasis$mesh) == "mesh.3D")
  {
    ndim = 3
    mydim = 3
  }else
  {
    stop('unknown mesh class.')
  }

  ##################### Checking parameters, sizes and conversion ################################

  # Preliminary consistency of optimization parameters
  if(lambda.selection.criterion == "grid")
  {
    optim = 0
  }else if(lambda.selection.criterion == "newton")
  {
    optim = 1
  }else if(lambda.selection.criterion == "newton_fd")
  {
    optim = 2
  }else
  {
    stop("'lambda.selection.criterion' must belong to the following list: 'none', 'grid', 'newton', 'newton_fd'.")
  }
  
  if(is.null(DOF.evaluation))
  {
    optim = c(optim,0)
  }else if(DOF.evaluation == 'stochastic')
  {
    optim = c(optim,1)
  }else if(DOF.evaluation == 'exact')
  {
    optim = c(optim,2)
  }else
  {
    stop("'DOF.evaluation' must be NULL, 'stochastic' or 'exact'.")
  }
  
  if(is.null(lambda.selection.lossfunction))
  {
    optim = c(optim,0)
  }else if(lambda.selection.lossfunction == 'GCV')
  {
    optim = c(optim,1)
  }else
  {
    stop("'lambda.selection.lossfunction' has to be 'GCV'.")
  }
  
  # OPTIMIZATION NOT IMPLEMENTED FOR GAM
  if(family != 'gaussian'& optim[1]!=0)
    stop("'lambda.selection.criterion' = 'grid' is the only method implemented for GAM problems")
  
  # --> General consistency rules
  if(optim[2]!=0 & optim[3]!=1)
  {
    warning("Dof are computed, setting 'lambda.selection.lossfunction' to 'GCV'")
    optim[3] = 1
  }
  if(optim[1]==1 & optim[2]!=2)
  {
    warning("This method needs evaluate DOF in an 'exact' way, selecting 'DOF.evaluation'='exact'")
    optim[2] = 2
  }
  if(!is.null(BC) & optim[1]==1)
  {
    warning("'newton' 'lambda.selection.criterion' can't be performed with non-NULL boundary conditions, using 'newton_fd' instead")
    optim[1] = 2
  }
  if((optim[1]==2 & optim[2]==0) || (optim[1]==0 & optim[2]==0 & optim[3]==1 & is.null(DOF.matrix)))
  {
    warning("This method needs evaluate DOF, selecting 'DOF.evaluation'='stochastic'")
    optim[2] = 1
  }
  if(optim[1]!=0 & optim[3]==0)
  {
    warning("An optimized method needs a loss function to perform the evaluation, selecting 'lambda.selection.lossfunction' as 'GCV'")
    optim[3] = 1
  }
  
  if(is.null(lambda) & optim[1] == 0)
  {
    warning("the lambda passed is NULL, passing to default optimized methods")
    optim = c(2,1,1)
  }
  

  if(any(lambda<=0))
  	stop("'lambda' can not be less than or equal to 0")

  # Search algorithm
  if(search=="naive")
  {
    search=1
  }else if(search=="tree")
  {
    search=2
  }else
  {
    stop("'search' must be either 'tree' or 'naive'.")
  }

  # If locations is null but bary.locations is not null, use the locations in bary.locations
  if(is.null(locations) & !is.null(bary.locations))
  {
    locations = bary.locations$locations
    locations = as.matrix(locations)
  }

  ## Converting to format for internal usage
  if(!is.null(locations))
    locations = as.matrix(locations)
  observations = as.matrix(observations)
  if(!is.null(covariates))
    covariates = as.matrix(covariates)
  if(!is.null(incidence_matrix))
    incidence_matrix = as.matrix(incidence_matrix)
  if(!is.null(BC))
  {
    BC$BC_indices = as.matrix(BC$BC_indices)
    BC$BC_values = as.matrix(BC$BC_values)
  }
  if(!is.null(lambda))
    lambda = as.matrix(lambda)
  if(!is.null(DOF.matrix))
    DOF.matrix = as.matrix(DOF.matrix)

  space_varying = checkSmoothingParameters(locations = locations, observations = observations, FEMbasis = FEMbasis,
    covariates = covariates, PDE_parameters = PDE_parameters, BC = BC,
    incidence_matrix = incidence_matrix, areal.data.avg = areal.data.avg,
    search = search, bary.locations = bary.locations,
    optim = optim, lambda = lambda, DOF.stochastic.realizations = DOF.stochastic.realizations, DOF.stochastic.seed = DOF.stochastic.seed,
    DOF.matrix = DOF.matrix, GCV.inflation.factor = GCV.inflation.factor, lambda.optimization.tolerance = lambda.optimization.tolerance)

  # If I have PDE non-sv case I need (constant) matrices as parameters
  if(!is.null(PDE_parameters) & space_varying == FALSE)
  {
    PDE_parameters$K = as.matrix(PDE_parameters$K)
    PDE_parameters$b = as.matrix(PDE_parameters$b)
    PDE_parameters$c = as.matrix(PDE_parameters$c)
  }


  checkSmoothingParametersSize(locations = locations, observations = observations, FEMbasis = FEMbasis,
    covariates = covariates, PDE_parameters = PDE_parameters, incidence_matrix = incidence_matrix,
    BC = BC, space_varying = space_varying, ndim = ndim, mydim = mydim,
    lambda = lambda, DOF.matrix = DOF.matrix)


  # Check whether the locations coincide with the mesh nodes (should be put after all the validations)
  if (!is.null(locations))
  {
    if(dim(locations)[1]==dim(FEMbasis$mesh$nodes)[1] & dim(locations)[2]==dim(FEMbasis$mesh$nodes)[2])
    {
      sum1=0
      sum2=0
      for (i in 1:nrow(locations))
      {
        sum1 = sum1 + abs(locations[i,1]-FEMbasis$mesh$nodes[i,1])
        sum2 = sum2 + abs(locations[i,2]-FEMbasis$mesh$nodes[i,2])
      }
      if (sum1==0 & sum2==0)
      {
        message("No search algorithm is used because the locations coincide with the nodes.")
        locations = NULL #In principle, R uses pass-by-value semantics in its function calls. So put ouside of checkSmoothingParameters function.
      }
    }
  }

  # FAMILY CHECK
  family_admit = c("binomial", "exponential", "gamma", "poisson", "gaussian")
  if(sum(family==family_admit)==0){
   stop("'family' parameter required.\nCheck if it is one of the following: binomial, exponential, gamma, poisson, gaussian")
  }



  ################## End checking parameters, sizes and conversion #############################
  if(family == "gaussian")
  {
    #----------------------------------------------------#
    ############# Standard Smooth method #################
    #----------------------------------------------------#
    if(class(FEMbasis$mesh) == 'mesh.2D' & is.null(PDE_parameters))
    {
      bigsol = NULL
      bigsol = CPP_smooth.FEM.basis(locations = locations, observations = observations, FEMbasis = FEMbasis,
        covariates = covariates, ndim = ndim, mydim = mydim, BC = BC,
        incidence_matrix = incidence_matrix, areal.data.avg = areal.data.avg,
        search = search, bary.locations = bary.locations,
        optim = optim, lambda = lambda, DOF.stochastic.realizations = DOF.stochastic.realizations, DOF.stochastic.seed = DOF.stochastic.seed,
        DOF.matrix = DOF.matrix, GCV.inflation.factor = GCV.inflation.factor, lambda.optimization.tolerance = lambda.optimization.tolerance)
      numnodes = nrow(FEMbasis$mesh$nodes)
    }else if(class(FEMbasis$mesh) == 'mesh.2D' & !is.null(PDE_parameters) & space_varying == FALSE)
    {
      bigsol = NULL
      bigsol = CPP_smooth.FEM.PDE.basis(locations = locations, observations = observations, FEMbasis = FEMbasis,
        covariates = covariates, PDE_parameters = PDE_parameters, ndim = ndim, mydim = mydim, BC = BC,
        incidence_matrix = incidence_matrix, areal.data.avg = areal.data.avg,
        search = search, bary.locations = bary.locations,
        optim = optim, lambda = lambda, DOF.stochastic.realizations = DOF.stochastic.realizations, DOF.stochastic.seed = DOF.stochastic.seed, 
        DOF.matrix = DOF.matrix, GCV.inflation.factor = GCV.inflation.factor, lambda.optimization.tolerance = lambda.optimization.tolerance)
      numnodes = nrow(FEMbasis$mesh$nodes)
    }else if(class(FEMbasis$mesh) == 'mesh.2D' & !is.null(PDE_parameters) & space_varying == TRUE)
    {
      bigsol = NULL
      bigsol = CPP_smooth.FEM.PDE.sv.basis(locations = locations, observations = observations, FEMbasis = FEMbasis,
        covariates=covariates, PDE_parameters = PDE_parameters, ndim = ndim, mydim = mydim, BC=BC,
        incidence_matrix=incidence_matrix, areal.data.avg = areal.data.avg,
        search=search, bary.locations = bary.locations,
        optim = optim, lambda = lambda, DOF.stochastic.realizations = DOF.stochastic.realizations, DOF.stochastic.seed = DOF.stochastic.seed,
        DOF.matrix = DOF.matrix, GCV.inflation.factor = GCV.inflation.factor, lambda.optimization.tolerance = lambda.optimization.tolerance)
      numnodes = nrow(FEMbasis$mesh$nodes)
    }else if(class(FEMbasis$mesh) == 'mesh.2.5D')
    {
      bigsol = NULL
      # if(!is.null(locations))
      #   stop("The option locations!=NULL for manifold domains is currently not implemented")
      bigsol = CPP_smooth.manifold.FEM.basis(locations = locations, observations = observations, FEMbasis = FEMbasis,
       covariates = covariates, ndim = ndim, mydim = mydim, BC = BC,
       incidence_matrix = incidence_matrix, areal.data.avg = areal.data.avg,
       search = search, bary.locations = bary.locations,
       optim = optim, lambda = lambda, DOF.stochastic.realizations = DOF.stochastic.realizations, DOF.stochastic.seed = DOF.stochastic.seed,
       DOF.matrix = DOF.matrix, GCV.inflation.factor = GCV.inflation.factor, lambda.optimization.tolerance = lambda.optimization.tolerance)
      numnodes = nrow(FEMbasis$mesh$nodes)
  	}else if(class(FEMbasis$mesh) == 'mesh.3D' & is.null(PDE_parameters))
    {
      bigsol = NULL
      bigsol = CPP_smooth.volume.FEM.basis(locations = locations, observations = observations, FEMbasis = FEMbasis,
        covariates = covariates, ndim = ndim, mydim = mydim, BC = BC,
        incidence_matrix = incidence_matrix, areal.data.avg = areal.data.avg,
        search = search, bary.locations = bary.locations,
        optim = optim, lambda = lambda, DOF.stochastic.realizations = DOF.stochastic.realizations, DOF.stochastic.seed = DOF.stochastic.seed,
        DOF.matrix = DOF.matrix, GCV.inflation.factor = GCV.inflation.factor, lambda.optimization.tolerance = lambda.optimization.tolerance)
      numnodes = nrow(FEMbasis$mesh$nodes)
    } else if(class(FEMbasis$mesh) == 'mesh.3D' & !is.null(PDE_parameters) & space_varying==FALSE)
    {
  	  bigsol = NULL
      bigsol = CPP_smooth.volume.FEM.PDE.basis(locations = locations, observations = observations, FEMbasis = FEMbasis,
        covariates = covariates, PDE_parameters=PDE_parameters, ndim = ndim, mydim = mydim, BC = BC,
        incidence_matrix = incidence_matrix, areal.data.avg = areal.data.avg,
        search = search, bary.locations = bary.locations,
        optim = optim, lambda = lambda, DOF.stochastic.realizations = DOF.stochastic.realizations, DOF.stochastic.seed = DOF.stochastic.seed,
        DOF.matrix = DOF.matrix, GCV.inflation.factor = GCV.inflation.factor, lambda.optimization.tolerance = lambda.optimization.tolerance)
      numnodes = nrow(FEMbasis$mesh$nodes)
  	} else if(class(FEMbasis$mesh) == 'mesh.3D' & !is.null(PDE_parameters) & space_varying==TRUE)
  	{
      bigsol = NULL
      bigsol = CPP_smooth.volume.FEM.PDE.sv.basis(locations = locations, observations = observations, FEMbasis = FEMbasis,
        covariates = covariates, PDE_parameters=PDE_parameters, ndim = ndim, mydim = mydim, BC = BC,
        incidence_matrix = incidence_matrix, areal.data.avg = areal.data.avg,
        search = search, bary.locations = bary.locations,
        optim = optim, lambda = lambda, DOF.stochastic.realizations = DOF.stochastic.realizations, DOF.stochastic.seed = DOF.stochastic.seed,
        DOF.matrix = DOF.matrix, GCV.inflation.factor = GCV.inflation.factor, lambda.optimization.tolerance = lambda.optimization.tolerance)
      numnodes = nrow(FEMbasis$mesh$nodes)
  	}
  } else
  {
    #----------------------------------------------------#
    ############# GAMs: FPIRLS algorithm #################
    #----------------------------------------------------#
    checkGAMParameters(observations = observations, max.steps.FPIRLS = max.steps.FPIRLS, mu0 = mu0, scale.param = scale.param, threshold.FPIRLS = threshold.FPIRLS, family = family)

    if(class(FEMbasis$mesh) == 'mesh.2D' & is.null(PDE_parameters))
    {
      bigsol = NULL
      bigsol = CPP_smooth.GAM.FEM(locations = locations, observations = observations, FEMbasis = FEMbasis,
        covariates = covariates, ndim = ndim, mydim = mydim, BC = BC,
        incidence_matrix = incidence_matrix, areal.data.avg = areal.data.avg,
        FAMILY=family, mu0 = mu0, max.steps.FPIRLS = max.steps.FPIRLS, scale.param = scale.param, threshold.FPIRLS = threshold.FPIRLS,
        search = search, bary.locations = bary.locations,
        optim = optim, lambda = lambda, DOF.stochastic.realizations = DOF.stochastic.realizations, DOF.stochastic.seed = DOF.stochastic.seed,
        DOF.matrix = DOF.matrix, GCV.inflation.factor = GCV.inflation.factor, lambda.optimization.tolerance = lambda.optimization.tolerance)
      numnodes = nrow(FEMbasis$mesh$nodes)
    }else if(class(FEMbasis$mesh) == 'mesh.2D' & !is.null(PDE_parameters) & space_varying == FALSE)
    {
        bigsol = NULL
        bigsol = CPP_smooth.GAM.FEM.PDE.basis(locations = locations, observations = observations, FEMbasis = FEMbasis,
          covariates = covariates, PDE_parameters = PDE_parameters, ndim = ndim, mydim = mydim, BC = BC,
          incidence_matrix = incidence_matrix, areal.data.avg = areal.data.avg,
          FAMILY = family, mu0 = mu0, max.steps.FPIRLS = max.steps.FPIRLS, scale.param = scale.param, threshold.FPIRLS = threshold.FPIRLS,
          search = search, bary.locations = bary.locations,
          optim = optim, lambda = lambda, DOF.stochastic.realizations = DOF.stochastic.realizations, DOF.stochastic.seed = DOF.stochastic.seed,
          DOF.matrix = DOF.matrix, GCV.inflation.factor = GCV.inflation.factor, lambda.optimization.tolerance = lambda.optimization.tolerance)
        numnodes = nrow(FEMbasis$mesh$nodes)
    }else if(class(FEMbasis$mesh) == 'mesh.2D' & !is.null(PDE_parameters) & space_varying == TRUE)
    {
      bigsol = NULL
      bigsol = CPP_smooth.GAM.FEM.PDE.sv.basis(locations = locations, observations = observations, FEMbasis = FEMbasis,
        covariates = covariates, PDE_parameters = PDE_parameters, ndim = ndim, mydim = mydim, BC = BC,
        incidence_matrix = incidence_matrix, areal.data.avg = areal.data.avg,
        FAMILY = family, mu0 = mu0, max.steps.FPIRLS = max.steps.FPIRLS, scale.param = scale.param, threshold.FPIRLS = threshold.FPIRLS,
        search = search, bary.locations = bary.locations,
        optim = optim, lambda = lambda, DOF.stochastic.realizations = DOF.stochastic.realizations, DOF.stochastic.seed = DOF.stochastic.seed,
        DOF.matrix = DOF.matrix, GCV.inflation.factor = GCV.inflation.factor, lambda.optimization.tolerance = lambda.optimization.tolerance)
      numnodes = nrow(FEMbasis$mesh$nodes)
    }else if(class(FEMbasis$mesh) == 'mesh.2.5D')
    {
      bigsol = NULL
      if(!is.null(locations))
        stop("The option locations!=NULL for manifold domains is currently not implemented")
      bigsol = CPP_smooth.manifold.GAM.FEM.basis(locations = locations, observations = observations, FEMbasis = FEMbasis,
        covariates = covariates, ndim = ndim, mydim = mydim, BC = BC,
        incidence_matrix = incidence_matrix, areal.data.avg = areal.data.avg,
        FAMILY = family, mu0 = mu0, max.steps.FPIRLS = max.steps.FPIRLS, scale.param = scale.param, threshold.FPIRLS = threshold.FPIRLS,
        search = search, bary.locations = bary.locations,
        optim = optim, lambda = lambda, DOF.stochastic.realizations = DOF.stochastic.realizations, DOF.stochastic.seed = DOF.stochastic.seed,
        DOF.matrix = DOF.matrix, GCV.inflation.factor = GCV.inflation.factor, lambda.optimization.tolerance = lambda.optimization.tolerance)
      numnodes = nrow(FEMbasis$mesh$nodes)
    }else if(class(FEMbasis$mesh) == 'mesh.3D')
    {
      bigsol = NULL
      bigsol = CPP_smooth.volume.GAM.FEM.basis(locations = locations, observations = observations, FEMbasis = FEMbasis,
        covariates = covariates, ndim = ndim, mydim = mydim, BC = BC,
        incidence_matrix = incidence_matrix, areal.data.avg = areal.data.avg,
        FAMILY = family, mu0 = mu0, max.steps.FPIRLS = max.steps.FPIRLS, scale.param = scale.param, threshold.FPIRLS = threshold.FPIRLS,
        search = search, bary.locations = bary.locations,
        optim = optim, lambda = lambda, DOF.stochastic.realizations = DOF.stochastic.realizations, DOF.stochastic.seed = DOF.stochastic.seed,
        DOF.matrix = DOF.matrix, GCV.inflation.factor = GCV.inflation.factor, lambda.optimization.tolerance = lambda.optimization.tolerance)
      numnodes = nrow(FEMbasis$mesh$nodes)
    }
  }

  # ---------- Solution -----------
  if(family != 'gaussian')
  {
    f = bigsol[[1]][1:numnodes,]
    g = bigsol[[1]][(numnodes+1):(2*numnodes),]

    dof = bigsol[[2]]
    GCV_ = bigsol[[3]]
    bestlambda = bigsol[[4]]+1

    if(!is.null(covariates))
    {
      beta = matrix(data=bigsol[[5]],nrow=ncol(covariates),ncol=length(lambda))
  	}
    else
    {
      beta = NULL
    }

    # Save information of Tree Mesh
    tree_mesh = list(
    treelev = bigsol[[6]][1],
    header_orig= bigsol[[7]],
    header_scale = bigsol[[8]],
    node_id = bigsol[[9]][,1],
    node_left_child = bigsol[[9]][,2],
    node_right_child = bigsol[[9]][,3],
    node_box= bigsol[[10]])

    # Reconstruct FEMbasis with tree mesh
    mesh.class= class(FEMbasis$mesh)
    if (is.null(FEMbasis$mesh$treelev))
    { #if doesn't exist the tree information
      FEMbasis$mesh = append(FEMbasis$mesh, tree_mesh)
    } #if already exist the tree information, don't append
    class(FEMbasis$mesh) = mesh.class

    # Save information of Barycenter
    if (is.null(bary.locations))
    {
        bary.locations = list(locations=locations, element_ids = bigsol[[11]], barycenters = bigsol[[12]])
    }
    class(bary.locations) = "bary.locations"

    # Make Functional objects object
    fit.FEM  = FEM(f, FEMbasis)
    PDEmisfit.FEM = FEM(g, FEMbasis)

    # Prepare return list
    reslist = NULL

    if(optim[3]==1)
    {
    	if(bestlambda == 1 || bestlambda == length(lambda))
    		warning("Your optimal 'GCV' is on the border of lambda sequence")
      stderr=sqrt(GCV_*(sum(!is.na(observations))-dof)/sum(!is.na(observations)))
      reslist=list(fit.FEM = fit.FEM, PDEmisfit.FEM = PDEmisfit.FEM,
              beta = beta, edf = dof, GCV = GCV_, stderr=stderr, bestlambda = bestlambda, bary.locations = bary.locations)
    }
    else
    {
      reslist=list(fit.FEM = fit.FEM, PDEmisfit.FEM = PDEmisfit.FEM, beta = beta, bary.locations = bary.locations)
    }

    # GAM outputs
    if(sum(family==c("binomial", "exponential", "gamma", "poisson")) == 1)
    {
      fn.eval = bigsol[[13]]
      J_minima = bigsol[[14]]
      variance.est=bigsol[[15]]
      if( variance.est[1]<0 ) variance.est = NULL
      reslist = c(reslist, list(fn.eval = fn.eval, J_minima = J_minima, variance.est = variance.est))
    }

    return(reslist)
  }
  else
  {
    if(!is.null(covariates))
    {
      if(optim[1]==0 & is.null(DOF.matrix) & optim[3]==0)
      {
        beta = matrix(data=bigsol[[15]],nrow=ncol(covariates),ncol=length(lambda))
      }
      else
      {
        beta = matrix(data=bigsol[[15]],nrow=ncol(covariates),ncol=1)
      }
    }
    else
    {
      beta = NULL
    }

    bestlambda=bigsol[[6]]
    if(optim[1]==0 & (bestlambda == 1 || bestlambda == length(lambda)))
            warning("Your optimal 'GCV' is on the border of lambda sequence")


    if (is.null(lambda.selection.lossfunction))
       { sd = -1 }
    else
       { sd = sqrt(bigsol[[4]])}

    solution = list(
      f = bigsol[[1]][1:numnodes,],
      g = bigsol[[1]][(numnodes+1):(2*numnodes),],
      z_hat = bigsol[[2]],
      beta = beta,
      rmse = bigsol[[3]],
      estimated_sd=sd
    )
    term = bigsol[[9]]
    ot = bigsol[[10]]

    if(term == 1)
    {
      termination = "reached tolerance"
    }else if(term == 2){
      termination = "reached max number iterations"
    }else{
      termination = "uninformative"
    }

    if(ot == 0)
    {
      optimization_type = "full optimization"
    }else if(ot == 1){
      optimization_type = "full DOF grid"
    }else{
      optimization_type = "uninformative"
    }
    
    optimization = list(
      lambda_solution = bigsol[[5]],
      lambda_position = bestlambda,
      GCV = bigsol[[7]],
      optimization_details = list(
          iterations = bigsol[[8]],
          termination = termination,
          optimization_type = optimization_type),
      dof = bigsol[[11]],
      lambda_vector = bigsol[[12]],
      GCV_vector = bigsol[[13]]
    )

    time = bigsol[[14]]

    # Save information of Tree Mesh
    tree_mesh = list(
      treelev = bigsol[[16]][1],
      header_orig= bigsol[[17]],
      header_scale = bigsol[[18]],
      node_id = bigsol[[19]][,1],
      node_left_child = bigsol[[19]][,2],
      node_right_child = bigsol[[19]][,3],
      node_box= bigsol[[20]])

    #Reconstruct FEMbasis with tree mesh
    mesh.class= class(FEMbasis$mesh)
    if (is.null(FEMbasis$mesh$treelev))
    { #if doesn't exist the tree information
      FEMbasis$mesh = append(FEMbasis$mesh, tree_mesh)
    } #if already exist the tree information, don't append
    class(FEMbasis$mesh) = mesh.class

    # Save information of Barycenter
    if (is.null(bary.locations))
    {
      bary.locations = list(locations=locations, element_ids = bigsol[[21]], barycenters = bigsol[[22]])
    }
    class(bary.locations) = "bary.locations"

    # Make Functional objects object
    fit.FEM  = FEM(solution$f, FEMbasis)
    PDEmisfit.FEM = FEM(solution$g, FEMbasis)

    reslist = list(fit.FEM = fit.FEM, PDEmisfit.FEM = PDEmisfit.FEM, solution = solution,
                optimization  = optimization, time = time, bary.locations = bary.locations)
    return(reslist)
  }
}
