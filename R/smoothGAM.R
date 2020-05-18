#' Spatial generalize GAM model with differential regularization
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
#' @param lambda A scalar or vector of smoothing parameters.
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
#' @param incidence_matrix A #regions-by-#triangles/tetrahedrons matrix where the element (i,j) equals 1 if the j-th 
#' triangle/tetrahedron is in the i-th region and 0 otherwise.
#' This is needed only for areal data. In case of pointwise data, this parameter is set to \code{NULL}.
#' @param BC A list with two vectors: 
#'  \code{BC_indices}, a vector with the indices in \code{nodes} of boundary nodes where a Dirichlet Boundary Condition should be applied;
#'  \code{BC_values}, a vector with the values that the spatial field must take at the nodes indicated in \code{BC_indices}.
#' @param GCV Boolean. If \code{TRUE} the following quantities are computed: the trace of the smoothing matrix, the estimated error standard deviation,  and 
#'        the Generalized Cross Validation criterion, for each value of the smoothing parameter specified in \code{lambda}.
#' @param GCVmethod This parameter is considered only when \code{GCV=TRUE}. It can be either "Exact" or "Stochastic". 
#' If set to "Exact" the algoritm performs an exact (but possibly slow) computation 
#' of the GCV index. If set to "Stochastic" the GCV is approximated by a stochastic algorithm.
#' @param nrealizations This parameter is considered only when \code{GCV=TRUE} and \code{GCVmethod = "Stochastic"}. 
#' It is a positive integer that represents the number of uniform random variables used in stochastic GCV computation.      
#' @param mu0 This parameter is a vector that set the starting point for FPIRLS algorithm. It represent an initial guess of the location parameter.
#' Default is set to observation for non binary distribution while equal to \code{0.5(observations + 0.5)} for binary data.
#' @param fam This parameter specify the distibution within exponential family used for GLM model.
#' The following distribution are implemented: "binomial", "exponential", "gamma", "poisson", "gaussian", "invgaussian".
#' The default link function for binomial is \code{logit} if you want either \code{probit} or \code{clogloc} set \code{fam = "probit"}, \code{fam = "cloglog"}.      
#' @param scale.param Dispersion parameter of the chosen distribution. This is only required for "gamma", "gaussian", "invgaussian".
#' User may specify the parameter as a positive real number. If the parameter is not supplied, it is estimated from data according to Wilhelm Sangalli 2016. 
#' @param tune Tuning parameter used for the estimation of GCV. Default value \code{tune = 1.8}.
#' It is advised to set it grather than 1 to avoid overfitting.
#' @param threshold This parameter is used for arresting algorithm iterations. Algorithm stops when two successive iterations lead to improvement in penalized log-likelihood smaller than threshold.
#' Default value \code{threshold = 0.0002020}.
#' @param max.steps This parameter is used to limit the maximum number of iteration.
#' Default value \code{max.steps=15}.
#' @return A list with the following variables:
#' \itemize{
#'    \item{\code{fit.FEM}}{A \code{FEM} object that represents the fitted spatial field.}
#'    \item{\code{PDEmisfit.FEM}}{A \code{FEM} object that represents the Laplacian of the estimated spatial field.}
#'    \item{\code{beta_hat}}{ If covariates is not \code{NULL}, a matrix with number of rows equal to the number of covariates and numer of columns equal to length of lambda.  The \code{j}th column represents the vector of regression coefficients when 
#' the smoothing parameter is equal to \code{lambda[j]}.}
#'    \item{\code{fn_hat}}{ A matrix with number of rows equal to number of locations and number of columns equal to length of lambda. Each column contain the evaluaton of the spatial field in the location points.}
#'    \item{\code{J_minima}}{A vector of the same length of lambda, containing the reached minima for each value of the smoothing parameter.}
#'    \item{\code{DoF}}{If GCV is \code{TRUE}, a scalar or vector with the trace of the smoothing matrix for each value of the smoothing parameter specified in \code{lambda}.}
#'    \item{\code{GCV}}{If GCV is \code{TRUE}, a  scalar or vector with the value of the GCV criterion for each value of the smoothing parameter specified in \code{lambda}.}
#' }
#' @description This function implements a spatial GAM model with differential regularization. 
#'  The regularizing term involves a Partial Differential Equation (PDE). In the simplest case the PDE involves only the 
#'  Laplacian of the spatial field, that induces an isotropic smoothing. When prior information about the anisotropy or 
#'  non-stationarity is available the PDE involves a general second order linear differential operator with possibly 
#'  space-varying coefficients. 
#'  The technique accurately handle data distributed over irregularly shaped domains. Moreover, various conditions 
#'  can be imposed at the domain boundaries.
#'  Non linear loglikelihood maximization is performed via Newton method according to FPIRLS algorithm.
#' @usage smooth.GAM.FEM= function(locations = NULL, observations, FEMbasis,
#'  lambda, covariates = NULL, BC = NULL, PDE_parameters=NULL, incidence_matrix = NULL, mu0=NULL,
#'  fam= "invgaussian",scale.param=NULL, GCV = FALSE, nrealizations=100, GCVmethod="Stochastic",
#'  tune=1.8, threshold=0.0002020, max.steps=15 )
#' @export        

#' @references 
#' \itemize{
#'    \item{Sangalli, L. M., Ramsay, J. O., Ramsay, T. O. (2013). Spatial spline regression models. 
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 75(4), 681-703.}
#'    \item{Matthieu Wilhelm & Laura M. Sangalli (2016). Generalized spatial regression with differential regularization. 
#'  Journal of Statistical Computation and Simulation, 86:13, 2497-2518.}
#' }
#' @examples
#' library(fdaPDE)

smooth.GAM.FEM = function(locations = NULL, observations, FEMbasis,
	lambda, covariates = NULL, BC = NULL, PDE_parameters=NULL, incidence_matrix = NULL, mu0=NULL,
	fam=c("binomial", "probit","cloglog", "exponential", "gamma", "poisson", "gaussian", "invgaussian"),
  scale.param=NULL, GCV = FALSE, nrealizations=100, GCVmethod="Stochastic",
	tune=1.8, threshold=0.0002020, max.steps=15, DOF_matrix=NULL, search = "tree", bary.locations = NULL)
{

  if(class(FEMbasis$mesh) == "mesh.2D"){
    ndim = 2
    mydim = 2
  }else if(class(FEMbasis$mesh) == "mesh.2.5D"){
    ndim = 3
    mydim = 2
  }else if(class(FEMbasis$mesh) == "mesh.3D"){
    ndim = 3
    mydim = 3
  }else{
    stop('Unknown mesh class')
  }
  ##################### Checking parameters, sizes and conversion ################################

  if(GCVmethod=="Stochastic")
    GCVMETHOD=2
  else if(GCVmethod=="Exact")
    GCVMETHOD=1
  else{
    stop("GCVmethod must be either Stochastic or Exact")
  }

  if(search=="naive")
    search=1
  else if(search=="tree")
    search=2
  else{
    stop("search must be either tree or naive.")
  }
  
  #!!!!! FAMILY CHECK  !!!!!
  family_admit = c("binomial", "probit", "cloglog", "exponential", "gamma", "poisson", "gaussian", "invgaussian")
  if(sum(fam==family_admit)==0 ){
     stop("'family' parameter required.\nCheck if it is one of the following: binomial, probit, cloglog, exponential, gamma, poisson, gaussian")
  }
  
  DOF=TRUE
  if(!is.null(DOF_matrix))
    DOF=FALSE
  
  #if locations is null but bary.locations is not null, use the locations in bary.locations
  if(is.null(locations) & !is.null(bary.locations)) {
    locations = bary.locations$locations
    locations = as.matrix(locations)
  }

  ## Converting to format for internal usage
  if(!is.null(locations))
    locations = as.matrix(locations)
    observations = as.matrix(observations)
    lambda = as.matrix(lambda)
  if(!is.null(covariates))
    covariates = as.matrix(covariates)
  if(!is.null(DOF_matrix))
    DOF_matrix = as.matrix(DOF_matrix)
  if(!is.null(incidence_matrix))
    incidence_matrix = as.matrix(incidence_matrix)
  if(!is.null(BC))
  {
    BC$BC_indices = as.matrix(BC$BC_indices)
    BC$BC_values = as.matrix(BC$BC_values)
  }
  
  space_varying=checkSmoothingParameters(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda, covariates=covariates, incidence_matrix=incidence_matrix, 
    BC=BC, GCV=GCV, PDE_parameters=PDE_parameters, GCVmethod=GCVMETHOD , nrealizations=nrealizations, search=search, bary.locations=bary.locations)

  # if I have PDE non-sv case I need (constant) matrices as parameters
  
  if(!is.null(PDE_parameters) & space_varying==FALSE) 
  {
    PDE_parameters$K = as.matrix(PDE_parameters$K)
    PDE_parameters$b = as.matrix(PDE_parameters$b)
    PDE_parameters$c = as.matrix(PDE_parameters$c)
  }

  checkSmoothingParametersSize(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda, covariates=covariates, incidence_matrix=incidence_matrix, 
    BC=BC, GCV=GCV, space_varying=space_varying, PDE_parameters=PDE_parameters, ndim=ndim, mydim=mydim)

  #checkGAMParameters(max.steps = max.steps, mu0 = mu0, observations.len = length(observations), scale.param = scale.param, tune = tune, threshold = threshold, fam = fam)
  
  ################## End checking parameters, sizes and conversion #############################

  if(class(FEMbasis$mesh) == 'mesh.2D' & is.null(PDE_parameters)){

    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.GAM.FEM(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda,
            covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
            BC=BC, GCV=GCV, GCVMETHOD=GCVMETHOD, nrealizations=nrealizations, FAMILY=fam,
            mu0 = mu0, max.steps=max.steps,
            scale.param=scale.param, tune=tune, threshold=threshold , DOF =DOF, DOF_matrix = DOF_matrix, search = search, bary.locations = bary.locations)

    numnodes = nrow(FEMbasis$mesh$nodes)

    } else if(class(FEMbasis$mesh) == 'mesh.2D' & !is.null(PDE_parameters) & space_varying==FALSE){
        bigsol = NULL
        print('C++ Code Execution')
        bigsol = CPP_smooth.GAM.FEM.PDE.basis(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda,
            PDE_parameters = PDE_parameters, 
            covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
            BC=BC, GCV=GCV, GCVMETHOD=GCVMETHOD, nrealizations=nrealizations, FAMILY=fam,
            mu0 = mu0, max.steps=max.steps,
            scale.param=scale.param, tune=tune,  threshold=threshold,  DOF =DOF, DOF_matrix = DOF_matrix, search = search, bary.locations = bary.locations )

        numnodes = nrow(FEMbasis$mesh$nodes)
   
   } else if(class(FEMbasis$mesh) == 'mesh.2D' & !is.null(PDE_parameters) & space_varying==TRUE){ 
    
    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.GAM.FEM.PDE.sv.basis(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda,
            PDE_parameters = PDE_parameters, 
            covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
            BC=BC, GCV=GCV, GCVMETHOD=GCVMETHOD, nrealizations=nrealizations, FAMILY=fam,
            mu0 = mu0, max.steps=max.steps,
            scale.param=scale.param, tune=tune, threshold=threshold,  DOF =DOF, DOF_matrix = DOF_matrix, search = search, bary.locations = bary.locations )
  
    numnodes = nrow(FEMbasis$mesh$nodes)
  } else if(class(FEMbasis$mesh) == 'mesh.2.5D'){
    
    bigsol = NULL  
    print('C++ Code Execution')
    if(!is.null(locations))
      stop("The option locations!=NULL for manifold domains is currently not implemented")
    bigsol = CPP_smooth.manifold.GAM.FEM.basis(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda,
            covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
            BC=BC, GCV=GCV, GCVMETHOD=GCVMETHOD, nrealizations=nrealizations, FAMILY=fam,
            mu0 = mu0, max.steps=max.steps, 
            scale.param=scale.param, tune=tune, threshold=threshold, DOF =DOF, DOF_matrix = DOF_matrix, search = search, bary.locations = bary.locations)
    
    numnodes = FEMbasis$mesh$nnodes
    
  }else if(class(FEMbasis$mesh) == 'mesh.3D'){
      
    bigsol = NULL  
    print('C++ Code Execution')
    bigsol = CPP_smooth.volume.GAM.FEM.basis(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda,
            covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
            BC=BC, GCV=GCV, GCVMETHOD=GCVMETHOD, nrealizations=nrealizations, FAMILY=fam,
            mu0 = mu0, max.steps=max.steps,
            scale.param=scale.param, tune=tune, threshold=threshold, DOF =DOF, DOF_matrix = DOF_matrix, search = search, bary.locations = bary.locations )
    
    numnodes = FEMbasis$mesh$nnodes
  }
    
  # Estimated functions (smooth fields)
  f = bigsol[[1]][1:numnodes,]
  g = bigsol[[1]][(numnodes+1):(2*numnodes),]# laplacian(f)
  # DoF estimate
  DoF = bigsol[[2]]
  computedGCV = bigsol[[3]]
  # Make Functional objects object
  fit.FEM  = FEM(f, FEMbasis)
  PDEmisfit.FEM = FEM(g, FEMbasis)
  
  bestlambda = bigsol[[4]]+1
 
  if(!is.null(covariates))
    beta_hat = matrix(data=bigsol[[5]],nrow=ncol(covariates),ncol=length(lambda))
  else
    beta_hat = NULL

  # Save information of Tree Mesh
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
  if (is.null(FEMbasis$mesh$treelev)) { #if doesn't exist the tree information
    FEMbasis$mesh = append(FEMbasis$mesh, tree_mesh)
  } #if already exist the tree information, don't append
  class(FEMbasis$mesh) = mesh.class  

  # Save information of Barycenter
  if (is.null(bary.locations)) {
      bary.locations = list(locations=locations, element_ids = bigsol[[11]], barycenters = bigsol[[12]])    
  }
  class(bary.locations) = "bary.locations"

  fn_hat = bigsol[[13]]
  J_minima = bigsol[[14]]

  scale.param.est = bigsol[[15]]

  reslist = NULL  
  # check if covariates is used
  if(!is.null(covariates)){
    reslist=list(fit.FEM = fit.FEM, PDEmisfit.FEM = PDEmisfit.FEM, beta_hat = beta_hat)
  }else{
    reslist=list(fit.FEM = fit.FEM, PDEmisfit.FEM = PDEmisfit.FEM)
  }

  # check if location is used
  if( !is.null(locations) ){
    reslist=c(reslist,list(fn_hat = fn_hat))
  }
  # check if GCV method is used  
  if(GCV == TRUE){
    reslist=c(reslist,list( J_minima = J_minima, DoF= DoF ,GCV = computedGCV), bestlambda = bestlambda)
    }else{
      reslist=c(reslist,list( J_minima = J_minima))
   }

   reslist=c(reslist, list( bary.locations = bary.locations))

  return(reslist)

  }
