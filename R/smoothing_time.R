#' @useDynLib fdaPDE
#' @import Matrix plot3D rgl plot3Drgl geometry
#' @importFrom grDevices heat.colors palette
#' @importFrom graphics plot segments points lines
NULL

#' Space-time regression with differential regularization

#' @param locations A matrix where each row specifies the spatial coordinates \code{x} and \code{y} (and \code{z} if ndim=3) of the corresponding observations in the vector \code{observations}.
#' This parameter can be \code{NULL}. In this case, if also the incidence matrix is \code{NULL} the spatial coordinates are assumed to coincide with the nodes of the \code{mesh}.
#' @param time_locations A vector containing the times of the corresponding observations in the vector \code{observations}. 
#' This parameter can be \code{NULL}. In this case the temporal locations are assumed to coincide with the nodes of the \code{time_mesh}.
#'
#' @param observations A matrix of #locations x #time_locations with the observed data values over the spatio-temporal domain.
#' The spatial locations of the observations can be specified with the \code{locations} argument.
#'
#' @param FEMbasis A \code{FEMbasis} object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param time_mesh A vector specifying the time mesh.
#' @param lambdaS A scalar or vector of smoothing parameters.
#' @param lambdaT A scalar or vector of smoothing parameters.
#' @param covariates A #observations-by-#covariates matrix where each row represents the covariates associated with the corresponding observed data value in \code{observations}.
#' @param PDE_parameters A list specifying the parameters of the PDE in the regularizing term. Default is NULL, i.e. regularization is by means of the Laplacian (stationary, isotropic case).
#'  If the PDE is elliptic it must contain: \code{K}, a 2-by-2 matrix of diffusion coefficients. This induces an anisotropic
#' smoothing with a preferential direction that corresponds to the first eigenvector of the diffusion matrix K; \code{b}, a vector of length 2 of advection coefficients. This induces a
#' smoothing only in the direction specified by the vector \code{b}; \code{c}, a scalar reaction coefficient. \code{c} induces a shrinkage of the surface to zero
#' If the PDE is space-varying it must contain: \code{K}, a function that for each spatial location in the spatial domain
#' (indicated by the vector of the 2 spatial coordinates) returns a 2-by-2 matrix of diffusion coefficients. This induces an anisotropic
#' smoothing with a local preferential direction that corresponds to the first eigenvector of the diffusion matrix K.The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
#' an array with dimensions 2-by-2-by-#points.\code{b}, a function that for each spatial location in the spatial domain returns
#' a vector of length 2 of transport coefficients. This induces a local smoothing only in the direction specified by the vector \code{b}. The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
#' a matrix with dimensions 2-by-#points; \code{c}, a function that for each spatial location in the spatial domain  returns a scalar reaction coefficient.
#' \code{c} induces a shrinkage of the surface to zero. The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
#' a vector with length #points; \code{u}, a function that for each spatial location in the spatial domain  returns a scalar reaction coefficient.
#' \code{u} induces a reaction effect. The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
#' a vector with length #points.
#' For 2.5D and 3D only the Laplacian is available (\code{PDE_parameters=NULL})

#' @param incidence_matrix A #regions-by-#triangles/tetrahedrons matrix where the element (i,j) equals 1 if the j-th triangle/tetrahedron is in the i-th region and 0 otherwise.
#' This is only for areal data. In case of pointwise data, this parameter is set to \code{NULL}.
#' @param BC A list with two vectors:
#'  \code{BC_indices}, a vector with the indices in \code{nodes} of boundary nodes where a Dirichlet Boundary Condition should be applied;
#'  \code{BC_values}, a vector with the values that the spatial field must take at the nodes indicated in \code{BC_indices}.
#' @param FLAG_MASS Boolean. This parameter is considerd only for separable problems i.e. when \code{FLAG_PARABOLIC==FALSE}. If \code{TRUE} the mass matrix in space and in time are used, if \code{FALSE} they are substituted with proper identity matrices.
#' @param FLAG_PARABOLIC Boolean. If \code{TRUE} the parabolic problem problem is selected, if \code{FALSE} the separable one.
#' @param IC Initial condition needed in case of parabolic problem i.e. when \code{FLAG_PARABOLIC==FALSE}.This parameter has to be set only for parabolic. If \code{FLAG_PARABOLIC=TRUE} and \code{IC=NULL} it is necessary to provide
#' also data at the initial time. IC will be estimated from it
#' @param GCV Boolean. If \code{TRUE} the following quantities are computed: the trace of the smoothing matrix, the estimated error standard deviation,  and
#'        the Generalized Cross Validation criterion, for each combination of the smoothing parameters specified in \code{lambdaS} and \code{lambdaT}.
#' @param GCVmethod either "Exact" or "Stochastic". If set to "Exact" perform an exact (but possibly slow) computation of the GCV index. If set to "Stochastic" approximate the GCV with a stochastic algorithm.
#'        This parameter is considered only when \code{GCV=TRUE}
#' @param nrealizations Number of random points used in the stochastic computation of the dofs
#' @param DOF_matrix Matrix of degrees of freedom. This parameter can be used if the DOF_matrix corresponding to \code{lambdaS} and \code{lambdaT} is available from precedent computation. This allows to save time
#' since the computation of the dof is the most expensive part of GCV.
#' @param search a flag to decide the search algorithm type (tree or naive or walking search algorithm).
#' @param bary.locations A list with three vectors:
#'  \code{locations}, location points which are same as the given locations options. (checks whether both locations are the same);
#'  \code{element ids}, a vector of element id of the points from the mesh where they are located;
#'  \code{barycenters}, a vector of barycenter of points from the located element.
#' @param GCV.inflation.factor Tuning parameter used for the estimation of GCV. Default value \code{GCV.inflation.factor = 1.8}.
#' It is advised to set it grather than 1 to avoid overfitting.
#' @param areal.data.avg Boolean. It involve the computation of Areal Data. If \code{TRUE} the areal data are averaged, otherwise not.
#' @return A list with the following variables:
#' \item{\code{fit.FEM.time}}{A \code{FEM.time} object that represents the fitted spatio-temporal field.}
#' \item{\code{PDEmisfit.FEM.time}}{A \code{FEM.time} object that represents the misfit of the penalized PDE.}
#' \item{\code{beta}}{If \code{covariates} is not \code{NULL}, a matrix with number of rows equal to the number of covariates and numer of columns equal to length of lambda.  The \code{j}th column represents the vector of regression coefficients when
#' the smoothing parameter is equal to \code{lambda[j]}.}
#' \item{\code{edf}}{If GCV is \code{TRUE}, a scalar or matrix with the trace of the smoothing matrix for each combination of the smoothing parameters specified in \code{lambdaS} and \code{lambdaT}.}
#' \item{\code{stderr}}{If GCV is \code{TRUE}, a scalar or matrix with the estimate of the standard deviation of the error for each combination of the smoothing parameters specified in \code{lambdaS} and \code{lambdaT}.}
#' \item{\code{GCV}}{If GCV is \code{TRUE}, a  scalar or matrix with the value of the GCV criterion for each combination of the smoothing parameters specified in \code{lambdaS} and \code{lambdaT}.}
#' \item{\code{bestlambda}}{If GCV is \code{TRUE}, a 2-elements vector with the indices of smoothing parameters returnig the lowest GCV}
#' \item{\code{ICestimated}}{If FLAG_PARABOLIC is \code{TRUE} and IC is \code{NULL}, a list containing a \code{FEM} object with the initial conditions, the value of the smoothing parameter lambda returning the lowest GCV and, in presence of covariates, the estimated beta coefficients}
#' \item{\code{bary.locations}}{A barycenter information of the given locations if the locations are not mesh nodes.}
#' @description Space-time regression  with differential regularization. Space-varying covariates can be included in the model. The technique accurately handle data distributed over irregularly shaped domains. Moreover, various conditions can be imposed at the domain boundaries.
#' @usage smooth.FEM.time(locations = NULL, time_locations=NULL, observations, FEMbasis, 
#'                       time_mesh=NULL, lambdaS, lambdaT = 1, covariates = NULL, 
#'                       PDE_parameters=NULL, incidence_matrix = NULL, BC = NULL, FLAG_MASS = FALSE, 
#'                       FLAG_PARABOLIC = FALSE, IC = NULL, GCV = FALSE,
#'                       GCVmethod = "Stochastic", nrealizations = 100, DOF_matrix=NULL, 
#'                       search = "tree", bary.locations = NULL, GCV.inflation.factor= 1, areal.data.avg = TRUE)
#' @export

#' @references Sangalli, L.M., Ramsay, J.O. & Ramsay, T.O., 2013. Spatial spline regression models. Journal of the Royal Statistical Society. Series B: Statistical Methodology, 75(4), pp. 681-703.
#' @examples
#' library(fdaPDE)
#' 
#' data(horseshoe2D)
#' boundary_nodes = horseshoe2D$boundary_nodes
#' boundary_segments = horseshoe2D$boundary_segments
#' locations = horseshoe2D$locations
#' time_locations = seq(0,1,length.out = 5)
#' 
#' mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)
#' 
#' space_time_locations = cbind(rep(time_locations,each=nrow(mesh$nodes)),
#'                              rep(mesh$nodes[,1],5),rep(mesh$nodes[,2],5))
#' 
#' FEMbasis = create.FEM.basis(mesh)
#' lambdaS = 10^-1
#' lambdaT = 10^-1
#' data = fs.test(space_time_locations[,2], 
#'                space_time_locations[,3])*cos(pi*space_time_locations[,1]) +
#'        rnorm(nrow(space_time_locations), sd = 0.5)
#' data = matrix(data, nrow = nrow(mesh$nodes), ncol = length(time_locations), byrow = TRUE)
#' 
#' solution = smooth.FEM.time(observations = data, time_locations = time_locations,
#'                            FEMbasis = FEMbasis, lambdaS = lambdaS, lambdaT = lambdaT)
#' plot(solution$fit.FEM)

smooth.FEM.time<-function(locations = NULL, time_locations=NULL, observations, FEMbasis, time_mesh=NULL, lambdaS, lambdaT = 1, covariates = NULL, PDE_parameters=NULL, incidence_matrix = NULL, BC = NULL, FLAG_MASS = FALSE, FLAG_PARABOLIC = FALSE, IC = NULL, GCV = FALSE, GCVmethod = "Stochastic", nrealizations = 100, DOF_matrix=NULL, search = "tree", bary.locations = NULL, GCV.inflation.factor= 1, areal.data.avg = TRUE)
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

  if(GCVmethod=="Stochastic"){
    GCVMETHOD=2
  }else if(GCVmethod=="Exact"){
    GCVMETHOD=1
  }else{
    stop("GCVmethod must be either Stochastic or Exact")
  }

  if(search=="naive"){
    search=1
  }else if(search=="tree"){
    search=2
  }else{
    stop("search must be either tree or naive.")
  }

  #if locations is null but bary.locations is not null, use the locations in bary.locations
  if(is.null(locations) & !is.null(bary.locations)) {
    locations = bary.locations$locations
    locations = as.matrix(locations)
  }
  
  DOF=TRUE
  if(!is.null(DOF_matrix))
      DOF=FALSE

  space_varying=checkSmoothingParameters_time(locations, time_locations, observations, FEMbasis, time_mesh, lambdaS, lambdaT, covariates, PDE_parameters, incidence_matrix, BC, FLAG_MASS, FLAG_PARABOLIC, IC, GCV, GCVMETHOD, nrealizations, DOF, DOF_matrix, search=search, bary.locations=bary.locations, GCV.inflation.factor = GCV.inflation.factor, areal.data.avg = areal.data.avg)

  ## Coverting to format for internal usage
  if(!is.null(locations))
    locations = as.matrix(locations)
  if(!is.null(time_locations))
    time_locations = as.matrix(time_locations)
  if(!is.null(time_mesh))
    time_mesh = as.matrix(time_mesh)
  observations = as.matrix(observations)
  lambdaS = as.matrix(lambdaS)
  lambdaT = as.matrix(lambdaT)
  if(!is.null(covariates))
    covariates = as.matrix(covariates)
  if(!is.null(DOF_matrix))
    DOF_matrix = as.matrix(DOF_matrix)
  if(!is.null(incidence_matrix))
    incidence_matrix = as.matrix(incidence_matrix)
  if(!is.null(IC))
    IC = as.matrix(IC)
  if(!is.null(BC))
  {
    BC$BC_indices = as.matrix(BC$BC_indices)
    BC$BC_values = as.matrix(BC$BC_values)
  }

  # if I have PDE non-sv case I need (constant) matrices as parameters

  if(!is.null(PDE_parameters) & space_varying==FALSE)
  {
    PDE_parameters$K = as.matrix(PDE_parameters$K)
    PDE_parameters$b = as.matrix(PDE_parameters$b)
    PDE_parameters$c = as.matrix(PDE_parameters$c)
  }


  checkSmoothingParametersSize_time(locations, time_locations, observations, FEMbasis, time_mesh, lambdaS, lambdaT, covariates, PDE_parameters, incidence_matrix, BC, FLAG_MASS, FLAG_PARABOLIC, IC, GCV, DOF, DOF_matrix, space_varying, ndim, mydim)
  observations<-as.vector(observations)

  if(is.null(time_locations))
  {
    if(FLAG_PARABOLIC && !is.null(IC))
      time_locations <- time_mesh[2:length(time_mesh)]
    else
      time_locations <- time_mesh
  }

  if(is.null(time_mesh))
  {
    if(FLAG_PARABOLIC && !is.null(IC))
      time_mesh <- rbind(2*time_locations[1]-time_locations[2],time_locations)
    else
      time_mesh<-time_locations
  }

  # Check whether the locations coincide with the mesh nodes (should be put after all the validations)
  if (!is.null(locations)) {
    if(dim(locations)[1]==dim(FEMbasis$mesh$nodes)[1] & dim(locations)[2]==dim(FEMbasis$mesh$nodes)[2]) {
      sum1=0
      sum2=0
      for (i in 1:nrow(locations)) {
      sum1 = sum1 + abs(locations[i,1]-FEMbasis$mesh$nodes[i,1])
      sum2 = sum2 + abs(locations[i,2]-FEMbasis$mesh$nodes[i,2])
      }
      if (sum1==0 & sum2==0) {
        message("No search algorithm is used because the locations coincide with the nodes.")
        locations = NULL #In principle, R uses pass-by-value semantics in its function calls. So put ouside of checkSmoothingParameters function.
      }
    }
  }

  ################## End checking parameters, sizes and conversion #############################
  if(class(FEMbasis$mesh) == 'mesh.2D' & is.null(PDE_parameters)){

    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.time(locations=locations, time_locations=time_locations, observations=observations, FEMbasis=FEMbasis,
                                  time_mesh=time_mesh, lambdaS=lambdaS, lambdaT=lambdaT, covariates=covariates, incidence_matrix=incidence_matrix,
                                  ndim=ndim, mydim=mydim, BC=BC, FLAG_MASS=FLAG_MASS, FLAG_PARABOLIC=FLAG_PARABOLIC, IC=IC, GCV=GCV,
                                  GCVMETHOD=GCVMETHOD, nrealizations=nrealizations,DOF=DOF,DOF_matrix=DOF_matrix, search=search, bary.locations=bary.locations, GCV.inflation.factor = GCV.inflation.factor, areal.data.avg = areal.data.avg)

  } else if(class(FEMbasis$mesh) == 'mesh.2D' & !is.null(PDE_parameters) & space_varying==FALSE){

    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.PDE.time(locations=locations, time_locations=time_locations, observations=observations, FEMbasis=FEMbasis,
                                      time_mesh=time_mesh, lambdaS=lambdaS, lambdaT=lambdaT, PDE_parameters=PDE_parameters,
                                      covariates=covariates, incidence_matrix=incidence_matrix,
                                      ndim=ndim, mydim=mydim, BC=BC, FLAG_MASS=FLAG_MASS, FLAG_PARABOLIC=FLAG_PARABOLIC, IC=IC, GCV=GCV,
                                      GCVMETHOD=GCVMETHOD, nrealizations=nrealizations,DOF=DOF,DOF_matrix=DOF_matrix, search=search, bary.locations=bary.locations, GCV.inflation.factor = GCV.inflation.factor, areal.data.avg = areal.data.avg)

  } else if(class(FEMbasis$mesh) == 'mesh.2D' & !is.null(PDE_parameters) & space_varying==TRUE){

    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.PDE.sv.time(locations=locations, time_locations=time_locations, observations=observations, FEMbasis=FEMbasis,
                                        time_mesh=time_mesh, lambdaS=lambdaS, lambdaT=lambdaT, PDE_parameters=PDE_parameters,
                                        covariates=covariates, incidence_matrix=incidence_matrix,
                                        ndim=ndim, mydim=mydim, BC=BC, FLAG_MASS=FLAG_MASS, FLAG_PARABOLIC=FLAG_PARABOLIC, IC=IC, GCV=GCV,
                                        GCVMETHOD=GCVMETHOD, nrealizations=nrealizations,DOF=DOF,DOF_matrix=DOF_matrix, search=search, bary.locations=bary.locations, GCV.inflation.factor = GCV.inflation.factor, areal.data.avg = areal.data.avg)

  }else if(class(FEMbasis$mesh) == 'mesh.2.5D'){

    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.manifold.FEM.time(locations=locations, time_locations=time_locations, observations=observations, FEMbasis=FEMbasis,
                                          time_mesh=time_mesh, lambdaS=lambdaS, lambdaT=lambdaT, covariates=covariates, incidence_matrix=incidence_matrix,
                                          ndim=ndim, mydim=mydim, BC=BC, FLAG_MASS=FLAG_MASS, FLAG_PARABOLIC=FLAG_PARABOLIC, IC=IC, GCV=GCV,
                                          GCVMETHOD=GCVMETHOD, nrealizations=nrealizations,DOF=DOF,DOF_matrix=DOF_matrix, search=search, bary.locations=bary.locations, GCV.inflation.factor = GCV.inflation.factor, areal.data.avg = areal.data.avg)

  }else if(class(FEMbasis$mesh) == 'mesh.3D'){

    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.volume.FEM.time(locations=locations, time_locations=time_locations, observations=observations, FEMbasis=FEMbasis,
                                        time_mesh=time_mesh, lambdaS=lambdaS, lambdaT=lambdaT, covariates=covariates, incidence_matrix=incidence_matrix,
                                        ndim=ndim, mydim=mydim, BC=BC, FLAG_MASS=FLAG_MASS, FLAG_PARABOLIC=FLAG_PARABOLIC, IC=IC, GCV=GCV,
                                        GCVMETHOD=GCVMETHOD, nrealizations=nrealizations,DOF=DOF,DOF_matrix=DOF_matrix, search=search, bary.locations=bary.locations, GCV.inflation.factor = GCV.inflation.factor, areal.data.avg = areal.data.avg)

  }

  N = nrow(FEMbasis$mesh$nodes)
  M = ifelse(FLAG_PARABOLIC,length(time_mesh)-1,length(time_mesh) + 2);
  if(is.null(IC) && FLAG_PARABOLIC)
    IC = bigsol[[13]]$coeff[,bigsol[[14]]]
  if(FLAG_PARABOLIC)
  {
    f = array(dim=c(length(IC)+M*N,length(lambdaS),length(lambdaT)))
    for (i in 1:length(lambdaS))
     for (j in 1:length(lambdaT))
       f[,i,j] = c(IC,bigsol[[1]][1:(N*M),i+(j-1)*length(lambdaS)])
  }
  else
    f = array(data=bigsol[[1]][1:(N*M),],dim = c(N*M,length(lambdaS),length(lambdaT)))
  if(FLAG_PARABOLIC)
  {
    g = array(dim=c(length(IC)+M*N,length(lambdaS),length(lambdaT)))
    for (i in 1:length(lambdaS))
      for (j in 1:length(lambdaT))
        g[,i,j] = c(rep(0,length(IC)),bigsol[[1]][(N*M+1):(2*N*M),i+(j-1)*length(lambdaS)])
  }
  else
    g = array(data=bigsol[[1]][(N*M+1):(2*N*M),],dim = c(N*M,length(lambdaS),length(lambdaT)))

  dof = bigsol[[2]]
  GCV_ = bigsol[[3]]
  bestlambda = bigsol[[4]]+1
  if(!is.null(covariates))
    beta = array(data=bigsol[[5]],dim=c(ncol(covariates),length(lambdaS),length(lambdaT)))
  else
    beta = NULL

  if(all(is.na(bigsol[[13]])))
    ICestimated = NULL
  else
    ICestimated = list(IC.FEM=bigsol[[13]],bestlambdaindex=bigsol[[14]],bestlambda=bigsol[[15]],beta=bigsol[[16]])

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

  # Make FEM.time objects
  fit.FEM.time  = FEM.time(f, time_mesh, FEMbasis, FLAG_PARABOLIC)
  PDEmisfit.FEM.time = FEM.time(g, time_mesh, FEMbasis, FLAG_PARABOLIC)

  # Prepare return list
  reslist = NULL

  if(GCV == TRUE)
  {
    stderr=sqrt(GCV_*(length(observations)-dof)/length(observations))
    reslist=list(fit.FEM.time = fit.FEM.time, PDEmisfit.FEM.time = PDEmisfit.FEM.time,
            beta = beta, edf = dof, GCV = GCV_, stderr=stderr, bestlambda = bestlambda, ICestimated=ICestimated, bary.locations = bary.locations)
  }else{
    reslist=list(fit.FEM.time = fit.FEM.time, PDEmisfit.FEM.time = PDEmisfit.FEM.time, beta = beta, ICestimated=ICestimated, bary.locations = bary.locations)
  }

  return(reslist)
}
