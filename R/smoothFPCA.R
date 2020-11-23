#' Smooth Functional Principal Component Analysis
#'
#' @param datamatrix A matrix of dimensions #samples-by-#locations with the observed data values over the domain
#' for each sample. The datamatrix needs to have zero mean.
#' If the \code{locations} argument is left \code{NULL} the datamatrix has to be dimensions #samples-by-#nodes where #nodes
#' is the number of nodes of the mesh in the FEMbasis. In this case, each observation is associated to the corresponding
#' node in the mesh.
#' If the data are observed only on a subset of the mesh nodes, fill with \code{NA} the values of the
#' \code{datamatrix} in correspondence of unobserved data.
#' @param locations A #observations-by-2 matrix in the 2D case and #observations-by-3 matrix in the 2.5D and 3D case, where
#' each row specifies the spatial coordinates \code{x} and \code{y} (and \code{z} in 2.5D and 3D) of the corresponding
#' observation in the \code{datamatrix}.
#' If the locations of the observations coincide with (or are a subset of) the nodes of the mesh in the \code{FEMbasis},
#' leave the parameter \code{locations = NULL} for a faster implementation.
#' @param FEMbasis A \code{FEMbasis} object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param lambda A scalar or vector of smoothing parameters.
#' @param nPC An integer specifying the number of Principal Components to compute.
#' @param validation A string specifying the type of validation to perform. If \code{lambda} is a vector, it has to
#' be specified as \code{"GCV"} or \code{"KFold"}. This parameter specify which method of cross-validation is used
#' to select the best parameter \code{lambda} among those values of the smoothing parameter specified in \code{lambda}
#' for each Principal Component.
#' @param NFolds This parameter is used only in case \code{validation = "KFold"}. It is an integer specifying
#' the number of folds to use if the KFold cross-validation method for the
#' selection of the best parameter \code{lambda} is chosen. Default value is 5.
#' @param GCVmethod This parameter is considered only when \code{validation = "GCV"}. It can be either "Exact" or
#' "Stochastic". If set to "Exact" the algoritm performs an exact (but possibly slow) computation
#' of the GCV index. If set to "Stochastic" the GCV is approximated by a stochastic algorithm.
#' @param nrealizations The number of realizations to be used in the stochastic algorithm for the estimation of GCV.
#' @param search a flag to decide the search algorithm type (tree or naive or walking search algorithm).
#' @param bary.locations A list with three vectors:
#'  \code{locations}, location points which are same as the given locations options. (checks whether both locations are the same);
#'  \code{element ids}, a vector of element id of the points from the mesh where they are located;
#'  \code{barycenters}, a vector of barycenter of points from the located element.
#' @return A list with the following variables:
#' \itemize{
#' \item{\code{loadings.FEM}}{A \code{FEM} object that represents the L^2-normalized functional loadings for each
#' Principal Component computed.}
#' \item{\code{scores}}{A #samples-by-#PrincipalComponents matrix that represents the unnormalized scores or PC vectors.}
#' \item{\code{lambda}}{A vector of length #PrincipalComponents with the values of the smoothing parameter \code{lambda}
#' chosen for that Principal Component.}
#' \item{\code{variance_explained}}{A vector of length #PrincipalComponents where each value represent the variance explained by that component.}
#' \item{\code{cumsum_percentage}}{A vector of length #PrincipalComponents containing the cumulative percentage of the variance explained by the first components.}
#' \item{\code{bary.locations}}{A barycenter information of the given locations if the locations are not mesh nodes.}
#' }
#' @description This function implements a smooth functional principal component analysis over a planar mesh,
#' a smooth manifold or a volume.
#' @usage FPCA.FEM(locations = NULL, datamatrix, FEMbasis, lambda, nPC = 1, validation = NULL,
#'                 NFolds = 5,GCVmethod = "Stochastic", nrealizations = 100, search = "tree",
#'                 bary.locations = NULL)
#' @references Lila, E., Aston, J.A.D.,  Sangalli, L.M., 2016a. Smooth Principal Component Analysis over two-dimensional
#' manifolds with an application to neuroimaging. Ann. Appl. Stat., 10(4), pp. 1854-1879.
#' @export
#' @examples
#' library(fdaPDE)
#'
#' ## Load the hub data
#' data(hub2.5D)
#' hub2.5D.nodes = hub2.5D$hub2.5D.nodes
#' hub2.5D.triangles = hub2.5D$hub2.5D.triangles
#'
#' mesh = create.mesh.2.5D(nodes = hub2.5D.nodes, triangles = hub2.5D.triangles)
#' ## Create the Finite Element basis
#' FEMbasis = create.FEM.basis(mesh)
#' ## Create a datamatrix
#' datamatrix = NULL
#' for(ii in 1:50){
#'   a1 = rnorm(1, mean = 1, sd = 1)
#'   a2 = rnorm(1, mean = 1, sd = 1)
#'   a3 = rnorm(1, mean = 1, sd = 1)
#'
#'   func_evaluation = numeric(nrow(mesh$nodes))
#'   for (i in 0:(nrow(mesh$nodes)-1)){
#'     func_evaluation[i+1] = a1* sin(2*pi*mesh$nodes[i+1,1]) +
#'                            a2* sin(2*pi*mesh$nodes[i+1,2]) +
#'                            a3*sin(2*pi*mesh$nodes[i+1,3]) + 1
#'   }
#'   data = func_evaluation + rnorm(nrow(mesh$nodes), mean = 0, sd = 0.5)
#'   datamatrix = rbind(datamatrix, data)
#' }
#' ## Compute the mean of the datamatrix and subtract it to the data
#' data_bar = colMeans(datamatrix)
#' data_demean = matrix(rep(data_bar,50), nrow=50, byrow=TRUE)
#'
#' datamatrix_demeaned = datamatrix - data_demean
#' ## Set the smoothing parameter lambda
#' lambda = 0.00375
#' ## Estimate the first 2 Principal Components
#' FPCA_solution = FPCA.FEM(datamatrix = datamatrix_demeaned,
#'                       FEMbasis = FEMbasis, lambda = lambda, nPC = 2)
#'
#' ## Plot the functional loadings of the estimated Principal Components
#' plot(FPCA_solution$loadings.FEM)

FPCA.FEM<-function(locations = NULL, datamatrix, FEMbasis, lambda, nPC = 1, validation = NULL, NFolds = 5,
                   GCVmethod = "Stochastic", nrealizations = 100, search = "tree", bary.locations = NULL)
{
  incidence_matrix=NULL # if areal fpca will be included in release, this should be put in the input

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

  if(GCVmethod=="Stochastic")
    GCVmethod=2
  else if(GCVmethod=="Exact")
    GCVmethod=1
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
##################### Checking parameters, sizes and conversion ################################
  #if locations is null but bary.locations is not null, use the locations in bary.locations
  if(is.null(locations) & !is.null(bary.locations)) {
    locations = bary.locations$locations
    locations = as.matrix(locations)
  }

  checkSmoothingParametersFPCA(locations=locations, datamatrix=datamatrix, FEMbasis=FEMbasis, incidence_matrix=incidence_matrix, lambda=lambda, nPC=nPC, validation=validation, NFolds=NFolds, GCVmethod=GCVmethod ,nrealizations=nrealizations,search=search, bary.locations=bary.locations)

  ## Coverting to format for internal usage
  if(!is.null(locations))
    locations = as.matrix(locations)
  datamatrix = as.matrix(datamatrix)
  if(!is.null(incidence_matrix))
	incidence_matrix = as.matrix(incidence_matrix)
  lambda = as.matrix(lambda)

  checkSmoothingParametersSizeFPCA(locations=locations, datamatrix=datamatrix, FEMbasis=FEMbasis, incidence_matrix=incidence_matrix, lambda=lambda, ndim=ndim, mydim=mydim, validation=validation, NFolds=NFolds)

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

  bigsol = NULL
  if(class(FEMbasis$mesh) == 'mesh.2D'){
	  
	  bigsol = CPP_smooth.FEM.FPCA(locations=locations, bary.locations=bary.locations, datamatrix=datamatrix, FEMbasis=FEMbasis, incidence_matrix=incidence_matrix,
                                 lambda=lambda, ndim=ndim, mydim=mydim, nPC=nPC, validation=validation, NFolds=NFolds, GCVmethod=GCVmethod, nrealizations=nrealizations, search=search)
	  numnodes = nrow(FEMbasis$mesh$nodes)
  } else if(class(FEMbasis$mesh) == 'mesh.2.5D'){
	  
	  bigsol = CPP_smooth.manifold.FEM.FPCA(locations=locations, bary.locations=bary.locations, datamatrix=datamatrix, FEMbasis=FEMbasis, incidence_matrix=incidence_matrix,
                                          lambda=lambda, ndim=ndim, mydim=mydim, nPC=nPC, validation=validation, NFolds=NFolds, GCVmethod=GCVmethod, nrealizations=nrealizations, search=search)
	  numnodes = nrow(FEMbasis$mesh$nodes)
  } else if(class(FEMbasis$mesh) == 'mesh.3D'){
	  
	  bigsol = CPP_smooth.volume.FEM.FPCA(locations=locations, bary.locations=bary.locations, datamatrix=datamatrix, FEMbasis=FEMbasis, incidence_matrix=incidence_matrix,
                                      lambda=lambda, ndim=ndim, mydim=mydim, nPC=nPC, validation=validation, NFolds=NFolds, GCVmethod=GCVmethod, nrealizations=nrealizations, search=search)
	  numnodes = nrow(FEMbasis$mesh$nodes)
  }

  loadings=bigsol[[1]]
   # Save information of Tree Mesh
    tree_mesh = list(
    treelev = bigsol[[7]][1],
    header_orig= bigsol[[8]],
    header_scale = bigsol[[9]],
    node_id = bigsol[[10]][,1],
    node_left_child = bigsol[[10]][,2],
    node_right_child = bigsol[[10]][,3],
    node_box= bigsol[[11]])


 # Reconstruct FEMbasis with tree mesh
  mesh.class= class(FEMbasis$mesh)
  if (is.null(FEMbasis$mesh$treelev)) { #if doesn't exist the tree information
    FEMbasis$mesh = append(FEMbasis$mesh, tree_mesh)
  } #if already exist the tree information, don't append
  class(FEMbasis$mesh) = mesh.class

  loadings.FEM=FEM(loadings,FEMbasis)

  scores=bigsol[[2]]

  lambda=bigsol[[3]]

  variance_explained=bigsol[[4]]

  cumsum_percentage=bigsol[[5]]

  var=bigsol[[6]]



  # Save information of Barycenter
  bary.locations = list(barycenters = bigsol[[12]], element_ids = bigsol[[13]])
  class(bary.locations) = "bary.locations"

    reslist=list(loadings.FEM=loadings.FEM, scores=scores, lambda=lambda, variance_explained=variance_explained, cumsum_percentage=cumsum_percentage, bary.locations = bary.locations)

  return(reslist)
}
