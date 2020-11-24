#' Evaluate a FEM object at a set of point locations
#'
#' @param FEM A \code{FEM} object to be evaluated.
#' @param locations A 2-columns (in 2D) or 3-columns (in 2.5D and 3D) matrix with the spatial locations where the
#' FEM object should be evaluated.
#' @param incidence_matrix In case of areal evaluations, the #regions-by-#elements incidence matrix defining the regions
#' where the FEM object should be evaluated.
#' @param search a flag to decide the search algorithm type (tree or naive or walking search algorithm).
#' @param bary.locations A list with three vectors:
#'  \code{locations}, location points which are same as the given locations options. (checks whether both locations are the same);
#'  \code{element ids}, a vector of element id of the points from the mesh where they are located;
#'  \code{barycenters}, a vector of barycenter of points from the located element.
#' @return
#' A vector or a matrix of numeric evaluations of the \code{FEM} object.
#' If the \code{FEM} object contains multiple finite element functions the output is a matrix, and
#' each row corresponds to the location (or areal region) where the evaluation has been taken, while each column
#' corresponds to the function evaluated.
#' @description It evaluates a FEM object at the specified set of locations or areal regions. The locations are used for
#' pointwise evaluations and incidence matrix for areal evaluations.
#' The locations and the incidence matrix cannot be both NULL or both provided.
#' @usage eval.FEM(FEM, locations = NULL, incidence_matrix = NULL, search = "tree", 
#'                 bary.locations = NULL)
#' @references
#' \itemize{
#'    \item{Sangalli, L. M., Ramsay, J. O., & Ramsay, T. O. (2013). Spatial spline regression models.
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 75(4), 681-703.}
#'    \item{Azzimonti, L., Sangalli, L. M., Secchi, P., Domanin, M., & Nobile, F. (2015). Blood flow velocity field estimation
#' via spatial regression with PDE penalization. Journal of the American Statistical Association, 110(511), 1057-1071.}
#' }
#' @examples
#' library(fdaPDE)
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
#'
#' ## Evaluate the finite element function in the location (1,0.5)
#' eval.FEM(FEMfunction, locations = matrix(c(1, 0.5), ncol = 2))
#'
#' ## Evaluate the mean of the finite element function over the fifth triangle of the mesh
#' incidence_matrix = matrix(0, ncol = nrow(mesh$triangles))
#' incidence_matrix[1,5] = 1
#' eval.FEM(FEMfunction, incidence_matrix = incidence_matrix)
#' @export

eval.FEM <- function(FEM, locations = NULL, incidence_matrix = NULL, search = "tree", bary.locations = NULL)
{
  ##################### Checking parameters, sizes and conversion ################################
  if (is.null(FEM))
    stop("FEM required;  is NULL.")
  if(class(FEM) != "FEM")
    stop("'FEM' is not of class 'FEM'")
  #if locations is null but bary.locations is not null, use the locations in bary.locations
  if(is.null(locations) && !is.null(bary.locations) && is.null(incidence_matrix)) {
    # print("'locations' and 'incidence_matrix' are NULL, evalutation performed in bary.locations")
    locations = bary.locations$locations
    locations = as.matrix(locations)
  }
  if (is.null(locations) && is.null(incidence_matrix))
    stop("'locations' NOR 'incidence_matrix' required;  both are NULL")
  if (!is.null(locations) && !is.null(incidence_matrix))
    stop("'locations' NOR 'incidence_matrix' required; both are given")

  # if(!is.null(locations))
  #  if(dim(locations)[1]==dim(FEM$FEMbasis$mesh$nodes)[1] & dim(locations)[2]==dim(FEM$FEMbasis$mesh$nodes)[2])
  #   warning("The locations matrix has the same dimensions as the mesh nodes. If you want to get the FEM object evaluation
  #           at the mesh nodes, use FEM$coeff instead")

  if(search == "naive" || search == 1)
    search=1
  else if(search == "tree" || search == 2)
    search=2
  else if(search == "walking" || search == 3)
    search=3

  if(class(FEM$FEMbasis$mesh)=='mesh.2.5D' && search ==3)
    stop("2.5D search must be either 'tree' or 'naive'")
  
  if (search != 1 && search != 2 && search != 3)
    stop("search must be either 'tree' or 'naive' or 'walking'")

  #Check the locations in 'bary.locations' and 'locations' are the same
  if(!is.null(bary.locations) && !is.null(locations))
  {
    flag=TRUE
    for (i in 1:nrow(locations)) {
      if (!(locations[i,1]==bary.locations$locations[i,1] & locations[i,2] == bary.locations$locations[i,2])) {
        flag = FALSE
        break
      }
    }
    if (flag == FALSE) {
      stop("Locations are not same as the one in barycenter information.")
    }
  }  # end of bary.locations

  if (is.null(locations))
    locations <- matrix(nrow = 0, ncol = 2)
  else
    incidence_matrix <- matrix(nrow = 0, ncol = 1)


  ################## End checking parameters, sizes and conversion #############################
  res <- NULL

  if(class(FEM$FEMbasis$mesh)=='mesh.2D'){
    ndim = 2
    mydim = 2
    res = CPP_eval.FEM(FEM, locations, incidence_matrix, TRUE, ndim, mydim, search, bary.locations)

  }else if(class(FEM$FEMbasis$mesh)=='mesh.2.5D'){
    ndim = 3
    mydim = 2
    res = CPP_eval.manifold.FEM(FEM, locations, incidence_matrix, TRUE, ndim, mydim, search, bary.locations)
  }else if(class(FEM$FEMbasis$mesh)=='mesh.3D'){
    ndim = 3
    mydim = 3
    res = CPP_eval.volume.FEM(FEM, locations, incidence_matrix, TRUE, ndim, mydim, search, bary.locations)
  }

  return(as.matrix(res))
}

#' Evaluate a FEM.time object at a set of point locations
#'
#' @param FEM.time A \code{FEM.time} object to be evaluated.
#' @param space.time.locations A 3-columns (in case of planar mesh) or 4-columns(in case of 2D manifold in a 3D space or a 3D volume) 
#' matrix with the time instants and spatial locations where the FEM.time object should be evaluated. 
#' The first column is for the time instants. If given, \code{locations}, \code{incidence_matrix} and \code{time.instants} must be NULL.
#' @param locations A 2-columns (in case of planar mesh) or 3-columns(in case of 2D manifold in a 3D space or a 3D volume) matrix with the spatial locations where the FEM.time object should be evaluated.
#' @param incidence_matrix In case of areal data, the #regions x #elements incidence matrix defining the regions.
#' @param time.instants A vector with the time instants where the FEM.time object should be evaluated.
#' @param lambdaS The index of the lambdaS choosen for the evaluation. 
#' @param lambdaT The index of the lambdaT choosen for the evaluation.
#' @param search a flag to decide the search algorithm type (tree or naive or walking search algorithm).
#' @param bary.locations A list with three vectors:
#'  \code{locations}, location points which are same as the given locations options. (checks whether both locations are the same);
#'  \code{element ids}, a vector of element id of the points from the mesh where they are located;
#'  \code{barycenters}, a vector of barycenter of points from the located element.
#' @return
#' A matrix of numeric evaluations of the \code{FEM.time} object. Each row indicates the location where 
#' the evaluation has been taken, the column indicates the function evaluated.
#' @description It evaluates a \code{FEM.time} object at the specified set of locations or regions. 
#' If \code{space.time.locations} is provided \code{locations}, \code{incidence_matrix} and 
#' \code{time.instants} must be NULL. Otherwise \code{time.instants} and one of \code{locations} and 
#' \code{incidence_matrix} must be given. In this case the evaluation is perform on the tensor grid
#' \code{time.instants}-by-\code{locations} (or \code{time.instants}-by-areal domains).
#' @usage eval.FEM.time(FEM.time, locations = NULL, time.instants = NULL, 
#'                      space.time.locations = NULL, incidence_matrix = NULL, lambdaS = 1, 
#'                      lambdaT = 1, search = "tree", bary.locations = NULL)
#' @references
#'  Devillers, O. et al. 2001. Walking in a Triangulation, Proceedings of the Seventeenth Annual Symposium on Computational Geometry
#' @export
#' @examples
#' library(fdaPDE)
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
#' time = 1:5
#' coeff = rep(fs.test(mesh$nodes[,1], mesh$nodes[,2]),5)*time
#' ## Create the FEM.time object
#' FEM_time_function = FEM.time(coeff=coeff, time_mesh=1:5, FEMbasis=FEMbasis, FLAG_PARABOLIC=TRUE)
#'
#' evaluations = eval.FEM.time(FEM_time_function, locations = matrix(c(-0.92,0), ncol=2), 
#'                             time.instants = time)

eval.FEM.time <- function(FEM.time, locations = NULL, time.instants = NULL, space.time.locations = NULL, 
                          incidence_matrix = NULL, lambdaS = 1, lambdaT = 1, search = "tree", bary.locations = NULL)
{
  if (is.null(FEM.time))
    stop("FEM.time required;  is NULL.")
  if(class(FEM.time) != "FEM.time")
    stop("'FEM.time' is not of class 'FEM.time'")
  
  # save flag for areal evaluation
  FLAG_AREAL_EVALUATION = FALSE
  if (!is.null(incidence_matrix)) 
    FLAG_AREAL_EVALUATION = TRUE
  
  # save flag for tensorization: if TRUE tensorization has to be performed, 
  # if FALSE the 'space.time.locations' has already have been constructed
  FLAG_TENSORIZE = TRUE
  if (!is.null(space.time.locations)) 
    FLAG_TENSORIZE = FALSE
  
  # check that when space.time.locations is provided the other fields are NULL
  if (!FLAG_TENSORIZE && FLAG_AREAL_EVALUATION) 
    stop("'incidence_matrix' and 'space.time.locations' both provided")
  if (!FLAG_TENSORIZE && !is.null(locations)) 
    stop("'locations' and 'space.time.locations' both provided")
  if (!FLAG_TENSORIZE && !is.null(time.instants)) 
    stop("'time.instants' and 'space.time.locations' both provided")
  
  if (FLAG_TENSORIZE) {
    if(is.null(time.instants))
      stop("'time.instants' must be given if space.time.locations is NULL")
    if (!is.vector(time.instants) && nrow(time.instants)!=1 && ncol(time.instants)!=1) 
      stop("'time.instants' must be a vector")
    time.instants = as.vector(time.instants)
    
    if (!FLAG_AREAL_EVALUATION) {
      incidence_matrix <- matrix(nrow = 0, ncol = 1)
      
      if(is.null(locations) && !is.null(bary.locations)) {
        # print("space.time.locations', 'locations' and 'incidence_matrix' are NULL, evalutation performed in bary.locations")
        locations = bary.locations$locations
        locations = as.matrix(locations)
      }
      if (is.null(locations))
        stop("'space.time.locations', 'locations' and 'incidence_matrix' are NULL")
      if (ncol(locations)!=2 && ncol(locations)!=3)
        stop("'locations' and 'bary.locations$locations' when provided must have 2 or 3 columns")
      
      #Check the locations in 'bary.locations' and 'locations' are the same
      if(!is.null(bary.locations) && !is.null(locations))
      {
        if (dim(locations)!=dim(bary.locations$locations))
          stop("'locations' and 'bary.locations$locations' are different")
        if (sum(abs(locations - bary.locations$locations))!=0)
          stop("'locations' and 'bary.locations$locations' are different")
        
        #Repeat 'bary.locations' to have the same size as 'space.time.locations' (needed to convert size for C++)
        time_len = length(FEM.time$mesh_time)
        bary.locations$locations = matrix( rep( t(bary.locations$locations) , time_len ) , ncol =  ncol(bary.locations$locations) , byrow = TRUE )
        bary.locations$element_ids = rep(bary.locations$element_ids, time_len)
        bary.locations$barycenters = matrix( rep( t(bary.locations$barycenters) , time_len ) , ncol =  ncol(bary.locations$barycenters) , byrow = TRUE ) 
      }  # end of bary.locations
        
      time_locations = rep(time.instants,each=nrow(locations))
      space.time.locations = cbind(rep(locations[,1],length(time.instants)),rep(locations[,2],length(time.instants)))
      if(ncol(locations)==3)
        space.time.locations=cbind(space.time.locations,rep(locations[,3],length(time.instants)))
    }else{
      space.time.locations <- matrix(nrow=0, ncol=1)
    }
  }else{
    if(dim(space.time.locations)[2]<3)
      stop("'space.time.locations' requires at least t,X,Y")
    # if(dim(space.time.locations)[1]==dim(FEM.time$FEMbasis$mesh$nodes)[1] & (dim(space.time.locations)[2]-1)==dim(FEM.time$FEMbasis$mesh$nodes)[2])
    #   warning("The space.time.locations matrix has the same dimensions as the mesh nodes. If you want to get the FEM.time object evaluation
    #       at the mesh nodes, use FEM.time$coeff instead")
    time_locations <- space.time.locations[,1]
    space.time.locations <- space.time.locations[,2:dim(space.time.locations)[2]]
    incidence_matrix <- matrix(nrow = 0, ncol = 1)
  }

  #conversion of search type
  if(search == "naive" || search == 1)
    search=1
  else if(search == "tree" || search == 2)
    search=2
  else if(search == "walking" || search == 3)
    search=3
  if (search != 1 & search != 2 & search != 3)
    stop("search must be either tree or naive or walking.")
  
  if(class(FEM.time$FEMbasis$mesh)=='mesh.2.5D' & search ==3)
  	stop("2.5D search must be either tree or naive.")


  if(dim(FEM.time$coeff)[2]>1||dim(FEM.time$coeff)[3]>1)
  {
    if(dim(FEM.time$coeff)[2]>1 && lambdaS==1)
      warning("the first value of lambdaS is being used")
    if(dim(FEM.time$coeff)[3]>1 && lambdaT==1)
      warning("the first value of lambdaT is being used")
    f = FEM.time(coeff=array(FEM.time$coeff[,lambdaS,lambdaT]),time_mesh=FEM.time$mesh_time,FEMbasis=FEM.time$FEMbasis,FLAG_PARABOLIC=FEM.time$FLAG_PARABOLIC)
  }
  else
    f = FEM.time

  res <- NULL
  if(class(FEM.time$FEMbasis$mesh)=='mesh.2D'){
    ndim = 2
    mydim = 2
    res = CPP_eval.FEM.time(f, space.time.locations, time_locations, incidence_matrix, FEM.time$FLAG_PARABOLIC, TRUE, ndim, mydim, search, bary.locations)
  }else if(class(FEM.time$FEMbasis$mesh)=='mesh.2.5D'){
    ndim = 3
    mydim = 2
    res = CPP_eval.manifold.FEM.time(f, space.time.locations, time_locations, incidence_matrix, FEM.time$FLAG_PARABOLIC, TRUE, ndim, mydim, search, bary.locations)
  }else if(class(FEM.time$FEMbasis$mesh)=='mesh.3D'){
    ndim = 3
    mydim = 3
    res = CPP_eval.volume.FEM.time(f, space.time.locations, time_locations, incidence_matrix, FEM.time$FLAG_PARABOLIC, TRUE, ndim, mydim, search, bary.locations)
  }

  return(as.matrix(res))
}


# R_eval.FEM <- function(FEM, locations)
# {
#   if (is.vector(locations))
#   {
#     locations = t(as.matrix(locations))
#   }
#
#   N = nrow(locations)
#
#   # Augment Xvec and Yvec by ones for computing barycentric coordinates
#   Pgpts = cbind(matrix(1,N,1),locations[,1],locations[,2])
#
#   # Get nodes and index
#
#   FEMbasis = FEM$FEMbasis
#   mesh = FEMbasis$mesh
#
#   nodes = mesh$nodes
#   triangles = mesh$triangles
#   coeff = FEM$coeff
#   nsurf = dim(coeff)[2]
#
#   FEMbasis = FEM$FEMbasis
#   order = FEMbasis$order
#   #nodeindex = params$nodeindex
#   detJ = FEMbasis$detJ
#
#   # 1st, 2nd, 3rd vertices of triangles
#
#   v1 = nodes[triangles[,1],]
#   v2 = nodes[triangles[,2],]
#   v3 = nodes[triangles[,3],]
#
#   if(order !=2 && order != 1)
#   {
#     stop('ORDER is neither 1 or 2.')
#   }
#
#   # Denominator of change of coordinates chsange matrix
#
#   modJ = FEMbasis$detJ
#   ones3 = matrix(1,3,1)
#   modJMat = modJ %*% t(ones3)
#
#   M1 = cbind(v2[,1]*v3[,2] - v3[,1]*v2[,2], v2[,2] - v3[,2], v3[,1] - v2[,1])/(modJMat)
#   M2 = cbind(v3[,1]*v1[,2] - v1[,1]*v3[,2], v3[,2] - v1[,2], v1[,1] - v3[,1])/(modJMat)
#   M3 = cbind(v1[,1]*v2[,2] - v2[,1]*v1[,2], v1[,2] - v2[,2], v2[,1] - v1[,1])/(modJMat)
#
#   ind = matrix(0,N,1)
#   for(i in 1:N)
#   {
#     ind[i] = R_insideIndex(mesh, as.numeric(locations[i,]))
#   }
#
#   evalmat = matrix(NA, nrow=N, ncol=nsurf)
#
#   for (isurf in 1:nsurf)
#   {
#     for(i in 1:N)
#     {
#       indi = ind[i]
#
#       if(!is.nan(indi))
#       {
#         baryc1 = (M1[indi,]*Pgpts[i,]) %*% ones3
#         baryc2 = (M2[indi,]*Pgpts[i,]) %*% ones3
#         baryc3 = (M3[indi,]*Pgpts[i,]) %*% ones3
#
#         if(order == 2)
#         {
#           c1 = coeff[triangles[indi,1],isurf]
#           c2 = coeff[triangles[indi,2],isurf]
#           c3 = coeff[triangles[indi,3],isurf]
#           c4 = coeff[triangles[indi,6],isurf]
#           c5 = coeff[triangles[indi,4],isurf]
#           c6 = coeff[triangles[indi,5],isurf]
#
#           fval =  c1*(2* baryc1^2 - baryc1) +
#             c2*(2* baryc2^2 - baryc2) +
#             c3*(2* baryc3^2 - baryc3) +
#             c4*(4* baryc1 * baryc2) +
#             c5*(4* baryc2 * baryc3) +
#             c6*(4* baryc3 * baryc1)
#           evalmat[i,isurf] = fval
#         }
#         else
#         {
#           c1 = coeff[triangles[indi,1],isurf]
#           c2 = coeff[triangles[indi,2],isurf]
#           c3 = coeff[triangles[indi,3],isurf]
#           fval = c1*baryc1 + c2*baryc2 + c3*baryc3
#           evalmat[i,isurf] = fval
#         }
#       }
#     }
#   }
#   return(evalmat)
# }
#
R_insideIndex = function (mesh, location)
{
  #  insideIndex returns the index of the triangle containing the point
  # (X,Y) if such a triangle exists, and NaN otherwise.
  #  TRICOEF may have already been calculated for efficiency,
  #  but if the function is called with four arguments, it is calculated.


  eps=2.2204e-016
  small = 10000*eps

  nodes = mesh$nodes
  triangles = mesh$triangles
  X = location[1]
  Y = location[2]

  ntri   = dim(triangles)[[1]]
  indtri   = matrix(1:ntri,ncol=1)

  #  compute coefficients for computing barycentric coordinates if needed

  tricoef = R_tricoefCal(mesh)

  #  compute barycentric coordinates
  r3 = X - nodes[triangles[,3],1]
  s3 = Y - nodes[triangles[,3],2]
  lam1 = ( tricoef[,4]*r3 - tricoef[,2]*s3)
  lam2 = (-tricoef[,3]*r3 + tricoef[,1]*s3)
  lam3 = 1 - lam1 - lam2

  #  test these coordinates for a triple that are all between 0 and 1
  int  = (-small <= lam1 & lam1 <= 1+small) &
    (-small <= lam2 & lam2 <= 1+small) &
    (-small <= lam3 & lam3 <= 1+small)

  #  return the index of this triple, or NaN if it doesn't exist
  indi = indtri[int]
  if (length(indi)<1)
  {
    ind = NA
  }else{
    ind = min(indi)
  }

  ind
}
