#' Create a FEM basis
#'
#' @param mesh A \code{mesh.2D}, \code{mesh.2.5D} or \code{mesh.3D} object representing the domain triangulation. See \link{create.mesh.2D}, \link{create.mesh.2.5D}, \link{create.mesh.3D}.
#' @param saveTree a flag to decide to save the tree mesh information in advance (default is FALSE)
#' @return A \code{FEMbasis} object. This contains the \code{mesh}, along with some additional quantities:
#' \itemize{
#'  \item{\code{order}}{Either "1" or "2" for the 2D and 2.5D case, and "1" for the 3D case.
#'  Order of the Finite Element basis.}
#'  \item{\code{nbasis}}{Scalar. The number of basis.}
#'  \item{\code{transf_coord}}{It takes value only in the 2D case. It is a list of 4 vectors: diff1x, diff1y, diff2x and diff2y.
#'  Each vector has length #triangles and encodes the information for the tranformation matrix that transforms the
#'  nodes of the reference triangle to the nodes of the i-th triangle.
#'  The tranformation matrix for the i-th triangle has the form [diff1x[i] diff2x[i]; diff1y[i] diff2y[i]].}
#'  \item{\code{detJ}}{It takes value only in the 2D case. A vector of length #triangles. The ith element contains
#'  the determinant of the transformation from the reference triangle to the nodes of the i-th triangle.
#'  Its value is also the double of the area of each triangle of the basis.}
#' }
#' @description Sets up a Finite Element basis. It requires a \code{mesh.2D}, \code{mesh.2.5D} or \code{mesh.3D} object,
#' as input.
#' The basis' functions are globally continuos functions, that are polynomials once restricted to a triangle in the mesh.
#' The current implementation includes linear finite elements (when \code{order = 1} in the input \code{mesh}) and
#' quadratic finite elements (when \code{order = 2} in the input \code{mesh}).
#' If saveTree flag is TRUE, it saves the tree mesh information in advance inside mesh object and can be used later on to save mesh construction time.
#' @usage create.FEM.basis(mesh, saveTree = FALSE)
#' @seealso \code{\link{create.mesh.2D}}, \code{\link{create.mesh.2.5D}},\code{\link{create.mesh.3D}}
#' @examples
#' ## Upload the quasicircle2D data
#' data(quasicircle2D)
#' boundary_nodes = quasicircle2D$boundary_nodes
#' boundary_segments = quasicircle2D$boundary_segments
#' locations = quasicircle2D$locations
#' data = quasicircle2D$data												
#'
#' ## Create the 2D mesh
#' mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)
#' ## Plot it
#' plot(mesh)
#' ## Create the basis
#' FEMbasis = create.FEM.basis(mesh)
#' ## Upload the hub2.5D data
#' data(hub2.5D)
#' hub2.5D.nodes = hub2.5D$hub2.5D.nodes
#' hub2.5D.triangles = hub2.5D$hub2.5D.triangles										
#'
#' ## Create the 2.5D mesh
#' mesh = create.mesh.2.5D(nodes = hub2.5D.nodes, triangles = hub2.5D.triangles)
#' ## Plot it
#' plot(mesh)
#' ## Create the basis
#' FEMbasis = create.FEM.basis(mesh)
#' @export

create.FEM.basis = function(mesh, saveTree = FALSE)
{
  if(class(mesh)!='mesh.2D' & class(mesh)!='mesh.2.5D' & class(mesh)!='mesh.3D')
    stop("Unknown mesh class")

  if (class(mesh)=="mesh.2D"){

    #  The number of basis functions corresponds to the number of vertices
    #  for order = 1, and to vertices plus edge midpoints for order = 2

    nbasis = dim(mesh$nodes)[[1]]
    eleProp = R_elementProperties(mesh)

    #eleProp = NULL
    #if(CPP_CODE == FALSE)
    #{
    #  eleProp = R_elementProperties(mesh)
    #}

  if (saveTree == TRUE) {
      ## Call C++ function
      orig_mesh = mesh
      mesh$triangles = mesh$triangles - 1
      mesh$edges = mesh$edges - 1
      mesh$neighbors[mesh$neighbors != -1] = mesh$neighbors[mesh$neighbors != -1] - 1
      myDim = 2
      nDim = 2

      storage.mode(mesh$nodes) <- "double"
      storage.mode(mesh$triangles) <- "integer"
      storage.mode(mesh$edges) <- "integer"
      storage.mode(mesh$neighbors) <- "integer"
      storage.mode(mesh$order) <- "integer"
      storage.mode(myDim) <- "integer"
      storage.mode(nDim) <- "integer"

      bigsol <- .Call("tree_mesh_construction", mesh, mesh$order, myDim, nDim, PACKAGE = "fdaPDE")
      tree_mesh = list(
      treelev = bigsol[[1]][1],
      header_orig= bigsol[[2]],
      header_scale = bigsol[[3]],
      node_id = bigsol[[4]][,1],
      node_left_child = bigsol[[4]][,2],
      node_right_child = bigsol[[4]][,3],
      node_box= bigsol[[5]])

      # Reconstruct FEMbasis with tree mesh
      mesh.class= class(orig_mesh)
      mesh = append(orig_mesh, tree_mesh)
      class(mesh) = mesh.class
  }

  FEMbasis = list(mesh = mesh, order = as.integer(mesh$order), nbasis = nbasis, detJ=eleProp$detJ, transf_coord = eleProp$transf_coord)
  class(FEMbasis) = "FEMbasis"

  FEMbasis
  } else if (class(mesh) == "mesh.2.5D" || class(mesh) == "mesh.3D"){
      if (saveTree == TRUE) {

        if (class(mesh) == "mesh.2.5D") {
          orig_mesh = mesh

          myDim = 2
          nDim = 3
          # C++ function for manifold works with vectors not with matrices
          mesh$triangles=c(t(mesh$triangles))
          mesh$nodes=c(t(mesh$nodes))
          # Indexes in C++ starts from 0, in R from 1, opportune transformation
          mesh$triangles=mesh$triangles-1

          storage.mode(mesh$order) <- "integer"
          storage.mode(mesh$nnodes) <- "integer"
          storage.mode(mesh$ntriangles) <- "integer"
          storage.mode(mesh$nodes) <- "double"
          storage.mode(mesh$triangles) <- "integer"
          storage.mode(myDim) <- "integer"
          storage.mode(nDim) <- "integer"

        } else if (class(mesh) == "mesh.3D") {
          orig_mesh = mesh

          myDim = 3
          nDim = 3
          # C++ function for volumetric works with vectors not with matrices
          mesh$tetrahedrons=c(t(mesh$tetrahedrons))
          mesh$nodes=c(t(mesh$nodes))
          # Indexes in C++ starts from 0, in R from 1, opportune transformation
          mesh$tetrahedrons=mesh$tetrahedrons-1

          storage.mode(mesh$order) <- "integer"
          storage.mode(mesh$nnodes) <- "integer"
          storage.mode(mesh$ntetrahedrons) <- "integer"
          storage.mode(mesh$nodes) <- "double"
          storage.mode(mesh$tetrahedrons) <- "integer"
          storage.mode(myDim) <- "integer"
          storage.mode(nDim) <- "integer"
        }


        bigsol <- .Call("tree_mesh_construction", mesh, mesh$order, myDim, nDim, PACKAGE = "fdaPDE")
        tree_mesh = list(
        treelev = bigsol[[1]][1],
        header_orig= bigsol[[2]],
        header_scale = bigsol[[3]],
        node_id = bigsol[[4]][,1],
        node_left_child = bigsol[[4]][,2],
        node_right_child = bigsol[[4]][,3],
        node_box= bigsol[[5]])

        # Reconstruct FEMbasis with tree mesh
        mesh.class= class(orig_mesh)
        mesh = append(orig_mesh, tree_mesh)
        class(mesh) = mesh.class
      }

      FEMbasis = list(mesh = mesh, order = as.integer(mesh$order), nbasis = mesh$nnodes)
      class(FEMbasis) = "FEMbasis"
      FEMbasis
  }
 }
#' Define a surface or spatial field by a Finite Element basis expansion
#'
#' @param coeff A vector or a matrix containing the coefficients for the Finite Element basis expansion. The number of rows
#' (or the vector's length) corresponds to the number of basis in \code{FEMbasis}.
#' The number of columns corresponds to the number of functions.
#' @param FEMbasis A \code{FEMbasis} object defining the Finite Element basis, created by \link{create.FEM.basis}.
#' @description This function defines a FEM object.
#' @usage FEM(coeff,FEMbasis)
#' @return An \code{FEM} object. This contains a list with components \code{coeff} and \code{FEMbasis}.
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
#' ## Plot it
#' plot(FEMfunction)
#' @export

FEM<-function(coeff,FEMbasis)
{
  if (is.null(coeff))
    stop("coeff required;  is NULL.")
  if (is.null(FEMbasis))
    stop("FEMbasis required;  is NULL.")
  if(class(FEMbasis) != "FEMbasis")
    stop("FEMbasis not of class 'FEMbasis'")
  coeff = as.matrix(coeff)
  if(nrow(coeff) != FEMbasis$nbasis)
    stop("Number of row of 'coeff' different from number of basis")

  fclass = NULL
  fclass = list(coeff=coeff, FEMbasis=FEMbasis)
  class(fclass)<-"FEM"
  return(fclass)
}

#' Define a spatio-temporal field by a Finite Element basis expansion
#'
#' @param coeff A vector or a matrix containing the coefficients for the spatio-temporal basis expansion. The number of rows
#' (or the vector's length) corresponds to the number of basis in \code{FEMbasis} times the number of knots in \code{time_mesh}.
#' @param FEMbasis A \code{FEMbasis} object defining the Finite Element basis, created by \link{create.FEM.basis}.
#' @param time_mesh A vector containing the b-splines knots for separable smoothing and the nodes for finite differences for parabolic smoothing
#' @param FLAG_PARABOLIC Boolean. If \code{TRUE} the coefficients are from parabolic smoothing, if \code{FALSE} the separable one.
#' @description This function defines a FEM.time object.
#' @usage FEM.time(coeff,time_mesh,FEMbasis,FLAG_PARABOLIC=FALSE)
#' @return A \code{FEM.time} object. This contains a list with components \code{coeff}, \code{mesh_time}, \code{FEMbasis} and \code{FLAG_PARABOLIC}.
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
#' coeff = rep(fs.test(mesh$nodes[,1], mesh$nodes[,2]),5)
#' time_mesh = seq(0,1,length.out = 5)
#' ## Create the FEM object
#' FEMfunction = FEM.time(coeff, time_mesh, FEMbasis, FLAG_PARABOLIC = TRUE)
#' ## Plot it at desired time
#' plot(FEMfunction,0.7)
#' @export

FEM.time<-function(coeff,time_mesh,FEMbasis,FLAG_PARABOLIC=FALSE)
{
  if(is.vector(coeff)){
    coeff = array(coeff, dim = c(length(coeff),1,1))
  }
  M = ifelse(FLAG_PARABOLIC,length(time_mesh),length(time_mesh)+2)
  if (is.null(coeff))
    stop("coeff required;  is NULL.")
  if (is.null(time_mesh))
    stop("time_mesh required;  is NULL.")
  if (is.null(FLAG_PARABOLIC))
    stop("FLAG_PARABOLIC required;  is NULL.")
  if (is.null(FEMbasis))
    stop("FEMbasis required;  is NULL.")
  if(class(FEMbasis) != "FEMbasis")
    stop("FEMbasis not of class 'FEMbasis'")
  if(dim(coeff)[1] != (FEMbasis$nbasis*M))
    stop("Number of row of 'coeff' different from number of basis")

  fclass = NULL
  fclass = list(coeff=coeff, mesh_time=time_mesh, FLAG_PARABOLIC=FLAG_PARABOLIC, FEMbasis=FEMbasis)
  class(fclass)<-"FEM.time"
  return(fclass)
}


R_elementProperties=function(mesh)
{
  nele = dim(mesh$triangles)[[1]]
  nodes = mesh$nodes
  triangles = mesh$triangles

  #detJ   = matrix(0,nele,1)      #  vector of determinant of transformations
  #metric = array(0,c(nele,2,2))  #  3-d array of metric matrices
  #transf = array(0,c(nele,2,2))

  transf_coord = NULL
  transf_coord$diff1x = nodes[triangles[,2],1] - nodes[triangles[,1],1]
  transf_coord$diff1y = nodes[triangles[,2],2] - nodes[triangles[,1],2]
  transf_coord$diff2x = nodes[triangles[,3],1] - nodes[triangles[,1],1]
  transf_coord$diff2y = nodes[triangles[,3],2] - nodes[triangles[,1],2]

  #  Jacobian or double of the area of triangle
  detJ = transf_coord$diff1x*transf_coord$diff2y - transf_coord$diff2x*transf_coord$diff1y

  #Too slow, computed only for stiff from diff1x,diff1y,..
  # for (i in 1:nele)
  # {
  #   #transf[i,,] = rbind(cbind(diff1x,diff2x),c(diff1y,diff2y))
  #   #  Compute controvariant transformation matrix OSS: This is (tranf)^(-T)
  #   Ael = matrix(c(diff2y, -diff1y, -diff2x,  diff1x),nrow=2,ncol=2,byrow=T)/detJ[i]
  #
  #   #  Compute metric matrix
  #   metric[i,,] = t(Ael)%*%Ael
  # }

  #FEStruct <- list(detJ=detJ, metric=metric, transf=transf)
  FEStruct <- list(detJ=detJ, transf_coord=transf_coord)
  return(FEStruct)
}