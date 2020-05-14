triangulate_native <- function(P, PB, PA, S, SB,H, TR, flags) {
  ## It is necessary to check for NAs and NaNs, as the triangulate C
  ## code crashes if fed with them
  
  P  <- as.matrix(P)
  PB <- as.integer(PB)
  PA <- as.matrix(PA)
  S  <- as.matrix(S)
  SB <- as.integer(SB)
  H  <- as.matrix(H)
  TR  <- as.matrix(TR)
  
  storage.mode(P)  <- "double"
  storage.mode(PA) <- "double"
  #storage.mode(PB) <- "integer"
  storage.mode(S)  <- "integer"
  #storage.mode(SB) <- "integer"
  storage.mode(H)  <- "double"
  storage.mode(TR) <- "integer"
  storage.mode(flags) <- 'character'
  ## Call the main routine
  out <- .Call("R_triangulate_native",
               t(P),
               PB,
               PA,
               t(S),
               SB,
               H,
               t(TR),
               flags,
               PACKAGE="fdaPDE")
  names(out) <- c("P", "PB", "PA", "T", "S", "SB", "E", "EB","TN", "VP", "VE", "VN", "VA")
  class(out) <- "triangulation"
  return(out)
}

#' Create a 2D triangular mesh
#' 
#' @param nodes A #nodes-by-2 matrix containing the x and y coordinates of the mesh nodes.
#' @param nodesattributes A matrix with #nodes rows containing nodes' attributes. 
#' These are passed unchanged to the output. If a node is added during the triangulation process or mesh refinement, its attributes are computed  
#' by linear interpolation using the attributes of neighboring nodes. This functionality is for instance used to compute the value 
#' of a Dirichlet boundary condition at boundary nodes added during the triangulation process.
#' @param segments A #segments-by-2 matrix. Each row contains the row's indices in \code{nodes} of the vertices where the segment starts from and ends to.
#' Segments are edges that are not splitted during the triangulation process. These are for instance used to define the boundaries
#' of the domain. If this is input is NULL, it generates a triangulation over the
#' convex hull of the points specified in \code{nodes}.
#' @param holes A #holes-by-2 matrix containing the x and y coordinates of a point internal to each hole of the mesh. These points are used to carve holes
#' in the triangulation, when the domain has holes.
#' @param triangles A #triangles-by-3 (when \code{order} = 1) or #triangles-by-6 (when \code{order} = 2) matrix.
#' This option is used when a triangulation is already available. It specifies the triangles giving the row's indices in \code{nodes} of the triangles' vertices and (when \code{nodes} = 2) also if the triangles' edges midpoints. The triangles' vertices and midpoints are ordered as described 
#' at \cr https://www.cs.cmu.edu/~quake/triangle.highorder.html.
#' In this case the function \code{create.mesh.2D} is used to produce a complete mesh.2D object. 
#' @param order Either '1' or '2'. It specifies wether each mesh triangle should be represented by 3 nodes (the triangle' vertices) or by 6 nodes (the triangle's vertices and midpoints). 
#' These are
#' respectively used for linear (order = 1) and quadratic (order = 2) Finite Elements. Default is \code{order} = 1.
#' @param verbosity This can be '0', '1' or '2'. It indicates the level of verbosity in the triangulation process. When \code{verbosity} = 0 no message is returned
#' during the triangulation. When \code{verbosity} = 2 the triangulation process is described step by step by displayed messages.
#' Default is \code{verbosity} = 0.
#' @description This function is a wrapper of the Triangle library (http://www.cs.cmu.edu/~quake/triangle.html). It can be used
#' to create a triangulation of the domain of interest starting from a list of points, to be used as triangles' vertices, and a list of segments, that define the domain boundary. The resulting
#' mesh is a Constrained Delaunay triangulation. This is constructed in a way to preserve segments provided in the input \code{segments} without splitting them. This imput can be used to define the boundaries
#' of the domain. If this imput is NULL, it generates a triangulation over the
#' convex hull of the points. 
#' It is also possible to create a mesh.2D from the nodes locations and the connectivity matrix.
#' @usage create.mesh.2D(nodes, nodesattributes = NA, segments = NA, holes = NA, 
#'                      triangles = NA, order = 1, verbosity = 0)
#' @seealso \code{\link{refine.mesh.2D}}, \code{\link{create.FEM.basis}}
#' @return An object of the class mesh.2D with the following output:
#' \itemize{
#' \item{\code{nodes}}{A #nodes-by-2 matrix containing the x and y coordinates of the mesh nodes.}
#' \item{\code{nodesmarkers}}{A vector of length #nodes, with entries either '1' or '0'. An entry '1' indicates that the corresponding node is a boundary node; an entry '0' indicates that the corresponding node is not a boundary node.}
#' \item{\code{nodesattributes}}{A matrix with #nodes rows containing nodes' attributes. 
#' These are passed unchanged from the input.}
#' \item{\code{triangles}}{A #triangles-by-3 (when \code{order} = 1) or #triangles-by-6 (when \code{order} = 2) matrix.
#' This option is used when a triangulation is already available. It specifies the triangles giving the indices in \code{nodes} of the triangles' vertices and (when \code{nodes} = 2) also if the triangles' edges midpoints. The triangles' vertices and midpoints are ordered as described 
#' at  \cr https://www.cs.cmu.edu/~quake/triangle.highorder.html.}
#' \item{\code{segmentsmarker}}{A vector of length #segments with entries either '1' or '0'. An entry '1' indicates that the corresponding element in \code{segments} is a boundary segment;  
#' an entry '0' indicates that the corresponding segment is not a boundary segment.}
#' \item{\code{edges}}{A #edges-by-2 matrix containing all the edges of the triangles in the output triangulation. Each row contains the row's indices in \code{nodes}, indicating the nodes where the edge starts from and ends to.}
#' \item{\code{edgesmarkers}}{A vector of lenght #edges with entries either '1' or '0'. An entry '1' indicates that the corresponding element in \code{edge} is a boundary edge;  
#' an entry '0' indicates that the corresponding edge is not a boundary edge.}
#' \item{\code{neighbors}}{A #triangles-by-3 matrix. Each row contains the indices of the three neighbouring triangles. An entry '-1' indicates that 
#' one edge of the triangle is a boundary edge.}
#' \item{\code{holes}}{A #holes-by-2 matrix containing the x and y coordinates of a point internal to each hole of the mesh. These points are used to carve holes
#' in the triangulation, when the domain has holes.}
#' \item{\code{order}}{Either '1' or '2'. It specifies wether each mesh triangle should be represented by 3 nodes (the triangle' vertices) or by 6 nodes (the triangle's vertices and midpoints). 
#' These are respectively used for linear (order = 1) and quadratic (order = 2) Finite Elements.}
#' }
#' @export
#' @examples 
#' library(fdaPDE)
#' 
#' ## Upload the quasicirle2D data
#' data(quasicircle2D)
#' boundary_nodes = quasicircle2D$boundary_nodes
#' boundary_segments = quasicircle2D$boundary_segments
#' locations = quasicircle2D$locations
#' data = quasicircle2D$data
#' 
#' ## Create mesh from boundary
#' ## if the domain is convex it is sufficient to call:
#' mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations))
#' plot(mesh)
#' 
#' ## if the domain is not convex, pass in addition the segments the compose the boundary:
#' mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)
#' 
#' ## Create mesh from data locations (without knowing the boundary)
#' mesh = create.mesh.2D(nodes = locations)
#' plot(mesh)
#' ## In this case the domain is the convex hull of the data locations. 
#' ## Do this only if you do not have any information about the shape of the domain of interest.  

create.mesh.2D <- function(nodes, nodesattributes = NA, segments = NA, holes = NA, triangles = NA, order = 1, verbosity = 0)
{ 
  ##########################
  ###   Input checking   ###
  ##########################
  
  # Triangle finds out which are on the border (see https://www.cs.cmu.edu/~quake/triangle.help.html)
  nodesmarkers = vector(mode = "integer", 0)
  segmentsmarkers = vector(mode = "integer", 0)
  
  nodes = as.matrix(nodes)
  if (ncol(nodes) != 2)
    stop("Matrix of nodes should have 2 columns")
  if (anyDuplicated(nodes))
    stop("Duplicated nodes")
  
  ## If attributes not specified, set them to a matrix with zero columns
  if (any(is.na(nodesattributes))) {
    nodesattributes <- matrix(0, nrow(nodes), 0)
  }else{
    nodesattributes <- as.matrix(nodesattributes)
    if (nrow(nodesattributes) != nrow(nodes))
      stop("Point attribute matrix \'nodesattributes\' does not have same number of rows the point matrix \'nodes\'")
  }
  
  ## If boundary nodes not specified, set them to 0
  #   if (any(is.na(nodesmarkers))) {
  #     nodesmarkers <- vector(mode = "integer", 0)
  #   }else{
  #     nodesmarkers = as.vector(nodesmarkers)
  #   }
  
  ## Deal with segments
  if (any(is.na(segments))) {
    segments <- matrix(0, 0, 2)
  } else {
    segments <- as.matrix(segments)
    if (ncol(segments) != 2) {
      stop("Matrix of segments should have 2 columns")
    }
  }
  
  ## If boundary segments not specified, set them to 0
  #   if (any(is.na(segmentsmarkers))) {
  #     segmentsmarkers <- vector(mode = "integer", 0)
  #   }else{
  #     segmentsmarkers = as.vector(segmentsmarkers)
  #   }
  
  ## If hole not specified, set it to empty matrix
  if (any(is.na(holes)))
    holes <- matrix(0, 0, 2)
  holes = as.matrix(holes)
  
  ## If triangles are not already specified
  if(any(is.na(triangles)))
    triangles = matrix(0,nrow = 0, ncol = 3)
  triangles = as.matrix(triangles)
  
  ## Set meshing parameters ##
  flags="ven"
  if(nrow(segments) == 0){
    flags = paste(flags,"c",sep = '')
  }
  
  if(nrow(segments)>0){
    flags = paste(flags,"p",sep = '')
  }
  
  #If order=2 add flag for second order nodes
  if(order == 2){
    flags = paste(flags,"o2",sep = '')
  }
  if(order < 1 || order >2){
    print('Order must be 1 or 2')
  }
  
  if(nrow(triangles) > 0){
    flags = paste(flags,"r",sep = '')
  }
  
  if (verbosity == 0) {
    flags = paste(flags,"Q",sep = '')
  }
  if (verbosity == 1) {
    flags = paste(flags,"V",sep = '')
  }
  if (verbosity == 2) {
    flags = paste(flags,"VV",sep = '')
  }
  
  out<-NULL
  #If triangles is null it makes the trianglulation
  #If triangle is not null it makes a refinement with no parameter, to compose the mesh object
  out <- triangulate_native(  
    nodes,
    nodesmarkers,
    nodesattributes,
    segments,
    segmentsmarkers,
    t(holes),
    triangles,
    flags
  )
  
  names(out)[1]<-"nodes"
  names(out)[2]<-"nodesmarkers"
  names(out)[3]<-"nodesattributes"
  names(out)[4]<-"triangles"
  names(out)[5]<-"segments"
  names(out)[6]<-"segmentsmarkers"
  names(out)[7]<-"edges"
  names(out)[8]<-"edgesmarkers"
  names(out)[9]<-"neighbors"
  
  out[13]<-NULL
  out[12]<-NULL
  out[11]<-NULL
  out[10]<-NULL
  
  out[[10]] = holes
  names(out)[10]<-"holes"
  out[[11]] = order
  names(out)[11]<-"order"
  
  
  class(out)<-"mesh.2D"
  
  return(out)
}

#' Refine a 2D triangular mesh
#' 
#' @param mesh A mesh.2D object representing the triangular mesh, created by \link{create.mesh.2D}.
#' @param minimum_angle A scalar specifying a minimun value for the triangles angles.
#' @param maximum_area A scalar specifying a maximum value for the triangles areas.
#' @param delaunay A boolean parameter indicating whether or not the output mesh should satisfy the Delaunay condition.
#' @param verbosity This can be '0', '1' or '2'. It indicates the level of verbosity in the triangulation process.
#' @description This function refines a Constrained Delaunay triangulation into a Conforming Delaunay triangulation. This is a wrapper of the Triangle library (http://www.cs.cmu.edu/~quake/triangle.html). It can be used to 
#' refine a mesh previously created with \link{create.mesh.2D}. The algorithm can add Steiner points (points through which the \code{segments} are splitted)
#' in order to meet the imposed refinement conditions.
#' @usage refine.mesh.2D(mesh, minimum_angle, maximum_area, delaunay, verbosity)
#' @seealso \code{\link{create.mesh.2D}}, \code{\link{create.FEM.basis}}
#' @return A mesh.2D object representing the refined triangular mesh,  with the following output:
#' \itemize{
#' \item{\code{nodes}}{A #nodes-by-2 matrix containing the x and y coordinates of the mesh nodes.}
#' \item{\code{nodesmarkers}}{A vector of length #nodes, with entries either '1' or '0'. An entry '1' indicates that the corresponding node is a boundary node; an entry '0' indicates that the corresponding node is not a boundary node.}
#' \item{\code{nodesattributes}}{nodesattributes A matrix with #nodes rows containing nodes' attributes. 
#' These are passed unchanged to the output. If a node is added during the triangulation process or mesh refinement, its attributes are computed  
#' by linear interpolation using the attributes of neighboring nodes. This functionality is for instance used to compute the value 
#' of a Dirichlet boundary condition at boundary nodes added during the triangulation process.}
#' \item{\code{triangles}}{A #triangles-by-3 (when \code{order} = 1) or #triangles-by-6 (when \code{order} = 2) matrix.}
#' \item{\code{edges}}{A #edges-by-2 matrix. Each row contains the row's indices of the nodes where the edge starts from and ends to.}
#' \item{\code{edgesmarkers}}{A vector of lenght #edges with entries either '1' or '0'. An entry '1' indicates that the corresponding element in \code{edge} is a boundary edge;  
#' an entry '0' indicates that the corresponding edge is not a boundary edge.}
#' \item{\code{neighbors}}{A #triangles-by-3 matrix. Each row contains the indices of the three neighbouring triangles. An entry '-1' indicates that 
#' one edge of the triangle is a boundary edge.}
#' \item{\code{holes}}{A #holes-by-2 matrix containing the x and y coordinates of a point internal to each hole of the mesh. These points are used to carve holes
#' in the triangulation, when the domain has holes.}
#' \item{\code{order}}{Either '1' or '2'. It specifies wether each mesh triangle should be represented by 3 nodes (the triangle' vertices) or by 6 nodes (the triangle's vertices and midpoints).
#' These are respectively used for linear (order = 1) and quadratic (order = 2) Finite Elements.}
#' }
#' @examples 
#' library(fdaPDE)
#' 
#' ## Upload the quasicircle2D data
#' data(quasicircle2D)
#' boundary_nodes = quasicircle2D$boundary_nodes
#' boundary_segments = quasicircle2D$boundary_segments
#' locations = quasicircle2D$locations
#' data = quasicircle2D$data
#' 
#' ## Create mesh from boundary:
#' mesh = create.mesh.2D(nodes = boundary_nodes, segments = boundary_segments)
#' plot(mesh)
#' ## Refine the mesh with the maximum area criterion:
#' finemesh = refine.mesh.2D(mesh = mesh, maximum_area = 0.1)
#' plot(finemesh)
#' ## Refine the mesh with the minimum angle criterion:
#' finemesh2 = refine.mesh.2D(mesh = mesh, minimum_angle = 30)
#' plot(finemesh2)
#' @export

refine.mesh.2D<-function(mesh, minimum_angle = NA, maximum_area = NA, delaunay = FALSE, verbosity = 0)
{ 
  if(class(mesh) !="mesh.2D")
    stop("Sorry, this function is implemented just for mesh.2D class ")
  
  flags="rpven" 
  
  if(!is.na(minimum_angle)){
    flags <- paste(flags, "q", sprintf("%.12f", minimum_angle), sep='')
  }
  
  if(!is.na(maximum_area)){
    flags <- paste(flags, "a", sprintf("%.12f", maximum_area), sep='')
  }
  
  if(delaunay){
    flags <- paste(flags, "D", sep='')
  }
  
  if(mesh$order==2){
    flags <- paste(flags, "o2", sep='')
  }
  
  if (verbosity == 0) {
    flags = paste(flags,"Q",sep = '')
  }
  if (verbosity == 1) {
    flags = paste(flags,"V",sep = '')
  }
  if (verbosity == 2) {
    flags = paste(flags,"VV",sep = '')
  }
  
  # Triangle finds out which are on the border (see https://www.cs.cmu.edu/~quake/triangle.help.html)
  mesh$nodesmarkers = vector(mode = "integer", 0)
  mesh$segmentsmarkers = vector(mode = "integer", 0)
  
  out<-NULL
  #If triangles is null it makes the trianglulation
  #If triangle is not null it makes a refinement with no parameter, to compose the mesh object
  out <- triangulate_native(
    mesh$nodes,
    mesh$nodesmarkers,
    mesh$nodesattributes,
    mesh$segments,
    mesh$segmentmarkers,
    t(mesh$holes),
    mesh$triangles,
    flags
  )
  
  names(out)[1]<-"nodes"
  names(out)[2]<-"nodesmarkers"
  names(out)[3]<-"nodesattributes"
  names(out)[4]<-"triangles"
  names(out)[5]<-"segments"
  names(out)[6]<-"segmentsmarkers"
  names(out)[7]<-"edges"
  names(out)[8]<-"edgesmarkers"
  names(out)[9]<-"neighbors"
  
  out[13]<-NULL
  out[12]<-NULL
  out[11]<-NULL
  out[10]<-NULL
  
  out[[10]] = mesh$holes
  names(out)[10]<-"holes"
  out[[11]] = mesh$order
  names(out)[11]<-"order"
  
  class(out)<-"mesh.2D"
  
  return(out)
}

#' Create a \code{mesh.2.5D} object from the nodes locations and the connectivty matrix 
#'
#' @param nodes A #nodes-by-3 matrix specifying the locations of each node.
#' @param triangles A #triangles-by-3 (when \code{order} = 1) or #triangles-by-6 (when \code{order} = 2) matrix, 
#' specifying the indices of the nodes in each triangle.
#' @param order Either '1' or '2'. It specifies wether each mesh triangle should be represented by 3 nodes 
#' (the triangle' vertices) or by 6 nodes (the triangle's vertices and midpoints). 
#' These are respectively used for linear (order = 1) and quadratic (order = 2) Finite Elements. 
#' Default is \code{order} = 1.
#' @return An object of the class \code{mesh.2.5D} with the following output:
#' \itemize{
#' \item{\code{nnodes}}{The #nodes in the mesh.}
#' \item{\code{ntriangles}}{The #triangles in the mesh.}
#' \item{\code{nodes}}{A #nodes-by-3 matrix containing the x,y and z coordinate for each point of the mesh.}
#' \item{\code{triangles}}{A #triangles-by-3 (when \code{order} = 1) or #triangles-by-6 (when \code{order} = 2) matrix, 
#' specifying the indices of the nodes in each triangle.}
#' \item{\code{order}}{Either '1' or '2'. It specifies wether each mesh triangle should be represented by 3 nodes 
#' (the triangle' vertices) or by 6 nodes (the triangle's vertices and midpoints). 
#' It is passed unchanged from the input.}
#' }
#' @export
#' @examples
#' library(fdaPDE)
#'
#' ## Upload the hub2.5D the data
#' data(hub2.5D)
#' hub2.5D.nodes = hub2.5D$hub2.5D.nodes
#' hub2.5D.triangles = hub2.5D$hub2.5D.triangles
#'
#' ## Create mesh from nodes and connectivity matrix:
#' mesh = create.mesh.2.5D(nodes = hub2.5D.nodes, triangles = hub2.5D.triangles)
#' plot(mesh)

create.mesh.2.5D<- function(nodes, triangles, order = 1)
{
  nnodes = dim(nodes)[1]
  
  ntriangles = dim(triangles)[1]
  
  if(order==1 & dim(triangles)[2]== 3){
    out = list(nnodes=nnodes, ntriangles=ntriangles, nodes=as.matrix(nodes), triangles = as.matrix(triangles), order=as.integer(order))
  }
  else if(order==2 & dim(triangles)[2]==6){ # triangles matrix contains both the true triangles and the midpoints ones
    out = list(nnodes=nnodes, ntriangles=ntriangles, nodes=as.matrix(nodes), triangles = as.matrix(triangles), order=as.integer(order))
  }
  else if(order==2 & dim(triangles)[2]==3){ # need to compute the midpoints
    print("You set order=2 but passed a matrix of triangles with just 3 columns. The midpoints for each edge will be computed.")
    out=second.order.mesh.2.5D(as.matrix(nodes),as.matrix(triangles))
  }
  else{
    stop("The number of columns of triangles matrix is not consistent with the order parameter")
  }
  
  class(out)<-"mesh.2.5D"
  
  return(out)
}

#' Project 3D points onto 2D 2.5D triangular mesh
#'
#' @param mesh A mesh.2.5D object representing the triangular mesh, created by \link{create.mesh.2.5D}.
#' @param locations 3D points to be projected onto 2.5D triangular mesh.
#' @description This function projects any 3D points onto 2.5D triangular mesh.
#' @return 3D points projected onto 2.5D triangluar mesh.
#' @export
#' @examples
#' library(fdaPDE)
#'
#' ## Upload the hub2.5D the data
#' data(hub2.5D)
#' hub2.5D.nodes = hub2.5D$hub2.5D.nodes
#' hub2.5D.triangles = hub2.5D$hub2.5D.triangles
#'
#' ## Create mesh
#' mesh = create.mesh.2.5D(nodes = hub2.5D.nodes, triangles = hub2.5D.triangles)
#' 
#' ## Create 3D points to be projected
#' x <- cos(seq(0,2*pi, length.out = 9))
#' y <- sin(seq(0,2*pi, length.out = 9))
#' z <- rep(0.5, 9)
#' locations = cbind(x,y,z)
#' 
#' ## Project the points on the mesh
#' loc = projection.points.2.5D(mesh, locations)

projection.points.2.5D<-function(mesh, locations) {
  if(class(mesh) !="mesh.2.5D")
  stop("Data projection is only available for 2.5D mesh ")

  if (mesh$order == 2)
    stop("Data projection is only available for order 1 ")

  # C++ function for 2.5D works with vectors not with matrices
  mesh$triangles=c(t(mesh$triangles))
  mesh$nodes=c(t(mesh$nodes))

  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  mesh$triangles=mesh$triangles-1

  # Imposing types, this is necessary for correct reading from C++
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(mesh$nnodes) <- "integer"
  storage.mode(mesh$ntriangles) <- "integer"
  storage.mode(mesh$nodes) <- "double"
  storage.mode(mesh$triangles) <- "integer"
  storage.mode(mesh$order) <- "integer"

  ## Call C++ function
  evalmat <- .Call("points_projection", mesh, locations, PACKAGE = "fdaPDE")

  #Returning the evaluation matrix
  return(evalmat)


}

second.order.mesh.2.5D<-function(nodes, triangles){
  
  toll=1e-5
  # T = matrix(mesh$triangles,nrow=mesh$ntriangles,ncol=3, byrow = TRUE)
  # V = matrix(mesh$nodes, nrow = mesh$nnodes, ncol= 3, byrow = TRUE)
  T=triangles
  V=nodes
  T <- cbind(T, matrix(0,nrow=nrow(T),ncol=3))
  nnodes=nrow(V)
  index=nrow(V)
  points = V[T[1,],]
  midpoints<-rbind((points[2,]+points[3,])/2,(points[1,]+points[3,])/2, (points[1,]+points[2,])/2);
  # if(!is.null(bc)){
  #   isBC<-c( any(bc==T[1,2]) & any(bc==T[1,3]),
  #            any(bc==T[1,1]) & any(bc==T[1,3]),
  #            any(bc==T[1,2]) & any(bc==T[1,1]))
  # }
  
  for (side in 1:3){
    point<-midpoints[side,]
    index<-index+1;
    V<-rbind(V,point)
    T[1,3+side]<-index;
    
    # if(!is.null(bc)&&isBC[side]==1){
    #   bc<-c(bc,index)
    # }
    
  }
  
  for (i in 2:nrow(T)){
    points = V[T[i,],]
    midpoints<-rbind((points[2,]+points[3,])/2,(points[1,]+points[3,])/2, (points[1,]+points[2,])/2);
    # if(!is.null(bc)){
    #   isBC<-c( any(bc==T[i,2]) & any(bc==T[i,3]),
    #            any(bc==T[i,1]) & any(bc==T[i,3]),
    #            any(bc==T[i,2]) & any(bc==T[i,1]))
    # }
    
    for (side in 1:3){
      point<-midpoints[side,]
      isthere<-apply(V[(nnodes+1):nrow(V),], 1, function(x) identical(as.vector(x), point))
      loc = which(isthere)
      if(length(loc)>0){
        loc = loc+nnodes
        T[i,3+side]<-loc[1]
      }else{
        index<-index+1;
        V<-rbind(V,point)
        T[i,3+side]<-index;
        # 
        #           if(!is.null(bc)&&isBC[side]==1){
        #             bc<-c(bc,index)
        #          }
      }
    }
  }
  #}
  # if(is.null(bc)){
  #   out = list(nnodes=nrow(V), ntriangles=nrow(T), nodes=V, triangles = T, order=2)
  #   class(out)<-"mesh.2.5D"
  #   return(out)
  #  }else{
  out = list(nnodes=nrow(V), ntriangles=nrow(T), nodes=V, triangles = T, order=2)
  #    class(out)<-"mesh.2.5D"
  #    retlist = list(mesh = out, bc_index=bc)
  return(out)
  
}

#' Create a \code{mesh.3D} object from the connectivty matrix and nodes locations
#'
#' @param nodes A #nodes-by-3 matrix specifying the locations of each node.
#' @param tetrahedrons A #tetrahedrons-by-4 matrix specifying the indices of the nodes in each tetrahedrons.
#' @param order Order of the Finite Element basis. Only order = 1 is currently implemented.
#' @return An object of the class \code{mesh.3D} with the following output:
#' \itemize{
#' \item{\code{nnodes}}{The #nodes in the mesh.}
#' \item{\code{ntetrahedrons}}{The #tetrahedrons in the mesh.}
#' \item{\code{nodes}}{A #nodes-by-3 matrix containing the x,y and z coordinate for each point of the mesh.}
#' \item{\code{tetrahedrons}}{A #tetrahedrons-by-4 matrix specifying the indices of the nodes in each tetrahedron of the mesh.}
#' \item{\code{order}}{It specifies the order of the associated Finite Element basis. When order = 1, each 
#' mesh tetrahedron is represented by 4 nodes (the tetrahedron vertices).}
#' }
#' @export
#' @examples
#' library(fdaPDE)
#'
#' ##Load the matrix nodes and tetrahedrons
#' data(sphere3Ddata)
#'
#' nodes = sphere3Ddata$nodes
#' tetrahedrons = sphere3Ddata$tetrahedrons
#'
#' ##Create the triangulated mesh from the connectivity matrix and nodes locations
#' mesh=create.mesh.3D(nodes,tetrahedrons)

create.mesh.3D<- function(nodes, tetrahedrons, order = 1)
{
  if(order!=1)
    stop("Sorry, for the moment only order 1 mesh creation is available. Please set order=1")
  
  nnodes = dim(nodes)[1]
  
  ntetrahedrons = dim(tetrahedrons)[1]
  
  if(dim(tetrahedrons)[2]!= 4*order){
    if (order==1)
      stop("The matrix 'tetrahedrons' has the wrong number of columns. See second.order.mesh(...)")
    stop("The matrix 'tetrahedrons' has wrong number of columns. Should be 4*order \n")
  }
  out = list(nnodes=nnodes, ntetrahedrons=ntetrahedrons, nodes=as.matrix(nodes), tetrahedrons = as.matrix(tetrahedrons), order=as.integer(order))
  
  class(out)<-"mesh.3D"
  
  return(out)
}
