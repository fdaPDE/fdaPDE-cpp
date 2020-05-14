#' Horseshoe domain
#'
#' The boundary and interior nodes and connectivity matrix of a triangular mesh of the horseshoe domain. 
#' This dataset can be used to create a \code{MESH.2D} object with the function \code{create.MESH.2D}.
#'
#' @name horseshoe2D
NULL


#' Quasicircle2D domain
#'
#' The boundary and interior nodes and connectivity matrix of a triangular mesh of a quasicircular domain, together 
#' with a non-stationary field observed over the nodes of the mesh.  
#' This dataset can be used to create a \code{MESH.2D} object with the function \code{create.MESH.2D} and to test
#' the smooth.FEM function.
#'
#' @name quasicircle2D
NULL


#' Quasicircle2Dareal domain
#'
#' The mesh of a quasicircular domain, together with a non-stationary field observed over seven circular subdomains and 
#' the incindence matrix defining the subdomains used by Azzimonti et. al 2015.   
#' This dataset can be used to test the smooth.FEM function for areal data. 
#'
#' @references Azzimonti, L., Sangalli, L. M., Secchi, P., Domanin, M., & Nobile, F. (2015). Blood flow velocity 
#' field estimation via spatial regression with PDE penalization. Journal of the American Statistical 
#' Association, 110(511), 1057-1071.
#' @name quasicircle2Dareal
NULL


#' Hub domain
#'
#' The nodes and connectivity matrix of a triangular mesh of a manifold representing a hub geometry. 
#' This dataset can be used to create a \code{MESH.2.5D} object with the function \code{create.MESH.2.5D}.
#'
#' @name hub2.5D
NULL


#' Sphere3Ddata
#'
#' A dataset with information about the connectivity matrix and the nodes locations of a sphere geometry. It containes:
#' \itemize{
#' 	\item nodes. A #nodes-by-3 matrix specifying the locations of each node.
#' 	\item tetrahedrons. A #tetrahedrons-by-4 matrix specifying the indices of the nodes in each tetrahedron.
#'	   }
#' This dataset can be used to create a \code{MESH.3D} object with the function \code{create.MESH.3D} 
#'
#' @name sphere3Ddata
NULL
