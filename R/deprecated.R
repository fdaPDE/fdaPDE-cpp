#' Deprecated Functions 
#' 
#' These functions are Deprecated in this release of fdaPDE, they will be 
#' marked as Defunct and removed in a future version. 
#' @name fdaPDE-deprecated


#' @param FEMbasis A \code{FEM} object representing the Finite Element basis. See \code{\link{create.FEM.basis}}.
#' @return A square matrix with the integrals of all the basis' functions pairwise products.
#' The dimension of the matrix is equal to the number of the nodes of the mesh.
#' @description Only executed when \code{smooth.FEM.basis} is run with the option  \code{CPP_CODE} = \code{FALSE}. It computes the mass matrix. The element (i,j) of this matrix contains the integral over the domain of the product between the ith and kth element 
#' of the Finite Element basis. As common practise in Finite Element Analysis, this quantities are computed iterating over all the mesh triangles. 
#' @rdname fdaPDE-deprecated
#' @export

R_mass=function(FEMbasis)
{
  .Deprecated("smooth.FEM",msg = 'This function is deprecated and will soon be removed. To perform smooth regression please refer to smooth.FEM')
  nodes = FEMbasis$mesh$nodes
  triangles = FEMbasis$mesh$triangles
  detJ = FEMbasis$detJ
  order = FEMbasis$order
  
  nele  = dim(triangles)[1]
  nnod  = dim(nodes)[1]
  
  if (order < 1 || order > 2){
    stop("ORDER not 1 or 2")
  }
  
  K0M = NULL
  
  if (order ==2)
  {   
    #  the integrals of products of basis functions for master element:
    
    K0M = matrix(c( 6, -1, -1, -4,  0,  0,
                    -1,  6, -1,  0, -4,  0,
                    -1, -1,  6,  0,  0, -4,
                    -4,  0,  0, 32, 16, 16,
                    0, -4,  0, 16, 32, 16,
                    0,  0, -4, 16, 16, 32), ncol=6, nrow=6, byrow=T)/360
  }else if (order == 1)
  {  
    #  the integrals of products of basis functions for master element:
    K0M = matrix( c( 2,  1,  1,
                     1,  2,  1,
                     1 , 1,  2), ncol=3, nrow=3, byrow=T) / 24
    
    
  }
  
  # assemble the mass matrix
  K0 = matrix(0,nrow=nnod,ncol=nnod)
  for (el in 1:nele)
  {  
    ind = triangles[el,]
    K0[ind,ind] = K0[ind,ind] + K0M * detJ[el]
  }
  
  K0
}


#' @param FEMbasis A \code{FEMbasis} object representing the basis; See \code{\link{create.FEM.basis}}.
#' @return A square matrix with the integrals of all the basis functions' gradients pairwise dot products.
#' The dimension of the matrix is equal to the number of the nodes of the mesh.
#' @description Only executed when \code{smooth.FEM.basis} is run with the option  \code{CPP_CODE} = \code{FALSE}. It computes the stifness matrix. The element (i,j) of this matrix contains the integral over the domain of the scalar product between the gradient of the ith and kth element 
#' of the Finite Element basis. As common practise in Finite Element Analysis, this quantities are computed iterating over all the mesh triangles. 
#' @rdname fdaPDE-deprecated
#' @export

R_stiff= function(FEMbasis)
{
  .Deprecated("smooth.FEM",msg = 'This function is deprecated and will soon be removed. To perform smooth regression please refer to smooth.FEM')
  
  nodes = FEMbasis$mesh$nodes
  triangles = FEMbasis$mesh$triangles
  
  nele  = dim(triangles)[[1]]
  nnod  = dim(nodes)[[1]]
  detJ     = FEMbasis$detJ
  order = FEMbasis$order
  metric = FEMbasis$metric
  
  
  
  KXX = NULL
  KXY = NULL
  KYY = NULL
  
  if (order < 1 || order > 2){
    stop("ORDER not 1 or 2")
  }
  
  if (order == 2)
  {   
    #  values of K1 for master elements
    KXX = matrix( c( 3, 1, 0, 0, 0,-4, 
                     1, 3, 0, 0, 0,-4,
                     0, 0, 0, 0, 0, 0,
                     0, 0, 0, 8,-8, 0,
                     0, 0, 0,-8, 8, 0,
                     -4,-4, 0, 0, 0, 8), ncol=6, nrow=6, byrow=T)/6
    
    KXY = matrix(c( 3, 0, 1, 0,-4, 0, 
                    1, 0,-1, 4, 0,-4, 
                    0, 0, 0, 0, 0, 0,
                    0, 0, 4, 4,-4,-4,
                    0, 0,-4,-4, 4, 4,
                    -4, 0, 0,-4, 4, 4), ncol=6, nrow=6, byrow=T)/6
    
    KYY = matrix( c( 3, 0, 1, 0,-4, 0,
                     0, 0, 0, 0, 0, 0,
                     1, 0, 3, 0,-4, 0,
                     0, 0, 0, 8, 0,-8,
                     -4, 0,-4, 0, 8, 0,
                     0, 0, 0,-8, 0, 8), ncol=6, nrow=6, byrow=T)/6
  }else if (order == 1)
  {   
    KXX = matrix( c(  1, -1,  0,
                      -1,  1,  0,
                      0,  0,  0), ncol=3, nrow=3, byrow=T) /2
    
    KXY = matrix( c(  1,  0, -1,
                      -1,  0,  1,
                      0,  0,  0), ncol=3, nrow=3, byrow=T) /2
    
    KYY = matrix( c(  1,  0, -1,
                      0,  0,  0,
                      -1,  0,  1), ncol=3, nrow=3, byrow=T) /2
  }
  
  #  assemble the stiffness matrix
  K1   = matrix(0,nrow=nnod,ncol=nnod)
  for (el in 1:nele)
  {
    ind    = triangles[el,]
    K1M = (metric[el,1,1]*KXX    + metric[el,1,2]*KXY +
             metric[el,2,1]*t(KXY) + metric[el,2,2]*KYY)
    K1[ind,ind] = K1[ind,ind] + K1M*detJ[el]
  }
  
  K1
}


#' @param observations A #observations vector with the observed data values over the domain. The locations of the observations can be specified with the \code{locations} argument.
#' Otherwise if only the vector of observations is given, these are consider to be located in the corresponding node in the table nodes of the mesh. In this last
#' case, an \code{NA} value in the observations vector indicates that there is no observation associated to the corresponding node.
#' @param locations A #observations-by-2 matrix where each row specifies the spatial coordinates of the corresponding observations in the vector \code{observations}.
#' @param FEMbasis A F\code{EMbasis} object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param lambda A scalar or vector of smoothing parameters.
#' @param covariates A #observations-by-#covariates matrix where each row represents the covariates associated with the corresponding observed data value in \code{observations}.
#' @param GCV Boolean. If \code{TRUE} the following quantities are computed: the trace of the smoothing matrix, the estimated error standard deviation,  and 
#'        the Generalized Cross Validation criterion, for each value of the smoothing parameter specified in \code{lambda}.
#' @return A list with the following quantities:
#'    \item{\code{fit.FEM}}{A \code{FEM} object that represents the fitted spatial field.}
#'    \item{\code{PDEmisfit.FEM}}{A \code{FEM} object that represents the Laplacian of the estimated spatial field.}
#'    \item{\code{beta}}{If covariates is not \code{NULL}, a vector of length #covariates with the regression coefficients associated with each covariate.}
#'    \item{\code{edf}}{If GCV is \code{TRUE}, a scalar or vector with the trace of the smoothing matrix for each value of the smoothing parameter specified in \code{lambda}.}
#'    \item{\code{stderr}}{If GCV is \code{TRUE}, a scalar or vector with the estimate of the standard deviation of the error for each value of the smoothing parameter specified in \code{lambda}.}
#'    \item{\code{GCV}}{If GCV is \code{TRUE}, a  scalar or vector with the value of the GCV criterion for each value of the smoothing parameter specified in \code{lambda}.}
#' @rdname fdaPDE-deprecated
#' @export


R_smooth.FEM.basis = function(locations, observations, FEMbasis, lambda, covariates = NULL, GCV)
{
  .Deprecated("smooth.FEM",msg = 'This function is deprecated and will soon be removed. To perform smooth regression please refer to smooth.FEM')
  
  # Stores the number of nodes of the mesh. This corresponds to the number of elements of the FE basis.
  numnodes = nrow(FEMbasis$mesh$nodes)
  if(!is.null(covariates))
  {
    covariates = as.matrix(covariates)
  }
  
  if(!is.null(locations))
  {
    locations <- as.matrix(locations)
  }
  
  #  ---------------------------------------------------------------
  # construct mass matrix K0 
  #  ---------------------------------------------------------------
  
  K0 = R_mass(FEMbasis)
  
  #  ---------------------------------------------------------------
  # construct stiffness matrix K1
  #  ---------------------------------------------------------------
  
  K1 = R_stiff(FEMbasis)
  
  
  #  ---------------------------------------------------------------
  # construct the penalty matrix P with ones on diagonal at data points
  #  ---------------------------------------------------------------
  
  #penalty = numeric(numnodes)
  #penalty[data[,1]] = 1
  #P = diag(penalty,nrow=numnodes)
  
  basismat = NULL
  Lmat = matrix(0,nrow = numnodes, ncol = numnodes)
  H = NULL
  Q = NULL
  
  if(!is.null(locations))
  {
    basismat = R_eval.FEM.basis(FEMbasis, locations)
  } 
  
  if(!is.null(covariates))
  {
    if(!is.null(locations))
    {
      H = covariates %*% ( solve( t(covariates) %*% covariates ) ) %*% t(covariates)
    }else{
      loc_nodes = (1:length(observations))[!is.na(observations)]
      H = covariates[loc_nodes,] %*% ( solve( t(covariates[loc_nodes,]) %*% covariates[loc_nodes,] ) ) %*% t(covariates[loc_nodes,])
    }
    Q=matrix(0,dim(H)[1], dim(H)[2])
    Q = diag(1,dim(H)[1])- H
  }
  
  if(!is.null(locations) && is.null(covariates))
  {
    Lmat = t(basismat) %*% basismat
  }
  if(!is.null(locations) && !is.null(covariates))
  {
    Lmat = t(basismat) %*% Q %*% basismat
  }
  if(is.null(locations) && is.null(covariates))
  {
    loc_nodes = (1:length(observations))[!is.na(observations)]
    Lmat[loc_nodes,loc_nodes]=diag(1,length(observations[loc_nodes]))
  }
  if(is.null(locations) && !is.null(covariates))
  {
    loc_nodes = (1:length(observations))[!is.na(observations)]
    Lmat[loc_nodes,loc_nodes] = Q
  }
  
  #  ---------------------------------------------------------------
  # construct vector b for system Ax=b
  #  ---------------------------------------------------------------
  
  b = matrix(numeric(numnodes*2),ncol=1)
  if(!is.null(locations) && is.null(covariates))
  {
    b[1:numnodes,] = t(basismat) %*% observations
  }
  if(!is.null(locations) && !is.null(covariates))
  {
    b[1:numnodes,] = t(basismat) %*% Q %*% observations
  }
  if(is.null(locations) && is.null(covariates))
  {
    loc_nodes = (1:length(observations))[!is.na(observations)]
    b[loc_nodes,] = observations[loc_nodes]
  }
  if(is.null(locations) && !is.null(covariates))
  {
    loc_nodes = (1:length(observations))[!is.na(observations)]
    b[loc_nodes,] = Q %*% observations[loc_nodes]
  }
  #print(b)
  
  #  ---------------------------------------------------------------
  # construct matrix A for system Ax=b.
  #  ---------------------------------------------------------------
  
  bigsol = list(solutions = matrix(nrow = 2*numnodes, ncol = length(lambda)), edf = vector(mode="numeric",length = length(lambda)))
  for(i in 1: length(lambda))
  {
    A  = rbind(cbind(Lmat, -lambda[i]*K1), cbind(K1, K0))
    #print(A)
    # solve system
    bigsol[[1]][,i] = solve(A,b)
    if(GCV == TRUE)
    {
      S = solve(Lmat + lambda[i] * K1 %*% solve(K0) %*% K1)
      if(is.null(locations))
      {
        loc_nodes = (1:length(observations))[!is.na(observations)]
        if(is.null(covariates))
        {
          bigsol[[2]][i]  = sum(diag(S[loc_nodes,loc_nodes]))
        }else{
          bigsol[[2]][i]  = ncol(covariates)+sum(diag(S[loc_nodes,loc_nodes]%*%Q))
        }
      }else{
        if(is.null(covariates))
        {
          bigsol[[2]][i]  = sum(diag(basismat%*%S[1:numnodes,1:numnodes]%*%t(basismat)))
        }else{
          bigsol[[2]][i]  = ncol(covariates)+sum(diag(basismat%*%S[1:numnodes,1:numnodes]%*%t(basismat)%*%Q))
        }
      }
    }else{
      bigsol[[2]][i]  = NA
    }
  }
  return(bigsol)
}


#' @param locations A 2-columns matrix with the spatial locations where the bases should be evaluated.
#' @param FEMbasis An \code{FEMbasis} object representing the Finite Element basis; See \link{create.FEM.basis}.
#' @param nderivs A vector of lenght 2 specifying the order of the partial derivatives of the bases to be evaluated. The vectors' entries can
#' be 0,1 or 2, where 0 indicates that only the basis functions, and not their derivatives, should be evaluated.
#' @return 
#' A matrix of basis function values. Each row indicates the location where the evaluation has been taken, the column indicates the 
#' basis function evaluated 
#' @description Only executed when the function \code{smooth.FEM.basis} is run with the option \code{CPP_CODE} = \code{FALSE}. It evaluates the Finite Element basis functions and their derivatives up to order 2 at the specified set of locations. 
#' This version of the function is implemented using only R code. It is called by \link{R_smooth.FEM.basis}.
#' @rdname fdaPDE-deprecated
#' @export

R_eval.FEM.basis <- function(FEMbasis, locations, nderivs = matrix(0,1,2))
{ 
  .Deprecated("eval.FEM")
  
  if(length(nderivs) != 2)
  {
    stop('NDERIVS not of length 2.')
  }
  
  if(sum(nderivs)>2)
  {
    stop('Maximum derivative order is greater than two.')
  }
  
  N = nrow(locations)
  nbasis = FEMbasis$nbasis
  
  # Augment Xvec and Yvec by ones for computing barycentric coordinates
  Pgpts = cbind(matrix(1,N,1),locations[,1],locations[,2])
  
  # Get nodes and index
  
  mesh = FEMbasis$mesh
  
  nodes = mesh$nodes
  triangles = mesh$triangles
  
  order = FEMbasis$order
  #nodeindex = params$nodeindex
  detJ = FEMbasis$detJ
  
  # 1st, 2nd, 3rd vertices of triangles
  
  v1 = nodes[triangles[,1],]
  v2 = nodes[triangles[,2],]
  v3 = nodes[triangles[,3],]
  
  if(order !=2 && order != 1)
  {
    stop('ORDER is neither 1 or 2.')
  }
  
  # Denominator of change of coordinates chsange matrix
  
  modJ = FEMbasis$detJ
  ones3 = matrix(1,3,1)
  modJMat = modJ %*% t(ones3)
  
  M1 = cbind(v2[,1]*v3[,2] - v3[,1]*v2[,2], v2[,2] - v3[,2], v3[,1] - v2[,1])/(modJMat)
  M2 = cbind(v3[,1]*v1[,2] - v1[,1]*v3[,2], v3[,2] - v1[,2], v1[,1] - v3[,1])/(modJMat)
  M3 = cbind(v1[,1]*v2[,2] - v2[,1]*v1[,2], v1[,2] - v2[,2], v2[,1] - v1[,1])/(modJMat)
  
  ind = matrix(0,N,1)
  for(i in 1:N)
  {
    ind[i] = R_insideIndex(mesh, locations[i,])
  }
  
  evalmat = matrix(0,N,nbasis)
  
  for(i in 1:N)
  {
    indi = ind[i]
    
    if(!is.nan(indi))
    {
      baryc1 = (M1[indi,]*Pgpts[i,]) %*% ones3
      baryc2 = (M2[indi,]*Pgpts[i,]) %*% ones3
      baryc3 = (M3[indi,]*Pgpts[i,]) %*% ones3
      
      if(order == 2)
      {
        if(sum(nderivs) == 0)
        {
          evalmat[i,triangles[indi,1]] = 2* baryc1^2 - baryc1
          evalmat[i,triangles[indi,2]] = 2* baryc2^2 - baryc2
          evalmat[i,triangles[indi,3]] = 2* baryc3^2 - baryc3
          evalmat[i,triangles[indi,6]] = 4* baryc1 * baryc2
          evalmat[i,triangles[indi,4]] = 4* baryc2 * baryc3
          evalmat[i,triangles[indi,5]] = 4* baryc3 * baryc1
        }
        else if(nderivs[1] == 1 && nderivs[2] == 0)
        {
          evalmat[i,triangles[indi,1]] = (4* baryc1 - 1) * M1[indi,2]
          evalmat[i,triangles[indi,2]] = (4* baryc2 - 1) * M2[indi,2]
          evalmat[i,triangles[indi,3]] = (4* baryc3 - 1) * M3[indi,2]
          evalmat[i,triangles[indi,6]] = (4* baryc2 ) * M1[indi,2] + 4*baryc1 * M2[indi,2]
          evalmat[i,triangles[indi,4]] = (4* baryc3 ) * M2[indi,2] + 4*baryc2 * M3[indi,2]
          evalmat[i,triangles[indi,5]] = (4* baryc1 ) * M3[indi,2] + 4*baryc3 * M1[indi,2]
        }
        else if(nderivs[1] == 0 && nderivs[2] == 1)
        {
          evalmat[i,triangles[indi,1]] = (4*baryc1 - 1)*M1[indi,3]
          evalmat[i,triangles[indi,2]] = (4*baryc2 - 1)*M2[indi,3]
          evalmat[i,triangles[indi,3]] = (4*baryc3 - 1)*M3[indi,3]
          evalmat[i,triangles[indi,6]] = 4*baryc2*M1[indi,3] + 4*baryc1*M2[indi,3]
          evalmat[i,triangles[indi,4]] = 4*baryc3*M2[indi,3] + 4*baryc2*M3[indi,3]
          evalmat[i,triangles[indi,5]] = 4*baryc1*M3[indi,3] + 4*baryc3*M1[indi,3]
        }
        else if(nderivs[1] == 1 && nderivs[2] == 1)
        {
          evalmat[i,triangles[indi,1]] = 4*M1[indi,2]%*%M1[indi,3];
          evalmat[i,triangles[indi,2]] = 4*M2[indi,2]%*%M2[indi,3];
          evalmat[i,triangles[indi,3]] = 4*M3[indi,2]%*%M3[indi,3];
          evalmat[i,triangles[indi,6]] = 4*M2[indi,2]%*%M1[indi,3] + 4*M2[indi,3]%*%M1[indi,2];
          evalmat[i,triangles[indi,4]] = 4*M3[indi,2]%*%M2[indi,3] + 4*M3[indi,3]%*%M2[indi,2];
          evalmat[i,triangles[indi,5]] = 4*M1[indi,2]%*%M3[indi,3] + 4*M1[indi,3]%*%M3[indi,2];
        }
        else if(nderivs[1] == 2 && nderivs[2] == 0)
        {
          evalmat[i,triangles[indi,1]] = 4*M1[indi,2]%*%M1[indi,2];
          evalmat[i,triangles[indi,2]] = 4*M2[indi,2]%*%M2[indi,2];
          evalmat[i,triangles[indi,3]] = 4*M3[indi,2]%*%M3[indi,2];
          evalmat[i,triangles[indi,6]] = 8*M2[indi,2]%*%M1[indi,2];
          evalmat[i,triangles[indi,4]] = 8*M3[indi,2]%*%M2[indi,2];
          evalmat[i,triangles[indi,5]] = 8*M1[indi,2]%*%M3[indi,2];
        }
        else if(nderivs[1] == 0 && nderivs[2] == 2)
        {
          evalmat[i,triangles[indi,1]] = 4*M1[indi,3]%*%M1[indi,3];
          evalmat[i,triangles[indi,2]] = 4*M2[indi,3]%*%M2[indi,3];
          evalmat[i,triangles[indi,3]] = 4*M3[indi,3]%*%M3[indi,3];
          evalmat[i,triangles[indi,6]] = 8*M2[indi,3]%*%M1[indi,3];
          evalmat[i,triangles[indi,4]] = 8*M3[indi,3]%*%M2[indi,3];
          evalmat[i,triangles[indi,5]] = 8*M1[indi,3]%*%M3[indi,3];
        }
        else
        {
          stop('Inadmissible derivative orders.')
        }
      }
      else
      {
        if(sum(nderivs) == 0)
        {
          evalmat[i,triangles[indi,1]] = baryc1;
          evalmat[i,triangles[indi,2]] = baryc2;
          evalmat[i,triangles[indi,3]] = baryc3;
        }
        else if(nderivs[1] == 1 && nderivs[2] == 0)
        {
          evalmat[i,triangles[indi,1]] = M1[indi[1],2];
          evalmat[i,triangles[indi,2]] = M1[indi[2],2];
          evalmat[i,triangles[indi,3]] = M1[indi[3],2]; 
        }
        else if(nderivs[1] == 0 && nderivs[2] == 1)
        {
          evalmat[i,triangles[indi,1]] = M1[indi[1],3];
          evalmat[i,triangles[indi,2]] = M1[indi[2],3];
          evalmat[i,triangles[indi,3]] = M1[indi[3],3]; 
        }
        else
        {
          stop('Inadmissible derivative orders.')
        }
      }
    }
  }
  return(evalmat)
}


#' @param FEM A \code{FEM} object to be evaluated
#' @param locations A 2-columns matrix with the spatial locations where the FEM object should be evaluated
#' @return 
#' A matrix of numeric evaluations of the \code{FEM} object. Each row indicates the location where the evaluation has been taken, the column indicates the 
#' function evaluated.
#' @description Only executed when the function \code{smooth.FEM.basis} is run with the option \code{CPP_CODE} = \code{FALSE}. It evaluates a FEM object at the specified set of locations. 
#' @rdname fdaPDE-deprecated
#' @export

R_eval.FEM <- function(FEM, locations)
{ 
  .Deprecated("eval.FEM")
  
  if (is.vector(locations))
  {
    locations = t(as.matrix(locations))
  }
  
  N = nrow(locations)
  
  # Augment Xvec and Yvec by ones for computing barycentric coordinates
  Pgpts = cbind(matrix(1,N,1),locations[,1],locations[,2])
  
  # Get nodes and index
  
  FEMbasis = FEM$FEMbasis
  mesh = FEMbasis$mesh
  
  nodes = mesh$nodes
  triangles = mesh$triangles
  coeff = FEM$coeff
  nsurf = dim(coeff)[2]
  
  FEMbasis = FEM$FEMbasis
  order = FEMbasis$order
  #nodeindex = params$nodeindex
  detJ = FEMbasis$detJ
  
  # 1st, 2nd, 3rd vertices of triangles
  
  v1 = nodes[triangles[,1],]
  v2 = nodes[triangles[,2],]
  v3 = nodes[triangles[,3],]
  
  if(order !=2 && order != 1)
  {
    stop('ORDER is neither 1 or 2.')
  }
  
  # Denominator of change of coordinates chsange matrix
  
  modJ = FEMbasis$detJ
  ones3 = matrix(1,3,1)
  modJMat = modJ %*% t(ones3)
  
  M1 = cbind(v2[,1]*v3[,2] - v3[,1]*v2[,2], v2[,2] - v3[,2], v3[,1] - v2[,1])/(modJMat)
  M2 = cbind(v3[,1]*v1[,2] - v1[,1]*v3[,2], v3[,2] - v1[,2], v1[,1] - v3[,1])/(modJMat)
  M3 = cbind(v1[,1]*v2[,2] - v2[,1]*v1[,2], v1[,2] - v2[,2], v2[,1] - v1[,1])/(modJMat)
  
  ind = matrix(0,N,1)
  for(i in 1:N)
  {
    ind[i] = R_insideIndex(mesh, as.numeric(locations[i,]))
  }
  
  evalmat = matrix(NA, nrow=N, ncol=nsurf)
  
  for (isurf in 1:nsurf)
  {
    for(i in 1:N)
    {
      indi = ind[i]
      
      if(!is.nan(indi))
      {
        baryc1 = (M1[indi,]*Pgpts[i,]) %*% ones3
        baryc2 = (M2[indi,]*Pgpts[i,]) %*% ones3
        baryc3 = (M3[indi,]*Pgpts[i,]) %*% ones3
        
        if(order == 2)
        {
          c1 = coeff[triangles[indi,1],isurf]
          c2 = coeff[triangles[indi,2],isurf]
          c3 = coeff[triangles[indi,3],isurf]
          c4 = coeff[triangles[indi,6],isurf]
          c5 = coeff[triangles[indi,4],isurf]
          c6 = coeff[triangles[indi,5],isurf]
          
          fval =  c1*(2* baryc1^2 - baryc1) +
            c2*(2* baryc2^2 - baryc2) +
            c3*(2* baryc3^2 - baryc3) +
            c4*(4* baryc1 * baryc2) +
            c5*(4* baryc2 * baryc3) +
            c6*(4* baryc3 * baryc1)
          evalmat[i,isurf] = fval
        }
        else
        {
          c1 = coeff[triangles[indi,1],isurf]
          c2 = coeff[triangles[indi,2],isurf]
          c3 = coeff[triangles[indi,3],isurf]
          fval = c1*baryc1 + c2*baryc2 + c3*baryc3
          evalmat[i,isurf] = fval
        }
      }
    }
  }
  return(evalmat)
}

R_tricoefCal = function(mesh)
{
  #  TRICOEFCAL compute the coefficient matrix TRICOEF
  #  required to test of a point is indside a triangle
  
  nodes = mesh$nodes
  triangles = mesh$triangles
  
  ntri   = dim(triangles)[[1]]
  
  #  compute coefficients for computing barycentric coordinates if
  #  needed
  
  tricoef = matrix(0, nrow=ntri, ncol=4)
  tricoef[,1] = nodes[triangles[,1],1]-nodes[triangles[,3],1]
  tricoef[,2] = nodes[triangles[,2],1]-nodes[triangles[,3],1]
  tricoef[,3] = nodes[triangles[,1],2]-nodes[triangles[,3],2]
  tricoef[,4] = nodes[triangles[,2],2]-nodes[triangles[,3],2]
  detT = matrix((tricoef[,1]*tricoef[,4] - tricoef[,2]*tricoef[,3]),ncol=1)
  tricoef = tricoef/(detT %*% matrix(1,nrow=1,ncol=4))
  
  return(tricoef)
}


#' @param observations A vector of length #observations with the observed data values over the domain. 
#' The locations of the observations can be specified with the \code{locations} argument. 
#' Otherwise if only the vector of observations is given, these are consider to be located in the corresponding node in the table
#' \code{nodes} of the mesh. In this last case, an \code{NA} value in the \code{observations} vector indicates that there is no observation associated to the corresponding
#'  node.
#' @param locations A #observations-by-2 matrix where each row specifies the spatial coordinates \code{x} and \code{y} of the corresponding observations in the vector \code{observations}.
#' This parameter can be \code{NULL}. In this case the spatial coordinates of the corresponding observations are assigned as specified in \code{observations}.
#' @param FEMbasis A \code{FEMbasis} object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param lambda A scalar or vector of smoothing parameters.
#' @param covariates A #observations-by-#covariates matrix where each row represents the covariates associated with the corresponding observed data value in \code{observations}.
#' @param BC A list with two vectors: 
#'  \code{BC_indices}, a vector with the indices in \code{nodes} of boundary nodes where a Dirichlet Boundary Condition should be applied;
#'  \code{BC_values}, a vector with the values that the spatial field must take at the nodes indicated in \code{BC_indices}.
#' @param GCV Boolean. If \code{TRUE} the following quantities are computed: the trace of the smoothing matrix, the estimated error standard deviation,  and 
#'        the Generalized Cross Validation criterion, for each value of the smoothing parameter specified in \code{lambda}.
#' @param CPP_CODE Boolean. If \code{TRUE} the computation relies on the C++ implementation of the algorithm. This usually ensures a much faster computation.
#' @return A list with the following variables:
#' \item{\code{fit.FEM}}{A \code{FEM} object that represents the fitted spatial field.}
#' \item{\code{PDEmisfit.FEM}}{A \code{FEM} object that represents the Laplacian of the estimated spatial field.}
#' \item{\code{beta}}{If covariates is not \code{NULL}, a matrix with number of rows equal to the number of covariates and numer of columns equal to length of lambda.  The \code{j}th column represents the vector of regression coefficients when 
#' the smoothing parameter is equal to \code{lambda[j]}.}
#' \item{\code{edf}}{If GCV is \code{TRUE}, a scalar or vector with the trace of the smoothing matrix for each value of the smoothing parameter specified in \code{lambda}.}
#' \item{\code{stderr}}{If GCV is \code{TRUE}, a scalar or vector with the estimate of the standard deviation of the error for each value of the smoothing parameter specified in \code{lambda}.}
#' \item{\code{GCV}}{If GCV is \code{TRUE}, a  scalar or vector with the value of the GCV criterion for each value of the smoothing parameter specified in \code{lambda}.}
#' @description This function implements a spatial regression model with differential regularization; isotropic and stationary case. In particular, the regularizing term involves the Laplacian of the spatial field. Space-varying covariates can be included in the model. The technique accurately handle data distributed over irregularly shaped domains. Moreover, various conditions can be imposed at the domain boundaries.
#' @rdname fdaPDE-deprecated
#' @export


smooth.FEM.basis<-function(locations = NULL, observations, FEMbasis, lambda, covariates = NULL, BC = NULL, GCV = FALSE, CPP_CODE = TRUE)
{
  .Deprecated("smooth.FEM", package = "fdaPDE")
  ans=smooth.FEM(locations=locations,observations=observations,FEMbasis=FEMbasis,lambda=lambda, covariates=covariates,BC=BC,GCV=GCV,GCVmethod = 'Exact')
  ans
}


#' @param observations A vector of length #observations with the observed data values over the domain. 
#' The locations of the observations can be specified with the \code{locations} argument. 
#' Otherwise if only the vector of observations is given, these are consider to be located in the corresponding node in the table
#' \code{nodes} of the mesh. In this last case, an \code{NA} value in the \code{observations} vector indicates that there is no observation associated to the corresponding
#'  node.
#' \code{NA} values are admissible to indicate that the node is not associated with any observed data value.
#' @param locations A #observations-by-2 matrix where each row specifies the spatial coordinates \code{x} and \code{y} of the corresponding observations in the vector \code{observations}.
#' This parameter can be \code{NULL}. In this case the spatial coordinates of the corresponding observations are assigned as specified in \code{observations}.
#' @param FEMbasis A \code{FEMbasis} object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param lambda A scalar or vector of smoothing parameters.
#' @param PDE_parameters A list specifying the parameters of the elliptic PDE in the regularizing term: \code{K}, a 2-by-2 matrix of diffusion coefficients. This induces an anisotropic 
#' smoothing with a preferential direction that corresponds to the first eigenvector of the diffusion matrix K; \code{b}, a vector of length 2 of advection coefficients. This induces a 
#' smoothing only in the direction specified by the vector \code{b}; \code{c}, a scalar reaction coefficient. \code{c} induces a shrinkage of the surface to zero
#' @param covariates A #observations-by-#covariates matrix where each row represents the covariates associated with the corresponding observed data value in \code{observations}.
#' @param BC A list with two vectors: 
#'  \code{BC_indices}, a vector with the indices in \code{nodes} of boundary nodes where a Dirichlet Boundary Condition should be applied;
#'  \code{BC_values}, a vector with the values that the spatial field must take at the nodes indicated in \code{BC_indices}.
#' @param GCV Boolean. If \code{TRUE} the following quantities are computed: the trace of the smoothing matrix, the estimated error standard deviation,  and 
#'        the Generalized Cross Validation criterion, for each value of the smoothing parameter specified in \code{lambda}.
#' @param CPP_CODE Boolean. If \code{TRUE} the computation relies on the C++ implementation of the algorithm. This usually ensures a much faster computation.
#' @return A list with the following variables:
#'          \item{\code{fit.FEM}}{A \code{FEM} object that represents the fitted spatial field.}
#'          \item{\code{PDEmisfit.FEM}}{A \code{FEM} object that represents the PDE misfit for the estimated spatial field.}
#'          \item{\code{beta}}{If covariates is not \code{NULL}, a matrix with number of rows equal to the number of covariates and numer of columns equal to length of lambda.  The \code{j}th column represents the vector of regression coefficients when 
#'          the smoothing parameter is equal to \code{lambda[j]}.}
#'          \item{\code{edf}}{If GCV is \code{TRUE}, a scalar or vector with the trace of the smoothing matrix for each value of the smoothing parameter specified in \code{lambda}.}
#'          \item{\code{stderr}}{If GCV is \code{TRUE}, a scalar or vector with the estimate of the standard deviation of the error for each value of the smoothing parameter specified in \code{lambda}.}
#'          \item{\code{GCV}}{If GCV is \code{TRUE}, a  scalar or vector with the value of the GCV criterion for each value of the smoothing parameter specified in \code{lambda}.}
#' @description This function implements a spatial regression model with differential regularization; anysotropic case. In particular, the regularizing term involves a second order elliptic PDE, that models the space-variation of the phenomenon. Space-varying covariates can be included in the model. The technique accurately handle data distributed over irregularly shaped domains. Moreover, various conditions can be imposed at the domain boundaries.
#' @rdname fdaPDE-deprecated
#' @export
smooth.FEM.PDE.basis<-function(locations = NULL, observations, FEMbasis, lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV = FALSE, CPP_CODE = TRUE)
{
  .Deprecated("smooth.FEM", package = "fdaPDE")
  ans=smooth.FEM(locations=locations,observations=observations,FEMbasis=FEMbasis,lambda=lambda,PDE_parameters=PDE_parameters,covariates=covariates,BC=BC,GCV=GCV,GCVmethod = 'Exact')
  ans
}



#' @param observations A vector of length #observations with the observed data values over the domain. 
#' The locations of the observations can be specified with the \code{locations} argument. 
#' Otherwise if only the vector of observations is given, these are consider to be located in the corresponding node in the table
#' \code{nodes} of the mesh. In this last case, an \code{NA} value in the \code{observations} vector indicates that there is no observation associated to the corresponding
#'  node.
#' @param locations A #observations-by-2 matrix where each row specifies the spatial coordinates \code{x} and \code{y} of the corresponding observations in the vector \code{observations}.
#' This parameter can be \code{NULL}. In this case the spatial coordinates of the corresponding observations are assigned as specified in \code{observations}.
#' @param FEMbasis A \code{FEMbasis} object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param lambda A scalar or vector of smoothing parameters.
#' @param PDE_parameters A list specifying the space-varying parameters of the elliptic PDE in the regularizing term: \code{K}, a function that for each spatial location in the spatial domain 
#' (indicated by the vector of the 2 spatial coordinates) returns a 2-by-2 matrix of diffusion coefficients. This induces an anisotropic 
#' smoothing with a local preferential direction that corresponds to the first eigenvector of the diffusion matrix K.The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
#' an array with dimensions 2-by-2-by-#points.\code{b}, a function that for each spatial location in the spatial domain returns 
#' a vector of length 2 of transport coefficients. This induces a local smoothing only in the direction specified by the vector \code{b}. The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
#' a matrix with dimensions 2-by-#points; \code{c}, a function that for each spatial location in the spatial domain  returns a scalar reaction coefficient.
#' \code{c} induces a shrinkage of the surface to zero. The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
#' a vector with length #points; \code{u}, a function that for each spatial location in the spatial domain  returns a scalar reaction coefficient.
#' \code{u} induces a reaction effect. The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
#' a vector with length #points.
#' @param covariates A #observations-by-#covariates matrix where each row represents the covariates associated with the corresponding observed data value in \code{observations}.
#' @param BC A list with two vectors: 
#'  \code{BC_indices}, a vector with the indices in \code{nodes} of boundary nodes where a Dirichlet Boundary Condition should be applied;
#'  \code{BC_values}, a vector with the values that the spatial field must take at the nodes indicated in \code{BC_indices}.
#' @param GCV Boolean. If \code{TRUE} the following quantities are computed: the trace of the smoothing matrix, the estimated error standard deviation,  and 
#'        the Generalized Cross Validation criterion, for each value of the smoothing parameter specified in \code{lambda}.
#' @param CPP_CODE Boolean. If \code{TRUE} the computation relies on the C++ implementation of the algorithm. This usually ensures a much faster computation.
#' @return A list with the following variables:
#'          \item{\code{fit.FEM}}{A \code{FEM} object that represents the fitted spatial field.}
#'          \item{\code{PDEmisfit.FEM}}{A \code{FEM} object that represents the PDE misfit for the estimated spatial field.}
#'          \item{\code{beta}}{If covariates is not \code{NULL}, a matrix with number of rows equal to the number of covariates and numer of columns equal to length of lambda.  The \code{j}th column represents the vector of regression coefficients when 
#'          the smoothing parameter is equal to \code{lambda[j]}.}
#'          \item{\code{edf}}{If GCV is \code{TRUE}, a scalar or vector with the trace of the smoothing matrix for each value of the smoothing parameter specified in \code{lambda}.}
#'          \item{\code{stderr}}{If GCV is \code{TRUE}, a scalar or vector with the estimate of the standard deviation of the error for each value of the smoothing parameter specified in \code{lambda}.}
#'          \item{\code{GCV}}{If GCV is \code{TRUE}, a  scalar or vector with the value of the GCV criterion for each value of the smoothing parameter specified in \code{lambda}.}
#' @description This function implements a spatial regression model with differential regularization; anysotropic and non-stationary case. In particular, the regularizing term involves a second order elliptic PDE with space-varying coefficients, that models the space-variation of the phenomenon. Space-varying covariates can be included in the model. The technique accurately handle data distributed over irregularly shaped domains. Moreover, various conditions can be imposed at the domain boundaries.
#' @rdname fdaPDE-deprecated
#' @export

smooth.FEM.PDE.sv.basis<-function(locations = NULL, observations, FEMbasis, lambda, PDE_parameters, covariates = NULL, BC = NULL, GCV = FALSE, CPP_CODE = TRUE)
{
  .Deprecated("smooth.FEM", package = "fdaPDE")
  ans=smooth.FEM(locations=locations,observations=observations,FEMbasis=FEMbasis,lambda=lambda,PDE_parameters=PDE_parameters,covariates=covariates,BC=BC,GCV=GCV,GCVmethod = 'Exact')
  ans
}


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
#' In this case the function \code{create.MESH.2D} is used to produce a complete MESH2D object. 
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
#' @usage create.MESH.2D(nodes, nodesattributes = NA, segments = NA, holes = NA, 
#'                      triangles = NA, order = 1, verbosity = 0)
#' @seealso \code{\link{refine.MESH.2D}}, \code{\link{create.FEM.basis}}
#' @return An object of the class MESH2D with the following output:
#' \item{\code{nodes}}{A #nodes-by-2 matrix containing the x and y coordinates of the mesh nodes.}
#' \item{\code{nodesmarkers}}{A vector of length #nodes, with entries either '1' or '0'. An entry '1' indicates that the corresponding node is a boundary node; an entry '0' indicates that the corresponding node is not a boundary node.}
#' \item{\code{nodesattributes}}{nodesattributes A matrix with #nodes rows containing nodes' attributes. 
#' These are passed unchanged to the output. If a node is added during the triangulation process or mesh refinement, its attributes are computed  
#' by linear interpolation using the attributes of neighboring nodes. This functionality is for instance used to compute the value 
#' of a Dirichlet boundary condition at boundary nodes added during the triangulation process.}
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
#' These are respectively used for linear (order = 1) and quadratic (order = 2) Finite Elements. Default is \code{order} = 1.}
#' @rdname fdaPDE-deprecated
#' @export

create.MESH.2D <- function(nodes, nodesattributes = NA, segments = NA, holes = NA, triangles = NA, order = 1, verbosity = 0)
{ 
  .Deprecated("create.mesh.2D")
  ans=create.mesh.2D(nodes=nodes,nodesattributes=nodesattributes,segments=segments,holes=holes,triangles=triangles,order=order,verbosity=verbosity)
  ans
}


#' @param mesh A MESH2D object representing the triangular mesh, created by \link{create.MESH.2D}.
#' @param minimum_angle A scalar specifying a minimun value for the triangles angles.
#' @param maximum_area A scalar specifying a maximum value for the triangles areas.
#' @param delaunay A boolean parameter indicating whether or not the output mesh should satisfy the Delaunay condition.
#' @param verbosity This can be '0', '1' or '2'. It indicates the level of verbosity in the triangulation process.
#' @description This function refines a Constrained Delaunay triangulation into a Conforming Delaunay triangulation. This is a wrapper of the Triangle library (http://www.cs.cmu.edu/~quake/triangle.html). It can be used to 
#' refine a mesh created previously with \link{create.MESH.2D}. The algorithm can add Steiner points (points through which the \code{segments} are splitted)
#' in order to meet the imposed refinement conditions.
#' @usage refine.MESH.2D(mesh, minimum_angle, maximum_area, delaunay, verbosity)
#' @seealso \code{\link{create.MESH.2D}}, \code{\link{create.FEM.basis}}
#' @return A MESH2D object representing the refined triangular mesh,  with the following output:
#' \item{\code{nodes}}{A #nodes-by-2 matrix containing the x and y coordinates of the mesh nodes.}
#' \item{\code{nodesmarkers}}{A vector of length #nodes, with entries either '1' or '0'. An entry '1' indicates that the corresponding node is a boundary node; an entry '0' indicates that the corresponding node is not a boundary node.}
#' \item{\code{nodesattributes}}{nodesattributes A matrix with #nodes rows containing nodes' attributes. 
#' These are passed unchanged to the output. If a node is added during the triangulation process or mesh refinement, its attributes are computed  
#' by linear interpolation using the attributes of neighboring nodes. This functionality is for instance used to compute the value 
#' of a Dirichlet boundary condition at boundary nodes added during the triangulation process.}
#' \item{\code{triangles}}{A #triangles-by-3 (when \code{order} = 1) or #triangles-by-6 (when \code{order} = 2) matrix.
#' This option is used when a triangulation is already available. It specifies the triangles giving the row's indices in \code{nodes} of the triangles' vertices and (when \code{nodes} = 2) also if the triangles' edges midpoints. The triangles' vertices and midpoints are ordered as described 
#' at \cr  https://www.cs.cmu.edu/~quake/triangle.highorder.html.}
#' \item{\code{edges}}{A #edges-by-2 matrix. Each row contains the row's indices of the nodes where the edge starts from and ends to.}
#' \item{\code{edgesmarkers}}{A vector of lenght #edges with entries either '1' or '0'. An entry '1' indicates that the corresponding element in \code{edge} is a boundary edge;  
#' an entry '0' indicates that the corresponding edge is not a boundary edge.}
#' \item{\code{neighbors}}{A #triangles-by-3 matrix. Each row contains the indices of the three neighbouring triangles. An entry '-1' indicates that 
#' one edge of the triangle is a boundary edge.}
#' \item{\code{holes}}{A #holes-by-2 matrix containing the x and y coordinates of a point internal to each hole of the mesh. These points are used to carve holes
#' in the triangulation, when the domain has holes.}
#' \item{\code{order}}{Either '1' or '2'. It specifies wether each mesh triangle should be represented by 3 nodes (the triangle' vertices) or by 6 nodes (the triangle's vertices and midpoints). 
#' These are respectively used for linear (order = 1) and quadratic (order = 2) Finite Elements. Default is \code{order} = 1.}
#' @export
#' @rdname fdaPDE-deprecated

refine.MESH.2D<-function(mesh, minimum_angle = NA, maximum_area = NA, delaunay = FALSE, verbosity = 0)
{ 
  .Deprecated("refine.mesh.2D")
  ans=refine.mesh.2D(mesh=mesh, minimum_angle=minimum_angle,maximum_area=maximum_area,delaunay=delaunay,verbosity=verbosity)
  ans
}


#' @param x A MESH2D object defining the triangular mesh, as generated by \code{create.Mesh.2D} or \code{refine.Mesh.2D}.
#' @param ... Arguments representing graphical options to be passed to \link[graphics]{par}.
#' @description Plot a mesh MESH2D object, generated by \code{create.MESH.2D} or \code{refine.MESH.2D}. Circles indicate the mesh nodes.
#' @usage \method{plot}{MESH2D}(x, ...)
#' @export
#' @rdname fdaPDE-deprecated

plot.MESH2D<-function(x, ...)
{
  .Deprecated("plot.mesh.2D")
  plot(x$nodes, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", ...)
  segments(x$nodes[x$edges[,1],1], x$nodes[x$edges[,1],2],
           x$nodes[x$edges[,2],1], x$nodes[x$edges[,2],2], ...)
  segments(x$nodes[x$segments[,1],1], x$nodes[x$segments[,1],2],
           x$nodes[x$segments[,2],1], x$nodes[x$segments[,2],2], col="red", ...)
}
