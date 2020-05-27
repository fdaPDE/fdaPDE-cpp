#' Nonparametric density estimation with differential regularization
#' 
#' @param data A matrix of dimensions #observations-by-ndim. Data are locations: each row corresponds to one point, 
#' the first column corresponds to the \code{x}-coordinates, the second column corresponds to the \code{y}-coordinates 
#' and, if ndim=3, the third column corresponds to the \code{z}-coordinates. 
#' @param FEMbasis A \code{FEMbasis} object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param lambda A scalar or vector of smoothing parameters. If it is a vector, the optimal smoothing parameter is choosen
#' with a \code{k}-fold cross-validation procedure based on the L2 norm.
#' @param fvec A vector of length #\code{nodes} of the mesh. It corresponds to the node values of the initial density function. 
#' If this is \code{NULL} the initial density is estimated thanks to a discretized heat diffusion 
#' process that starts from the empirical density of the data. Default is \code{NULL}.
#' N.B. This vector cannot be the constant vector of zeros since the algortihm works with the log(f).
#' @param heatStep Real specifying the time step for the discretized heat diffusionn process.
#' @param heatIter Integer specifying the number of iteriations to perform the discretized heat diffusion process.
#' @param stepProposals A scalar or a vector containing the step parameters useful for the descent algotihm. If there is a 
#' vector of parameters, the biggest one such that the functional decreases at each iteration is choosen. If it is \code{NULL}
#' the following vector \code{c(0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 1e-7, 1e-8, 1e-9)} is proposed. Default is \code{NULL}.
#' N.B. If the program does not receive a right parameter, it abort the R session. Try a smaller parameter.
#' @param tol1 A scalar specifying the tolerance to use for the termination criterion based on the percentage difference
#' between two consecutive iterations of the minimization algorithm of the loss function, the log-likelihood and the
#' penalization. Default is 1e-5.
#' @param tol2 A scalar specifying the tolerance to use for the termination criterion based on the norm of the gradient 
#' of the functional to be minimized (the true minimum is such that this norm is zero). The default does not use this 
#' criterion. Default is 0.
#' @param print A boolean that is \code{TRUE} if the user wants the value of the functional, of the loglikelihood and of the
#' penalization term printed on console at each iteration of the descent algorithm. Default is \code{FALSE}.
#' N.B. We suggest to let it \code{FALSE} if \code{preprocess_method} is 'RightCV' or 'SimplifiedCV'.
#' @param nfolds An integer specifying the number of folds used in cross validation techinque to find the best \code{lambda} parameter.
#' If there is only one \code{lambda} it can be \code{NULL}. Default is \code{NULL}.
#' @param nsimulations An integer specifying the number of iterations used in the optimization algorithms. Default value is 500.
#' @param step_method String. This parameter specifies which step method use in the descent algorithm. 
#' If it is \code{Fixed_Step}, the step is constant during all the algorithm and it is choosen according to \code{stepProposals};
#' if it is \code{Backtracking_Method}, the step is computed at each iteration according to the backtracking method; finally
#' if it is \code{Wolfe_Method}, the step is computed at each iteration according to the Wolfe method. Default is \code{Fixed_Step}.
#' @param direction_method String. This parameter specifies which descent direction use in the descent algorithm. 
#' If it is \code{Gradient}, the direction is the one given by the gradient descent method (the opposite to the gradient of
#' the functional); if instead it is \code{BFGS} the direction is the one given by the BFGS method
#' (Broyden Fletcher Goldfarb and Shanno, a Quasi-Newton method). Default is \code{BFGS}.
#' @param preprocess_method String. This parameter specifies the k fold cross validation technique to use, if there is more
#' than one smoothing parameter \code{lambda} (otherwise it should be \code{NULL}). If it is \code{RightCV} the usual k fold 
#' cross validation method is performed. If it is \code{SimplifiedCV} a simplified version is performed. 
#' In the latter case the number of smoothing parameters \code{lambda} must be equal to the number of folds \code{nfolds}.
#' Default is \code{NULL}.
#' @param search An integer specifying the search algorithm to use. It is either 1 (Naive search algorithm) or 2 (Tree search algorithm).
#' The default is 2.
#' @return A list with the following variables:
#' \item{\code{FEMbasis}}{Given FEMbasis with tree informations.}
#' \item{\code{g}}{A vector of length #\code{nodes} that represents the value of the g-function estimated for each \code{node} of the mesh.
#' The density is the exponential of this function.}
#' \item{\code{f_init}}{A #\code{nodes}-by-#\code{lambda} parameters matrix. Each column contains the node values of the initial
#' density used for the lambda given by the column.}
#' \item{\code{lambda}}{A scalar representing the optimal smoothing parameter selected via k fold cross validation, if in the 
#' input there is a vector of parameters; the scalar given in input otherwise.}
#' \item{\code{data}}{A matrix of dimensions #observations-by-ndim containing the data used in the algorithm. They are the 
#' same given in input if the domain is 2D pr 3D; they are the original data projected on the mesh if the domain is 2.5D.}
#' \item{\code{CV_err}}{A vector of length \code{nfolds} containing the cross validation errors obtained in each fold, if 
#' \code{preprocess_method} is either \code{RightCV} or \code{SimplifiedCV}.}
#' @description This function implements a nonparametric density estimation method with differential regularization 
#' (given by the square root of the L2 norm of the laplacian of the density function), when points are located over a 
#' planar mesh. The computation relies only on the C++ implementation of the algorithm.
#' @usage DE.FEM(data, FEMbasis, lambda, fvec=NULL, heatStep=0.1, heatIter=500, 
#'               stepProposals=NULL,tol1=1e-4, tol2=0, print=FALSE, nfolds=NULL, 
#'               nsimulations=500, step_method="Fixed_Step", direction_method="BFGS", 
#'               preprocess_method="NoCrossValidation", search = 2)
#' @export
#' @examples
#' library(fdaPDE)
#' 
#' ## Create a 2D mesh over a squared domain
#' Xbound <- seq(-3, 3, length.out = 10)
#' Ybound <- seq(-3, 3, length.out = 10)
#' grid_XY <- expand.grid(Xbound, Ybound)
#' Bounds <- grid_XY[(grid_XY$Var1 %in% c(-3, 3)) | (grid_XY$Var2 %in% c(-3, 3)), ]
#' mesh <- create.mesh.2D(nodes = Bounds, order = 1)
#' mesh <- refine.mesh.2D(mesh, maximum_area = 0.2)
#' FEMbasis <- create.FEM.basis(mesh)
#' 
#' ## Generate data
#' n <- 50
#' set.seed(10)
#' data_x <- rnorm(n)
#' data_y <- rnorm(n)
#' data <- cbind(data_x, data_y)
#' 
#' plot(mesh)
#' points(data, col="red", pch=19, cex=0.5)
#' 
#' ## Density Estimation
#' lambda = 0.1
#' sol <- DE.FEM(data = data, FEMbasis = FEMbasis, lambda = lambda, fvec=NULL, heatStep=0.1,
#'                   heatIter=500, stepProposals=NULL, tol1=1e-4, tol2=0, print=FALSE, 
#'                   nfolds=NULL, nsimulations=500,step_method = "Fixed_Step", 
#'                   direction_method = "BFGS",preprocess_method="NoCrossValidation", 
#'                   search = 2)
#' 
#' ## Visualization 
#' n = 100
#' X <- seq(-3, 3, length.out = n)
#' Y<- seq(-3, 3, length.out = n)
#' grid <- expand.grid(X, Y)
#'
#' evaluation <- eval.FEM(FEM(FEMbasis, coeff = sol$g), locations = grid)
#' evaluation <- exp(evaluation)
#' eval <- matrix(evaluation, n, n)
#' 
#' image2D(x = X, y = Y, z = eval, col = heat.colors(100), xlab = "x", ylab = "y", 
#'         contour = list(drawlabels = FALSE), main = "Estimated density")


DE.FEM <- function(data, FEMbasis, lambda, fvec=NULL, heatStep=0.1, heatIter=500, stepProposals=NULL,
                  tol1=1e-4, tol2=0, print=FALSE, nfolds=NULL, nsimulations=500, step_method="Fixed_Step",
                  direction_method="BFGS", preprocess_method="NoCrossValidation", search = 2) 
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
  
  ###################### Checking parameters, sizes and conversion #################################
  checkParametersDE(data, FEMbasis, lambda, step_method, direction_method, preprocess_method, tol1, tol2, nfolds, nsimulations, heatStep, heatIter, search) 
  
  ## Coverting to format for internal usage
  data = as.matrix(data)
  lambda = as.vector(lambda)
  if(!is.null(fvec))
    fvec = as.vector(fvec)
  if(!is.null(stepProposals))
    stepProposals = as.vector(stepProposals)
  
  checkParametersSizeDE(data, FEMbasis, ndim, fvec, preprocess_method, nfolds) 
  ###################### End checking parameters, sizes and conversion #############################
  
  
  ###################### C++ Code Execution #########################################################
  bigsol = NULL
  if(class(FEMbasis$mesh) == 'mesh.2D'){	  
    print('C++ Code Execution')
    bigsol = CPP_FEM.DE(data, FEMbasis, lambda, fvec, heatStep, heatIter, ndim, mydim, step_method, direction_method, preprocess_method,
                        stepProposals, tol1, tol2, print, nfolds, nsimulations, search)
    
  } else if(class(FEMbasis$mesh) == 'mesh.2.5D'){
    print('C++ Code Execution')
    bigsol = CPP_FEM.manifold.DE(data, FEMbasis, lambda, fvec, heatStep, heatIter, ndim, mydim, step_method, direction_method, preprocess_method,
                                 stepProposals, tol1, tol2, print, nfolds, nsimulations, search)
    
  } else if(class(FEMbasis$mesh) == 'mesh.3D'){
    bigsol = CPP_FEM.volume.DE(data, FEMbasis, lambda, fvec, heatStep, heatIter, ndim, mydim, step_method, direction_method, preprocess_method,
                               stepProposals, tol1, tol2, print, nfolds, nsimulations, search)
  }
  
  ###################### Collect Results ############################################################  
  
  g = bigsol[[1]]
  f_init = bigsol[[2]]
  lambda = bigsol[[3]]
  data = bigsol[[4]]
  CV_err = bigsol[[5]]
  
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
  
  reslist = list(FEMbasis = FEMbasis, g = g, f_init = f_init, lambda = lambda, data = data, CV_err = CV_err)
  return(reslist)
}
