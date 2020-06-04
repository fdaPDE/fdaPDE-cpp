CPP_FEM.DE_init <- function(data, FEMbasis, lambda, fvec, heatStep, heatIter, ndim, mydim, step_method, direction_method, preprocess_method,
                       stepProposals, tol1, tol2, print, nfolds, nsimulations, search, init, nFolds)
{
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  
  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

  ## Set proper type for correct C++ reading
  storage.mode(data) <- "double"                        
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(fvec) <- "double"
  storage.mode(heatStep) <- "double"
  heatIter <- as.integer(heatIter)
  storage.mode(heatIter) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(stepProposals) <- "double"
  storage.mode(tol1) <- "double"
  storage.mode(tol2) <- "double"
  storage.mode(print) <- "logical"
  nfolds <- as.integer(nfolds)
  storage.mode(nfolds) <- "integer"
  nsimulations <- as.integer(nsimulations)
  storage.mode(nsimulations) <- "integer"
  storage.mode(search) <- "integer"
  init <- as.character(init)
  storage.mode(init) <- "character"
  nFolds <- as.integer(nFolds)
  storage.mode(nFolds) <- "integer"
  
  ## Call C++ function
  bigsol <- .Call("Density_Initialization", data, FEMbasis$mesh, FEMbasis$order, mydim, ndim, fvec, heatStep, heatIter, lambda,
                  nfolds, nsimulations, stepProposals, tol1, tol2, print, 
                  search, init, nFolds, PACKAGE = "fdaPDE")
  
  return(bigsol)
}


CPP_FEM.manifold.DE_init <- function(data, FEMbasis, lambda, fvec, heatStep, heatIter, ndim, mydim, step_method, direction_method, preprocess_method,
                            stepProposals, tol1, tol2, print, nfolds, nsimulations, search, init, nFolds)
{
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  
  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  
  ## Set proper type for correct C++ reading
  storage.mode(data) <- "double"                        
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(fvec) <- "double"
  storage.mode(heatStep) <- "double"
  heatIter <- as.integer(heatIter)
  storage.mode(heatIter) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(stepProposals) <- "double"
  storage.mode(tol1) <- "double"
  storage.mode(tol2) <- "double"
  storage.mode(print) <- "logical"
  nfolds <- as.integer(nfolds)
  storage.mode(nfolds) <- "integer"
  nsimulations <- as.integer(nsimulations)
  storage.mode(nsimulations) <- "integer"
  storage.mode(search) <- "integer"
  init <- as.character(init)
  storage.mode(init) <- "character"
  nFolds <- as.integer(nFolds)
  storage.mode(nFolds) <- "integer"
  
  ## Call C++ function
  bigsol <- .Call("Density_Initialization", data, FEMbasis$mesh, FEMbasis$order, mydim, ndim, fvec, heatStep, heatIter, lambda,
                  nfolds, nsimulations, stepProposals, tol1, tol2, print,
                  search, init, nFolds, PACKAGE = "fdaPDE")
  
  return(bigsol)
}


CPP_FEM.volume.DE_init <- function(data, FEMbasis, lambda, fvec, heatStep, heatIter, ndim, mydim, step_method, direction_method, preprocess_method, 
                              stepProposals, tol1, tol2, print, nfolds, nsimulations, search, init, nFolds)
{
  
  FEMbasis$mesh$tetrahedrons=c(t(FEMbasis$mesh$tetrahedrons))
  FEMbasis$mesh$nodes=c(t(FEMbasis$mesh$nodes))
  
  #riporto in R lo shift degli indici
  FEMbasis$mesh$tetrahedrons=FEMbasis$mesh$tetrahedrons-1
  
  
  ## Set proper type for correct C++ reading
  storage.mode(data) <- "double"                        
  storage.mode(FEMbasis$mesh$order) <- "integer"
  storage.mode(FEMbasis$mesh$nnodes) <- "integer"
  storage.mode(FEMbasis$mesh$ntetrahedrons) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$tetrahedrons) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(fvec) <- "double"
  storage.mode(heatStep) <- "double"
  heatIter <- as.integer(heatIter)
  storage.mode(heatIter) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(stepProposals) <- "double"
  storage.mode(tol1) <- "double"
  storage.mode(tol2) <- "double"
  storage.mode(print) <- "logical"
  storage.mode(nfolds) <- "integer"
  storage.mode(nsimulations) <- "integer"
  storage.mode(search) <- "integer"
  init <- as.character(init)
  storage.mode(init) <- "character"
  nFolds <- as.integer(nFolds)
  storage.mode(nFolds) <- "integer"
  
  
  bigsol <- .Call("Density_Initialization", data, FEMbasis$mesh, FEMbasis$order, mydim, ndim, fvec, heatStep, heatIter,
                  lambda, nfolds, nsimulations, stepProposals, tol1, tol2, print,
                  search, init, nFolds,
                  PACKAGE = "fdaPDE")
  
  return(bigsol)
}