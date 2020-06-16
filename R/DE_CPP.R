CPP_FEM.DE <- function(data, FEMbasis, lambda, fvec, heatStep, heatIter, ndim, mydim, step_method, direction_method, preprocess_method,
                       stepProposals, tol1, tol2, print, nfolds, nsimulations, search)
{
  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

  if(is.null(stepProposals))
    stepProposals = c(0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 1e-7, 1e-8, 1e-9)

  if(is.null(preprocess_method))
    preprocess_method = ""

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
  step_method <- as.character(step_method)
  storage.mode(step_method) <- "character"
  direction_method <- as.character(direction_method)
  storage.mode(direction_method) <- "character"
  preprocess_method <- as.character(preprocess_method)
  storage.mode(preprocess_method) <- "character"
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

  ## Call C++ function
  bigsol <- .Call("Density_Estimation", data, FEMbasis$mesh, FEMbasis$order, mydim, ndim, fvec,
                  heatStep, heatIter, lambda, nfolds, nsimulations, stepProposals, tol1, tol2,
                  print, step_method, direction_method, preprocess_method, search,
                  PACKAGE = "fdaPDE")

  ## Reset them correctly
  # FEMbasis$mesh$triangles = FEMbasis$mesh$triangles + 1
  # FEMbasis$mesh$edges = FEMbasis$mesh$edges + 1
  # FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] + 1
  return(bigsol)
}


CPP_FEM.manifold.DE <- function(data, FEMbasis, lambda, fvec, heatStep, heatIter, ndim, mydim, step_method, direction_method, preprocess_method,
                                stepProposals, tol1, tol2, print, nfolds, nsimulations, search)
{
  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1


  if(is.null(stepProposals))
    stepProposals = c(0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 1e-7, 1e-8, 1e-9)

  if(is.null(preprocess_method))
    preprocess_method = ""

  ## Set proper type for correct C++ reading
  data <- as.matrix(data)
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
  step_method <- as.character(step_method)
  storage.mode(step_method) <- "character"
  direction_method <- as.character(direction_method)
  storage.mode(direction_method) <- "character"
  preprocess_method <- as.character(preprocess_method)
  storage.mode(preprocess_method) <- "character"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(stepProposals) <- "double"
  storage.mode(tol1) <- "double"
  storage.mode(tol2) <- "double"
  storage.mode(print) <- "logical"
  storage.mode(nfolds) <- "integer"
  storage.mode(nsimulations) <- "integer"
  storage.mode(search) <- "integer"

  bigsol <- .Call("Density_Estimation", data, FEMbasis$mesh, FEMbasis$order, mydim, ndim, fvec,
                  heatStep, heatIter, lambda, nfolds, nsimulations, stepProposals, tol1, tol2,
                  print, step_method, direction_method, preprocess_method, search,
                  PACKAGE = "fdaPDE")

  return(bigsol)
}


CPP_FEM.volume.DE <- function(data, FEMbasis, lambda, fvec, heatStep, heatIter, ndim, mydim, step_method, direction_method, preprocess_method,
                              stepProposals, tol1, tol2, print, nfolds, nsimulations, search)
{

  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$tetrahedrons = FEMbasis$mesh$tetrahedrons - 1
  FEMbasis$mesh$faces = FEMbasis$mesh$faces - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

  if(is.null(stepProposals))
    stepProposals = c(0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 1e-7, 1e-8, 1e-9)

  if(is.null(preprocess_method))
    preprocess_method = ""

  ## Set proper type for correct C++ reading
  storage.mode(data) <- "double"
  storage.mode(FEMbasis$mesh$order) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$faces) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$mesh$tetrahedrons) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(fvec) <- "double"
  storage.mode(heatStep) <- "double"
  heatIter <- as.integer(heatIter)
  storage.mode(heatIter) <- "integer"
  step_method <- as.character(step_method)
  storage.mode(step_method) <- "character"
  direction_method <- as.character(direction_method)
  storage.mode(direction_method) <- "character"
  preprocess_method <- as.character(preprocess_method)
  storage.mode(preprocess_method) <- "character"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(stepProposals) <- "double"
  storage.mode(tol1) <- "double"
  storage.mode(tol2) <- "double"
  storage.mode(print) <- "logical"
  storage.mode(nfolds) <- "integer"
  storage.mode(nsimulations) <- "integer"
  storage.mode(search) <- "integer"

  bigsol <- .Call("Density_Estimation", data, FEMbasis$mesh, FEMbasis$order, mydim, ndim, fvec,
                  heatStep, heatIter, lambda, nfolds, nsimulations, stepProposals, tol1, tol2,
                  print, step_method, direction_method, preprocess_method, search,
                  PACKAGE = "fdaPDE")

  return(bigsol)
}
