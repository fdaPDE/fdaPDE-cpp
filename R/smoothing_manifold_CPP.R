CPP_smooth.manifold.FEM.basis<-function(locations, bary.locations, observations, FEMbasis, lambda, covariates = NULL, incidence_matrix = NULL, ndim, mydim, BC = NULL, GCV, GCVMETHOD = 2, nrealizations = 100, DOF=TRUE, DOF_matrix=NULL, search, GCV.inflation.factor = GCV.inflation.factor, areal.data.avg = TRUE)
{

  # C++ function for manifold works with vectors not with matrices

  FEMbasis$mesh$triangles=c(t(FEMbasis$mesh$triangles))
  FEMbasis$mesh$nodes=c(t(FEMbasis$mesh$nodes))

  # Indexes in C++ starts from 0, in R from 1, opporGCV.inflation.factor transformation

  FEMbasis$mesh$triangles=FEMbasis$mesh$triangles-1


  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(DOF_matrix))
  {
    DOF_matrix<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = ndim)
  }

  if(is.null(incidence_matrix))
  {
    incidence_matrix<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(BC$BC_indices))
  {
    BC$BC_indices<-vector(length=0)
  }else
  {
    BC$BC_indices<-as.vector(BC$BC_indices)-1

  }

  if(is.null(BC$BC_values))
  {
    BC$BC_values<-vector(length=0)
  }else
  {
    BC$BC_values<-as.vector(BC$BC_values)
  }

  ## Set propr type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  data <- as.vector(observations)
  storage.mode(observations) <- "double"
  storage.mode(FEMbasis$mesh$order) <- "integer"
  storage.mode(FEMbasis$mesh$nnodes) <- "integer"
  storage.mode(FEMbasis$mesh$ntriangles) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values)  <- "double"
  GCV <- as.integer(GCV)
  storage.mode(GCV) <- "integer"
  DOF <- as.integer(DOF)
  storage.mode(DOF) <- "integer"

  storage.mode(nrealizations) <- "integer"
  storage.mode(GCVMETHOD) <- "integer"
  storage.mode(search) <- "integer"

  storage.mode(GCV.inflation.factor) <- "double"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <-"integer"

  ## Call C++ function
  bigsol <- .Call("regression_Laplace", locations, bary.locations, data, FEMbasis$mesh, FEMbasis$mesh$order, mydim, ndim, lambda, covariates,
                  incidence_matrix, BC$BC_indices, BC$BC_values, GCV, GCVMETHOD, nrealizations, DOF, DOF_matrix, search, GCV.inflation.factor, areal.data.avg, PACKAGE = "fdaPDE")

  return(bigsol)
}

CPP_eval.manifold.FEM = function(FEM, locations, incidence_matrix, redundancy, ndim, mydim, search, bary.locations)
{
  FEMbasis = FEM$FEMbasis

  # C++ function for manifold works with vectors not with matrices

  FEMbasis$mesh$triangles=c(t(FEMbasis$mesh$triangles))
  FEMbasis$mesh$nodes=c(t(FEMbasis$mesh$nodes))

  #NEW: riporto shift indici in R
  FEMbasis$mesh$triangles=FEMbasis$mesh$triangles-1

  # Imposing types, this is necessary for correct reading from C++
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  coeff <- as.matrix(FEM$coeff)
  storage.mode(coeff) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(locations) <- "double"
  storage.mode(redundancy) <- "integer"
  storage.mode(search) <- "integer"

  if(!is.null(bary.locations))
  {
    storage.mode(bary.locations$element_ids) <- "integer"
    element_ids <- as.matrix(bary.locations$element_ids)
    storage.mode(bary.locations$barycenters) <- "double"
    barycenters <- as.matrix(bary.locations$barycenters)
  }

  # if (search == 1) { #use Naive search
  #   print('This is Naive Search')
  # } else if (search == 2)  { #use Tree search (default)
  #   print('This is Tree Search')
  # }

  #Calling the C++ function "eval_FEM_fd" in RPDE_interface.cpp
  evalmat = matrix(0,max(nrow(locations),nrow(incidence_matrix)),ncol(coeff))
  for (i in 1:ncol(coeff)){
    evalmat[,i] <- .Call("eval_FEM_fd", FEMbasis$mesh, locations, incidence_matrix, coeff[,i],
                         FEMbasis$order, redundancy, mydim, ndim, search, bary.locations, PACKAGE = "fdaPDE")
  }

  #Returning the evaluation matrix
  evalmat
}
