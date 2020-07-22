CPP_smooth.GAM.FEM<-function(locations, bary.locations ,observations, FEMbasis, lambda,
                                  covariates=NULL, incidence_matrix=NULL, ndim, mydim,
                                  BC=NULL, GCV, GCVMETHOD = 2, nrealizations=100, mu0=NULL, max.steps.FPIRLS=15, FAMILY, 
                                  scale.param=NULL, GCV.inflation.factor=1.8, threshold.FPIRLS=0.0002020, DOF, DOF_matrix=NULL, search, areal.data.avg = FALSE )
{
  # Indexes in C++ starts from 0, in R from 1, opporGCV.inflation.factor transformation
  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  max.steps.FPIRLS = max.steps.FPIRLS - 1
  
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
    locations<-matrix(nrow = 0, ncol = 2)
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

  if(is.null(mu0))
  {
    mu0<-matrix(nrow = 0, ncol = 1)
  }
 
  if(is.null(scale.param))
  {
    scale.param<- -1
  }

  
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <-"double"
  
  # Type conversion: form integer to bool(logical)
  GCV <- as.integer(GCV)
  storage.mode(GCV) <-"integer"
  DOF <- as.integer(DOF)
  storage.mode(DOF) <-"integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <-"integer"


  storage.mode(nrealizations) <- "integer"
  storage.mode(GCVMETHOD) <- "integer"
  storage.mode(FAMILY) <- "character"
  #storage.mode(FAMILY) <- "integer"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(mu0) <- "double"
  storage.mode(scale.param) <- "double"
  storage.mode(threshold.FPIRLS) <- "double"
  storage.mode(search) <- "integer"

  ## Call C++ function
  bigsol <- .Call("gam_Laplace", locations, bary.locations, observations, FEMbasis$mesh, FEMbasis$order,
                 mydim, ndim, lambda, covariates, incidence_matrix, BC$BC_indices, BC$BC_values,
                 GCV, GCVMETHOD, nrealizations, FAMILY, max.steps.FPIRLS, threshold.FPIRLS, GCV.inflation.factor, mu0, scale.param, DOF, DOF_matrix, search, areal.data.avg, PACKAGE = "fdaPDE")
  
  
  return(bigsol)
}


CPP_smooth.GAM.FEM.PDE.basis<-function(locations, bary.locations ,observations, FEMbasis, lambda, PDE_parameters,
                                  covariates=NULL, incidence_matrix=NULL, ndim, mydim,
                                  BC=NULL, GCV, GCVMETHOD = 2, nrealizations=100, mu0=NULL, max.steps.FPIRLS=15, FAMILY, 
                                  scale.param=NULL, GCV.inflation.factor=1.8, threshold.FPIRLS=0.0004, DOF, DOF_matrix=NULL, search, areal.data.avg = FALSE )
{
  # Indexes in C++ starts from 0, in R from 1, opporGCV.inflation.factor transformation
  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  max.steps.FPIRLS = max.steps.FPIRLS - 1
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
    locations<-matrix(nrow = 0, ncol = 2)
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

  if(is.null(mu0))
  {
    mu0<-matrix(nrow = 0, ncol = 1)
  }
  if(is.null(scale.param))
  {
    scale.param<- -1
  }

  threshold.FPIRLS<-0.0004
  
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <-"double"
  
  # Type conversion: form integer to bool(logical)
  GCV <- as.integer(GCV)
  storage.mode(GCV) <-"integer"
  DOF <- as.integer(DOF)
  storage.mode(DOF) <-"integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <-"integer"

  storage.mode(nrealizations) <- "integer"
  storage.mode(GCVMETHOD) <- "integer"
  storage.mode(FAMILY) <- "character"
  #storage.mode(FAMILY) <- "integer"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(mu0) <- "double"

  storage.mode(scale.param) <- "double"
  storage.mode(threshold.FPIRLS) <- "double"

  storage.mode(PDE_parameters$K) <- "double"
  storage.mode(PDE_parameters$b) <- "double"
  storage.mode(PDE_parameters$c) <- "double"
  storage.mode(search) <- "integer"
  ## Call C++ function
  bigsol <- .Call("gam_PDE", locations, bary.locations, observations, FEMbasis$mesh, FEMbasis$order,
                 mydim, ndim, lambda, PDE_parameters$K, PDE_parameters$b, PDE_parameters$c, covariates, incidence_matrix, BC$BC_indices, BC$BC_values,
                 GCV, GCVMETHOD, nrealizations, FAMILY, max.steps.FPIRLS, threshold.FPIRLS, GCV.inflation.factor, mu0, scale.param, DOF, DOF_matrix, search, areal.data.avg,PACKAGE = "fdaPDE")
  # GCV became DOF variable in C++, the motivation is if GCV = TRUE we need to compute DoF otherwise no. 
  
  return(bigsol)
}

CPP_smooth.GAM.FEM.PDE.sv.basis<-function(locations, bary.locations, observations, FEMbasis, lambda, PDE_parameters,
                                  covariates=NULL, incidence_matrix=NULL, ndim, mydim,
                                  BC=NULL, GCV, GCVMETHOD = 2, nrealizations=100, mu0=NULL, max.steps.FPIRLS=15, FAMILY, 
                                  scale.param=NULL, GCV.inflation.factor=1.8, threshold.FPIRLS=0.0004, DOF, DOF_matrix=NULL, search, areal.data.avg = FALSE )
{
  # Indexes in C++ starts from 0, in R from 1, opporGCV.inflation.factor transformation
  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  max.steps.FPIRLS = max.steps.FPIRLS - 1
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
    locations<-matrix(nrow = 0, ncol = 2)
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
    BC$BC_indices<-as.vector(BC$BC_indices) -1
  }
  
  if(is.null(BC$BC_values))
  {
    BC$BC_values<-vector(length=0)
  }else
  {
    BC$BC_values<-as.vector(BC$BC_values) 
  }

  if(is.null(mu0))
  {
    mu0<-matrix(nrow = 0, ncol = 1)
  }
  if(is.null(scale.param))
  {
    scale.param<- -1
  }

  PDE_param_eval = NULL
  points_eval = matrix(CPP_get_evaluations_points(mesh = FEMbasis$mesh, order = FEMbasis$order),ncol = 2)
  PDE_param_eval$K = (PDE_parameters$K)(points_eval)
  PDE_param_eval$b = (PDE_parameters$b)(points_eval)
  PDE_param_eval$c = (PDE_parameters$c)(points_eval)
  PDE_param_eval$u = (PDE_parameters$u)(points_eval)

  
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <-"double"
  
  # Type conversion: form integer to bool(logical)
  GCV <- as.integer(GCV)
  storage.mode(GCV) <-"integer"
  DOF <- as.integer(DOF)
  storage.mode(DOF) <-"integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <-"integer"

  storage.mode(nrealizations) <- "integer"
  storage.mode(GCVMETHOD) <- "integer"
  storage.mode(FAMILY) <- "character"
  #storage.mode(FAMILY) <- "integer"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(mu0) <- "double"
  storage.mode(scale.param) <- "double"
  storage.mode(threshold.FPIRLS) <- "double"

  storage.mode(PDE_param_eval$K) <- "double"
  storage.mode(PDE_param_eval$b) <- "double"
  storage.mode(PDE_param_eval$c) <- "double"
  storage.mode(PDE_param_eval$u) <- "double"

  storage.mode(search) <- "integer"
  ## Call C++ function
  bigsol <- .Call("gam_PDE_space_varying", locations, bary.locations, observations, FEMbasis$mesh, FEMbasis$order,
                 mydim, ndim, lambda, PDE_param_eval$K, PDE_param_eval$b, PDE_param_eval$c, PDE_param_eval$u, covariates, incidence_matrix, BC$BC_indices, BC$BC_values,
                 GCV, GCVMETHOD, nrealizations, FAMILY, max.steps.FPIRLS, threshold.FPIRLS, GCV.inflation.factor, mu0, scale.param, DOF, DOF_matrix, search, areal.data.avg,PACKAGE = "fdaPDE")
  # GCV became DOF variable in C++, the motivation is if GCV = TRUE we need to compute DoF otherwise no. 
  
  return(bigsol)
}

CPP_smooth.manifold.GAM.FEM.basis<-function(locations, bary.locations, observations, FEMbasis, lambda,
                                  covariates=NULL, incidence_matrix=NULL, ndim, mydim,
                                  BC=NULL, GCV, GCVMETHOD = 2, nrealizations=100, mu0=NULL, max.steps.FPIRLS=15, FAMILY, 
                                  scale.param=NULL, GCV.inflation.factor=1.8, threshold.FPIRLS=0.0004, DOF, DOF_matrix=NULL, search, areal.data.avg = FALSE )
{
  
  # C++ function for manifold works with vectors not with matrices
  
  FEMbasis$mesh$triangles=c(t(FEMbasis$mesh$triangles))
  FEMbasis$mesh$nodes=c(t(FEMbasis$mesh$nodes))
  
  # Indexes in C++ starts from 0, in R from 1, opporGCV.inflation.factor transformation
 
  FEMbasis$mesh$triangles=FEMbasis$mesh$triangles-1

  max.steps.FPIRLS = max.steps.FPIRLS - 1
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
    locations<-matrix(nrow = 0, ncol = 2)
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

  if(is.null(mu0))
  {
    mu0<-matrix(nrow = 0, ncol = 1)
  }
 
  if(is.null(scale.param))
  {
    scale.param<- -1
  }

  
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(FEMbasis$mesh$order) <- "integer"
  storage.mode(FEMbasis$mesh$nnodes) <- "integer"
  storage.mode(FEMbasis$mesh$ntriangles) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <-"double"
  
  # Type conversion: form integer to bool(logical)
  GCV <- as.integer(GCV)
  storage.mode(GCV) <-"integer"
  DOF <- as.integer(DOF)
  storage.mode(DOF) <-"integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <-"integer"
  
  storage.mode(nrealizations) <- "integer"
  storage.mode(GCVMETHOD) <- "integer"
  storage.mode(FAMILY) <- "character"
  #storage.mode(FAMILY) <- "integer"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(mu0) <- "double"
  storage.mode(scale.param) <- "double"
  storage.mode(threshold.FPIRLS) <- "double"
  storage.mode(search) <- "integer"
  ## Call C++ function
  bigsol <- .Call("gam_Laplace", locations, bary.locations, observations, FEMbasis$mesh, FEMbasis$mesh$order,
                 mydim, ndim, lambda, covariates, incidence_matrix, BC$BC_indices, BC$BC_values,
                 GCV, GCVMETHOD, nrealizations, FAMILY, max.steps.FPIRLS, threshold.FPIRLS, GCV.inflation.factor, mu0, scale.param, DOF, DOF_matrix, search, areal.data.avg,PACKAGE = "fdaPDE")
  # GCV became DOF variable in C++, the motivation is if GCV = TRUE we need to compute DoF otherwise no. 
  
  return(bigsol)
}

CPP_smooth.volume.GAM.FEM.basis<-function(locations, bary.locations, observations, FEMbasis, lambda,
                                  covariates=NULL, incidence_matrix=NULL, ndim, mydim,
                                  BC=NULL, GCV, GCVMETHOD = 2, nrealizations=100, mu0=NULL, max.steps.FPIRLS=15, FAMILY, 
                                  scale.param=NULL, GCV.inflation.factor=1.8, threshold.FPIRLS=0.0004, DOF, DOF_matrix=NULL, search, areal.data.avg = FALSE )
{
  
  # C++ function for volumetric works with vectors not with matrices
  
  FEMbasis$mesh$tetrahedrons=c(t(FEMbasis$mesh$tetrahedrons))
  FEMbasis$mesh$nodes=c(t(FEMbasis$mesh$nodes))
  
  # Indexes in C++ starts from 0, in R from 1, opporGCV.inflation.factor transformation

  FEMbasis$mesh$tetrahedrons=FEMbasis$mesh$tetrahedrons-1

  max.steps.FPIRLS = max.steps.FPIRLS - 1
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
    locations<-matrix(nrow = 0, ncol = 2)
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

  if(is.null(mu0))
  {
    mu0<-matrix(nrow = 0, ncol = 1)
  }
 
  if(is.null(scale.param))
  {
    scale.param<- -1
  }

  
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(FEMbasis$mesh$order) <- "integer"
  storage.mode(FEMbasis$mesh$nnodes) <- "integer"
  storage.mode(FEMbasis$mesh$ntetrahedrons) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$tetrahedrons) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <-"double"
  
  # Type conversion: form integer to bool(logical)
  GCV <- as.integer(GCV)
  storage.mode(GCV) <-"integer"
  DOF <- as.integer(DOF)
  storage.mode(DOF) <-"integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <-"integer"


  storage.mode(nrealizations) <- "integer"
  storage.mode(GCVMETHOD) <- "integer"
  storage.mode(FAMILY) <- "character"

  #storage.mode(FAMILY) <- "integer"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(mu0) <- "double"
  storage.mode(scale.param) <- "double"
  storage.mode(threshold.FPIRLS) <- "double"
  storage.mode(search) <- "integer"
  ## Call C++ function
  bigsol <- .Call("gam_Laplace", locations, bary.locations, observations, FEMbasis$mesh, FEMbasis$mesh$order,
                 mydim, ndim, lambda, covariates, incidence_matrix, BC$BC_indices, BC$BC_values,
                 GCV, GCVMETHOD, nrealizations, FAMILY, max.steps.FPIRLS, threshold.FPIRLS, GCV.inflation.factor, mu0, scale.param, DOF, DOF_matrix, search, areal.data.avg, PACKAGE = "fdaPDE")
  # GCV became DOF variable in C++, the motivation is if GCV = TRUE we need to compute DoF otherwise no. 
  
  return(bigsol)
}


