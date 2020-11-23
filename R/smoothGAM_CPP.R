CPP_smooth.GAM.FEM<-function(locations, observations, FEMbasis, covariates = NULL, ndim, mydim, BC = NULL, incidence_matrix = NULL, areal.data.avg = FALSE, FAMILY, mu0 = NULL, max.steps.FPIRLS = 15, scale.param = NULL, threshold.FPIRLS = 0.0002020, search, bary.locations, optim, lambda = NULL, DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0, DOF.matrix = NULL, GCV.inflation.factor = 1.8, lambda.optimization.tolerance = 0.05)
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

  if(is.null(DOF.matrix))
  {
    DOF.matrix<-matrix(nrow = 0, ncol = 1)
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
  
  if(is.null(lambda))
  {
    lambda<-vector(length=0)
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
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <-"double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <-"integer"
  storage.mode(FAMILY) <- "character"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(mu0) <- "double"
  storage.mode(scale.param) <- "double"
  storage.mode(threshold.FPIRLS) <- "double"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"  
  storage.mode(lambda) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"

  ## Call C++ function
  bigsol <- .Call("gam_Laplace", locations, bary.locations, observations, FEMbasis$mesh, FEMbasis$order,
                 mydim, ndim, covariates, BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg,
                 FAMILY, max.steps.FPIRLS, threshold.FPIRLS, mu0, scale.param, search,
                 optim, lambda, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix, GCV.inflation.factor, lambda.optimization.tolerance, PACKAGE = "fdaPDE")
  
  return(bigsol)
}


CPP_smooth.GAM.FEM.PDE.basis<-function(locations, observations, FEMbasis, covariates = NULL, PDE_parameters, ndim, mydim, BC = NULL, incidence_matrix = NULL, areal.data.avg = FALSE, FAMILY, mu0 = NULL, max.steps.FPIRLS = 15, scale.param = NULL, threshold.FPIRLS = 0.0004, search, bary.locations, optim, lambda = NULL, DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0, DOF.matrix = NULL, GCV.inflation.factor = 1.8, lambda.optimization.tolerance = 0.05)
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

  if(is.null(DOF.matrix))
  {
    DOF.matrix<-matrix(nrow = 0, ncol = 1)
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
  
  if(is.null(lambda))
  {
    lambda<-vector(length=0)
  }
  

  threshold.FPIRLS<-0.0004 #??
  
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
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer" 
  storage.mode(PDE_parameters$K) <- "double"
  storage.mode(PDE_parameters$b) <- "double"
  storage.mode(PDE_parameters$c) <- "double"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <-"double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <-"integer"
  storage.mode(FAMILY) <- "character"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(mu0) <- "double"
  storage.mode(scale.param) <- "double"
  storage.mode(threshold.FPIRLS) <- "double"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"  
  storage.mode(lambda) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"
  
  ## Call C++ function
  bigsol <- .Call("gam_PDE", locations, bary.locations, observations, FEMbasis$mesh, FEMbasis$order,
                 mydim, ndim, PDE_parameters$K, PDE_parameters$b, PDE_parameters$c, covariates, BC$BC_indices, BC$BC_values,
                 incidence_matrix, areal.data.avg, FAMILY, max.steps.FPIRLS, threshold.FPIRLS, mu0, scale.param, search,
                 optim, lambda, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix, GCV.inflation.factor, lambda.optimization.tolerance, PACKAGE = "fdaPDE")
  
  return(bigsol)
}

CPP_smooth.GAM.FEM.PDE.sv.basis<-function(locations, observations, FEMbasis, covariates = NULL, PDE_parameters, ndim, mydim, BC = NULL, incidence_matrix = NULL, areal.data.avg = FALSE, FAMILY, mu0 = NULL, max.steps.FPIRLS = 15, scale.param = NULL, threshold.FPIRLS = 0.0004, search, bary.locations, optim, lambda = NULL, DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0, DOF.matrix = NULL, GCV.inflation.factor = 1.8, lambda.optimization.tolerance = 0.05)
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
  
  if(is.null(DOF.matrix))
  {
    DOF.matrix<-matrix(nrow = 0, ncol = 1)
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
  
  if(is.null(lambda))
  {
    lambda<-vector(length=0)
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
  storage.mode(PDE_param_eval$K) <- "double"
  storage.mode(PDE_param_eval$b) <- "double"
  storage.mode(PDE_param_eval$c) <- "double"
  storage.mode(PDE_param_eval$u) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <-"double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <-"integer"
  storage.mode(FAMILY) <- "character"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(mu0) <- "double"
  storage.mode(scale.param) <- "double"
  storage.mode(threshold.FPIRLS) <- "double"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"  
  storage.mode(lambda) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"

  ## Call C++ function
  bigsol <- .Call("gam_PDE_space_varying", locations, bary.locations, observations, FEMbasis$mesh, FEMbasis$order,
                 mydim, ndim, PDE_param_eval$K, PDE_param_eval$b, PDE_param_eval$c, PDE_param_eval$u, covariates, BC$BC_indices, BC$BC_values,
                 incidence_matrix, areal.data.avg, FAMILY, max.steps.FPIRLS, threshold.FPIRLS, mu0, scale.param,
                 search, optim, lambda, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix, GCV.inflation.factor, lambda.optimization.tolerance, PACKAGE = "fdaPDE")
  
  return(bigsol)
}

CPP_smooth.manifold.GAM.FEM.basis<-function(locations, observations, FEMbasis, covariates = NULL, ndim, mydim, BC = NULL, incidence_matrix = NULL, areal.data.avg = FALSE, FAMILY, mu0 = NULL, max.steps.FPIRLS = 15, scale.param = NULL, threshold.FPIRLS = 0.0004, search, bary.locations, optim, lambda = NULL, DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0, DOF.matrix = NULL, GCV.inflation.factor = 1.8, lambda.optimization.tolerance = 0.05)
{
  
  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  

  max.steps.FPIRLS = max.steps.FPIRLS - 1
  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(DOF.matrix))
  {
    DOF.matrix<-matrix(nrow = 0, ncol = 1)
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
  
  if(is.null(lambda))
  {
    lambda<-vector(length=0)
  }
  
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(FEMbasis$mesh$order) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"

  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <-"double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <-"integer"
  storage.mode(FAMILY) <- "character"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(mu0) <- "double"
  storage.mode(scale.param) <- "double"
  storage.mode(threshold.FPIRLS) <- "double"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"  
  storage.mode(lambda) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"
  
  ## Call C++ function
  bigsol <- .Call("gam_Laplace", locations, bary.locations, observations, FEMbasis$mesh, FEMbasis$mesh$order,
                mydim, ndim, covariates, BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg,
                FAMILY, max.steps.FPIRLS, threshold.FPIRLS, mu0, scale.param, search, 
                optim, lambda, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix, GCV.inflation.factor, lambda.optimization.tolerance, PACKAGE = "fdaPDE")
  
  return(bigsol)
}

CPP_smooth.volume.GAM.FEM.basis<-function(locations, observations, FEMbasis, covariates = NULL, ndim, mydim, BC = NULL, incidence_matrix = NULL, areal.data.avg = FALSE, FAMILY, mu0 = NULL, max.steps.FPIRLS = 15, scale.param = NULL, threshold.FPIRLS = 0.0004, search, bary.locations, optim, lambda = NULL, DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0, DOF.matrix = NULL, GCV.inflation.factor = 1.8, lambda.optimization.tolerance = 0.05)
{
  
  # Indexes in C++ starts from 0, in R from 1

  FEMbasis$mesh$tetrahedrons = FEMbasis$mesh$tetrahedrons - 1
  FEMbasis$mesh$faces = FEMbasis$mesh$faces - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  

  max.steps.FPIRLS = max.steps.FPIRLS - 1
  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(DOF.matrix))
  {
    DOF.matrix<-matrix(nrow = 0, ncol = 1)
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
  
  if(is.null(lambda))
  {
    lambda<-vector(length=0)
  }
  
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(FEMbasis$mesh$order) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$mesh$faces) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$tetrahedrons) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <-"double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <-"integer"
  storage.mode(FAMILY) <- "character"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(mu0) <- "double"
  storage.mode(scale.param) <- "double"
  storage.mode(threshold.FPIRLS) <- "double"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"  
  storage.mode(lambda) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"
  
  ## Call C++ function
  bigsol <- .Call("gam_Laplace", locations, bary.locations, observations, FEMbasis$mesh, FEMbasis$mesh$order,
                 mydim, ndim, covariates, BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg,
                 FAMILY, max.steps.FPIRLS, threshold.FPIRLS, mu0, scale.param, search,  
                 optim, lambda, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix, GCV.inflation.factor, lambda.optimization.tolerance, PACKAGE = "fdaPDE")
  
  return(bigsol)
}
