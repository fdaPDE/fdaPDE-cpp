CPP_smooth.GAM.FEM.time <- function(locations, time_locations, observations, FEMbasis, time_mesh,
                                    covariates = NULL, ndim, mydim, BC = NULL, incidence_matrix = NULL,
                                    areal.data.avg = TRUE, FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold, IC, FAMILY, 
                                    mu0 = NULL, max.steps.FPIRLS = 15, scale.param = NULL, threshold.FPIRLS = 0.0002020,
                                    search, bary.locations, optim, lambdaS = NULL, lambdaT = NULL,
                                    DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0,
                                    DOF.matrix = NULL, GCV.inflation.factor = 1,
                                    lambda.optimization.tolerance = 0.05) 
{
  # Indexes in C++ starts from 0, in R from 1
  FEMbasis$mesh$triangles <- FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges <- FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] <-
    FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

  if (is.null(covariates)) {
    covariates <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(DOF.matrix)) {
    DOF.matrix <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(locations)) {
    locations <- matrix(nrow = 0, ncol = 2)
  }

  if (is.null(incidence_matrix)) {
    incidence_matrix <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(IC)) {
    IC <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(BC$BC_indices)) {
    BC$BC_indices <- vector(length = 0)
  } else {
    BC$BC_indices <- as.vector(BC$BC_indices) - 1
  }

  if (is.null(BC$BC_values)) {
    BC$BC_values <- vector(length = 0)
  } else {
    BC$BC_values <- as.vector(BC$BC_values)
  }
  if (is.null(mu0)) {
    mu0 <- matrix(nrow = 0, ncol = 1)
  }

  if (is.null(scale.param)) {
    scale.param <- -1
  }

  if (is.null(lambdaS)) {
    lambdaS <- vector(length = 0)
  } else {
    lambdaS <- as.vector(lambdaS)
  }

  if (is.null(lambdaT)) {
    lambdaT <- vector(length = 0)
  } else {
    lambdaT <- as.vector(lambdaT)
  }

  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  time_mesh <- as.matrix(time_mesh)
  storage.mode(time_mesh) <- "double"
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
  storage.mode(BC$BC_values) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <- "integer"
  storage.mode(FAMILY) <- "character"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(mu0) <- "double"
  storage.mode(scale.param) <- "double"
  storage.mode(threshold.FPIRLS) <- "double"
  FLAG_MASS <- as.integer(FLAG_MASS)
  storage.mode(FLAG_MASS) <- "integer"
  FLAG_PARABOLIC <- as.integer(FLAG_PARABOLIC)
  storage.mode(FLAG_PARABOLIC) <- "integer"
  FLAG_ITERATIVE<-as.integer(FLAG_ITERATIVE)
  storage.mode(FLAG_ITERATIVE)<-"integer"
  storage.mode(max.steps) <- "integer"
  storage.mode(threshold) <- "double"
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"
  storage.mode(lambdaS) <- "double"
  storage.mode(lambdaT) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"

  ## IC estimation for parabolic smoothing from the first column of observations
  ICsol <- NA
  # empty dof matrix
  DOF.matrix_IC <- matrix(nrow = 0, ncol = 1)
  if (nrow(IC) == 0 && FLAG_PARABOLIC) {
    NobsIC <- length(observations) %/% nrow(time_locations)
    notNAIC <- which(!is.na(observations[1:NobsIC]))
    observationsIC <- observations[notNAIC]

    if (nrow(locations) == 0) {
      locationsIC <- locations
    } else {
      locationsIC <- as.matrix(locations[notNAIC, ])
      storage.mode(locationsIC) <- "double"
    }

    if (nrow(covariates) == 0) {
      covariatesIC <- covariates
    } else {
      covariatesIC <- covariates[notNAIC, ]
      covariatesIC <- as.matrix(covariatesIC)
      storage.mode(covariatesIC) <- "double"
    }

    ## set of lambdas for GCV in IC estimation
    lambdaSIC <- lambdaS
    lambdaSIC <- as.matrix(lambdaSIC)
    storage.mode(lambdaSIC) <- "double"
    ## call the smoothing function with initial observations to estimates the IC

    ICsol <- .Call(
      "gam_Laplace", locationsIC, bary.locations, observationsIC,
      FEMbasis$mesh, FEMbasis$order, mydim, ndim, covariatesIC, BC$BC_indices,
      BC$BC_values, incidence_matrix, areal.data.avg, FAMILY,
      max.steps.FPIRLS, threshold.FPIRLS, mu0, scale.param, search,
      optim, lambdaSIC, DOF.stochastic.realizations,
      DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor,
      lambda.optimization.tolerance,
      PACKAGE = "fdaPDE"
    )

    if (nrow(covariates) != 0) {
      betaIC = matrix(data=ICsol[[5]],nrow=ncol(covariatesIC),ncol=length(lambdaSIC))
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC
      # estimation
      covariates <- covariates[(NobsIC + 1):nrow(covariates), ]
      covariates <- as.matrix(covariates)
    } else {
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC
      # estimation
      betaIC <- NULL
    }
    ## return a FEM object containing IC estimates with best lambda and best lambda index
    ICsol <-
      list(
        IC.FEM = FEM(ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ], FEMbasis),
        bestlambdaindex = ICsol[[4]] + 1,
        bestlambda = lambdaSIC[ICsol[[4]] + 1], beta = betaIC, fn.eval=ICsol[[13]]
      )
    time_locations <- time_locations[2:nrow(time_locations)]
    observations <- observations[(NobsIC + 1):length(observations)]
  }
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"

  M <- ifelse(FLAG_PARABOLIC, length(time_mesh) - 1, length(time_mesh) + 2)
  BC$BC_indices <- rep((0:(M - 1)) * nrow(FEMbasis$mesh$nodes),
    each = length(BC$BC_indices)
  ) + rep(BC$BC_indices, M)
  BC$BC_values <- rep(BC$BC_values, M)
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"


  ## Call C++ function
  bigsol <- .Call(
    "gam_Laplace_time", locations, bary.locations, time_locations,
    observations, FEMbasis$mesh, time_mesh, FEMbasis$order, mydim, ndim,
    covariates, BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg,
    FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold, IC, FAMILY, max.steps.FPIRLS, threshold.FPIRLS,
    mu0, scale.param, search, optim, lambdaS, lambdaT,
    DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix,
    GCV.inflation.factor, lambda.optimization.tolerance,
    PACKAGE = "fdaPDE"
  )
  return(c(bigsol, ICsol))
}

CPP_smooth.manifold.GAM.FEM.time <- function(locations, time_locations, observations, FEMbasis, time_mesh,
                                    covariates = NULL, ndim, mydim, BC = NULL, incidence_matrix = NULL,
                                    areal.data.avg = TRUE, FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold, IC, FAMILY, 
                                    mu0 = NULL, max.steps.FPIRLS = 15, scale.param = NULL, threshold.FPIRLS = 0.0002020,
                                    search, bary.locations, optim, lambdaS = NULL, lambdaT = NULL,
                                    DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0,
                                    DOF.matrix = NULL, GCV.inflation.factor = 1,
                                    lambda.optimization.tolerance = 0.05) 
{
  # Indexes in C++ starts from 0, in R from 1
  FEMbasis$mesh$triangles <- FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges <- FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] <-
    FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  
  if (is.null(covariates)) {
    covariates <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(DOF.matrix)) {
    DOF.matrix <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(locations)) {
    locations <- matrix(nrow = 0, ncol = 2)
  }
  
  if (is.null(incidence_matrix)) {
    incidence_matrix <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(IC)) {
    IC <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(BC$BC_indices)) {
    BC$BC_indices <- vector(length = 0)
  } else {
    BC$BC_indices <- as.vector(BC$BC_indices) - 1
  }
  
  if (is.null(BC$BC_values)) {
    BC$BC_values <- vector(length = 0)
  } else {
    BC$BC_values <- as.vector(BC$BC_values)
  }
  if (is.null(mu0)) {
    mu0 <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(scale.param)) {
    scale.param <- -1
  }
  
  if (is.null(lambdaS)) {
    lambdaS <- vector(length = 0)
  } else {
    lambdaS <- as.vector(lambdaS)
  }
  
  if (is.null(lambdaT)) {
    lambdaT <- vector(length = 0)
  } else {
    lambdaT <- as.vector(lambdaT)
  }
  
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  time_mesh <- as.matrix(time_mesh)
  storage.mode(time_mesh) <- "double"
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
  storage.mode(BC$BC_values) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <- "integer"
  storage.mode(FAMILY) <- "character"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(mu0) <- "double"
  storage.mode(scale.param) <- "double"
  storage.mode(threshold.FPIRLS) <- "double"
  FLAG_MASS <- as.integer(FLAG_MASS)
  storage.mode(FLAG_MASS) <- "integer"
  FLAG_PARABOLIC <- as.integer(FLAG_PARABOLIC)
  storage.mode(FLAG_PARABOLIC) <- "integer"
  FLAG_ITERATIVE<-as.integer(FLAG_ITERATIVE)
  storage.mode(FLAG_ITERATIVE)<-"integer"
  storage.mode(max.steps) <- "integer"
  storage.mode(threshold) <- "double"
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"
  storage.mode(lambdaS) <- "double"
  storage.mode(lambdaT) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"
  
  ## IC estimation for parabolic smoothing from the first column of observations
  ICsol <- NA
  # empty dof matrix
  DOF.matrix_IC <- matrix(nrow = 0, ncol = 1)
  if (nrow(IC) == 0 && FLAG_PARABOLIC) {
    NobsIC <- length(observations) %/% nrow(time_locations)
    notNAIC <- which(!is.na(observations[1:NobsIC]))
    observationsIC <- observations[notNAIC]
    
    if (nrow(locations) == 0) {
      locationsIC <- locations
    } else {
      locationsIC <- as.matrix(locations[notNAIC, ])
      storage.mode(locationsIC) <- "double"
    }
    
    if (nrow(covariates) == 0) {
      covariatesIC <- covariates
    } else {
      covariatesIC <- covariates[notNAIC, ]
      covariatesIC <- as.matrix(covariatesIC)
      storage.mode(covariatesIC) <- "double"
    }
    
    ## set of lambdas for GCV in IC estimation
    lambdaSIC <- lambdaS
    lambdaSIC <- as.matrix(lambdaSIC)
    storage.mode(lambdaSIC) <- "double"
    ## call the smoothing function with initial observations to estimates the IC
    
    ICsol <- .Call(
      "gam_Laplace", locationsIC, bary.locations, observationsIC,
      FEMbasis$mesh, FEMbasis$order, mydim, ndim, covariatesIC, BC$BC_indices,
      BC$BC_values, incidence_matrix, areal.data.avg, FAMILY,
      max.steps.FPIRLS, threshold.FPIRLS, mu0, scale.param, search,
      optim, lambdaSIC, DOF.stochastic.realizations,
      DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor,
      lambda.optimization.tolerance,
      PACKAGE = "fdaPDE"
    )
    
    if (nrow(covariates) != 0) {
      betaIC = matrix(data=ICsol[[5]],nrow=ncol(covariatesIC),ncol=length(lambdaSIC))
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC
      # estimation
      covariates <- covariates[(NobsIC + 1):nrow(covariates), ]
      covariates <- as.matrix(covariates)
    } else {
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC
      # estimation
      betaIC <- NULL
    }
    ## return a FEM object containing IC estimates with best lambda and best lambda index
    ICsol <-
      list(
        IC.FEM = FEM(ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ], FEMbasis),
        bestlambdaindex = ICsol[[4]] + 1,
        bestlambda = lambdaSIC[ICsol[[4]] + 1], beta = betaIC, fn.eval=ICsol[[13]]
      )
    time_locations <- time_locations[2:nrow(time_locations)]
    observations <- observations[(NobsIC + 1):length(observations)]
  }
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  
  M <- ifelse(FLAG_PARABOLIC, length(time_mesh) - 1, length(time_mesh) + 2)
  BC$BC_indices <- rep((0:(M - 1)) * nrow(FEMbasis$mesh$nodes),
                       each = length(BC$BC_indices)
  ) + rep(BC$BC_indices, M)
  BC$BC_values <- rep(BC$BC_values, M)
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"
  
  
  ## Call C++ function
  bigsol <- .Call(
    "gam_Laplace_time", locations, bary.locations, time_locations,
    observations, FEMbasis$mesh, time_mesh, FEMbasis$order, mydim, ndim,
    covariates, BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg,
    FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold, IC, FAMILY, max.steps.FPIRLS, threshold.FPIRLS,
    mu0, scale.param, search, optim, lambdaS, lambdaT,
    DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix,
    GCV.inflation.factor, lambda.optimization.tolerance,
    PACKAGE = "fdaPDE"
  )
  return(c(bigsol, ICsol))
}

CPP_smooth.volume.GAM.FEM.time <- function(locations, time_locations, observations, FEMbasis, time_mesh,
                                             covariates = NULL, ndim, mydim, BC = NULL, incidence_matrix = NULL,
                                             areal.data.avg = TRUE, FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold, IC, FAMILY, 
                                             mu0 = NULL, max.steps.FPIRLS = 15, scale.param = NULL, threshold.FPIRLS = 0.0002020,
                                             search, bary.locations, optim, lambdaS = NULL, lambdaT = NULL,
                                             DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0,
                                             DOF.matrix = NULL, GCV.inflation.factor = 1,
                                             lambda.optimization.tolerance = 0.05) 
{
  # Indexes in C++ starts from 0, in R from 1
  FEMbasis$mesh$tetrahedrons = FEMbasis$mesh$tetrahedrons - 1
  FEMbasis$mesh$faces = FEMbasis$mesh$faces - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = 
    FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  
  if (is.null(covariates)) {
    covariates <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(DOF.matrix)) {
    DOF.matrix <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(locations)) {
    locations <- matrix(nrow = 0, ncol = 2)
  }
  
  if (is.null(incidence_matrix)) {
    incidence_matrix <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(IC)) {
    IC <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(BC$BC_indices)) {
    BC$BC_indices <- vector(length = 0)
  } else {
    BC$BC_indices <- as.vector(BC$BC_indices) - 1
  }
  
  if (is.null(BC$BC_values)) {
    BC$BC_values <- vector(length = 0)
  } else {
    BC$BC_values <- as.vector(BC$BC_values)
  }
  if (is.null(mu0)) {
    mu0 <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(scale.param)) {
    scale.param <- -1
  }
  
  if (is.null(lambdaS)) {
    lambdaS <- vector(length = 0)
  } else {
    lambdaS <- as.vector(lambdaS)
  }
  
  if (is.null(lambdaT)) {
    lambdaT <- vector(length = 0)
  } else {
    lambdaT <- as.vector(lambdaT)
  }
  
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  time_mesh <- as.matrix(time_mesh)
  storage.mode(time_mesh) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$tetrahedrons) <- "integer"
  storage.mode(FEMbasis$mesh$faces) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <- "integer"
  storage.mode(FAMILY) <- "character"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(mu0) <- "double"
  storage.mode(scale.param) <- "double"
  storage.mode(threshold.FPIRLS) <- "double"
  FLAG_MASS <- as.integer(FLAG_MASS)
  storage.mode(FLAG_MASS) <- "integer"
  FLAG_PARABOLIC <- as.integer(FLAG_PARABOLIC)
  storage.mode(FLAG_PARABOLIC) <- "integer"
  FLAG_ITERATIVE<-as.integer(FLAG_ITERATIVE)
  storage.mode(FLAG_ITERATIVE)<-"integer"
  storage.mode(max.steps) <- "integer"
  storage.mode(threshold) <- "double"
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"
  storage.mode(lambdaS) <- "double"
  storage.mode(lambdaT) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"
  
  ## IC estimation for parabolic smoothing from the first column of observations
  ICsol <- NA
  # empty dof matrix
  DOF.matrix_IC <- matrix(nrow = 0, ncol = 1)
  if (nrow(IC) == 0 && FLAG_PARABOLIC) {
    NobsIC <- length(observations) %/% nrow(time_locations)
    notNAIC <- which(!is.na(observations[1:NobsIC]))
    observationsIC <- observations[notNAIC]
    
    if (nrow(locations) == 0) {
      locationsIC <- locations
    } else {
      locationsIC <- as.matrix(locations[notNAIC, ])
      storage.mode(locationsIC) <- "double"
    }
    
    if (nrow(covariates) == 0) {
      covariatesIC <- covariates
    } else {
      covariatesIC <- covariates[notNAIC, ]
      covariatesIC <- as.matrix(covariatesIC)
      storage.mode(covariatesIC) <- "double"
    }
    
    ## set of lambdas for GCV in IC estimation
    lambdaSIC <- lambdaS
    lambdaSIC <- as.matrix(lambdaSIC)
    storage.mode(lambdaSIC) <- "double"
    ## call the smoothing function with initial observations to estimates the IC
    
    ICsol <- .Call(
      "gam_Laplace", locationsIC, bary.locations, observationsIC,
      FEMbasis$mesh, FEMbasis$order, mydim, ndim, covariatesIC, BC$BC_indices,
      BC$BC_values, incidence_matrix, areal.data.avg, FAMILY,
      max.steps.FPIRLS, threshold.FPIRLS, mu0, scale.param, search,
      optim, lambdaSIC, DOF.stochastic.realizations,
      DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor,
      lambda.optimization.tolerance,
      PACKAGE = "fdaPDE"
    )
    
    if (nrow(covariates) != 0) {
      betaIC = matrix(data=ICsol[[5]],nrow=ncol(covariatesIC),ncol=length(lambdaSIC))
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC
      # estimation
      covariates <- covariates[(NobsIC + 1):nrow(covariates), ]
      covariates <- as.matrix(covariates)
    } else {
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC
      # estimation
      betaIC <- NULL
    }
    ## return a FEM object containing IC estimates with best lambda and best lambda index
    ICsol <-
      list(
        IC.FEM = FEM(ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ], FEMbasis),
        bestlambdaindex = ICsol[[4]] + 1,
        bestlambda = lambdaSIC[ICsol[[4]] + 1], beta = betaIC, fn.eval=ICsol[[13]]
      )
    time_locations <- time_locations[2:nrow(time_locations)]
    observations <- observations[(NobsIC + 1):length(observations)]
  }
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  
  M <- ifelse(FLAG_PARABOLIC, length(time_mesh) - 1, length(time_mesh) + 2)
  BC$BC_indices <- rep((0:(M - 1)) * nrow(FEMbasis$mesh$nodes),
                       each = length(BC$BC_indices)
  ) + rep(BC$BC_indices, M)
  BC$BC_values <- rep(BC$BC_values, M)
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"
  
  
  ## Call C++ function
  bigsol <- .Call(
    "gam_Laplace_time", locations, bary.locations, time_locations,
    observations, FEMbasis$mesh, time_mesh, FEMbasis$order, mydim, ndim,
    covariates, BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg,
    FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold, IC, FAMILY, max.steps.FPIRLS, threshold.FPIRLS,
    mu0, scale.param, search, optim, lambdaS, lambdaT,
    DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix,
    GCV.inflation.factor, lambda.optimization.tolerance,
    PACKAGE = "fdaPDE"
  )
  return(c(bigsol, ICsol))
}

CPP_smooth.graph.GAM.FEM.time <- function(locations, time_locations, observations, FEMbasis, time_mesh,
                                             covariates = NULL, ndim, mydim, BC = NULL, incidence_matrix = NULL,
                                             areal.data.avg = TRUE, FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold, IC, FAMILY, 
                                             mu0 = NULL, max.steps.FPIRLS = 15, scale.param = NULL, threshold.FPIRLS = 0.0002020,
                                             search, bary.locations, optim, lambdaS = NULL, lambdaT = NULL,
                                             DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0,
                                             DOF.matrix = NULL, GCV.inflation.factor = 1,
                                             lambda.optimization.tolerance = 0.05) 
{
  # Indexes in C++ starts from 0, in R from 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  num_sides = 2*dim(FEMbasis$mesh$edges)[1] 
  for(i in 1:num_sides){
    if( dim(FEMbasis$mesh$neighbors[[i]] )[1] > 0)
      FEMbasis$mesh$neighbors[[i]] = FEMbasis$mesh$neighbors[[i]] - 1
  }
  
  
  if (is.null(covariates)) {
    covariates <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(DOF.matrix)) {
    DOF.matrix <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(locations)) {
    locations <- matrix(nrow = 0, ncol = 2)
  }
  
  if (is.null(incidence_matrix)) {
    incidence_matrix <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(IC)) {
    IC <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(BC$BC_indices)) {
    BC$BC_indices <- vector(length = 0)
  } else {
    BC$BC_indices <- as.vector(BC$BC_indices) - 1
  }
  
  if (is.null(BC$BC_values)) {
    BC$BC_values <- vector(length = 0)
  } else {
    BC$BC_values <- as.vector(BC$BC_values)
  }
  if (is.null(mu0)) {
    mu0 <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(scale.param)) {
    scale.param <- -1
  }
  
  if (is.null(lambdaS)) {
    lambdaS <- vector(length = 0)
  } else {
    lambdaS <- as.vector(lambdaS)
  }
  
  if (is.null(lambdaT)) {
    lambdaT <- vector(length = 0)
  } else {
    lambdaT <- as.vector(lambdaT)
  }
  
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  time_mesh <- as.matrix(time_mesh)
  storage.mode(time_mesh) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  for(i in 1:num_sides)
    storage.mode(FEMbasis$mesh$neighbors[[i]]) <- "integer" 
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <- "integer"
  storage.mode(FAMILY) <- "character"
  storage.mode(max.steps.FPIRLS) <- "integer"
  storage.mode(mu0) <- "double"
  storage.mode(scale.param) <- "double"
  storage.mode(threshold.FPIRLS) <- "double"
  FLAG_MASS <- as.integer(FLAG_MASS)
  storage.mode(FLAG_MASS) <- "integer"
  FLAG_PARABOLIC <- as.integer(FLAG_PARABOLIC)
  storage.mode(FLAG_PARABOLIC) <- "integer"
  FLAG_ITERATIVE<-as.integer(FLAG_ITERATIVE)
  storage.mode(FLAG_ITERATIVE)<-"integer"
  storage.mode(max.steps) <- "integer"
  storage.mode(threshold) <- "double"
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  storage.mode(search) <- "integer"
  storage.mode(optim) <- "integer"
  storage.mode(lambdaS) <- "double"
  storage.mode(lambdaT) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"
  
  ## IC estimation for parabolic smoothing from the first column of observations
  ICsol <- NA
  # empty dof matrix
  DOF.matrix_IC <- matrix(nrow = 0, ncol = 1)
  if (nrow(IC) == 0 && FLAG_PARABOLIC) {
    NobsIC <- length(observations) %/% nrow(time_locations)
    notNAIC <- which(!is.na(observations[1:NobsIC]))
    observationsIC <- observations[notNAIC]
    
    if (nrow(locations) == 0) {
      locationsIC <- locations
    } else {
      locationsIC <- as.matrix(locations[notNAIC, ])
      storage.mode(locationsIC) <- "double"
    }
    
    if (nrow(covariates) == 0) {
      covariatesIC <- covariates
    } else {
      covariatesIC <- covariates[notNAIC, ]
      covariatesIC <- as.matrix(covariatesIC)
      storage.mode(covariatesIC) <- "double"
    }
    
    ## set of lambdas for GCV in IC estimation
    lambdaSIC <- lambdaS
    lambdaSIC <- as.matrix(lambdaSIC)
    storage.mode(lambdaSIC) <- "double"
    ## call the smoothing function with initial observations to estimates the IC
    
    ICsol <- .Call(
      "gam_Laplace", locationsIC, bary.locations, observationsIC,
      FEMbasis$mesh, FEMbasis$order, mydim, ndim, covariatesIC, BC$BC_indices,
      BC$BC_values, incidence_matrix, areal.data.avg, FAMILY,
      max.steps.FPIRLS, threshold.FPIRLS, mu0, scale.param, search,
      optim, lambdaSIC, DOF.stochastic.realizations,
      DOF.stochastic.seed, DOF.matrix_IC, GCV.inflation.factor,
      lambda.optimization.tolerance,
      PACKAGE = "fdaPDE"
    )
    
    if (nrow(covariates) != 0) {
      betaIC = matrix(data=ICsol[[5]],nrow=ncol(covariatesIC),ncol=length(lambdaSIC))
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC
      # estimation
      covariates <- covariates[(NobsIC + 1):nrow(covariates), ]
      covariates <- as.matrix(covariates)
    } else {
      IC <- ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ICsol[[4]] + 1] # best IC
      # estimation
      betaIC <- NULL
    }
    ## return a FEM object containing IC estimates with best lambda and best lambda index
    ICsol <-
      list(
        IC.FEM = FEM(ICsol[[1]][1:nrow(FEMbasis$mesh$nodes), ], FEMbasis),
        bestlambdaindex = ICsol[[4]] + 1,
        bestlambda = lambdaSIC[ICsol[[4]] + 1], beta = betaIC, fn.eval=ICsol[[13]]
      )
    time_locations <- time_locations[2:nrow(time_locations)]
    observations <- observations[(NobsIC + 1):length(observations)]
  }
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  
  M <- ifelse(FLAG_PARABOLIC, length(time_mesh) - 1, length(time_mesh) + 2)
  BC$BC_indices <- rep((0:(M - 1)) * nrow(FEMbasis$mesh$nodes),
                       each = length(BC$BC_indices)
  ) + rep(BC$BC_indices, M)
  BC$BC_values <- rep(BC$BC_values, M)
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"
  
  
  ## Call C++ function
  bigsol <- .Call(
    "gam_Laplace_time", locations, bary.locations, time_locations,
    observations, FEMbasis$mesh, time_mesh, FEMbasis$order, mydim, ndim,
    covariates, BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg,
    FLAG_MASS, FLAG_PARABOLIC, FLAG_ITERATIVE, max.steps, threshold, IC, FAMILY, max.steps.FPIRLS, threshold.FPIRLS,
    mu0, scale.param, search, optim, lambdaS, lambdaT,
    DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix,
    GCV.inflation.factor, lambda.optimization.tolerance,
    PACKAGE = "fdaPDE"
  )
  return(c(bigsol, ICsol))
}
