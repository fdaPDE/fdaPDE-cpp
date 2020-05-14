CPP_smooth.volume.FEM.time<-function(locations, bary.locations, time_locations, observations, FEMbasis, time_mesh, lambdaS, lambdaT, covariates = NULL, incidence_matrix = NULL, ndim, mydim, BC = NULL, FLAG_MASS, FLAG_PARABOLIC, IC, GCV,GCVMETHOD = 2, nrealizations = 100, DOF=TRUE, DOF_matrix=NULL, search)
{

  # C++ function for volumetric works with vectors not with matrices

  FEMbasis$mesh$tetrahedrons=c(t(FEMbasis$mesh$tetrahedrons))
  FEMbasis$mesh$nodes=c(t(FEMbasis$mesh$nodes))

  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$tetrahedrons=FEMbasis$mesh$tetrahedrons-1

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

  if(is.null(IC))
  {
    IC<-matrix(nrow = 0, ncol = 1)
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
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  time_mesh <- as.matrix(time_mesh)
  storage.mode(time_mesh) <- "double"
  observations <- as.vector(observations)
  storage.mode(observations) <- "double"
  storage.mode(FEMbasis$mesh$order) <- "integer"
  storage.mode(FEMbasis$mesh$nnodes) <- "integer"
  storage.mode(FEMbasis$mesh$ntetrahedrons) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$tetrahedrons) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(lambdaS) <- "double"
  storage.mode(lambdaT) <- "double"
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values)  <- "double"
  GCV <- as.integer(GCV)
  storage.mode(GCV) <-"integer"
  DOF <- as.integer(DOF)
  storage.mode(DOF) <- "integer"

  FLAG_MASS <- as.integer(FLAG_MASS)
  storage.mode(FLAG_MASS) <-"integer"

  FLAG_PARABOLIC <- as.integer(FLAG_PARABOLIC)
  storage.mode(FLAG_PARABOLIC) <-"integer"

  storage.mode(nrealizations) <- "integer"
  storage.mode(GCVMETHOD) <- "integer"
  storage.mode(search) <- "integer"

  ## IC estimation for parabolic smoothing from the first column of observations
  ICsol=NA
  if(nrow(IC)==0 && FLAG_PARABOLIC)
  {
    NobsIC = length(observations)%/%nrow(time_locations)

    if(nrow(covariates)==0)
    covariatesIC = covariates
    else
    {
      covariatesIC = covariates[1:NobsIC,]
      covariatesIC = as.matrix(covariatesIC)
      storage.mode(covariatesIC) <- "double"
    }

    ## set of lambdas for GCV in IC estimation
    lambdaSIC <- 10^seq(-7,3,0.1)
    lambdaSIC <- as.matrix(lambdaSIC)
    storage.mode(lambdaSIC) <- "double"
    ## call the smoothing function with initial observations to estimates the IC
    ICsol <- .Call("regression_Laplace", locations, bary.locations, observations[1:NobsIC],
      FEMbasis$mesh, FEMbasis$order, mydim, ndim, lambdaSIC, covariatesIC,
      incidence_matrix, BC$BC_indices, BC$BC_values,
      T, as.integer(1), nrealizations, T, DOF_matrix, search, PACKAGE = "fdaPDE")

    ## shifting the lambdas interval if the best lambda is the smaller one and retry smoothing
    if((ICsol[[4]][1]+1)==1)
    {
      lambdaSIC <- 10^seq(-9,-7,0.1)
      lambdaSIC <- as.matrix(lambdaSIC)
      storage.mode(lambdaSIC) <- "double"
      ICsol <- .Call("regression_Laplace", locations, bary.locations, observations[1:NobsIC],
        FEMbasis$mesh, FEMbasis$order, mydim, ndim, lambdaSIC, covariatesIC,
        incidence_matrix, BC$BC_indices, BC$BC_values,
        T, as.integer(1), nrealizations, T, DOF_matrix, search, PACKAGE = "fdaPDE")
    }
    else
    {
      ## shifting the lambdas interval if the best lambda is the higher one and retry smoothing
      if((ICsol[[4]][1]+1)==length(lambdaSIC))
      {
        lambdaSIC <- 10^seq(3,5,0.1)
        lambdaSIC <- as.matrix(lambdaSIC)
        storage.mode(lambdaSIC) <- "double"
        ICsol <- .Call("regression_Laplace", locations, bary.locations, observations[1:NobsIC],
          FEMbasis$mesh, FEMbasis$order, mydim, ndim, lambdaSIC, covariatesIC,
          incidence_matrix, BC$BC_indices, BC$BC_values,
          T, as.integer(1), nrealizations, T, DOF_matrix, search, PACKAGE = "fdaPDE")
      }
    }

    if(nrow(covariates)!=0)
    {
      betaIC = ICsol[[5]]
      IC = ICsol[[1]][1:FEMbasis$mesh$nnodes,ICsol[[4]][1]+1] ## best IC estimation
      covariates=covariates[(NobsIC+1):nrow(covariates),]
      covariates <- as.matrix(covariates)
    }
    else
    {
      IC = ICsol[[1]][1:FEMbasis$mesh$nnodes,ICsol[[4]][1]+1] ## best IC estimation
      betaIC = NULL
    }
    ## return a FEM object containing IC estimates with best lambda and best lambda index
    ICsol = list(IC.FEM=FEM(ICsol[[1]][1:FEMbasis$mesh$nnodes,],FEMbasis),bestlambdaindex=ICsol[[4]][1]+1,bestlambda=lambdaSIC[ICsol[[4]][1]+1],beta=betaIC)
    time_locations=time_locations[2:nrow(time_locations)]
    observations = observations[(NobsIC+1):length(observations)]
  }
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"

  M = ifelse(FLAG_PARABOLIC,length(time_mesh)-1,length(time_mesh) + 2);
  BC$BC_indices = rep((0:(M-1))*FEMbasis$mesh$nnodes,each=length(BC$BC_indices)) + rep(BC$BC_indices,M)
  BC$BC_values = rep(BC$BC_values,M)
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <-"double"
  storage.mode(search) <- "integer"

  ## Call C++ function
  bigsol <- .Call("regression_Laplace_time", locations, bary.locations, time_locations, observations, FEMbasis$mesh, time_mesh, FEMbasis$mesh$order, mydim, ndim, lambdaS, lambdaT, covariates,
    incidence_matrix, BC$BC_indices, BC$BC_values, FLAG_MASS, FLAG_PARABOLIC,
    IC, GCV, GCVMETHOD, nrealizations, DOF, DOF_matrix, search, PACKAGE = "fdaPDE")

  return(c(bigsol,ICsol))
}

CPP_eval.volume.FEM.time = function(FEM.time, locations, time_locations, incidence_matrix, FLAG_PARABOLIC, redundancy, ndim, mydim, search, bary.locations)
{
  FEMbasis = FEM.time$FEMbasis

  # C++ function for volumetric works with vectors not with matrices

  FEMbasis$mesh$tetrahedrons=c(t(FEMbasis$mesh$tetrahedrons))
  FEMbasis$mesh$nodes=c(t(FEMbasis$mesh$nodes))

  #NEW: riporto shift indici in R
  FEMbasis$mesh$tetrahedrons=FEMbasis$mesh$tetrahedrons-1

  # Imposing types, this is necessary for correct reading from C++
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$tetrahedrons) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  storage.mode(FEM.time$mesh_time) <- "double"
  coeff <- as.matrix(FEM.time$coeff)
  storage.mode(coeff) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(locations) <- "double"
  storage.mode(redundancy) <- "integer"
  storage.mode(FLAG_PARABOLIC) <- "integer"
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
    evalmat[,i] <- .Call("eval_FEM_time", FEMbasis$mesh, FEM.time$mesh_time, locations, time_locations, incidence_matrix, coeff[,i],
     FEMbasis$order, redundancy, FLAG_PARABOLIC, mydim, ndim, search, bary.locations, PACKAGE = "fdaPDE")
  }

  #Returning the evaluation matrix
  evalmat
}
