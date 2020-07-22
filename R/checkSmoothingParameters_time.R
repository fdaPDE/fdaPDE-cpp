checkSmoothingParameters_time<-function(locations = NULL, time_locations=NULL, observations, FEMbasis, time_mesh=NULL, lambdaS, lambdaT = 1, covariates = NULL, PDE_parameters=NULL, incidence_matrix = NULL, BC = NULL, FLAG_MASS = FALSE, FLAG_PARABOLIC = FALSE, IC = NULL, GCV = FALSE,GCVmethod = 2,nrealizations = 100, DOF=NULL, DOF_matrix=NULL, search, bary.locations=bary.locations, GCV.inflation.factor = 1, areal.data.avg = TRUE)
{
  #################### Parameter Check #########################

  if(is.null(FEMbasis))
    stop("FEMbasis required;  is NULL.")
  if(class(FEMbasis)!= "FEMbasis")
    stop("'FEMbasis' is not class 'FEMbasis'")

  if(class(FEMbasis$mesh)!='mesh.2D' & class(FEMbasis$mesh) != "mesh.2.5D" & class(FEMbasis$mesh) != "mesh.3D")
    stop('Unknown mesh class')

  if((class(FEMbasis$mesh) == "mesh.2.5D" || class(FEMbasis$mesh) == "mesh.3D") & !is.null(PDE_parameters) )
    stop('For mesh classes different from mesh.2D, anysotropic regularization is not yet implemented.
         Use Laplacian regularization instead')



  if(!is.null(locations))
  {
    if(any(is.na(locations)))
      stop("Missing values not admitted in 'locations'.")
  }

  if (is.null(time_locations) && is.null(time_mesh))
    stop("at least one between 'time_locations' and 'time_mesh' is required;  both NULL")

  if (any(is.na(time_mesh)))
    stop("Missing values not admitted in 'time_mesh'")

  if (is.null(observations))
    stop("observations required;  is NULL.")

  if (is.null(lambdaS))
    stop("lambda required;  is NULL.")

  if (is.null(lambdaT))
    stop("lambda required;  is NULL.")

  if (!is.null(locations) && !is.null(incidence_matrix))
    stop("Both 'locations' and 'incidence_matrix' are given. In case of pointwise data, set 'incidence_matrix to NULL. In case of areal data, set 'locations' to NULL.")

  if (any(incidence_matrix!=0 & incidence_matrix!=1))
    stop("Value different than 0 or 1 in 'incidence_matrix'.")

  if(!is.null(BC))
  {
    if (is.null(BC$BC_indices))
      stop("'BC_indices' required in BC;  is NULL.")
    if (is.null(BC$BC_values))
      stop("'BC_indices' required in BC;  is NULL.")
  }

  if (is.null(GCV))
    stop("GCV required;  is NULL.")
  if(!is.logical(GCV))
    stop("'GCV' is not logical")

  if (is.null(FLAG_MASS))
    stop("FLAG_MASS required;  is NULL.")
  if(!is.logical(FLAG_MASS))
    stop("'FLAG_MASS' is not logical")

  if (is.null(FLAG_PARABOLIC))
    stop("FLAG_PARABOLIC required;  is NULL.")
  if(!is.logical(FLAG_PARABOLIC))
    stop("'FLAG_PARABOLIC' is not logical")

  if(!is.null(PDE_parameters))
  {
    if (is.null(PDE_parameters$K))
      stop("'K' required in PDE_parameters;  is NULL.")
    if (is.null(PDE_parameters$b))
      stop("'b' required in PDE_parameters;  is NULL.")
    if (is.null(PDE_parameters$c))
      stop("'c' required in PDE_parameters;  is NULL.")
  }

  #Check the locations in 'bary.locations' and 'locations' are the same
 if(!is.null(bary.locations) & !is.null(locations))
 {
   flag=TRUE
   for (i in 1:nrow(locations)) {
     if (!(locations[i,1]==bary.locations$locations[i,1] & locations[i,2] == bary.locations$locations[i,2])) {
       flag = FALSE
       break
     }
   }

   if (flag == FALSE) {
     stop("Locations are not same as the one in barycenter information.")
   }
 }  # end of bary.locations
 
  space_varying=FALSE

  if(!is.null(PDE_parameters$u)){

    space_varying=TRUE

    message("Smoothing: anysotropic and non-stationary case")

    if(!is.function(PDE_parameters$K))
      stop("'K' in 'PDE_parameters' is not a function")
    if(!is.function(PDE_parameters$b))
      stop("'b' in 'PDE_parameters' is not a function")
    if(!is.function(PDE_parameters$c))
      stop("'c' in 'PDE_parameters' is not a function")
    if(!is.function(PDE_parameters$u))
      stop("'u' in 'PDE_parameters' is not a function")

  }
  else if(!is.null(PDE_parameters)){
    message("Smoothing: anysotropic and stationary case")
  }

  if(is.null(PDE_parameters))
    message("Smoothing: isotropic and stationary case")


  if (GCVmethod != 1 && GCVmethod != 2)
    stop("GCVmethod must be either Exact or Stochastic")

  if( !is.numeric(nrealizations) || nrealizations < 1)
    stop("nrealizations must be a positive integer")

  # check GCV.inflation.factor
  if (is.null(GCV.inflation.factor)){ 
    stop("'GCV.inflation.factor' required;  is NULL.")
  }else if( !is.numeric(GCV.inflation.factor) || GCV.inflation.factor < 0){
      stop("'GCV.inflation.factor' must be a real positive")
  }

   if (is.null(areal.data.avg))
    stop("'areal.data.avg' required;  is NULL.")
  if(!is.logical(areal.data.avg))
    stop("'areal.data.avg' is not logical")


  ans=space_varying

  ans
}

checkSmoothingParametersSize_time<-function(locations = NULL, time_locations=NULL, observations, FEMbasis, time_mesh=NULL, lambdaS, lambdaT = 1, covariates = NULL, PDE_parameters=NULL, incidence_matrix = NULL, BC = NULL, FLAG_MASS = FALSE, FLAG_PARABOLIC = FALSE, IC = NULL, GCV = FALSE, DOF=FALSE, DOF_matrix=NULL, space_varying, ndim, mydim)
{
  #################### Parameter Check #########################
  if(ncol(observations) < 1)
    stop("'observations' must contain at least one element")
  if(nrow(observations) < 1)
    stop("'observations' must contain at least one element")

  if(!is.null(time_locations)){
    if(ncol(time_locations) != 1)
      stop("'time_locations' must be a column vector")
    if (ncol(observations) != nrow(time_locations))
      stop("'observations' must be a #locations x #time_locations matrix")
  }else{
    if (ncol(observations) != nrow(time_mesh))
      stop("'observations' must be a #locations x #time_locations matrix")
  }

  if(!is.null(time_mesh))
    if(ncol(time_mesh) != 1)
      stop("'time_locations' must be a column vector")

  if(!is.null(time_mesh) && length(time_mesh)==1)
      stop("'time_mesh' must be of length bigger than 1. For only space problems use smooth.FEM")

  if(!is.null(time_locations) && length(time_locations)==1)
      stop("'time_locations' must be of length bigger than 1. For only space problems use smooth.FEM")


  if(is.null(locations) && is.null(incidence_matrix))
  {
    if(!is.null(time_locations))
      if(ifelse(class(FEMbasis$mesh) == "mesh.2D", nrow(FEMbasis$mesh$nodes),FEMbasis$mesh$nnodes) != nrow(observations) || ncol(observations) != nrow(time_locations))
        stop("'locations' and 'observations' have incompatible size;")

    if(is.null(time_locations))
    {
      if(FLAG_PARABOLIC && !is.null(IC))
      {
        if(ifelse(class(FEMbasis$mesh) == "mesh.2D", nrow(FEMbasis$mesh$nodes),FEMbasis$mesh$nnodes) != nrow(observations) || (nrow(time_mesh)-1) != ncol(observations))
          stop("'locations' and 'observations' have incompatible size;")
      }
      else
      {
        if(ifelse(class(FEMbasis$mesh) == "mesh.2D", nrow(FEMbasis$mesh$nodes),FEMbasis$mesh$nnodes) != nrow(observations) || nrow(time_mesh) != ncol(observations))
          stop("'locations' and 'observations' have incompatible size;")
      }
    }
  }

  if(!is.null(locations))
  {
    if(ncol(locations) != ndim)
      stop("'locations' must be a ndim-columns matrix;")

    if(!is.null(time_locations) && is.null(incidence_matrix))
    {
      if(nrow(locations) != nrow(observations))
        stop("'locations' and 'observations' have incompatible size;")

      if(nrow(time_locations) != ncol(observations))
        stop("'time_locations' and 'observations' have incompatible size;")
    }
    if(is.null(time_locations) && is.null(incidence_matrix))
    {
      if(FLAG_PARABOLIC)
      {
        if(nrow(locations) != nrow(observations))
          stop("'locations' and 'observations' have incompatible size;")
        if(!is.null(IC))
          if((nrow(time_mesh)-1) != ncol(observations))
            stop("'time_mesh' and 'observations' have incompatible size;")
      }
      if(!FLAG_PARABOLIC)
      {
        if(nrow(locations) != nrow(observations))
          stop("'locations' and 'observations' have incompatible size;")

        if(nrow(time_mesh) != ncol(observations))
          stop("'time_mesh' and 'observations' have incompatible size;")
      }
    }
    if(dim(locations)[1]==dim(FEMbasis$mesh$nodes)[1] & dim(locations)[2]==dim(FEMbasis$mesh$nodes)[2])
      warning("The locations matrix has the same dimensions as the mesh nodes. If the locations you are using are the
              mesh nodes, set locations=NULL instead")
  }


  if(ncol(lambdaS) != 1)
    stop("'lambdaS' must be a column vector")
  if(nrow(lambdaS) < 1)
    stop("'lambdaS' must contain at least one element")

  if(ncol(lambdaT) != 1)
    stop("'lambdaT' must be a column vector")
  if(nrow(lambdaT) < 1)
    stop("'lambdaT' must contain at least one element")

  if(nrow(lambdaS)+nrow(lambdaT)>2 && !GCV)
    warning("different values of lambdaS or lambdaT have been passed but GCV=FALSE,
              if you want to compute the GCV please set GCV=TRUE")

  if(!is.null(covariates))
  {
    if(nrow(covariates) != nrow(observations)*ncol(observations))
      stop("'covariates' and 'observations' have incompatible size;")
  }

  if (!is.null(incidence_matrix))
  {
    if (nrow(incidence_matrix) != nrow(observations))
      stop("'incidence_matrix' and 'observations' have incompatible size;")
    if (class(FEMbasis$mesh) == 'mesh.2D' && ncol(incidence_matrix) != nrow(FEMbasis$mesh$triangles))
      stop("'incidence_matrix' must be a ntriangles-columns matrix;")
    else if (class(FEMbasis$mesh) == 'mesh.2.5D' && ncol(incidence_matrix) != FEMbasis$mesh$ntriangles)
      stop("'incidence_matrix' must be a ntriangles-columns matrix;")
    else if (class(FEMbasis$mesh) == 'mesh.3D' && ncol(incidence_matrix) != FEMbasis$mesh$ntetrahedrons)
      stop("'incidence_matrix' must be a ntetrahedrons-columns matrix;")
  }

  if(!is.null(BC))
  {
    if(ncol(BC$BC_indices) != 1)
      stop("'BC_indices' must be a column vector")
    if(ncol(BC$BC_values) != 1)
      stop("'BC_values' must be a column vector")
    if(nrow(BC$BC_indices)>ifelse(class(FEMbasis$mesh)=='mesh.2D',nrow(FEMbasis$mesh$nodes),FEMbasis$mesh$nnodes))
      stop("'BC_indices' longer than the mesh")
    if(nrow(BC$BC_indices) != nrow(BC$BC_values))
      stop("'BC_indices' and 'BC_values' have incompatible size;")

    N=ifelse(class(FEMbasis$mesh) == "mesh.2D",nrow(FEMbasis$mesh$nodes),FEMbasis$mesh$nnodes)
    #M=ifelse(FLAG_PARABOLIC==TRUE,length(time_mesh)-1,length(time_mesh)+2)
    if(any(BC$BC_indices)<0)
      stop("'BC_indices' elements must be non negative")
    if(any(BC$BC_indices)>ifelse(class(FEMbasis$mesh)=='mesh.2D',nrow(FEMbasis$mesh$nodes),FEMbasis$mesh$nnodes))
      stop("At least one index in 'BC_indices' larger then the number of 'nodes' in the mesh")
  }

  if (FLAG_PARABOLIC==TRUE && is.null(IC))
  {
    if (!is.null(time_mesh))
    {
      if(ncol(observations)!=nrow(time_mesh))
        stop("IC is required for parabolic smoothing")
      else
        message("IC is required for parabolic smoothing, will be estimated from the first column of 'observations'")
    }
    else
      message("IC is required for parabolic smoothing, will be estimated from the first column of 'observations'")
  }
  if(!is.null(IC))
  {
    if(ncol(IC) != 1)
      stop("'IC' must be a column vector")
    if(length(IC) != FEMbasis$nbasis)
      stop("'IC' must be a vector of length 'FEMbasis$nbasis'")
  }

  if(!is.null(PDE_parameters) & space_varying==FALSE)
  {
    if(!all.equal(dim(PDE_parameters$K), c(2,2)))
      stop("'K' in 'PDE_parameters must be a 2x2 matrix")
    if(!all.equal(dim(PDE_parameters$b), c(2,1)))
      stop("'b' in 'PDE_parameters must be a column vector of size 2")
    if(!all.equal(dim(PDE_parameters$c), c(1,1)))
      stop("'c' in 'PDE_parameters must be a double")
  }

  if(!is.null(PDE_parameters) & space_varying==TRUE)
  {

    n_test_points = min(nrow(FEMbasis$mesh$nodes), 5)
    test_points = FEMbasis$mesh$nodes[1:n_test_points, ]

    try_K_func = PDE_parameters$K(test_points)
    try_b_func = PDE_parameters$b(test_points)
    try_c_func = PDE_parameters$c(test_points)
    try_u_func = PDE_parameters$u(test_points)

    if(!is.numeric(try_K_func))
      stop("Test on function 'K' in 'PDE_parameters' not passed; output is not numeric")
    if(!all.equal(dim(try_K_func), c(2,2,n_test_points)) )
      stop("Test on function 'K' in 'PDE_parameters' not passed; wrong size of the output")

    if(!is.numeric(try_b_func))
      stop("Test on function 'b' in 'PDE_parameters' not passed; output is not numeric")
    if(!all.equal(dim(try_b_func), c(2,n_test_points)))
      stop("Test on function 'b' in 'PDE_parameters' not passed; wrong size of the output")

    if(!is.numeric(try_c_func))
      stop("Test on function 'c' in 'PDE_parameters' not passed; output is not numeric")
    if(length(try_c_func) != n_test_points)
      stop("Test on function 'c' in 'PDE_parameters' not passed; wrong size of the output")

    if(!is.numeric(try_u_func))
      stop("Test on function 'u' in 'PDE_parameters' not passed; output is not numeric")
    if(length(try_u_func) != n_test_points)
      stop("Test on function 'u' in 'PDE_parameters' not passed; wrong size of the output")
  }

  if(!is.null(DOF_matrix))
  {
    if(GCV==FALSE)
      warning("GCV=FALSE. DOF_matrix is passed but GCV is not computed")
    if(nrow(DOF_matrix)!=length(lambdaS))
      stop("The number of rows of DOF_matrix is different from the number of lambdaS")
    if(ncol(DOF_matrix)!=length(lambdaT))
      stop("The number of columns of DOF_matrix is different from the number of lambdaT")
  }

}
