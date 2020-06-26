checkSmoothingParametersFPCA<-function(locations = NULL, datamatrix, FEMbasis, incidence_matrix = NULL, lambda, nPC, validation, NFolds, GCVmethod = 2, nrealizations = 100, search, bary.locations=bary.locations)
{
  #################### Parameter Check #########################
  if(!is.null(locations))
  {
    if(any(is.na(locations)))
      stop("Missing values not admitted in 'locations'.")
    if(any(is.na(datamatrix)))
      stop("Missing values not admitted in 'datamatrix' when 'locations' are specified.")
  # if (search == 1) { #use Naive search
    #   print('This is Naive Search')
    # } else if (search == 2)  { #use Tree search (default)
    #   print('This is Tree Search')
    # }

  } # end of locations
  if (is.null(datamatrix))
    stop("observations required;  is NULL.")
  if (is.null(FEMbasis))
    stop("FEMbasis required;  is NULL.")
  if(class(FEMbasis)!= "FEMbasis")
    stop("'FEMbasis' is not class 'FEMbasis'")

  if (!is.null(locations) && !is.null(incidence_matrix))
    stop("Both 'locations' and 'incidence_matrix' are given. In case of pointwise data, set 'incidence_matrix to NULL. In case of areal data, set 'locations' to NULL.")

  if (any(incidence_matrix!=0 & incidence_matrix!=1))
    stop("Value different than 0 or 1 in 'incidence_matrix'.")


  if (is.null(lambda))
    stop("lambda required;  is NULL.")
  if(is.null(nPC))
    stop("nPC required; is NULL.")
  if(!is.null(validation)){
    if(validation!="GCV" && validation!="KFold")
   	stop("'validation' needs to be 'GCV' or 'KFold'")
    if(validation=="KFold" && is.null(NFolds))
   	stop("NFolds is required if 'validation' is 'KFold'")
   }
   if (GCVmethod != 1 && GCVmethod != 2)
    stop("GCVmethod must be either 1(exact calculation) or 2(stochastic estimation)")

  if( !is.numeric(nrealizations) || nrealizations < 1)
    stop("nrealizations must be a positive integer")

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

}

checkSmoothingParametersSizeFPCA<-function(locations = NULL, datamatrix, FEMbasis, incidence_matrix, lambda, ndim, mydim, validation, NFolds)
{
  #################### Parameter Check #########################
  if(nrow(datamatrix) < 1)
    stop("'datamatrix' must contain at least one element")
  if(is.null(locations))
  {
    if(class(FEMbasis$mesh) == "mesh.2D"){
    	if(ncol(datamatrix) > nrow(FEMbasis$mesh$nodes))
     	 stop("Size of 'datamatrix' is larger then the size of 'nodes' in the mesh")
    }else if(class(FEMbasis$mesh) == "mesh.2.5D" || class(FEMbasis$mesh) == "mesh.3D"){
    	if(ncol(datamatrix) > nrow(FEMbasis$mesh$nodes))
     	 stop("Size of 'datamatrix' is larger then the size of 'nodes' in the mesh")
    }
  }
  if(!is.null(locations))
  {
    if(ncol(locations) != ndim)
      stop("'locations' and the mesh points have incompatible size;")
    if(nrow(locations) != ncol(datamatrix))
      stop("'locations' and 'datamatrix' have incompatible size;")
    if(dim(locations)[1]==dim(FEMbasis$mesh$nodes)[1] & dim(locations)[2]==dim(FEMbasis$mesh$nodes)[2])
      warning("The locations matrix has the same dimensions as the mesh nodes. If the locations you are using are the
              mesh nodes, set locations=NULL instead")

  }
  if (!is.null(incidence_matrix))
  {
    if (nrow(incidence_matrix) != ncol(datamatrix))
      stop("'incidence_matrix' and 'datamatrix' have incompatible size;")
    if (class(FEMbasis$mesh) == 'mesh.2D' && ncol(incidence_matrix) != nrow(FEMbasis$mesh$triangles))
      stop("'incidence_matrix' must be a ntriangles-columns matrix;")
    else if (class(FEMbasis$mesh) == 'mesh.2.5D' && ncol(incidence_matrix) != nrow(FEMbasis$mesh$triangles))
      stop("'incidence_matrix' must be a ntriangles-columns matrix;")
    else if (class(FEMbasis$mesh) == 'mesh.3D' && ncol(incidence_matrix) != nrow(FEMbasis$mesh$tetrahedrons))
      stop("'incidence_matrix' must be a ntetrahedrons-columns matrix;")
  }
  if(ncol(lambda) != 1)
    stop("'lambda' must be a column vector")
  if(nrow(lambda) < 1)
    stop("'lambda' must contain at least one element")
  if(nrow(lambda)>1 && is.null(validation))
    stop("If 'lambda' contains more than one element, 'validation' needs to be specified as 'GCV' or 'KFold'")
}
