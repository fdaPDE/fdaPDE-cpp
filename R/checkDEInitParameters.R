checkParametersDE_init <- function(data, FEMbasis, lambda, heatStep, heatIter, init, search) 
{
  #################### Parameter Check #########################
  if (is.null(data)) 
    stop("'data' required;  is NULL.")
  else{
    if(any(is.na(data)))
      stop("Missing values not admitted in 'data'.")
  }
  
  if (is.null(FEMbasis)) 
    stop("'FEMbasis' required;  is NULL.")
  if(class(FEMbasis)!= "FEMbasis")
    stop("'FEMbasis' is not class 'FEMbasis'.")
  
  if(init=="Heat"){
    if (is.null(lambda))  
      stop("'lambda' required if init='Heat'; is NULL.")
    else{
      for(i in 1:length(lambda)){
        if(lambda[i]<=0)
          stop("'lambda' has to have positive members.")
      }
    }
  }
  
  if(!is.numeric(search))
    stop("'search' needs to be an integer.")
  
  if(search != 1 && search != 2)
    stop("'search' needs to be an integer either equal to 1 or equal to 2.")

  if (is.null(init)) 
    stop("'init' is required;  is NULL.")
  else{
    if(init!="Heat" && init!="CV")
      stop("'init' needs to be either 'Heat' or 'CV'.")
  }
  
  if(init=="CV" && length(lambda)>1)
    stop("The initialization procedure via cross-validation is only for a given lamnda.")
}


checkParametersSizeDE_init <- function(data, ndim) 
{
  if(nrow(data) < 1)
    stop("'data' must contain at least one element.")
  if(ncol(data) != ndim)                                       
    stop("'data' and the mesh points have incompatible size.")
}
