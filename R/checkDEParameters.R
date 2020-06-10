checkParametersDE <- function(data, FEMbasis, lambda, step_method, direction_method, preprocess_method, tol1, tol2, nfolds, nsimulations, heatStep, heatIter, search) 
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
  
  if (is.null(lambda))  
    stop("'lambda' required;  is NULL.")
  else{
    for(i in 1:length(lambda)){
      if(lambda[i]<=0)
        stop("'lambda' has to have positive members.")
    }
  }
  
  if (is.null(step_method)) 
    stop("'step_method' is required;  is NULL.")
  else{
    if(step_method!="Fixed_Step" && step_method!="Backtracking_Method" && step_method!="Wolfe_Method")
      stop("'step_method' needs to be either 'Fixed_Step' or 'Backtarcking_Method' or 'Wolfe_Method'.")
  }
  
  if (is.null(direction_method)) 
    stop("'direction_method' is required;  is NULL.")
  else{
    if(direction_method!="Gradient" && direction_method!="BFGS")
      stop("'direction_method' needs to be either 'Gradient' or 'BFGS'.")
  }

  if(length(lambda)>1 && preprocess_method!="RightCV" && preprocess_method!="SimplifiedCV")
    stop("'preprocess_method' needs to be either 'RightCV' or 'SimplifiedCV' if there are more than one smoothing parameters 'lambda'.")
  
  if(length(lambda)==1 && preprocess_method!="NoCrossValidation")
    stop("'preprocess_method' needs to be 'NoCrossValidation' if there is only one smoothing parameter 'lambda'.")
  
  if(preprocess_method=="SimplifiedCV" && length(lambda)!=nfolds)
    stop("'SimplifiedCV' requires the number of lambdas equal to the number of folds.")
  
  if(tol1 < 0 || tol2 < 0)
    stop("Tolerances 'tol1' and 'tol2' needs to be non negative numbers")
  
  if(length(lambda) > 1 && (!is.numeric(nfolds) || floor(nfolds)<=1))
    stop("'nfolds' needs to be an integer greater or equal than two.")

  if(!is.numeric(nsimulations) || nsimulations<1)
    stop("'nrealizations' needs to be a positive integer.")
  
  if(!is.numeric(heatStep) || heatStep<0 || heatStep>1)
    stop("'heatStep' needs to be a positive real number not greater than 1.")
  
  if(!is.numeric(heatIter) || heatIter<1)
    stop("'heatIter' needs to be a positive integer.")
  
  if(!is.numeric(search))
    stop("'search' needs to be an integer.")
  
  if(search != 1 && search != 2)
    stop("'search' needs to be an integer either equal to 1 or equal to 2.")
  
}


checkParametersSizeDE <- function(data, FEMbasis, ndim, fvec, preprocess_method, nfolds) 
{
  if(nrow(data) < 1)
    stop("'data' must contain at least one element.")
  if(preprocess_method!="NoCrossValidation" && nrow(data) < floor(nfolds))
    stop("The number of folds needs to be less than the number of data.")
  if(ncol(data) != ndim)                                       
    stop("'data' and the mesh points have incompatible size.")
  
  if(!is.null(fvec)){
    if(length(fvec) != nrow(FEMbasis$mesh$nodes))
      stop("The length of fvec has to be equal to the number of mesh nodes")
  }
 
   if(preprocess_method!="NoCrossValidation")
     if (nrow(data)*(floor(nfolds)-1)/floor(nfolds) < 30)
       stop("The training set needs to have at least 30 data: increase the number of folds.")
  
}
