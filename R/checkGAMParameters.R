checkGAMParameters<-function(max.steps, mu0, observations.len, scale.param, tune, threshold, fam)
{
	  #################### Parameter Check #########################
	# Check max.steps 
	if(!all.equal(max.steps, as.integer(max.steps)) || max.steps <= 0 )
		stop("max.steps must be a positive integer.")
	
	# Check mu0 
	if(!is.null(mu0))
  	{
    	if(any(is.na(mu0)))
      		stop("Missing values not admitted in 'mu0'.")

      	if(ncol(mu0) != 1)
    		stop("'mu0' must be a column vector")
  		if(nrow(mu0) < 1)
  			stop("'mu0' must contain at least one element")
  		
  		if(length(mu0) != observations.len )
    		stop(" 'mu0' and 'observations' must have equal length")

  		fam_positive = c("exponential", "gamma", "poisson", "cloglog")
  		if(sum(fam==fam_positive)==1){
  			if(any(mu0)<0)
  				stop("mu0 must be composed by real positive number for your distribution")
  		}
  	}
	
	# check scale.param
	if(!is.null(scale.param)){
		if( !is.numeric(scale.param) || scale.param <= 0 )
    		stop("The dispersion parameter of your distribution ('scale.param') must be a positive real number")
	}

	# check tune
	if (is.null(tune)){ 
		stop("'tune' required;  is NULL.")
	}else if( !is.numeric(tune) || tune < 0){
    	stop("'tune' must be a real positive")
	}
	

	# check threshold
	if (is.null(threshold)){ 
		stop("'threshold' required;  is NULL.")
	}else if( !is.numeric(threshold) || threshold <= 0){
    	stop("'threshold' must be a real positive")
	}


}