checkGAMParameters<-function(observations, max.steps.FPIRLS, mu0, scale.param, threshold.FPIRLS, family)
{
  observations.len = length(observations)
	  #################### Parameter Check #########################
	# Check max.steps.FPIRLS 
	if(!all.equal(max.steps.FPIRLS, as.integer(max.steps.FPIRLS)) || max.steps.FPIRLS <= 0 )
		stop("'max.steps.FPIRLS' must be a positive integer.")
	
  #check observations
  if( family == "binomial"){
      if(length(c(which(observations==0),which(observations==1))) != observations.len )
        stop("'observations' must be composed by 0 and 1 values for your distribution.")
  }
  if( length(c(which(observations==0),which(observations==1))) == observations.len && family != "binomial"){
    stop("'observations' have only 0 and 1 values. Your distribution must be 'bernulli'.")
  }
  if(family == "gamma"){
      if(any(observations<=0))
          stop("'observations' must be composed by real positive number for your distribution.")
  }
  if(family == "exponential"){
      if(any(observations<0))
          stop("'observations' must be composed by real positive number for your distribution.")
  }
  if(family == "poisson"){
      if(any(observations<0))
          stop("'observations' must be composed by real positive number for your distribution.")
      if(any(observations != floor(observations)))
          stop("'observations' must be composed by natural positive number for your distribution.")
  }
  if( all(observations == floor(observations)) && sum(family==c("binomial","poisson")) == 0 ){
      stop("'observations' is composed by natural positive number. Select a discrete distribution.")
  }



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

  		family_positive = c("exponential", "gamma", "poisson")
  		if(sum(family==family_positive)==1){
  			if(any(mu0<0))
  				stop("mu0 must be composed by real positive number for your distribution")
  		}
  	}
	
	# check scale.param
	if(!is.null(scale.param)){
		if( !is.numeric(scale.param) || scale.param <= 0 )
    		stop("The dispersion parameter of your distribution ('scale.param') must be a positive real number")
	}	

	# check threshold.FPIRLS
	if (is.null(threshold.FPIRLS)){ 
		stop("'threshold.FPIRLS' required;  is NULL.")
	}else if( !is.numeric(threshold.FPIRLS) || threshold.FPIRLS <= 0){
    	stop("'threshold.FPIRLS' must be a real positive")
	}

  # observations check
  if(any(is.na(observations)))
    stop("Missing values not admitted in 'observations' in FPIRLS method")

}