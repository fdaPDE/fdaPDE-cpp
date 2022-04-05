template <unsigned int N>
template <typename... Args>
std::pair<SVector<N>, double> findMinimum(const ScalarField<N>& objective_,
					  const Args&... args){
  // can be set true by some extension to cause a foced stop at any point of the execution
  bool customStop = false; 
  customStop |= Extension::executeInitOptimization(*this, objective, args...);

  // number of iterations for each dimension
  array<int,N> numberIteration;
  // prepare data structure to execute the N-dimensional scan
  SVector<N> gridPoint;
  
  // init helper data structures
  for(size_t i = 0; i < N; ++i){
    numberIteration[i] = std::ceil((domainLimits[i].second - domainLimits[i].first)/steps[i]) + 1;
    gridPoint[i] = domainLimits[i].first;
  }

  // initial minimum
  double minimum = objective(gridPoint);
  SVector<N> minPoint = gridPoint;
  
  // counter to keep track of iterations executed
  // this is a vector of N elements started at [0, 0, ..., 0]. The algorithms always increments the vector at
  // position 0 and when it reaches the maximum number of iterations for that dimension it is restarted from 0 incrementing
  // recursively the next positions following the same schema. Loop stops when currentIteration = numberIteration
  array<int, N> currentIteration = {};
  
  // search for minimum
  while(currentIteration[N-1] < numberIteration[N-1] && !customStop){ // loop until the last dimension has not been scan entirely
    customStop |= Extension::executeInitIteration(*this, objective, args...);

    // evaluate function at current point
    gridPoint[0] += steps[0];
    double x = objective(gridPoint);

    // update minimum
    if(x < minimum){
      minimum = x;
      minPoint = gridPoint;
    }

    // increment counter at first dimension
    currentIteration[0]++; 
    
    // trigger increment in next elements of the counter if necessary
    for(int i = 0; currentIteration[i] == numberIteration[i] && i < N-1; ++i){ // dimension i has been spanned
      currentIteration[i] = 0;               // reset dimension i
      gridPoint[i] = domainLimits[i].first;
      
      currentIteration[i+1]++;               // increment counter for next dimension
      gridPoint[i+1] += steps[i+1];
    }
    customStop |= Extension::executeEndIteration(*this, objective, args...);
  }

  customStop |= Extension::executeEndOptimization(*this, objective, args...);
  return std::pair<SVector<N>,double>(minPoint, minimum);
}
