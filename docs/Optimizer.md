The optimization module provides a set of general routines for optimizing a generic [ScalarField](ScalarField.md) \( f : \mathbb{R}^N \rightarrow \mathbb{R} \).

# Optimizer

> :fontawesome-solid-file-code: core/OPT/Optimizer.h &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; :fontawesome-solid-sitemap: Extends: -

Base interface for the whole optimization module. 

``` c++
template <unsigned int N> class Optimizer { ... };
```

!!! info
    The template parameter N denotes the dimension of the space where to search for the optimum point. This information should be known at compile time in order to define the data structure required during the optimization.

!!! caution "Developer's advice"
    Any class meant to work as an optimizer **must** derive this class or one of its direct child classes.
	
## Methods offered by the interface

``` c++
virtual std::pair<SVector<N>, double> findMinimum();
```

!!! quote ""
	
    The main entry point for the optimization routine. A call to this method should produce as output a `std::pair<SVector<N>, double>` where the first component is the point where the minimum of the objective is reached while the second component is the actual minimum value found. 
	
	The algorithm used for producing the result is implementation dependent.

# IterativeOptimizer <a name="IterativeOptimizer"></a> 

> :fontawesome-solid-file-code: core/OPT/IterativeOptimizer.h &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; :fontawesome-solid-sitemap: Extends: [Optimizer](Optimizer.md)

Template abstract class representing a general iterative optimizer. An iterative optimizer is any optimization routine which can be schematized as follow:

``` c++ linenums="1"
let x[0] the starting point
let k = 0 
while( some stopping condition is not met ){
	// update the current solution
	x[k+1] = x[k] + a[k]*d[k];
	k++;
}

return x[k];
```

Even if algorithms can drammatically differ in the computation of the update step this general footprint is still mantained. For this reason a rich family of optimization algorithms fall under this class.

The `IterativeOptimizer` interface should be used to inform the user that an iterative optimization procedure is implemented under the hood.

``` c++
template <unsigned int N> class IterativeOptimizer : public Optimizer<N> { ... };
```

!!! info
    The `IterativeOptimizer` class does not implement the general schema of an iterative optimizer, whose implementation is left to the derived classes. Forcing any deriving class to follow a too rigid schema could introduce useless complications in the implementation of the optimization procedure.
	
	This class introduce anyway the possibility to control in some way the flow of execution of an iterative optimizer. For example it would be possible to monitor during the execution of the algorithm itself some quantities of interest for the particular problem at hand and force the stop of the procedure on the base of some custom stopping criterion bypassing the default one.
	
    The way this mechanism is reached is by extending any concrete implementation of the `IterativeOptimizer` class and overloading the wanted methods exposed by `IterativeOptimizer` itself.

!!! warning "Developer's advice"
    Is up to the derived classes of `IterativeOptimizer` to implement properly the stated mechanism. 
	
	When you want to give the possibility to execute a (possibly) custom action of a (possibly) deriving class you should insert i.e. `this->init()` for executing a custom initialization.
	See the semantic of `IterativeOptimizer`' methods to see what kind of customizations can be added and decide the proper place where to execute the action.


## Methods

``` c++
virtual void init();
```

!!! quote ""
	
	This method should be called once before entering the iterative loop.
	
``` c++
virtual void preStep();
```

!!! quote ""
	
	This method should be inserted in the iterative loop, possibly as first loop instruction. In any case insert it before the update step.

``` c++
virtual void postStep();
```

!!! quote ""
	
    This method should be inserted in the iterative loop. Insert it after the update step. A good place is immediately after the error update.


``` c++
virtual bool stopCondition();
```

!!! quote ""
	
	This method should be called in a `while` of `for` statement togheter with the default termination criteria. The good practice is to write something like
	
	```c++
	while(!default_stopping_condition() && !this->stopCondition()) { ... }
	```

	It must return false if the custom stopping condition is not met.

``` c++
virtual void finalize();
```

!!! quote ""
	
	This method should be called once outside the iterative loop, possibly right before the return statement.

As an example consider the implementation of the optimization routine for the [GradientDescent](GradientDescentOptimizer.md)


``` c++ linenums="1" hl_lines="4 13 14 23 29"
// gradient descent optimization routine
template <unsigned int N>
std::pair<SVector<N>, double> GradientDescentOptimizer<N>::findMinimum(){
  this->init();             // execute custom action
  
  // algorithm initialization
  x_old = x0;
  unsigned int numIteration = 0;
  
  // standard termination criteria based on l^2 norm of the gradient
  error = objective.derive()(x_old).squaredNorm();
  
  while (numIteration < maxIteration && error > tolerance && !this->stopCondition()){
    this->preStep();        // execute custom action
    
    // compute exact gradient
    gradientExact = objective.derive()(x_old);
    // update step    
    x_new = x_old - step*gradientExact;
    // error update: standard termination criteria based on l^2 norm of the gradient
    error = gradientExact.squaredNorm();

    this->postStep();       // execute custom action
    // prepare next iteration    
    x_old = x_new;    
    numIteration++;
  }
  
  this->finalize();         // execute custom action
  return std::pair<SVector<N>, double>(x_old, objective(x_old));
}
```

!!! tip
    The `IterativeOptimizer` class offers a `std::unordered_map<std::string, std::list<double>> controllerData` which can be used **by a customization** to store or record values needed to perform custom actions. Use the `init()` method to initialize any field you might require before the actual optimization starts.
