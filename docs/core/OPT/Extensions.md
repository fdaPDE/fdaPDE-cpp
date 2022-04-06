# Extensions

This page contains more detailed informations about available extensions and how to extending an already existing optimizer by writing a new extension.

## Available Extensions

> :fontawesome-solid-file-code: core/OPT/extensions/... &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; :fontawesome-solid-sitemap: Namespace: fdaPDE::core::OPT

This is the list of currently available extensions:

* **BacktrackingAdaptiveStep**

	> :fontawesome-solid-file-code: core/OPT/extensions/BacktrackingAdaptiveStep.h &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; :fontawesome-solid-sitemap: Extends: -


	This extension is for enabling the backtracking line search method to adapt the step size of an iterative optimizer during execution. 
	
	For a function \( f : \mathbb{R}^N \rightarrow \mathbb{R} \) whose gradient is known, given an iterative optimization scheme:
	$$ x_{k+1} = x_k + \lambda_k p_k $$
	Let \(\alpha, \beta \) and \( \gamma\) three fixed parameters. Given \( x_k \) the current solution of the iterative method and \( p_k \) the computed update direction, the backtraking method computes \( \lambda_k \) before computing the solution update \( x_{k+1} \) by producing a sequence \( \alpha_k \) using the following scheme:
	
	* set \( \alpha_0 = \alpha \) and \( i = 1 \)
	* let \( \alpha_i = \beta \alpha_{i-1} \)
	* if \( f(x_k) - f( x - \alpha_i \nabla f(x_k) ) + \gamma \alpha_i \langle \nabla f(x_k), p_k \rangle < 0 \) stop and set \( \lambda_k = \alpha_i \), otherwise increment \(i\) by one and repeat from the previous point.

	!!! info
		The extension injects the computation in the `InitIteration` checkpoint of the iterative procedure. 
	
``` c++
BacktrackingAdaptiveStep(double alpha_, double beta_, double gamma_);
```

!!! quote ""
	Constructor.

    | Args                                                                                          | Description                                                 |
    |------------------------------------------------------------------------------------------|-------------------------------------------------------------|
    | `double alpha_` | The \( \alpha \) value in the backtracking method (i.e. the starting point for the line search algorithm). |
    | `double beta_` | The \( \beta \) value in the backtracking method. Should be a value in between 0 and 1. |
    | `double gamma_` | The \( \gamma \) value in the backtracking method. Should be a value in between 0 and 1. |
	
	!!! note
		The class provides also a default constructor which sets the three parameters to some default value. Parameter selection is really probelm specific and default initialization should be avoided.

* **Summarize**

	> :fontawesome-solid-file-code: core/OPT/extensions/Summarize.h &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; :fontawesome-solid-sitemap: Extends: -
	
	An extension to report a small result of the optimization procedure.

	!!! info
		The extension injects the computation in the `EndOptimization` checkpoint of the iterative procedure. 

## Write a new extension

!!! info
	When the code finds a class intended to be an extension passed as argument to the optimizer's entry point it looks inside its definition to search for a well defined set of methods which will be then called in a well defined position. All happens at compile-time with no run-time overhead.
	
You can refer to the following snippet as basic footprint for future extensions:

``` c++ linenums="1"
#ifndef __MYEXT_H__
#define __MYEXT_H__

...

class MyExt {
private:
	... some internal state of the extension ...

public: 
	MyExt() = default;
	
	template <typename Optimizer, typename Objective>
    bool initOptimization(Optimizer& opt, Objective& obj){
	   /* this will be called at the start of the optimization */
       return false;
    }
	
	template <typename Optimizer, typename Objective>
    bool initIteration(Optimizer& opt, Objective& obj){
	   /* this will be called at the start of each iteration */
       return false;
    }
	
	template <typename Optimizer, typename Objective>
    bool endIteration(Optimizer& opt, Objective& obj){
	   /* this will be called at the end of each iteration */
       return false;
    }

    template <typename Optimizer, typename Objective>
    bool endOptimization(Optimizer& opt, Objective& obj){
	   /* this will be called at the end of the optimization */
       return false;
    }
	
	... other methods ...
};

#endif // __MYEXT_H__

```

The return value of each method is used to stop the optimization at any point (i.e. you can make one method to return `true` if you want to stop the optimization procedure). 

!!! note "Developer's advice"
	In case the extension has effect in just one single phase of the optimization scheme you can avoid (and this is recommended) to define all the other checkpoint methods.


## Make an optimizer open to extensions

An optimizer which supports the extension mechanism must have a variadic templated signature for its entry point as follow:

``` c++
template <typename... Args>
std::pair<SVector<N>, double> findMinimum(const ScalarField<N>& objective, const SVector<N>& x0,
										  const Args&... args);
```

Another condition to make the implementation in the position to process possible extensions passed as input is to insert in the implementation of `findMinimum()` a well defined set of static members (provided by `core/OPT/extensions/Extension.h`) which will activate the extension machinery.

You can refer to the following snippet as basic footprint for an extensible optimizer:

``` c++ linenums="1"
template <unsigned int N>
template <typename... Args>
std::pair<SVector<N>, double> findMinimum(const ScalarField<N>& objective, const SVector<N>& x0,
										  const Args&... args){
	bool customStop = false;
	customStop |= Extension::executeInitOptimization(*this, objective, args...);
	
	... initialization ...
	
	while(... && !customStop){
		customStop |= Extension::executeInitIteration(*this, objective, args...);
		
		... optimization logic ...
	
		customStop |= Extension::executeEndIteration (*this, objective, args...);
	}
	
	Extension::executeEndOptimization(*this, objective, args...);
	return ...;
}
```

!!! note
	The position in which checkpoints are inserted is the one commonly adopted in the whole module. You can change the position according to your needs even if is recommended to mantain this scheme in order to have a consistent logic across different optimizers.
