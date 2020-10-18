#ifndef __OPTIMIZATON_METHODS_FACTORY_H__
#define __OPTIMIZATON_METHODS_FACTORY_H__

// HEADERS
#include "Newton.h"
#include "../../Global_Utilities/Include/Make_Unique.h"
#include <memory>

//! A Factory class: A class for the choice of the cross-validation method to use for the selection of the parameter lambda.
/* \tparam Function the function used to create the shared pointer to the optimization method
 * \tparam Tuple type of gradient emplyed
 * \tparam Hessian type of hessian employed
 * \tparam EvaluationType type of the evaluated function
*/
template<typename Tuple, typename Hessian, typename EvaluationType>
class Opt_method_factory
{
	public:
        	//! A method that takes as parameter a string and builds a pointer to the right object
		/*!
		 \param validation a string code to decide which pointer to create
		 \param F a function from which to create the pointer
		 \return a pointer to the validated optimization method
		*/
        	static std::unique_ptr<Opt_methods<Tuple,Hessian,EvaluationType>> create_Opt_method(const std::string & validation, Function_Wrapper<Tuple, Real, Tuple, Hessian,EvaluationType> & F)
                {
                	if(validation=="newton")
                                return make_unique<Newton_ex<Tuple, Hessian, EvaluationType>>(F);
                	if(validation=="newton_fd")
                                return make_unique<Newton_fd<Tuple, Hessian, EvaluationType>>(F);
			else // default is fd
			{
				Rprintf("Method not found, using Newton_fd");
				return make_unique<Newton_fd<Tuple, Hessian, EvaluationType>>(F);
			}
        	}
};

#endif
