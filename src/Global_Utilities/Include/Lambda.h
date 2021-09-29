#ifndef __LAMBDA_H__
#define __LAMBDA_H__

#include "../../FdaPDE.h"

namespace lambda
{
	template<UInt size>
	using type = typename std::conditional<size==1, Real, VectorXr>::type;
	
	type<2> make_pair(Real lambdaS, Real lambdaT);
	type<2> make_pair(lambda::type<2> lambda, Real lambdaT); //degenerate case useful in parabolic case

	template<UInt size>
	typename std::enable_if<size==1, type<1>>::type
	init(Real value) {return value;}

	template<UInt size>
	typename std::enable_if<size==2, type<2>>::type
	init(Real value) {return make_pair(value, value);}
}

#endif
