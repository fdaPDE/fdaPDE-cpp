#ifndef __LAMBDA_H__
#define __LAMBDA_H__

#include "../../FdaPDE.h"

namespace lambda
{
	template<UInt size>
	using type = typename std::conditional<size==1, Real, VectorXr>::type;
	
	type<2> make_pair(Real lambdaS, Real lambdaT);
	type<2> make_pair(lambda::type<2> lambda, Real lambdaT); //degenerate case useful in parabolic case
}

#endif
