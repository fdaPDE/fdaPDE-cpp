#include "../Include/Lambda.h"

lambda::type<2> lambda::make_pair(Real lambdaS, Real lambdaT)
{
	return (VectorXr(2,1) << lambdaS, lambdaT).finished();
}

lambda::type<2> lambda::make_pair(lambda::type<2> lambda, Real lambdaT)
{
	return lambda;
}
