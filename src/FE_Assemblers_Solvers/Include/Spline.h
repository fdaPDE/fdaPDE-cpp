#ifndef __SPLINE_H__
#define __SPLINE_H__

#include <iostream>

template<typename Integrator, UInt DEGREE, UInt ORDER_DERIVATIVE>
class Spline
{
    public:

        //! A Constructor that initializes the vector of knots of the spline (including multiple knots)
        Spline(const std::vector<Real>& t_instants)
        {
            UInt n_time_instants = t_instants.size();
            //std::cout << n_time_instants << std::endl;

            for (UInt i = 0; i < DEGREE; ++i)
                knots_.push_back(t_instants[0]);

            for (UInt i = 0; i < n_time_instants; ++i)
                knots_.push_back(t_instants[i]);

            for (UInt i = 0; i < DEGREE; ++i)
                knots_.push_back(t_instants[n_time_instants-1]);
        }

        Spline(const Real *t_instants, const UInt n_time_instants)
        {
	    // UInt n_time_instants = t_instants.size();
            //std::cout << n_time_instants << std::endl;

            for (UInt i = 0; i < DEGREE; ++i)
                knots_.push_back(t_instants[0]);

            for (UInt i = 0; i < n_time_instants; ++i)
                knots_.push_back(t_instants[i]);

            for (UInt i = 0; i < DEGREE; ++i)
                knots_.push_back(t_instants[n_time_instants-1]);
        }


        //! Method that return the number of knots of the spline
        UInt num_knots() const
        {
            return knots_.size();
        }

        //! Method that returns the i-th node of the spline
        const Real& getKnot(UInt i) const
        {
            return knots_[i];
        }

        //! Method that returns the degree of the spline
//        Real getDegree()
//        {
//            return DEGREE;
//        }

//        //! Method that prints the degree of the spline
//        void printDegree()
//        {
//            std::cout << "Degree of the spline: " << DEGREE << std::endl;
//        }


        //! Method that computes the value of the i-th basis function in point u
        Real BasisFunction(UInt degree, UInt i, Real u) const 
        {
            if(degree == 0)
            {
                if( (u >= knots_[i] && u < knots_[i+1]) || ((u == knots_[knots_.size()-1]) && (i == knots_.size()-DEGREE-2)))
                    return 1;
                else
                    return 0;
            }
            else
            {
                if((knots_[i+degree]-knots_[i]) == 0)
                    return (knots_[i+degree+1]-u)/(knots_[i+degree+1]-knots_[i+1])*BasisFunction(degree-1, i+1, u);
                else if((knots_[i+degree+1]-knots_[i+1]) == 0)
                    return (u-knots_[i])/(knots_[i+degree]-knots_[i])*BasisFunction(degree-1, i, u);
                else
                    return (u-knots_[i])/(knots_[i+degree]-knots_[i])*BasisFunction(degree-1, i, u) +
                       (knots_[i+degree+1]-u)/(knots_[i+degree+1]-knots_[i+1])*BasisFunction(degree-1, i+1, u);
            }
        }

        //! Method that computes the value of the orderDerivative-th derivative of the i-th basis function in point u
        Real BasisFunctionDerivative(UInt degree, UInt orderDerivative, UInt i, Real u) const
        {
            if(degree == 0)
	            return 0;
            else
				if(orderDerivative == 0) return BasisFunction(degree, i, u);
                else if(orderDerivative == 1)
                    if((knots_[i+degree]-knots_[i]) == 0) return -degree/(knots_[i+degree+1]-knots_[i+1])*BasisFunction(degree-1, i+1, u);
                    else if((knots_[i+degree+1]-knots_[i+1]) == 0) return degree/(knots_[i+degree]-knots_[i])*BasisFunction(degree-1, i, u);
                    else
                        return degree/(knots_[i+degree]-knots_[i])*BasisFunction(degree-1, i, u) -
                            degree/(knots_[i+degree+1]-knots_[i+1])*BasisFunction(degree-1, i+1, u);

                 else //if(orderDerivative == 2)
                    if((knots_[i+degree]-knots_[i]) == 0) return -degree/(knots_[i+degree+1]-knots_[i+1])*BasisFunctionDerivative(degree-1, orderDerivative-1, i+1, u);
                    else if((knots_[i+degree+1]-knots_[i+1]) == 0) return degree/(knots_[i+degree]-knots_[i])*BasisFunctionDerivative(degree-1, orderDerivative-1, i, u);
                    else
                        return degree/(knots_[i+degree]-knots_[i])*BasisFunctionDerivative(degree-1, orderDerivative-1, i, u) -
                            degree/(knots_[i+degree+1]-knots_[i+1])*BasisFunctionDerivative(degree-1, orderDerivative-1, i+1, u);
         }

         Real time_mass_impl(UInt i, UInt j, Real u) const {
           return BasisFunctionDerivative(DEGREE, ORDER_DERIVATIVE, i, u) * BasisFunctionDerivative(DEGREE, ORDER_DERIVATIVE, j, u);
         }

    private:
        std::vector<Real> knots_;
};


#endif
