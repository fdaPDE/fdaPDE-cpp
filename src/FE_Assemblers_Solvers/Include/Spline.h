#ifndef __SPLINE_H__
#define __SPLINE_H__

#include "Integration.h"

template<UInt DEGREE, UInt ORDER_DERIVATIVE>
class Spline
{
    static_assert(DEGREE>=0 && ORDER_DERIVATIVE>=0,
                    "ERROR! TRYING TO INSTANTIATE SPLINE WITH NEGATIVE DEGREE/ORDER OF DERIVATIVE! See Spline.h");

    public:
        using Integrator=IntegratorGaussP5;

        //! A Constructor that initializes the vector of knots of the spline (including multiple knots)
        Spline(const Real *t_instants, const UInt n_time_instants)
        {
            knots_.reserve(2*DEGREE+n_time_instants);

            for (UInt i = 0; i < DEGREE; ++i)
                knots_.push_back(t_instants[0]);

            for (UInt i = 0; i < n_time_instants; ++i)
                knots_.push_back(t_instants[i]);

            for (UInt i = 0; i < DEGREE; ++i)
                knots_.push_back(t_instants[n_time_instants-1]);
        }


        Spline(const std::vector<Real>& t_instants) : 
            Spline(t_instants.data(), t_instants.size()) {}

        //! Method that return the number of knots of the spline
        UInt num_knots() const {return knots_.size();}

        //! Method that returns the i-th node of the spline
        const Real& getKnot(UInt i) const {return knots_[i];}

        Real BasisFunction(UInt i, Real u) const {return BasisFunction_impl(DEGREE, i, u);}

        Real BasisFunctionDerivative(UInt i, Real u) const {return BasisFunctionDerivative_impl(DEGREE, ORDER_DERIVATIVE, i, u);}

        Real time_mass_impl(UInt i, UInt j, Real u) const {return BasisFunctionDerivative(i, u) * BasisFunctionDerivative(j, u);}
        
    private:
        std::vector<Real> knots_;

        //! Method that computes the value of the i-th basis function in point u
        Real BasisFunction_impl(UInt degree, UInt i, Real u) const 
        {
            if(degree == 0) 
                return u >= knots_[i] && u < knots_[i+1] || u == knots_.back() && i == knots_.size()-DEGREE-2;

            else if(knots_[i+degree] == knots_[i])
                return (knots_[i+degree+1]-u)/(knots_[i+degree+1]-knots_[i+1])*BasisFunction_impl(degree-1, i+1, u);

            else if(knots_[i+degree+1] == knots_[i+1])
                return (u-knots_[i])/(knots_[i+degree]-knots_[i])*BasisFunction_impl(degree-1, i, u);

            else
                return (u-knots_[i])/(knots_[i+degree]-knots_[i])*BasisFunction_impl(degree-1, i, u) +
                       (knots_[i+degree+1]-u)/(knots_[i+degree+1]-knots_[i+1])*BasisFunction_impl(degree-1, i+1, u);
            
        }

        //! Method that computes the value of the orderDerivative-th derivative of the i-th basis function in point u
        Real BasisFunctionDerivative_impl(UInt degree, UInt orderDerivative, UInt i, Real u) const
        {
            if(degree == 0)
                return 0;
            else if(orderDerivative == 0) 
                return BasisFunction_impl(degree, i, u);
            else if(orderDerivative == 1)
                    if(knots_[i+degree] == knots_[i]) 
                        return -degree/(knots_[i+degree+1]-knots_[i+1])*BasisFunction_impl(degree-1, i+1, u);
                    else if(knots_[i+degree+1] == knots_[i+1]) 
                        return degree/(knots_[i+degree]-knots_[i])*BasisFunction_impl(degree-1, i, u);
                    else
                        return degree/(knots_[i+degree]-knots_[i])*BasisFunction_impl(degree-1, i, u) -
                            degree/(knots_[i+degree+1]-knots_[i+1])*BasisFunction_impl(degree-1, i+1, u);
            else
                if(knots_[i+degree] == knots_[i]) 
                    return -degree/(knots_[i+degree+1]-knots_[i+1])*BasisFunctionDerivative_impl(degree-1, orderDerivative-1, i+1, u);
                else if(knots_[i+degree+1] == knots_[i+1]) 
                    return degree/(knots_[i+degree]-knots_[i])*BasisFunctionDerivative_impl(degree-1, orderDerivative-1, i, u);
                else
                        return degree/(knots_[i+degree]-knots_[i])*BasisFunctionDerivative_impl(degree-1, orderDerivative-1, i, u) -
                            degree/(knots_[i+degree+1]-knots_[i+1])*BasisFunctionDerivative_impl(degree-1, orderDerivative-1, i+1, u);
         }


};


#endif
