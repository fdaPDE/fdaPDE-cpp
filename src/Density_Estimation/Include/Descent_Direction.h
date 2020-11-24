#ifndef __DESCENT_DIRECTION_H__
#define __DESCENT_DIRECTION_H__

#include "../../Global_Utilities/Include/Make_Unique.h"

// This file contains the direction search technique useful for the optimization algorithm of the Density Estimation problem

/*! @brief An abstract class for computing the descent direction. The father is pure virtual; the right direction
is computed inside the children according to a proper method chosen.
*/
template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class DirectionBase{
  protected:
    // to give generality if you want to add other children
    const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& funcProblem_;

  public:
    //! A constructor
    DirectionBase(const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& fp): funcProblem_(fp){};
    //! A destructor.
    virtual ~DirectionBase(){};
    //! A pure virtual clone method.
    virtual std::unique_ptr<DirectionBase<Integrator_noPoly, ORDER, mydim, ndim>> clone() const = 0;
    //! A pure virtual method to compute the descent direction.
    virtual VectorXr computeDirection(const VectorXr& g, const VectorXr& grad) = 0;
    //! A pure virtual method to reset all the old parameters.
    virtual void resetParameters() = 0;
};


//! @brief A class for computing the gradient descent direction.
template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class DirectionGradient : public DirectionBase<Integrator_noPoly, ORDER, mydim, ndim>{
  public:
    //! A delegating constructor.
    DirectionGradient(const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& fp):
    DirectionBase<Integrator_noPoly, ORDER, mydim, ndim>(fp){};
    //! Clone method overridden.
    std::unique_ptr<DirectionBase<Integrator_noPoly, ORDER, mydim, ndim>> clone() const override;
    //! A method to compute the gradient descent direction.
    VectorXr computeDirection(const VectorXr& g, const VectorXr& grad) override;
    //! A method to reset all the old parameters. In the gradient method they aren't.
    void resetParameters() override {};

};


//! @brief A class for computing the BFGS descent direction.
template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class DirectionBFGS : public DirectionBase<Integrator_noPoly, ORDER, mydim, ndim>{
  private:
    MatrixXr HInit_, HOld_;
    VectorXr gOld_, gradOld_;
    bool updateH_;

  public:
    //! A constructor.
    DirectionBFGS(const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& fp, UInt k):
    DirectionBase<Integrator_noPoly, ORDER, mydim, ndim>(fp), HInit_(MatrixXr::Identity(k, k)), HOld_(MatrixXr::Identity(k, k)), updateH_(false){};
    //! A copy constructor: it just creates a DirectionBFGS object with the same features of rhs but it initializes the matrices.
    DirectionBFGS(const DirectionBFGS<Integrator_noPoly, ORDER, mydim, ndim>& rhs);
    //! Clone method overridden.
    std::unique_ptr<DirectionBase<Integrator_noPoly, ORDER, mydim, ndim>> clone() const override;
    //! A method to compute the BFGS descent direction.
    VectorXr computeDirection(const VectorXr& g, const VectorXr& grad) override;
    //! A method to reset all the old parameters.
    void resetParameters() override;

};

#include "Descent_Direction_imp.h"
#endif
