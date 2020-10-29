#ifndef __FE_DENSITY_ESTIMATION_H__
#define __FE_DENSITY_ESTIMATION_H__

#include "Preprocess_Phase.h"
#include "Preprocess_Factory.h"

// This file is useful to perform the Density Estimation problem

/*! @brief A class to perform the whole density estimation problem.
*/
template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class FEDE{
  private:
    // A member to acess data problem methods
    const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dataProblem_;
    // A member to acess functional methods
    const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& funcProblem_;
    // A member to do the minimization phase
    std::shared_ptr<MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>> minAlgo_;
    // A member to do the preprocess phase
    std::unique_ptr<Preprocess<Integrator_noPoly, ORDER, mydim, ndim>> preprocess_;
    // A member to store the final density estimated
    VectorXr gcoeff_;
    // A member to store the initial densities selected
    std::vector<const VectorXr*> fInit_;
    // A member to save the best lambda
    Real bestLambda_;

    std::vector<Real> CV_errors_;

  public:
    //! A costructor
    FEDE(const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dp,
      const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& fp,
      std::shared_ptr<MinimizationAlgorithm<Integrator_noPoly, ORDER, mydim, ndim>> ma, const std::string& p);

    //! A method to perform the whole density estimation task.
    void apply();

    // Getters
    //! A method returning the estimated density coefficients.
    inline VectorXr getDensity_g() const {return gcoeff_;}
    //! A method returning initial densities.
    inline std::vector<const VectorXr*> getInitialDensity() const {return fInit_;}
    //! A method returning the smmothing parameter selected.
    inline Real getBestLambda() const {return bestLambda_;}

    // to delete
    inline std::vector<Real> getCvError() const {return CV_errors_;}


};

#include "FE_Density_Estimation_imp.h"

#endif
