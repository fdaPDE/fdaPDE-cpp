#ifndef __INITIAL_CONDITION_ESTIMATOR_H__
#define __INITIAL_CONDITION_ESTIMATOR_H__

// core imports
#include "../core/FEM/PDE.h"
#include "../core/FEM/operators/Laplacian.h"
using fdaPDE::core::FEM::PDE;
using fdaPDE::core::FEM::Laplacian;
#include "../core/OPT/optimizers/GridOptimizer.h"
using fdaPDE::core::OPT::GridOptimizer;
// calibration module imports
#include "../calibration/GCV.h"
using fdaPDE::calibration::GCV;
using fdaPDE::calibration::ExactEDF;
// models module imports
#include "../models/ModelTraits.h"
using fdaPDE::models::is_generalized;
#include "../models/regression/SRPDE.h"
#include "../models/regression/GSRPDE.h"
using fdaPDE::models::SpaceTimeParabolicTag;
using fdaPDE::models::SRPDE;
using fdaPDE::models::GSRPDE;

namespace fdaPDE{
namespace preprocess {

  // trait to select the type of space-only model to use for estimation of initial condition
  template <typename Model, typename PDE>
  struct ICEstimator_internal_solver{
    typedef typename std::decay<Model>::type Model_;
    using type = typename std::conditional<
      !is_generalized<Model_>::value,
      // STRPDE model
      SRPDE <PDE, model_traits<Model_>::sampling>,
      // generalized STRPDE model
      GSRPDE<PDE, typename fdaPDE::models::SpaceOnlyTag, model_traits<Model_>::sampling,
	     model_traits<Model_>::solver, typename model_traits<Model_>::DistributionType>
      >::type;
  };
  
  // for a space-time regression model builds an estimation of the initial condition from the data at time step 0
  // the estimation is obtained selecting the best spatial field estimate obtained from an SRPDE model via GCV optimization
  template <typename Model>
  class InitialConditionEstimator {
    static_assert(std::is_same<typename model_traits<Model>::RegularizationType, SpaceTimeParabolicTag>::value,
		  "Initial condition estimation is for parabolic regularization only");
  private:
    const Model& model_; // model to initialize
    DMatrix<double> estimate_;
  public:
    // constructor
    InitialConditionEstimator(const Model& model) : model_(model) {};
    // builds the initial condition estimate
    void apply(const std::vector<SVector<1>>& lambdas){
      // extract data at time step 0
      std::size_t n = model_.n_locs();
      BlockFrame<double, int> df = model_.data()(0, n-1).extract();

      // cast space-time differential operator df/dt + Lf = u to space-only Lf = u
      auto L = model_.pde().bilinearForm().template remove_operator<dT>();
      // prepare regularization term
      typedef PDE<Model::M, Model::N, Model::K, decltype(L), DMatrix<double>> PDE_;
      PDE_ problem(model_.domain());
      DMatrix<double> u = model_.pde().forcingData().col(0);
      problem.setBilinearForm(L);
      problem.setForcing(u);
      problem.init(); // init PDE object
      
      // define solver for initial condition estimation
      typename ICEstimator_internal_solver<Model, decltype(problem)>::type solver(problem, model_.locs());
      solver.setData(df); // impose data
      solver.init();
      
      // find optimal smoothing parameter
      GCV<decltype(solver), ExactEDF<decltype(solver)>> GCV(solver);
      GridOptimizer<1> opt;

      ScalarField<1, decltype(GCV)> obj(GCV);
      opt.findMinimum(obj, lambdas); // optimize gcv field
      SVector<1> best_lambda = opt.getSolution();

      // fit model with optimal lambda
      solver.setLambda(best_lambda);
      solver.solve();
      // store initial condition estimate
      estimate_ = solver.f();
      return;
    }

    // getters
    const DMatrix<double>& get() const { return estimate_; }
  };
  
}}
#endif // __INITIAL_CONDITION_ESTIMATOR_H__
