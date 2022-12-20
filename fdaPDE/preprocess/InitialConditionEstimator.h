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
using fdaPDE::calibration::ExactGCV;
// models module imports
#include "../models/ModelTraits.h"
using fdaPDE::models::is_sampling_pointwise_at_mesh;
#include "../models/regression/SRPDE.h"
using fdaPDE::models::SpaceTimeParabolicTag;
using fdaPDE::models::SRPDE;

namespace fdaPDE{
namespace preprocess {

  // tries to build an estimation of the initial condition from the data at time step 0
  // the estimation is obtained selecting the best (according to GCV) spatial field estimate obtained from
  // an SRPDE model under simple laplacian regularization
  template <typename Model>
  class InitialConditionEstimator {
    static_assert(std::is_same<typename model_traits<Model>::RegularizationType, SpaceTimeParabolicTag>::value,
		  "Initial condition estimation is for parabolic regularization only");
  private:
    const Model& model_;
  public:
    InitialConditionEstimator(const Model& model) : model_(model) {};

    DMatrix<double> apply(const std::vector<SVector<1>>& lambdas){
      // extract data at time step 0
      std::size_t n = model_.n_locs();
      BlockFrame<double, int> df = model_.data()(0, n-1).extract();

      // define simple laplacian regularization term
      auto L = Laplacian();
      typedef PDE<Model::M, Model::N, Model::K, decltype(L), DMatrix<double>> PDE_;
      PDE_ problem(model_.domain());
      DMatrix<double> u = DMatrix<double>::Zero(problem.domain().elements()*problem.integrator().numNodes(), 1);
      problem.setBilinearForm(L);
      problem.setForcing(u);

      problem.init();
      // define SRPDE model for initial condition estimation
      SRPDE<PDE_, model_traits<Model>::sampling> m(problem, model_.locs());

      // impose data
      m.setData(df);
      
      // define stochastic GCV calibrator
      ExactGCV<decltype(m)> calibrator(m);
      GridOptimizer<1> gcv_optimizer;
      
      // find optimal smoothing parameter
      SVector<1> opt = calibrator.approx(gcv_optimizer, 0.01, lambdas);
      // check if lambda is on border of domain

      // fit model with optimal lambda
      m.setLambda(opt[0]);
      m.solve();
      
      // return initial condition estimate
      return m.f();
    }
    
  };
  
}}
#endif // __INITIAL_CONDITION_ESTIMATOR_H__
