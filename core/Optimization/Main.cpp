#include "GridOptimizer.h"
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <math.h>
#include <ostream>
#include <utility>

#include "Function.h"
#include "NewtonOptimizer.h"
#include "ExactNewtonOptimizer.h"
#include "GradientDescentOptimizer.h"
#include "BFGSOptimizer.h"

//#include "SVector.h"


int main() {

  // define target to optimize
  // std::function<double(std::array<double,1>)> fun = [](std::array<double,1> x) -> double { return 2*x[0] - std::sqrt(2)*std::pow(x[0], 2) + 4*std::pow(x[0], 3); };

  // // perform a 1D grid optimization

  // // domain of optimization
  // std::array<std::pair<double, double>, 1> domain = {std::pair<double, double>(1,5)};
  // std::array<double,1> lambda = {0.01};

  // GridOptimizer<1> opt(domain, lambda);
  // std::pair<array<double, 1>,double> min = opt.findMinimum(fun);

  // std::cout << min.second << std::endl;


  // // define target to optimize
  // std::function<double(std::array<double,2>)> g = [](std::array<double,2> x) -> double { return 2*std::pow(x[0],2) - 2*std::pow(x[1],2); };

  // // perform a 2D grid optimization

  // // domain of optimization
  // std::array<std::pair<double, double>, 2> domain2D = {std::pair<double, double>(1,5), std::pair<double, double>(1,5)};
  // std::array<double,2> lambda2D = {0.01, 0.01};

  // GridOptimizer<2> opt2D(domain2D, lambda2D);
  // std::pair<array<double, 2>,double> min_g = opt2D.findMinimum(g);

  // std::cout << "2D optimization" << std::endl;
  
  // std::cout << min_g.first[0] << std::endl;
  // std::cout << min_g.first[1] << std::endl;
  // std::cout << min_g.second << std::endl;

  // std::function<double(Point<2>)> g = [](Point<2> x) -> double { return 2*std::pow(x[0],2) - 2*std::pow(x[1],2)*x[0]; };

  // std::cout << g(Point<2>(4,1)) << std::endl;

  // ScalarField<2> fun(g);

  // std::cout << fun.evaluateAtPoint(Point<2>(4,1)) << std::endl;
  
  // Point<2> grad = fun.getGradientApprox(Point<2>(2,1), 0.001);
  
  
  // std::cout << "evaluation of gradient at point" << std::endl;
  // std::cout << grad << std::endl;

  // Eigen::Matrix<double, 2, 2> hessian = fun.getHessianApprox(Point<2>(2,1), 0.001);

  // std::cout << "hessian at point" << std::endl;
  // std::cout << hessian << std::endl;

  // // define differentiable scalar field
  // std::function<Point<2>(Point<2>)> dg = [](Point<2> x) -> Point<2> { return Point<2>(4*x[0] - 2*std::pow(x[1],2), -4*x[1]*x[0]); };


  // Point<2> x(1,1);
  // std::cout << "printo punto" << std::endl;
  // std::cout << dg(x) << std::endl;

  // DifferentiableScalarField<2> ddfun(g, dg);

  // std::cout << "printo punto da scalar field" << std::endl;
  // std::cout << ddfun.derive()(x) << std::endl;

  // newton optimization
  // std::function<double(Point<2>)> g_newton = [](Point<2> x) -> double { return 2*std::pow(x[0],2) + 2*std::pow(x[1],2); };

  // // perform a 2D grid optimization

  // std::array<std::pair<double, double>, 2> domain = {std::pair<double, double>(-1,1), std::pair<double, double>(-1,1)};
  // double lambda = 0.01;
  
  // // domain of optimization
  // NewtonOptimizer<2> opt_newton2D(domain, lambda);
  // ScalarField<2> newton_objective(g_newton);

  // // add check to see if initial point is out of domain
  // std::pair<Point<2>,double> min_g = opt_newton2D.findMinimum(Point<2>(1,1), 1000, 0.001, newton_objective);

  // std::cout << "newton ottimo" << std::endl;
  // std::cout << min_g.second << std::endl;
  // std::cout << "point" << std::endl;
  // std::cout << min_g.first << std::endl;
    
  std::cout << "*******************" << std::endl;

  //std::function<double(SVector<2>)> g = [](SVector<2> x) -> double { return 2*std::pow(x[0],2) - 2*std::pow(x[1],2)*x[0]; };

  //std::cout << g(SVector<2>({4,1})) << std::endl;

  //ScalarField<2> fun(g);

  //std::cout << fun.evaluateAtPoint(SVector<2>({4,1})) << std::endl;
  
  //SVector<2> grad = fun.getGradientApprox(SVector<2>({2,1}), 0.001);
  
  
  // std::cout << "evaluation of gradient at point" << std::endl;
  // std::cout << grad << std::endl;

  // Eigen::Matrix<double, 2, 2> hessian = fun.getHessianApprox(SVector<2>({2,1}), 0.001);

  // std::cout << "hessian at point" << std::endl;
  // std::cout << hessian << std::endl;

  // define differentiable scalar field
  //std::function<SVector<2>(SVector<2>)> dg = [](SVector<2> x) -> SVector<2> { return SVector<2>({4*x[0] - 2*std::pow(x[1],2), -4*x[1]*x[0]}); };


  // SVector<2> x({1,1});
  // std::cout << "printo punto" << std::endl;
  // std::cout << dg(x) << std::endl;

  // std::cout << "printo punto da scalar field" << std::endl;
  // std::cout << ddfun.derive()(x) << std::endl;

  // std::cout << "#########################" << std::endl;

  // newton optimization
  std::function<double(SVector<2>)> g_newton = [](SVector<2> x) -> double { return 2*std::pow(x[0],2) + x[0] + 2*std::pow(x[1],2); };

  // perform a 2D grid optimization

  //std::array<std::pair<double, double>, 2> domain = {std::pair<double, double>(-1,1), std::pair<double, double>(-1,1)};
  double lambda = 0.01;
  
  // domain of optimization
  NewtonOptimizer<2> opt_newton2D(lambda);
  ScalarField<2> newton_objective(g_newton);

  // add check to see if initial point is out of domain
  std::pair<SVector<2>,double> min_g = opt_newton2D.findMinimum(SVector<2>(1,1), 1000, 0.001, newton_objective);

  std::cout << "newton ottimo" << std::endl;
  std::cout << min_g.second << std::endl;
  std::cout << "point" << std::endl;
  std::cout << min_g.first << std::endl;

  // define differentiable scalar field
  std::function<SVector<2>(SVector<2>)> dg = [](SVector<2> x) -> SVector<2> { return SVector<2>({4*x[0] + 1, 4*x[1]}); };
  std::function<SMatrix<2>(SVector<2>)> ddg = [](SVector<2> x) -> SMatrix<2> { return SMatrix<2>({{4, 0},
												  {0, 4}});};

  TwiceDifferentiableScalarField<2> dddfun(g_newton, dg, ddg);

  ExactNewtonOptimizer<2> opt_newton2Dexact(lambda);

  // add check to see if initial point is out of domain
  std::pair<SVector<2>,double> min_g_exact = opt_newton2Dexact.findMinimumExact(SVector<2>(1,1), 1000, 0.001, dddfun);

  std::cout << "exact newton ottimo" << std::endl;
  std::cout << min_g_exact.second << std::endl;
  std::cout << "point" << std::endl;
  std::cout << min_g_exact.first << std::endl;

  // std::cout << "!!!!! prova 1D !!!!!" << std::endl;

  // std::function<double(SVector<1>)> f1D = [](SVector<1> x) -> double { return 2*std::pow(x[0], 4) - 32*std::sqrt(x[0]) + 1; };
  
  // NewtonOptimizer<1> opt_newton1D(lambda);
  // ScalarField<1> field_f(f1D);

  // std::pair<SVector<1>,double> min_f = opt_newton1D.findMinimum(SVector<1>(1), 1000, 0.001, field_f);

  // std::cout << "newton ottimo" << std::endl;
  // std::cout << min_f.second << std::endl;
  // std::cout << "point" << std::endl;
  // std::cout << min_f.first << std::endl;

  std::cout << "----------------- gradient descent" << std::endl;

  DifferentiableScalarField<2> ddfun(g_newton, dg);
  GradientDescentOptimizer<2> gradientOptim(0.001, SVector<2>(1,1), 5000, 0.001, ddfun);

  //gradientOptim.init = [&gradientOptim]() mutable -> void { gradientOptim.setControllerData("it", {}); };

  //gradientOptim.beginIteration = [&gradientOptim]() mutable -> void {
  //  gradientOptim.setControllerData("it", gradientOptim.getError());
  //};
  
  // inject this controller in the optimization process  
  std::pair<SVector<2>,double> min_grad = gradientOptim.findMinimum();

  std::cout << "gradient fixed step" << std::endl;
  std::cout << min_grad.second << std::endl;
  std::cout << "point" << std::endl;
  std::cout << min_grad.first << std::endl;

  std::cout << "----------------- BFGS" << std::endl;

  BFGSOptimizer<2> bfgsOptim(0.001, SVector<2>(1,1), 1000, 0.001, ddfun);
  
  std::pair<SVector<2>,double> min_bfgs = bfgsOptim.findMinimum();

  std::cout << "BFGS fixed step" << std::endl;
  std::cout << min_bfgs.second << std::endl;
  std::cout << "point" << std::endl;
  std::cout << min_bfgs.first << std::endl;
  
  std::cout << "~~~~~~~~~~~~~~~~~ Grid" << std::endl;
  // perform a 2D grid optimization

  // domain of optimization
  std::array<std::pair<double, double>, 2> domain2D = {std::pair<double, double>(-1,1), std::pair<double, double>(-1,1)};
  std::array<double,2> lambda2D = {0.001, 0.001};

  GridOptimizer<2> opt2D(domain2D, lambda2D, newton_objective);
  std::pair<SVector<2>,double> min_grid = opt2D.findMinimum();

  std::cout << "2D optimization" << std::endl;
  
  std::cout << min_grid.second << std::endl;
  std::cout << min_grid.first << std::endl;
  
  return 0;
}
