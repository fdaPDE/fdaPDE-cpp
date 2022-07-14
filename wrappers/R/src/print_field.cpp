#include <Rcpp.h>
#include <fdaPDE/Core.h>
using fdaPDE::core::ScalarField;

using namespace Rcpp;

// just a function to print a ScalarField
// [[Rcpp::export]]
List print_field() {

  // define a scalar field: x*e^{-x^2 - y^2} + (x^2 + y^2)/20
  ScalarField<2> field([](SVector<2> x) -> double {
    return x[0]*std::exp(- std::pow(x[0],2) - std::pow(x[1], 2)) + (std::pow(x[0],2) + std::pow(x[1], 2))/20;
  });

  NumericVector result;
  NumericVector x_coord;
  NumericVector y_coord;
  
  for(double x = -2.0; x <= 2.0; x+=0.1){
    for(double y = -2.0; y <= 2.0; y+=0.1){
      double eval = field(SVector<2>(x,y));
      result.push_back(eval);
      x_coord.push_back(x);
      y_coord.push_back(y);
    }
  }
  return List::create(Named("x") = x_coord, Named("y") = y_coord, Named("field") = result);
}
