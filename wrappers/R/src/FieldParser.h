#ifndef __FIELD_PARSER_H__
#define __FIELD_PARSER_H__

#include <mpParser.h> // muparserx main header
#include <memory>
#include <fdaPDE/Core.h>
using fdaPDE::core::ScalarField;

// a simple template class to wrap the result of muparserx into a ScalarField object
template <unsigned int N>
class FieldParser {
private:
  // wrap by pointers to let FieldParser copy constructible and avoid SEGFAULT
  std::unique_ptr<mup::ParserX> p_; // parser instance
  std::unique_ptr<mup::Value> x_vect_; // memory reserved for parsed variable
  std::string expr_; // field's expression
  
public:
  // construcotr
  FieldParser() = default;
  FieldParser(std::string expr)
    : p_(std::make_unique<mup::ParserX>(mup::pckALL_NON_COMPLEX)),
      x_vect_(std::make_unique<mup::Value>(N,0)), expr_(expr) {
    // init parser
    p_->DefineVar("x", mup::Variable(x_vect_.get()));
    p_->SetExpr(expr_); // parse the field expression
  }

  // returns a ScalarField wrapping the parsed field. Observe that the FieldParsed must be in the same scope
  // of the ScalarField to avoid SEGFAULT for the whole life of the ScalarField
  ScalarField<N> get() {
    return ScalarField<N>
      ([this](const SVector<N>& x) mutable -> double {
	// update the memory observed by the parsed with the given point x
	for(std::size_t i = 0; i < N; ++i)
	  x_vect_->At(i) = x[i];
	
	// evaluate the field
	return p_->Eval().GetFloat();
      });
  }

  // ask directly for an evaluation of the field, preserve callable interface
  double operator()(const SVector<N>& x) {
    // update the memory observed by the parsed with the given point x
    for(std::size_t i = 0; i < N; ++i)
      x_vect_->At(i) = x[i];
	
    // evaluate the field
    return p_->Eval().GetFloat();    
  }
  
};

#endif // __FIELD_PARSER_H__
