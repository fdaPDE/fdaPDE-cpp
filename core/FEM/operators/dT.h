#ifndef __DT_H__
#define __DT_H__

#include "BilinearFormExpressions.h"

// class representing the time derivative operator. This operator just acts as a tag in a bilinear form expression
// to activate the time loop in the assembler phase. dT() is detected by the is_symmetric trait
template <unsigned int L = 0> class dT : public BilinearFormExpr<dT<L>> {
public:
  std::tuple<dT<L>> getTypeList() const { return std::make_tuple(*this); }
};

#endif // __DT_H__
