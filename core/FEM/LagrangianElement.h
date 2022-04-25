#ifndef __LAGRANGIAN_ELEMENT__
#define __LAGRANGIAN_ELEMENT__

#include "../MESH/Element.h"
#include "LagrangianBasis.h"

using fdaPDE::core::MESH::Element;


// a class representing a Lagrangian finite element
template <unsigned int N, unsigned int M>
class LagrangianElement{

  Element<N,M> domain;
  //LagrangianBasis<N,ORDER> basis;
};


#endif // __LAGRANGIAN_ELEMENT__
