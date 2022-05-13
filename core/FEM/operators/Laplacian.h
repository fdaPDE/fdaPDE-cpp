#ifndef __LAPLACIAN_H__
#define __LAPLACIAN_H__

#include "../utils/fields/VectorField.h"
using fdaPDE::core::InnerProduct;
using fdaPDE::core::VectorField;
#include "../MESH/Element.h"
using fdaPDE::core::MESH::Element;

#include "BilinearFormExpressions.h"

// class representing the laplacian operator (isotropic diffusion)
struct Laplacian : public BilinearFormExpr<Laplacian>{

  // provide the discretization for the laplacian operator. In particular this method implements a custom quadrature rule
  // for approximating the (i,j)-th element of the stiffness matrix \int_e [\Nabla phi_i * \Nabla phi_j]
  // integrate() will be called by Integrator as a result of the expression template expansion of the problem's bilinear form

  // basis: any type compliant with a functional basis behaviour. See LagrangianBasis for an example
  // e: the element where we are integrating
  // i,j: indexes of the stiffness matrix element we are computing
  // quadrature_point: the point where to evaluate the integrand
  template <unsigned int N, int M, unsigned int ORDER, typename B>
  double integrate(const B& basis, const Element<ORDER, N>& e, int i , int j, const SVector<M>& quadrature_point) const{
    // express gradient of f in terms of gradients of basis functions over reference element.
    // This entails to compute (J^{-1})^T * \Nabla phi_i. In the following we assume basis[i] = phi_i

    Eigen::Matrix<double, N, ORDER> invJ = e.getInvBaryMatrix().transpose();
    // Given \Nabla phi_i premultiply it by (J^{-1})^T = invJ.
    // NOTE: we assume "basis" to provide functions already defined on the reference element
    VectorField<N> NablaPhi_i = invJ * basis[i].gradient();  
    VectorField<N> NablaPhi_j = invJ * basis[j].gradient();
    
    // for laplacian: \nabla phi_i * \nabla * phi_j
    return NablaPhi_i.dot(NablaPhi_j)(quadrature_point);
  }
};

#endif // __LAPLACIAN_H__
