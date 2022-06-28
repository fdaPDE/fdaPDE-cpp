#ifndef __SMW_H__
#define __SMW_H__

#include "../utils/Symbols.h"
#include <Eigen/QR>

// A general implementation of a system solver using the Sherman–Morrison–Woodbury formula for matrix inversion

// Consider a linear system Mx = b. The direct solution of this system would require the inversion of M, i.e. x = M^{-1}*b.
// If we are able to decompose matrix M as (A + UCV) we can take advantage of the SMW formula to compute
//
//        (A + UCV)^{-1} = A^{-1} - A^{-1}*U*(C^{-1} + V*A^{-1}*U)^{-1}*V*A^{-1}
//
// In case A is a sparse matrix and C a small dense matrix, computing M^{-1} using the above decomposition is much more efficient
// than computing M^{-1} directly

struct SMW{
  // solves linear system (A + U*C^{-1}*V)x = b, observe that we assume to supply the already computed inversion of the dense matrix C.
  // Indeed in some cases the inverse of C can be directly computed with known quantites without performing any inversion. In the case
  // this is not possible, you can find the inverse of C and then call this .solve() method
  static DVector solve(SpMatrix& A, SpMatrix& U, DMatrix& invC, SpMatrix V, DVector b){
    // bring A in compressed format
    A.makeCompressed();

    // compute sparseLU factorization of matrix A and store it for fast reuse
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > SparseLU;
    SparseLU.analyzePattern(A);    // Compute the ordering permutation vector from the structural pattern of A
    SparseLU.factorize(A);         // compute LU factorization of matrix A
    
    DVector y = SparseLU.solve(b);  // compute y = A^{-1}b
    DMatrix Y = SparseLU.solve(U);  // compute Y = A^{-1}U
    DMatrix G = invC + V*Y;         // compute dense matrix G = C^{-1} + V*A^{-1}*U = C^{-1} + V*y

    // solve dense linear system (C^{-1} + V*A^{-1}*U)z = V*A^{-1}*b <-> z = (C^{-1} + V*Y)^{-1}*V*y = G^{-1}*V*y
    Eigen::ColPivHouseholderQR<DMatrix> DenseQR;
    DenseQR.compute(G);
    DVector t = DenseQR.solve(V*y);
 
    // compute v = A^{-1}*U*t = A^{-1}*U*(C^{-1} + V*A^{-1}*U)^{-1}*V*A^{-1}*b by solving linear system A*v = U*t
    DVector v = SparseLU.solve(U*t);

    // return system solution
    return y - v;
  }
  
};

#endif // _SMW_H__
