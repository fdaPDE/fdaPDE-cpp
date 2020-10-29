#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "../../FdaPDE.h"

//!  A Linear System QR solver class
/*!
 * This class gives offers a standard interface to the QR resolutor for dense matrices
*/
class QR{
	public:
	static void solve(const MatrixXr & A, const VectorXr & b,VectorXr &x){x=A.householderQr().solve(b);};
};

//!  A Linear System LU Partial Pivoting solver class
/*!
 * This class gives offers a standard interface to the LU Partial Pivoting resolutor for dense matrices.
 * OBS: The matrix should be invertible.
*/
class LUPV{
	public:
	static void solve(MatrixXr const & A, VectorXr const & b,VectorXr &x){x=A.partialPivLu().solve(b);};
};

//!  A Linear System LDLT solver class
/*!
 * This class gives offers a standard interface to the LDLT resolutor for dense matrices.
 * OBS: The matrix should be symmetric and SDP.
*/
class Symmetric{
	public:
	static void solve(MatrixXr const & A, VectorXr const & b,VectorXr &x){x=A.ldlt().solve(b);};
};

//!  A Linear System Cholesky solver class
/*!
 * This class gives offers a standard interface to the Cholesky resolutor for dense matrices.
 * OBS: The matrix should be symmetric and SDP, faster and more stable than others.
*/
class Cholesky{
	public:
	static void solve(MatrixXr const & A, VectorXr const & b,VectorXr &x){x=A.ldlt().solve(b);};
};

//!  A Linear System LU sparse solver class
/*!
 * This class gives offers a standard interface to the LU resolutor for sparse matrices.
*/
class SpLU{
	public:
	static void solve(SpMat const & A, VectorXr const & b, VectorXr &x )
	{
		Eigen::SparseLU<SpMat> solver;
		solver.compute(A);
		if(solver.info()!=Eigen::Success){
		//std::cerr<<"Decomposition failed!"<<std::endl;
		}
		x=solver.solve(b);
		if(solver.info()!=Eigen::Success)
		{
		//std::cerr<<"solving failed!"<<std::endl;
		}
	};
};

//!  A Linear System QR sparse solver class
/*!
 * This class gives offers a standard interface to the QR resolutor for sparse matrices.
*/
class SpQR{
	public:
	static void solve(SpMat const & A, VectorXr const & b, VectorXr &x )
	{
		Eigen::SparseQR<SpMat,Eigen::COLAMDOrdering<int> > solver;
		solver.compute(A);
		if(solver.info()!=Eigen::Success){
		//std::cerr<<"Decomposition failed!"<<std::endl;
		}
		x=solver.solve(b);
		if(solver.info()!=Eigen::Success)
		{
		//std::cerr<<"solving failed!"<<std::endl;
		}
	};
};

//!  A Linear System Cholesky sparse solver class
/*!
 * This class gives offers a standard interface to the Cholesky resolutor for sparse matrices.
*/
class SpCholesky{
	public:
	static void solve(SpMat const & A, VectorXr const & b, VectorXr &x )
	{
		Eigen::SimplicialLDLT<SpMat> solver;
		solver.compute(A);
		if(solver.info()!=Eigen::Success){
		//std::cerr<<"Decomposition failed!"<<std::endl;
		}
		x=solver.solve(b);
		if(solver.info()!=Eigen::Success)
		{
		//std::cerr<<"solving failed!"<<std::endl;
		}
	};
};

//!  A Linear System Conjugate Gradient sparse solver class
/*!
 * This class gives offers a standard interface to the Conjugate Gradient resolutor for sparse matrices.
*/
class SpConjGrad{
	public:
	static void solve(SpMat const & A, VectorXr const & b, VectorXr &x )
	{
		Eigen::ConjugateGradient<SpMat> solver;
		solver.compute(A);
		if(solver.info()!=Eigen::Success){
		//std::cerr<<"Decomposition failed!"<<std::endl;
		}
		x=solver.solve(b);
		if(solver.info()!=Eigen::Success)
		{
		//std::cerr<<"solving failed!"<<std::endl;
		}
	};
};

//!  A Linear System BiConjugate Gradient stabilized sparse solver class
/*!
 * This class gives offers a standard interface to the BiConjugate Gradient stabilized resolutor for sparse matrices.
*/

class BiCGSTAB{
	public:
	static void solve(SpMat const & A, VectorXr const & b, VectorXr &x )
	{
		Eigen::BiCGSTAB<SpMat> solver;
		solver.compute(A);
		if(solver.info()!=Eigen::Success){
		//std::cerr<<"Decomposition failed!"<<std::endl;
		}
		x=solver.solve(b);
		if(solver.info()!=Eigen::Success)
		{
		//std::cerr<<"solving failed!"<<std::endl;
		}
	};
};

//!  A Linear System BiConjugate Gradient stabilized with Incomplete LUT preconditioner sparse solver class
/*!
 * This class gives offers a standard interface to the BiConjugate Gradient stabilized BiConjugate Gradient stabilized with Incomplete LUT preconditioner resolutor for sparse matrices.
*/

class BiCGSTABILUT{
	public:
	static void solve(SpMat const & A, VectorXr const & b, VectorXr &x )
	{
		Eigen::BiCGSTAB<SpMat,Eigen::IncompleteLUT<Real>> solver;
		solver.compute(A);
		if(solver.info()!=Eigen::Success){
		//std::cerr<<"Decomposition failed!"<<std::endl;
		}
		x=solver.solve(b);
		if(solver.info()!=Eigen::Success)
		{
		//std::cerr<<"solving failed!"<<std::endl;
		}
	};
};

#endif
