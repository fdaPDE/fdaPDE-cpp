#include "../Include/Auxiliary_Optimizer.h"

//! Utility method for recycle boundary conditions correction
/*!
 \param mat the matrix on which to perform the correction, passed by reference
 \param bc_idxp pointer of boundary condition indices
 \note version for sparse matrices
*/
void AuxiliaryOptimizer::bc_utility(MatrixXr & mat, const std::vector<UInt> * bc_idxp)
{
        /*
        	In the case of no boundary conditions, nbc_indices=0
        	In the case of boundary conditions, we refer to the number of matrix
        	rows as this function is used both in the monolithic and iterative
        	cases. In the iterative case, the number of boundary conditions is
        	(length of time mesh) times the number of spaces. We apply only the
        	first (number of spaces) boundary conditions to the reduced matrix.
	*/
        UInt nbc_indices = (bc_idxp->size() < mat.rows()) ? bc_idxp->size() : mat.rows();
        if(nbc_indices!=0) // Add boundary conditions
        {
                Real pen = 10e20;
                for(UInt i=0; i<nbc_indices; i++)
                {
                        UInt id = (*bc_idxp)[i];
                        mat(id,id) = pen;
                }

        }
}

//! Utility method for recycle boundary conditions correction
/*!
 \param mat the matrix on which to perform the correction, passed by reference
 \param bc_idxp pointer of boundary condition indices
 \note version for full matrices
*/
void AuxiliaryOptimizer::bc_utility(SpMat & mat, const std::vector<UInt> * bc_idxp)
{
	//see comment on the function overload above
        UInt nbc_indices = (bc_idxp->size() < mat.rows()) ? bc_idxp->size() : mat.rows();
        if(nbc_indices!=0) // Add boundary conditions
        {
                Real pen = 10e20;
                for(UInt i=0; i<nbc_indices; i++)
                {
                        UInt id = (*bc_idxp)[i];
                        mat.coeffRef(id,id) = pen;
                }

        }
}

//! Utility method to compute matrix E in locationbynodes pointwise
/*!
 \param E the matrix to fill, passed by reference
 \param kp pointer to identiy locations
 \param Qp pointer to Q projection matrix
 \param nr number of nodes
 \param s number of observations
*/
void AuxiliaryOptimizer::set_E_ln_W_ptw(MatrixXr & E, const std::vector<UInt> * kp, const MatrixXr * Qp, UInt nr, UInt s)
{
        E = MatrixXr::Zero(nr, s);

        for (UInt i = 0; i < s ; i++)
                for (int j = 0; j < s; j++)
                        E.coeffRef((*kp)[i], j) += (*Qp).coeff(i, j);
}

//! Utility method to compute matrix E in not-locationbynodes pointwise
/*!
 \param E the matrix to fill, passed by reference
 \param psi_tp pointer to the transpose of Psi matrix
 \param Qp pointer to Q projection matrix
*/
void AuxiliaryOptimizer::set_E_lnn_W_ptw(MatrixXr & E, const SpMat * psi_tp, const MatrixXr * Qp)
{
        E = ((*psi_tp)*(*Qp));
}

//! Utility method to compute matrix E in areal setting, with regression
/*!
 \param E the matrix to fill, passed by reference
 \param psi_tp pointer to the transpose of Psi matrix
 \param Qp pointer to Q projection matrix
 \param Ap pointer to areal vector
*/
void AuxiliaryOptimizer::set_E_W_a(MatrixXr & E, const SpMat * psi_tp, const MatrixXr * Qp, const VectorXr * Ap)
{
        E = ((*psi_tp)*(*Ap).asDiagonal()*(*Qp));
}

//! Utility method to compute matrix E in areal setting, without regression
/*!
 \param E the matrix to fill, passed by reference
 \param psi_tp pointer to the transpose of Psi matrix
 \param Ap pointer to areal vector
*/
void AuxiliaryOptimizer::set_E_nW_a(MatrixXr & E, const SpMat * psi_tp, const VectorXr * Ap)
{
        E = ((*psi_tp)*(*Ap).asDiagonal());
}
