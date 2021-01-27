#ifndef __DATA_PROBLEM_H__
#define __DATA_PROBLEM_H__

#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <utility>
#include <numeric>
//#include <omp.h>
#include "../../FdaPDE.h"
#include "DE_Data.h"
#include "../../FE_Assemblers_Solvers/Include/Projection.h"
#include "../../FE_Assemblers_Solvers/Include/Integration.h"

// This file contains data informations for the Density Estimation problem

//! @brief A class to store common data for the problem.
template<UInt ORDER, UInt mydim, UInt ndim>
class DataProblem{
private:
    using Integrator = typename DensityIntegratorHelper::Integrator<mydim>;
    static constexpr UInt EL_NNODES = how_many_nodes(ORDER,mydim);
    DEData<ndim> deData_;
    MeshHandler<ORDER, mydim, ndim> mesh_;
    SpMat R0_, R1_, GlobalPsi_;
    MatrixXr P_;
    Eigen::Matrix<Real, Integrator::NNODES, EL_NNODES> PsiQuad_;

    //! A method to compute the finite element matrices.
    void fillFEMatrices();
    //! A method to compute the matrix which evaluates the basis function at the quadrature EL_NNODES.
    void fillPsiQuad();

public:
    //! A constructor: it delegates DEData and MeshHandler costructors.
    DataProblem(SEXP Rdata, SEXP Rorder, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter, SEXP Rlambda, SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals,
      SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP Rsearch, SEXP Rmesh);

    //! A method to compute the integral of a function.
    Real FEintegrate(const VectorXr& f) const {return (R0_*f).sum();}
    //! A method to compute the integral of the square of a function.
    Real FEintegrate_square(const VectorXr& f) const {return f.dot(R0_*f);}
    //! A method to compute the integral of the exponential of a function.
    Real FEintegrate_exponential(const VectorXr& g) const;
    //! A method to compute the matrix which evaluates the basis function at the data points.
    SpMat computePsi(const std::vector<UInt>& indices) const;

    // Getters
	//! A method to access the data. It calls the same method of DEData class.
    const std::vector<Point<ndim> >& data() const {return deData_.data();}
    //! A method returning a datum. It calls the same method of DEData class.
    const Point<ndim>& data(UInt i) const {return deData_.data(i);}

    //! A method returning the number of observations. It calls the same method of DEData class.
    UInt dataSize() const {return deData_.dataSize();}
		//! A method returning the the input order. It calls the same method of DEData class.
    UInt getOrder() const {return deData_.getOrder();}
		//! A method returning the initial coefficients for the density. It calls the same method of DEData class.
    VectorXr getFvec() const {return deData_.getFvec();}
		//! A method returning a bool which says if there is a user's initial density. It calls the same method of DEData class.
    bool isFvecEmpty() const {return deData_.isFvecEmpty();}
    //! A method returning the heat diffusion process alpha parameter. It calls the same method of DEData class.
    Real getHeatStep() const {return deData_.getHeatStep();}
    //! A method returning the number of iterations for the heat diffusion process. It calls the same method of DEData class.
    UInt getHeatIter() const {return deData_.getHeatIter();}
		//! A method returning the penalization parameters. It calls the same method of DEData class.
    Real getLambda(UInt i) const {return deData_.getLambda(i);}
		//! A method returning the number of lambdas. It calls the same method of DEData class.
    UInt getNlambda()  const {return deData_.getNlambda();}
		//! A method returning the number of folds for CV. It calls the same method of DEData class.
    UInt getNfolds()  const {return deData_.getNfolds();}
		//! A method returning the number of iterations to use in the optimization algorithm. It calls the same method of DEData class.
    UInt getNsimulations() const {return deData_.getNsimulations();}
		//! A method returning the number of parameters for fixed step methods. It calls the same method of DEData class.
    UInt getNstepProposals() const {return deData_.getNstepProposals();}
		//! A method returning a parameter for fixed step methods. It calls the same method of DEData class.
    Real getStepProposals(UInt i) const {return deData_.getStepProposals(i);}
    //! A method returning the tolerance for optimization algorithm first termination criteria. It calls the same method of DEData class.
    Real getTol1() const {return deData_.getTol1();}
    //! A method returning the tolerance for optimization algorithm second termination criteria. It calls the same method of DEData class.
    Real getTol2() const {return deData_.getTol2();}
    //! A method returning the boolean print member. It calls the same method of DEData class.
    bool Print() const {return deData_.Print();}
    //! A method returning the integer that specifies the search algorithm type.
    UInt getSearch() const {return deData_.getSearch();}

    //getter for mesh
    //! A method returning the mesh.
    const MeshHandler<ORDER, mydim, ndim>& getMesh() const {return mesh_;}
    //getter for specific mesh features
    //! A method returning the number of mesh EL_NNODES. It calls the same method of MeshHandler class.
    UInt getNumNodes() const {return mesh_.num_nodes();}
    //! A method returning the number of mesh elements. It calls the same method of MeshHandler class.
    UInt getNumElements() const {return mesh_.num_elements();}
    //! A method returning a node. It calls the same method of MeshHandler class.
    Point<ndim> getPoint(Id id) const {return mesh_.getPoint(id);}
    //! A method returning an element. It calls the same method of MeshHandler class.
    Element<EL_NNODES,mydim,ndim> getElement(Id id) const {return mesh_.getElement(id);}
    //! A method returning the element in which the point in input is located. It calls the same method of MeshHandler class.
    Element<EL_NNODES,mydim,ndim> findLocation(const Point<ndim>& point) const {return mesh_.findLocation(point);}

    //getter for matrices
    //! A method returning the P matrix.
    MatrixXr getP() const {return P_;}
    //! A method returning the PsiQuad_ matrix.
    const Eigen::Matrix<Real, Integrator::NNODES, EL_NNODES>& getPsiQuad() const {return PsiQuad_;}
    //! A method returning the GlobalPsi_ matrix.
    SpMat getGlobalPsi() const {return GlobalPsi_;}
};


#include "Data_Problem_imp.h"

#endif
