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

// This file contains data informations for the Density Estimation problem

//! @brief A class to store common data for the problem.
template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
class DataProblem{
  private:
    DEData<ndim> deData_;
    MeshHandler<ORDER, mydim, ndim> mesh_;
    SpMat R0_, R1_, GlobalPsi_;
    MatrixXr P_, PsiQuad_;
    static constexpr UInt Nodes = (mydim==2) ? 3*ORDER : 6*ORDER-2;

    //! A method to compute the finite element matrices.
    void fillFEMatrices();
    //! A method to compute the matrix which evaluates the basis function at the quadrature nodes.
    void fillPsiQuad();

  public:
    //! A constructor: it delegates DEData and MeshHandler costructors.
    DataProblem(SEXP Rdata, SEXP Rorder, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter, SEXP Rlambda, SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals,
      SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP Rsearch, SEXP Rmesh);

    //! A method to compute the integral of a function.
    inline Real FEintegrate(const VectorXr& f) const {return (R0_*f).sum();}
    //! A method to compute the integral of the square of a function.
    inline Real FEintegrate_square(const VectorXr& f) const {return f.dot(R0_*f);}
    //! A method to compute the integral of the exponential of a function.
    Real FEintegrate_exponential(const VectorXr& g) const;
    //! A method to compute the matrix which evaluates the basis function at the data points.
    SpMat computePsi(const std::vector<UInt>& indices) const;

    // Getters
		//! A method returning the data. It calls the same method of DEData class.
		inline std::vector<Point<ndim> > getData() const {return deData_.getData();}
    //! A method returning a datum. It calls the same method of DEData class.
    inline Point<ndim> getDatum(UInt i) const {return deData_.getDatum(i);}
    //! A method returning the number of observations. It calls the same method of DEData class.
		inline UInt getNumberofData() const {return deData_.getNumberofData();}
		//! A method returning the the input order. It calls the same method of DEData class.
		inline UInt getOrder() const {return deData_.getOrder();}
		//! A method returning the initial coefficients for the density. It calls the same method of DEData class.
		inline VectorXr getFvec() const {return deData_.getFvec();}
		//! A method returning a bool which says if there is a user's initial density. It calls the same method of DEData class.
		inline bool isFvecEmpty() const {return deData_.isFvecEmpty();}
    //! A method returning the heat diffusion process alpha parameter. It calls the same method of DEData class.
    inline Real getHeatStep() const {return deData_.getHeatStep();}
    //! A method returning the number of iterations for the heat diffusion process. It calls the same method of DEData class.
    inline UInt getHeatIter() const {return deData_.getHeatIter();}
		//! A method returning the penalization parameters. It calls the same method of DEData class.
		inline Real getLambda(UInt i) const {return deData_.getLambda(i);}
		//! A method returning the number of lambdas. It calls the same method of DEData class.
		inline UInt getNlambda()  const {return deData_.getNlambda();}
		//! A method returning the number of folds for CV. It calls the same method of DEData class.
		inline UInt getNfolds()  const {return deData_.getNfolds();}
		//! A method returning the number of iterations to use in the optimization algorithm. It calls the same method of DEData class.
		inline UInt getNsimulations() const {return deData_.getNsimulations();}
		//! A method returning the number of parameters for fixed step methods. It calls the same method of DEData class.
		inline UInt getNstepProposals() const {return deData_.getNstepProposals();}
		//! A method returning a parameter for fixed step methods. It calls the same method of DEData class.
		inline Real getStepProposals(UInt i) const {return deData_.getStepProposals(i);}
    //! A method returning the tolerance for optimization algorithm first termination criteria. It calls the same method of DEData class.
    inline Real getTol1() const {return deData_.getTol1();}
    //! A method returning the tolerance for optimization algorithm second termination criteria. It calls the same method of DEData class.
    inline Real getTol2() const {return deData_.getTol2();}
    //! A method returning the boolean print member. It calls the same method of DEData class.
    inline bool Print() const {return deData_.Print();}
    //! A method returning the integer that specifies the search algorithm type.
    inline UInt getSearch() const {return deData_.getSearch();}

    //getter for mesh
    //! A method returning the mesh.
    inline const MeshHandler<ORDER, mydim, ndim>& getMesh() const {return mesh_;}
    //getter for specific mesh features
    //! A method returning the number of mesh nodes. It calls the same method of MeshHandler class.
    inline UInt getNumNodes() const {return mesh_.num_nodes();}
    //! A method returning the number of mesh elements. It calls the same method of MeshHandler class.
    inline UInt getNumElements() const {return mesh_.num_elements();}
    //! A method returning a node. It calls the same method of MeshHandler class.
    inline Point<ndim> getPoint(Id id) const {return mesh_.getPoint(id);}
    //! A method returning an element. It calls the same method of MeshHandler class.
    inline Element<Nodes,mydim,ndim> getElement(Id id) const {return mesh_.getElement(id);}
    //! A method returning the element in which he point in input is located by using a naive search. It calls the same method of MeshHandler class.
    inline Element<Nodes,mydim,ndim> findLocationNaive(Point<ndim> point) const {return mesh_.findLocationNaive(point);}
    //! A method returning the element in which he point in input is located by using a tree search. It calls the same method of MeshHandler class.
    inline Element<Nodes,mydim,ndim> findLocationTree(Point<ndim> point) const {return mesh_.findLocationTree(point);}

    //getter for matrices
    //! A method returning the P matrix.
    inline MatrixXr getP() const {return P_;}
    //! A method returning the PsiQuad_ matrix.
    inline MatrixXr getPsiQuad() const {return PsiQuad_;}
    //! A method returning the GlobalPsi_ matrix.
    inline SpMat getGlobalPsi() const {return GlobalPsi_;}
};


#include "Data_Problem_imp.h"

#endif
