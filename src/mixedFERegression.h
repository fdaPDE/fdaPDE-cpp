#ifndef __MIXEDFEREGRESSION_HPP__
#define __MIXEDFEREGRESSION_HPP__

#include "fdaPDE.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "mesh.h"
#include "param_functors.h"
#include "regressionData.h"
#include "solver.h"
#include "integratePsi.h"
#include "kronecker_product.h"
#include <memory>


/*! A base class for the smooth regression.
*/
template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
class MixedFERegressionBase
{
	protected:

	const MeshHandler<ORDER, mydim, ndim> &mesh_;
	const std::vector<Real> mesh_time_;
	const UInt N_; //!< Number of spatial basis functions.
	const UInt M_;

	const InputHandler& regressionData_;

	// For only space problems
	//  system matrix= 	|psi^T * A *psi | lambda R1^T  |   +  |psi^T * A * (-H) * psi |  O |   =  matrixNoCov + matrixOnlyCov
	//	                |     R1        | R0	      |      |         O             |  O |

	//For space time problems
	// Separable case:
	//  system matrix= 	| B^T * Ak *B + lambdaT*Ptk |  -lambdaS*R1k^T  |   +  |B^T * Ak * (-H) * B |  O |   =  matrixNoCov + matrixOnlyCov
	//	                |      -lambdaS*R1k^T       |  -lambdaS*R0k	   |      |         O          |  O |

	// Parabolic case:
	//  system matrix= 	|          B^T * Ak *B           | -lambdaS*(R1k^T+lambdaT*LR0k)  |   +  |B^T * Ak * (-H) * B |  O |   =  matrixNoCov + matrixOnlyCov
	//	                | -lambdaS*(R1k^T+lambdaT*LR0k)  |        -lambdaS*R0k	          |      |         O          |  O |

	SpMat matrixNoCov_;	//!< System matrix without

	SpMat R1_;	 //!< R1 matrix of the model
	SpMat R0_;	 //!< Mass matrix in space
	SpMat psi_;  //!< Psi matrix of the model
	MatrixXr R_; //!< R1 ^T * R0^-1 * R1


	SpMat Ptk_; 	//!< kron(Pt,IN) (separable version)
	SpMat LR0k_; 	//!< kron(L,R0) (parabolic version)
	VectorXr A_; 	//!< A_.asDiagonal() areal matrix


	MatrixXr U_;	//!< psi^T * W or psi^T * A * W padded with zeros, needed for Woodbury decomposition
	MatrixXr V_;   //!< W^T*psi, if pointwise data is U^T, needed for Woodbury decomposition

	MatrixXr barycenters_; //!< barycenter information
	VectorXi element_ids_; //!< elements id information

	Eigen::SparseLU<SpMat> matrixNoCovdec_; //!< Stores the factorization of matrixNoCov_
	Eigen::PartialPivLU<MatrixXr> Gdec_;	//!< Stores factorization of G =  C + [V * matrixNoCov^-1 * U]
	Eigen::PartialPivLU<MatrixXr> WTW_;	//!< Stores the factorization of W^T * W
	bool isWTWfactorized_ = false; 
	bool isRcomputed_ = false;
	Eigen::SparseLU<SpMat> R0dec_; //!< Stores the factorization of R0_


	VectorXr rhs_ft_correction_;	//!< right hand side correction for the forcing term:
	VectorXr rhs_ic_correction_;	//!< Initial condition correction (parabolic case)
	VectorXr _rightHandSide;      //!< A Eigen::VectorXr: Stores the system right hand side.
	MatrixXv _solution; 		//!< A Eigen::MatrixXv: Stores the system solution.
	MatrixXr _dof;       //!< A Eigen::MatrixXr storing the computed dofs
	MatrixXr _GCV;	 //!< A Eigen::MatrixXr storing the computed GCV
	UInt bestLambdaS_=0;	//!< Stores the index of the best lambdaS according to GCV
	UInt bestLambdaT_=0;	//!< Stores the index of the best lambdaT according to GCV
	Real _bestGCV=10e20;	//!< Stores the value of the best GCV
	MatrixXv _beta;		//!< A Eigen::MatrixXv storing the computed beta coefficients

	//Flag to avoid the computation of R0,R1,Psi_
	bool isAComputed = false;
	bool isPsiComputed = false;
	bool isR0Computed = false;
	bool isR1Computed = false;

	bool isSpaceVarying = false; //!< used to distinguish whether to use the forcing term u in apply() or not

	//! A member function computing the Psi matrix
	void setPsi();
	//! A method computing the no-covariates version of the system matrix
	void buildMatrixNoCov(const SpMat& NWblock,  const SpMat& SWblock,  const SpMat& SEblock);
	//! A function that given a vector u, performs Q*u efficiently
	MatrixXr LeftMultiplybyQ(const MatrixXr& u);
	//! A function which adds Dirichlet boundary conditions before solving the system ( Remark: BC for areal data are not implemented!)
	void addDirichletBC();
	//! A method which takes care of missing values setting to 0 the corresponding rows of B_
	void addNA();
 	//! A member function which builds the A vector containing the areas of the regions in case of areal data
	void setA();
	//! A member function returning the system right hand data
	void getRightHandData(VectorXr& rightHandData);
	//! A method which builds all the matrices needed for assembling matrixNoCov_
	void buildSpaceTimeMatrices();
	//! A method computing dofs in case of exact GCV, it is called by computeDegreesOfFreedom
	void computeDegreesOfFreedomExact(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT);
	//! A method computing dofs in case of stochastic GCV, it is called by computeDegreesOfFreedom
	void computeDegreesOfFreedomStochastic(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT);
	//! A method computing GCV from the dofs
	void computeGeneralizedCrossValidation(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT);

  	//! A function to factorize the system, using Woodbury decomposition when there are covariates
	void system_factorize();
	//! A function which solves the factorized system
	template<typename Derived>
	MatrixXr system_solve(const Eigen::MatrixBase<Derived>&);

	public:
	//!A Constructor.
	MixedFERegressionBase(const MeshHandler<ORDER,mydim,ndim>& mesh, const InputHandler& regressionData):
			mesh_(mesh), N_(mesh.num_nodes()), M_(1), regressionData_(regressionData), _dof(regressionData.getDOF_matrix()){};
	MixedFERegressionBase(const MeshHandler<ORDER,mydim,ndim>& mesh, const std::vector<Real>& mesh_time, const InputHandler& regressionData):
  		mesh_(mesh), mesh_time_(mesh_time), N_(mesh.num_nodes()), M_(regressionData.getFlagParabolic()? mesh_time.size()-1 : mesh_time.size()+SPLINE_DEGREE-1), regressionData_(regressionData), _dof(regressionData.getDOF_matrix()){};

	//! The function solving the system, used by the children classes. Saves the result in _solution
	/*!
	    \param oper an operator, which is the Stiffness operator in case of Laplacian regularization
	    \param u the forcing term, will be used only in case of anysotropic nonstationary regression
	*/
	//! A method which builds all the space matrices
	template<typename A>
	void apply(EOExpr<A> oper,const ForcingTerm & u);

	//! A member function computing the dofs for external calls
	//template<typename A>
	//void computeDegreesOfFreedom(EOExpr<A> oper);

	//! A method computing the dofs
	void computeDegreesOfFreedom(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT);
	//! A method that set WTW flag to false, in order to recompute the matrix WTW.
	inline void recomputeWTW(){ this->isWTWfactorized_ = false;}
	//! A function returning the computed barycenters of the locationss
	inline MatrixXr const & getBarycenters() const{return barycenters_;};
	//! A function returning the element ids of the locations
	inline VectorXi const & getElementIds() const{return element_ids_;};
	//! A inline member that returns a VectorXr, returns the whole solution_.
	inline MatrixXv const & getSolution() const{return _solution;}
	//! A function returning the computed dofs of the model
	inline MatrixXr const & getDOF() const{return _dof;}
	//! A method returning the computed GCV of the model
	inline MatrixXr const & getGCV() const{return _GCV;}
	//! A method returning the computed beta coefficients of the model
	inline MatrixXv const & getBeta() const{return _beta;}
	//! A method returning the index of the best lambdaS according to GCV
	inline UInt getBestLambdaS(){return bestLambdaS_;}
	//! A method returning the index of the best lambdaT according to GCV
	inline UInt getBestLambdaT(){return bestLambdaT_;}
	//! A method returning the psi matrix
	inline SpMat const getPsi()const{return psi_;}
	//! A method returning the R0 matrix
	inline SpMat const getR0()const{return R0_;}
	//! A method returning the R1 matrix
	inline SpMat const getR1()const{return R1_;}
	


};

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
class MixedFERegression : public MixedFERegressionBase<InputHandler,IntegratorSpace,ORDER,IntegratorTime,SPLINE_DEGREE,ORDER_DERIVATIVE,mydim,ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER,mydim,ndim>& mesh, const InputHandler& regressionData):
			MixedFERegressionBase<InputHandler,IntegratorSpace,ORDER,IntegratorTime,SPLINE_DEGREE,ORDER_DERIVATIVE,mydim,ndim>(mesh, regressionData){};
	MixedFERegression(const MeshHandler<ORDER,mydim,ndim>& mesh, const std::vector<Real>& mesh_time, const InputHandler& regressionData):
			MixedFERegressionBase<InputHandler,IntegratorSpace,ORDER,IntegratorTime,SPLINE_DEGREE,ORDER_DERIVATIVE,mydim,ndim>(mesh, mesh_time, regressionData){};

	void apply()
	{
		#ifdef R_VERSION_
			Rprintf("Option not implemented! \n");
		#endif
		std::cout << "Option not implemented! \n";
	}

};


//! A class for the construction of the temporal matrices needed for the parabolic case
template<typename InputHandler, typename Integrator, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE>
class MixedSplineRegression
{
	private:
		const std::vector<Real>& mesh_time_;
		const InputHandler& regressionData_;

		SpMat phi_;   //! Matrix of the evaluations of the spline basis functions in the time locations
		SpMat Pt_;
		SpMat timeMass_; //! Mass matrix in time

	public:
		MixedSplineRegression(const std::vector<Real>& mesh_time, const InputHandler& regressionData):mesh_time_(mesh_time), regressionData_(regressionData){};

    void setPhi();
		void setTimeMass();
    void smoothSecondDerivative();

		inline SpMat const & getPt() const { return Pt_; }
		inline SpMat const & getPhi() const { return phi_; }
		inline SpMat const & getTimeMass() const { return timeMass_; }

};

//! A class for the construction of the temporal matrices needed for the separable case
template<typename InputHandler>
class MixedFDRegression
{
	private:
		const std::vector<Real>& mesh_time_;
		const InputHandler& regressionData_;

		SpMat derOpL_; //!matrix associated with derivation in time

	public:
		MixedFDRegression(const std::vector<Real>& mesh_time, const InputHandler& regressionData):mesh_time_(mesh_time), regressionData_(regressionData){};

    void setDerOperator(); //! sets derOpL_
		inline SpMat const & getDerOpL() const { return derOpL_; }

};

#include "mixedFERegression_imp.h"

#endif
