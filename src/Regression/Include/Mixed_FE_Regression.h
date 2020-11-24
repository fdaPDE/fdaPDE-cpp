#ifndef __MIXED_FE_REGRESSION_H__
#define __MIXED_FE_REGRESSION_H__

#include <memory>
#include <type_traits>
#include "../../FdaPDE.h"
#include "../../FE_Assemblers_Solvers/Include/Finite_Element.h"
#include "../../FE_Assemblers_Solvers/Include/Spline.h"
#include "../../FE_Assemblers_Solvers/Include/Kronecker_Product.h"
#include "../../FE_Assemblers_Solvers/Include/Matrix_Assembler.h"
#include "../../FE_Assemblers_Solvers/Include/Param_Functors.h"
#include "../../FE_Assemblers_Solvers/Include/Solver.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../Lambda_Optimization/Include/Optimization_Data.h"
#include "Regression_Data.h"

/*! A base class for the smooth regression.
*/
template<typename InputHandler>
class MixedFERegressionBase
{
	protected:
		const std::vector<Real> mesh_time_;
		const UInt N_; 			//!< Number of spatial basis functions.
		const UInt M_;

		const InputHandler & regressionData_;
                OptimizationData & optimizationData_; //!<COnst reference to OptimizationData class
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

		SpMat 		matrixNoCov_;	//!< System matrix without
		SpMat 		DMat_;
		SpMat 		R1_;		//!< R1 matrix of the model
		SpMat 		R0_;	 	//!< Mass matrix in space
		SpMat 		R0_lambda;
		SpMat 		R1_lambda;
		SpMat 		psi_;  		//!< Psi matrix of the model
		SpMat 		psi_t_;  	//!< Psi ^T matrix of the model
		SpMat 		Ptk_; 		//!< kron(Pt,IN) (separable version)
		SpMat 		LR0k_; 		//!< kron(L,R0) (parabolic version)
		MatrixXr 	R_; 		//!< R1 ^T * R0^-1 * R1
		MatrixXr 	H_; 		//!< The hat matrix of the regression
		MatrixXr	Q_; 		//!< Identity - H, projects onto the orthogonal subspace
		VectorXr 	A_; 		//!< A_.asDiagonal() areal matrix
		MatrixXr 	U_;		//!< psi^T * W or psi^T * A * W padded with zeros, needed for Woodbury decomposition
		MatrixXr 	V_;  		//!< W^T*psi, if pointwise data is U^T, needed for Woodbury decomposition
		MatrixXr 	barycenters_; 	//!< barycenter information
		VectorXi 	element_ids_; 	//!< elements id information

		// Factorizations
		Eigen::SparseLU<SpMat> matrixNoCovdec_; //!< Stores the factorization of matrixNoCov_
		//std::unique_ptr<Eigen::PartialPivLU<MatrixXr>>  matrixNoCovdec_{new Eigen::PartialPivLU<MatrixXr>}; //!< Stores the factorization of matrixNoCov_
		Eigen::PartialPivLU<MatrixXr> Gdec_;	//!< Stores factorization of G =  C + [V * matrixNoCov^-1 * U]

		Eigen::PartialPivLU<MatrixXr> WTW_;	//!< Stores the factorization of W^T * W
		bool isWTWfactorized_ = false;
		bool isRcomputed_ = false;
		Eigen::SparseLU<SpMat> R0dec_; 		//!< Stores the factorization of R0_

		VectorXr rhs_ft_correction_;	//!< right hand side correction for the forcing term:
		VectorXr rhs_ic_correction_;	//!< Initial condition correction (parabolic case)
		VectorXr _rightHandSide;      	//!< A Eigen::VectorXr: Stores the system right hand side.
		MatrixXv _solution; 		//!< A Eigen::MatrixXv: Stores the system solution.
		MatrixXr _dof;      		//!< A Eigen::MatrixXr storing the computed dofs
		MatrixXr _GCV;			//!< A Eigen::MatrixXr storing the computed GCV
		MatrixXv _beta;			//!< A Eigen::MatrixXv storing the computed beta coefficients

		//Flag to avoid the computation of R0, R1, Psi_ onece already performed
		bool isAComputed   = false;
		bool isPsiComputed = false;
		bool isR0Computed  = false;
		bool isR1Computed  = false;

		bool isSpaceVarying = false; //!< used to distinguish whether to use the forcing term u in apply() or not
		bool isGAMData;

	        // -- SETTERS --
		template<UInt ORDER, UInt mydim, UInt ndim>
	    void setPsi(const MeshHandler<ORDER, mydim, ndim> & mesh_);
		//! A method computing the no-covariates version of the system matrix
		void buildMatrixNoCov(const SpMat & NWblock, const SpMat & SWblock,  const SpMat & SEblock);

		//! A function which adds Dirichlet boundary conditions before solving the system ( Remark: BC for areal data are not implemented!)
		void addDirichletBC();
		//! A function which adds Dirichlet boundary conditions only to MatrixnoCov( Remark: BC for areal data are not implemented!)
		void addDirichletBC_matrix();
		//! A method which takes care of missing values setting to 0 the corresponding rows of B_
		void addNA();
	 	//! A member function which builds the A vector containing the areas of the regions in case of areal data
	    template<UInt ORDER, UInt mydim, UInt ndim>
		void setA(const MeshHandler<ORDER, mydim, ndim> & mesh_);
		//! A member function which sets psi_t_
		void setpsi_t_(void);
	        //! A member function which builds DMat, to be changed in apply for the temporal case
		void setDMat(void);
		//! A member function which builds the Q matrix
		void setQ(void);
		//! A member function which builds the H matrix
		void setH(void);
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

		// -- BUILD SYSTEM --
		 //! Spatial version
		void buildSystemMatrix(Real lambda);
		//! Space-time version
		void buildSystemMatrix(Real lambdaS, Real lambdaT);

		// -- FACTORIZER --
	  	//! A function to factorize the system, using Woodbury decomposition when there are covariates
		void system_factorize();

		// -- SOLVER --
		//! A function which solves the factorized system
		template<typename Derived>
		MatrixXr system_solve(const Eigen::MatrixBase<Derived>&);

	public:
		//!A Constructor.
		MixedFERegressionBase( const InputHandler & regressionData, OptimizationData & optimizationData,  UInt nnodes_) :
			N_(nnodes_), M_(1), regressionData_(regressionData), optimizationData_(optimizationData), _dof(optimizationData.get_DOF_matrix()){isGAMData = regressionData.getisGAM();};

		MixedFERegressionBase(const std::vector<Real> & mesh_time, const InputHandler & regressionData, OptimizationData & optimizationData, UInt nnodes_, UInt spline_degree) :
			mesh_time_(mesh_time), N_(nnodes_), M_(regressionData.getFlagParabolic() ? mesh_time.size()-1 : mesh_time.size()+spline_degree-1),
			regressionData_(regressionData), optimizationData_(optimizationData), _dof(optimizationData.get_DOF_matrix()){isGAMData = regressionData.getisGAM();};

		//! A member function computing the dofs for external calls
		//template<typename A>
		//void computeDegreesOfFreedom(EOExpr<A> oper);

		// -- UTILITIES --
		//! A method computing the dofs
		void computeDegreesOfFreedom(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT);
		//! A method that set WTW flag to false, in order to recompute the matrix WTW.
		inline void recomputeWTW(void){ this->isWTWfactorized_ = false;}

		// -- GETTERS --
		//! A function returning the computed barycenters of the locationss
		inline MatrixXr const & getBarycenters(void) const {return barycenters_;}; //returns a const reference as in rergressionData
		//! A function returning the element ids of the locations
		inline VectorXi const & getElementIds(void) const {return element_ids_;};
		//! A inline member that returns a VectorXr, returns the whole solution_.
		inline MatrixXv const & getSolution(void) const {return _solution;}
		//! A function returning the computed dofs of the model
		inline MatrixXr const & getDOF(void) const { if (optimizationData_.get_DOF_matrix().rows()!=0 && optimizationData_.get_DOF_matrix().cols()!=0)
                                                             	return  optimizationData_.get_DOF_matrix();
							     else return this->_dof;}
		//! A method returning the computed GCV of the model
		inline MatrixXr const & getGCV(void) const {return _GCV;}
		//! A method returning the computed beta coefficients of the model
		inline MatrixXv const & getBeta(void) const {return _beta;}
		//! A method returning the psi matrix
		inline const SpMat * getpsi_(void) const {return &this->psi_;}
		//! A method returning the psi matrix transposed
		inline const SpMat * getpsi_t_(void) const {return &this->psi_t_;}
		//! A method returning the R0 matrix
		inline const SpMat * getR0_(void) const {return &this->R0_;}
		//! A method returning the R1 matrix
		inline const SpMat * getR1_(void) const {return &this->R1_;}
		//! A method returning the DMat matrix, da implementare la DMat
		inline const SpMat * getDMat_(void) const {return &this->DMat_;}
		//! A method returning the Q_ matrix -> da impementare la Q
		inline const MatrixXr *	getQ_(void) const {return &this->Q_;}
		//! A method returning the H_ matrix da implementare la H
		inline const MatrixXr *	getH_(void) const {return &this->H_;}
		//! A method returning the A_ matrix
		inline const VectorXr *	getA_(void) const {return &this->A_;}
		//! A method returning the rhs
		inline const VectorXr *	getrhs_(void) const {return &this->_rightHandSide;}
		//! A method returning the forcing term
		inline const VectorXr *	getu_(void) const {return &this->rhs_ft_correction_;}
		//! A method returning the number of nodes of the mesh
		inline UInt getnnodes_(void) const {return this->N_;}
		inline bool isSV(void) const {return this->isSpaceVarying;}

		//! A function that given a vector u, performs Q*u efficiently
		MatrixXr LeftMultiplybyQ(const MatrixXr & u);

		// -- APPLY --
		//! The function solving the system, used by the children classes. Saves the result in _solution
		/*!
		    \param oper an operator, which is the Stiffness operator in case of Laplacian regularization
		    \param u the forcing term, will be used only in case of anysotropic nonstationary regression
		*/
		//! A method which builds all the space matrices
		template<UInt ORDER, UInt mydim, UInt ndim, typename A>
		void preapply(EOExpr<A> oper, const ForcingTerm & u, const MeshHandler<ORDER, mydim, ndim> & mesh_ );

		MatrixXv apply(void);
		MatrixXr apply_to_b(const MatrixXr & b);
};

//----------------------------------------------------------------------------//

template<typename InputHandler>
class MixedFERegression : public MixedFERegressionBase<InputHandler>
{
	public:
		MixedFERegression(const InputHandler & regressionData,  OptimizationData & optimizationData, UInt nnodes_):
			MixedFERegressionBase<InputHandler>(regressionData, optimizationData, nnodes_) {};
		MixedFERegression(const std::vector<Real> & mesh_time, const InputHandler & regressionData,  OptimizationData & optimizationData, UInt nnodes_, UInt spline_degree):
			MixedFERegressionBase<InputHandler>(mesh_time, regressionData, optimizationData, nnodes_, spline_degree) {};

		void apply(void)
		{
			Rprintf("Option not implemented!\n");
		}
};

//! A class for the construction of the temporal matrices needed for the parabolic case
template<typename InputHandler>
class MixedSplineRegression
{
	private:
		const std::vector<Real> & mesh_time_;
		const InputHandler & regressionData_;

		SpMat phi_;   //!< Matrix of the evaluations of the spline basis functions in the time locations
		SpMat Pt_;
		SpMat timeMass_; //!< Mass matrix in time

	public:
		static constexpr UInt SPLINE_DEGREE=3;
		static constexpr UInt ORDER_DERIVATIVE=2;

		MixedSplineRegression(const std::vector<Real> & mesh_time, const InputHandler & regressionData):
			mesh_time_(mesh_time), regressionData_(regressionData) {};

		void setPhi(void);
		void setTimeMass(void);
		void smoothSecondDerivative(void);

		inline const SpMat & getPt(void) const {return Pt_;}
		inline const SpMat & getPhi(void) const {return phi_;}
		inline const SpMat & getTimeMass(void) const {return timeMass_;}

};

//! A class for the construction of the temporal matrices needed for the separable case
template<typename InputHandler>
class MixedFDRegression
{
	private:
		const std::vector<Real> & mesh_time_;
		const InputHandler & regressionData_;

		SpMat derOpL_; //!< matrix associated with derivation in time

	public:
		MixedFDRegression(const std::vector<Real> & mesh_time, const InputHandler & regressionData):
			mesh_time_(mesh_time), regressionData_(regressionData) {};

    	void setDerOperator(void); //!< sets derOpL_
		inline const SpMat & getDerOpL(void) const {return derOpL_;}

};

#include "Mixed_FE_Regression_imp.h"

#endif
