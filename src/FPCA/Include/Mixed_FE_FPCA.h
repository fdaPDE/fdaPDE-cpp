#ifndef __MIXED_FE_FPCA_H__
#define __MIXED_FE_FPCA_H__

#include "../../FdaPDE.h"
#include "../../FE_Assemblers_Solvers/Include/Finite_Element.h"
#include "../../FE_Assemblers_Solvers/Include/Matrix_Assembler.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../FE_Assemblers_Solvers/Include/Param_Functors.h"
#include "../../FE_Assemblers_Solvers/Include/Solver.h"
#include "FPCA_Data.h"
#include "FPCA_Object.h"
#include "../../Global_Utilities/Include/Solver_Definitions.h"
#include <memory>

//! A virtual Base class for the implementation of the SF-PCA algorithm: it contains some methods useful for the construction and the resolution of the linear system that has to be done iteratively in the algorithm.

class MixedFEFPCABase
{
protected:
	//!A FPCAData object containing the data passed from R.
	const FPCAData& fpcaData_;
	std::vector<coeff> tripletsData_;

	//!A Eigen::SparseMatrix<Real> containing the basis function.
	//SpMat psi_;

	MatrixXr barycenters_; //barycenter information
	VectorXi element_ids_; //elements id information

	bool isRcomputed_;
	MatrixXr R_; //R1 ^T * R0^-1 * R1

	//SpMat NWblock_;
	//!A Eigen::SparseMatrix<Real> containing the NorthWest block of the linear system.
	SpMat DMat_;
	//!A Eigen::SparseMatrix<Real> containing the anti-diagonal block of the linear system.
	SpMat AMat_;
	//!A Eigen::SparseMatrix<Real> containing the SouthEast block of the linear system.
	SpMat MMat_;

	//!A Eigen::SparseMatrix<Real> containing the basis function.
	SpMat Psi_;

	//! A Eigen::VectorXr: Stores the area/volume of each region
	VectorXr Delta_; //Delta_.asDiagonal() = diag(|D_1|,...,|D_N|)

	//!A Eigen::SparseMatrix<Real> containing the matrix of the linear system
	SpMat coeffmatrix_;       //!A Eigen::VectorXr: Stores the system right hand side.
	VectorXr b_;			  //!A Eigen::VectorXr: Stores the system solution.
	std::vector<VectorXr> solution_;

	Sparse_LU sparseSolver_;
	//!A Real : Stores the variance of the edf computation using the stochastic method.
	std::vector<Real> var_;

	UInt nnodes_;
	
	//!A Eigen::VectorXr : Stores the final scores computed for each PC.
	std::vector<VectorXr> scores_mat_;
	//!A Eigen::VectorXr : Stores the final loadings computed for each PC.
	std::vector<VectorXr> loadings_mat_;
	//!A std::vector<Real> : Stores the lambda selected by the cross-validation method computed for each PC.
	std::vector<Real> lambda_PC_;
	//!A std::vector<Real> : Stores the variance explained by each PC.
	std::vector<Real> variance_explained_;
	//!A std::vector<Real> : Stores the cumulative percentage explained by each PC.
	std::vector<Real> cumsum_percentage_;
	//!A Eigen::MatrixXr : Stores the datamatrix and its updates after each PC is computed.
	MatrixXr datamatrixResiduals_ ;

	//! A method for the computation of Delta_
	template<UInt ORDER, UInt mydim, UInt ndim>
	void computeDelta(const MeshHandler<ORDER, mydim, ndim> & mesh);
	//! A method for the computation of the basis function matrix Psi_.
	template<UInt ORDER, UInt mydim, UInt ndim>
	void computeBasisEvaluations(const MeshHandler<ORDER, mydim, ndim> & mesh);
	//! A method for the assembling of the matrix of the system.
	void buildCoeffMatrix(const SpMat& DMat,  const SpMat& AMat,  const SpMat& MMat);
	//! A method for the computation of the NW block of the system.
	void computeDataMatrix(SpMat& DMat);
	//! A method for the computation of the NW block of the system when there are NAs in the datamatrix.
	void computeDataMatrixByIndices(SpMat& DMat);
	//! A method for the computation of the right hand side of the system.
	void computeRightHandData(VectorXr& rightHandData,FPCAObject& FPCAinput);
	//! A method for the computation of the variance explained by the PCs.
	void computeVarianceExplained();
	//! A method for the computation of the cumulative percentage of the variance explained by the PCs.
	void computeCumulativePercentageExplained();
	//! A method for the computation of the iterations of the SF-PCA algorithm.
	void computeIterations(MatrixXr & datamatrixResiduals_,FPCAObject & FPCAinput, UInt lambda_index, UInt nnodes);
	//! A method for the initialization of all the parameters used in the iteration of the SF-PCA algorithm.
	void SetAndFixParameters();

public:
	//!A Constructor.
	MixedFEFPCABase(const FPCAData& fpcaData): fpcaData_(fpcaData),isRcomputed_(false) {};

	//! A method for the initialization of all the parameters used in the iteration of the SF-PCA algorithm.
	template<UInt ORDER, UInt mydim, UInt ndim>
	void SetAndFixParameters(const MeshHandler<ORDER, mydim, ndim> & mesh);

	//!A destructor.
	virtual ~MixedFEFPCABase(){};

	//!A virtual method that will be called for performing the algorithm. This method is specified in the children classes.
	virtual  void apply()=0;

	//! A inline member that returns a VectorXr, returns the whole solution_.
	inline std::vector<VectorXr> const & getSolution() const{return solution_;};
	//! A method returning a reference to the scores matrix.
	inline std::vector<VectorXr> const & getScoresMat() const {return scores_mat_;}
	//! A method returning a reference to the loadings matrix.
	inline std::vector<VectorXr> const & getLoadingsMat() const {return loadings_mat_;}
	//! A method returning a reference to the vector of lambda taken for each PC.
	inline std::vector<Real> const & getLambdaPC() const {return lambda_PC_;}
	//! A method returning a reference to the vector of variance explained for each PC.
	inline std::vector<Real> const & getVarianceExplained() const {return variance_explained_;}
	//! A method returning a reference to the vector of the percentage explained cumulatively by the first N PC.
	inline std::vector<Real> const & getCumulativePercentage() const {return cumsum_percentage_;}
	//! A method returning the vector of the variance of the stochastic computation of the edf.
	inline std::vector<Real> const & getVar() const{return var_;};
	//! A function returning the computed barycenters of the locationss
	inline MatrixXr const & getBarycenters() const{return barycenters_;};
	//! A function returning the element ids of the locations
	inline VectorXi const & getElementIds() const{return element_ids_;};
};


//! A class for the implementation of the SF-PCA algorithm without the use of a cross-validation method for the selection of the parameter lambda
class MixedFEFPCA : public MixedFEFPCABase
{
public:
	//!A constructor.
	MixedFEFPCA(const FPCAData& fpcaData):MixedFEFPCABase(fpcaData){};

	//!A Destructor.
	virtual ~MixedFEFPCA(){};

	//!A specification of the virtual method for performing the SF-PCA algorithm with no cross-validation method for the choice of the parameter lambda.
	void apply();
};

//! A class for the implementation of the SF-PCA algorithm with the GCV used as cross-validation method for the selection of the parameter lambda
class MixedFEFPCAGCV: public MixedFEFPCABase
{
protected:
	std::vector<VectorXr> loadings_lambda_;
	std::vector<VectorXr> scores_lambda_;
	//! A std::vector<Real>: Stores the degrees of freedom for each lambda
	std::vector<Real> dof_;
	//! A std::vector<Real>: Stores the GCV computed for each lambda
	std::vector<Real> GCV_;
	//! A method for the selection of the exact or stochastic way to compute the degrees of freedom
	void computeDegreesOfFreedom(UInt output_index, Real lambda);
	//! A method for the exact computation of the degres of freedom for a specific lambda
	void computeDegreesOfFreedomExact(UInt output_index, Real lambda);
	//! A method for the stochastic computation of the degres of freedom for a specific lambda
	void computeDegreesOfFreedomStochastic(UInt output_index, Real lambda);

	//! A method for the computation of the iterations of the SF-PCA algorithm
	void computeIterationsGCV(MatrixXr &datamatrixResiduals_, UInt nnodes, UInt np);
	//! A method for the computation of the GCV
	void computeGCV(FPCAObject& FPCAinput,UInt output_index);

public:
	//!A Constructor.
	MixedFEFPCAGCV(const FPCAData& fpcaData):MixedFEFPCABase(fpcaData){};

	//!A Destructor.
	virtual ~MixedFEFPCAGCV(){};

	//!A specification of the virtual method for performing the SF-PCA algorithm with GCV cross-validation method for the choice of the parameter lambda.
	void apply();
	//!A method returning the degrees of freedom of the problem.
	inline std::vector<Real> const & getDOF() const{return dof_;};

};

//! A class for the implementation of the SF-PCA algorithm with the K-Fold used as cross-validation method for the selection of the parameter lamb
class MixedFEFPCAKFold : public MixedFEFPCABase
{
protected:
	std::vector<Real> KFold_;
	//! A UInt specifying the number of folds to use in the algorithm
	UInt nFolds;
	//!A method for computing and validating the K-folds
	void computeKFolds(MatrixXr & datamatrixResiduals_, UInt lambda_index, UInt nnodes,UInt nFolds);
public:
	//!A Constructor.
	MixedFEFPCAKFold(const FPCAData& fpcaData):MixedFEFPCABase(fpcaData){};

	//!A Destructor.
	virtual ~MixedFEFPCAKFold(){};

	//!A specification of the virtual method for performing the SF-PCA algorithm with the K-Fold cross-validation method for the choice of the parameter lambda.
	void apply();
};

#include "Mixed_FE_FPCA_imp.h"

#endif
