#ifndef __FPCA_OBJECT_H__
#define __FPCA_OBJECT_H__

#include "../../FdaPDE.h"
#include "../../Mesh/Include/Mesh_Objects.h"
#include "../../FE_Assemblers_Solvers/Include/Param_Functors.h"

//!  A class for handling the update of the loadings, the scores and the observation data at each iteration
/*!
 * This class job is to handle the update of the loadings and the scores at each
 * iteration during the SF-PCA algorithm.
*/

class  FPCAObject
{
	private:


		//Loadings and scores estimation
		VectorXr scores_;
		VectorXr loadings_;

		//Regression data
		VectorXr ObservationData_;

	public:

		FPCAObject(){};

		//!A constructor of the FPCAObject. As parameter it only takes the reference to the datamatrix X. In the constructor, the SVD decomposition of the matrix is performed and the loadings and the scores vector are initialized.
		explicit FPCAObject(const MatrixXr& datamatrix_);

		//!A method for updating the Scores vector given the datamatrix
		void setScores(const MatrixXr& datamatrix_);
		//!A method for updating the ObservationData vector given the datamatrix
		void setObservationData(const MatrixXr& datamatrix_);
		//void setObservationData(const MatrixXr& datamatrix_, const SpMat& psi_);
		//!A method for updating the Loadings vector given the solution of the linear system and the basis function matrix Psi
		void setLoadingsPsi(UInt nnodes, const VectorXr& f_sol, const SpMat& psi);
		//!A method for updating the Loadings vector given the solution of the linear system and the indices of the non NA values
		void setLoadings(UInt nnodes, const VectorXr& f_sol, const std::vector<UInt>& obs_indices);
		//!A method for finalize the Loadings vector
		void finalizeLoadings(const std::vector<UInt>& obs_indices, UInt nlocations);



		//!A method for printing the Scores vector on a specified output
		void printScores(std::ostream & out) const;
		//!A method for printing the Loadings vector on a specified output
		void printLoadings(std::ostream & out) const;
		//!A method for printing the ObservationData on a specified output
		void printObservationData(std::ostream & out) const;
		//void printDatamatrix(std::ostream & out) const;


		//! A method returning a reference to the scores vector
		inline VectorXr const & getScores() const {return scores_;}
		//! A method returning a reference to the loadings vector
		inline VectorXr const & getLoadings() const {return loadings_;}
		//! A method returning a reference to the observation data vector
		inline VectorXr const & getObservationData() const {return ObservationData_;}
};

#endif
