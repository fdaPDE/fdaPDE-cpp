#ifndef __FPCA_DATA_H__
#define __FPCA_DATA_H__

#include "../../FdaPDE.h"
#include "../../Mesh/Include/Mesh_Objects.h"
#include "../../FE_Assemblers_Solvers/Include/Param_Functors.h"

//!  An IO handler class for objects passed from R
/*!
 * This class, given the data from R, convert them in a C++ format, offering a
 * series of method for their access, so isolating the more the possible the specific
 * code for R/C++ data conversion.
*/
class  FPCAData{
	private:

		const RNumericMatrix locations_;

		//barycenter information
		MatrixXr barycenters_; //barycenter information
		VectorXi element_ids_; //elements id information
		bool locations_by_barycenter_;

		//Design matrix
		MatrixXr datamatrix_;

		UInt order_;

		//Areal data
		MatrixXi incidenceMatrix_;
		//lambda
		std::vector<Real> lambda_;

		//Number of Principal Components
		UInt nPC_;

		//Number of Folds for KFold
		UInt nFolds_;

		//Parameters for better GCV timings
		UInt GCVmethod_;
		UInt nrealizations_;      // Number of relizations for the stochastic estimation of GCV


		std::vector<UInt> observations_indices_;
		UInt n_;
		UInt p_;

		UInt nRegions_;

		bool locations_by_nodes_;
		UInt search_;


		void setDatamatrix(SEXP Rdatamatrix);
		void setBaryLocations(SEXP RbaryLocations);
		void setNrealizations(SEXP Rnrealizations);
		void setIncidenceMatrix(SEXP RincidenceMatrix);

	public:
		//! A basic version of the constructor.

		/*!
			It initializes the object storing the R given objects. This is the simplest of the two possible interfaces with R
			\param Rdatamatrix an R-matrix containing the datamatrix.

			\param Rsearch an R-integer to decide the search algorithm type (tree or naive or walking search algorithm).

			\param Rlocations an R-matrix containing the location of the observations.

			\param Rorder an R-integer containing the order of the approximating basis.

			\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions in the model with areal data

			\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.

			\param RnPC an R-integer specifying the number of principal components to compute.

			\param RnFolds an R-integer specifying the number of folds to use if K-Fold cross validation method is chosen.

			\param RGCVmethod an R-integer specifying if the GCV computation has to be exact(if = 1) or stochastic (if = 2).

			\param Rnrealizations an R-integer specifying the number of realizations to use when computing the GCV stochastically.

		*/

		// FPCAData(){};

		explicit FPCAData(SEXP Rlocations, SEXP RbaryLocations, SEXP Rdatamatrix, SEXP Rorder, SEXP RincidenceMatrix,
		SEXP Rlambda, SEXP RnPC, SEXP RnFolds, SEXP RGCVmethod, SEXP Rnrealizations, SEXP Rsearch);


		explicit FPCAData(Real* locations, UInt n_locations, UInt ndim, MatrixXr& datamatrix,
		UInt order, MatrixXi& incidenceMatrix, std::vector<Real> lambda, UInt nPC, UInt nFolds, UInt search);


		void printDatamatrix(std::ostream & out) const;
		void printLocations(std::ostream & out) const;
		void printIncidenceMatrix(std::ostream & out) const;

		//! A method returning the locations of the observations
		template<UInt ndim>
		inline Point<ndim> getLocations(UInt i) const {return Point<ndim>(i, locations_);}
		//! A method returning TRUE if the observations are located in the nodes of the mesh or FALSE otherwise
		inline bool isLocationsByNodes() const {return locations_by_nodes_;}

		//void newDatamatrix(const VectorXr& scores_,const VectorXr& loadings_);

		//! A method returning a reference to the observations vector
		inline MatrixXr const & getDatamatrix() const {return datamatrix_;}

		//! A method returning a reference to the incidence matrix
		inline MatrixXi const & getIncidenceMatrix() const {return incidenceMatrix_;}

		//! A method returning the number of observations
		inline UInt const getNumberofObservations() const {return datamatrix_.cols();}
		//! A method returning the observations indices
		inline std::vector<UInt> const & getObservationsIndices() const {return observations_indices_;}

		//! A method returning the number of regions
		inline UInt const getNumberOfRegions() const {return nRegions_;}

		//! A method returning the number of Principal Components to compute
		inline UInt const getNPC() const {return nPC_;}

		//! A method returning the the penalization term
		inline std::vector<Real> const & getLambda() const {return lambda_;}
		//! A method returning the input order
		inline UInt const getOrder() const {return order_;}
		//! A method returning the input search
		inline UInt const getSearch() const {return search_;}
		//! A method returning the input nFolds
		inline UInt const getNFolds() const {return nFolds_;}
		//! A method returning the method that should be used to compute the GCV:
		//! 1: exact calculation
		//! 2: stochastic estimation
		inline UInt const & getGCVmethod() const {return GCVmethod_;}
		//! A method returning the number of vectors to use to stochastically estimate the edf
		inline UInt const & getNrealizations() const {return nrealizations_;}
		inline MatrixXr const & getBarycenters() const {return barycenters_;}
		inline VectorXi const & getElementIds() const {return element_ids_;}
		inline Real const & getBarycenter(int i, int j) const {return barycenters_(i,j);}
		inline UInt const & getElementId(Id i) const {return element_ids_(i);}

		inline bool isLocationsByBarycenter() const {return locations_by_barycenter_;}


};

#endif
