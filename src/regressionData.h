#ifndef __REGRESSIONDATA_HPP__
#define __REGRESSIONDATA_HPP__

#include "fdaPDE.h"
#include "mesh_objects.h"
#include "param_functors.h"

//!  An IO handler class for objects passed from R
/*!
 * This class, given the data from R, convert them in a C++ format, offering a
 * series of method for their access, so isolating the more the possible the specific
 * code for R/C++ data conversion.
*/
class  RegressionData{
	protected:
		// Design matrix pointer and dimensions
		std::vector<Point> locations_;

		VectorXr observations_;
		
		bool locations_by_nodes_; // if location is on the mesh nodes or not

		UInt nRegions_; //for areal data

		VectorXr WeightsMatrix_; //Weighted regression

		std::vector<Real> lambdaS_; //space penalization

		bool DOF_; // do we need to compute DoF?


	private:

		//Vector of the time locations
		std::vector<Real> time_locations_;

		std::vector<UInt> observations_indices_;
		std::vector<UInt> observations_na_;

		//barycenter information
		VectorXi element_ids_; //elements id information
		MatrixXr barycenters_; //barycenter information
		bool locations_by_barycenter_;

		//Design matrix
		MatrixXr covariates_;
		UInt n_;
		UInt p_;

		//Areal data
		MatrixXi incidenceMatrix_;

		//Other parameters
		UInt order_;

		std::vector<Real> lambdaT_;	//time penalization
		UInt GCVmethod_;
		UInt nrealizations_;      // Number of relizations for the stochastic estimation of GCV

		std::vector<Real> bc_values_;
		std::vector<UInt> bc_indices_;

		VectorXr ic_; //initial conditions
		bool flag_mass_;	//mass penalization, only for separable version (flag_parabolic_==FALSE)
		bool flag_parabolic_;

		//bool inputType;
		bool GCV_;
		MatrixXr dof_matrix_;

		bool flag_SpaceTime_; // TRUE if space time smoothing

		UInt search_; // search algorith type

		#ifdef R_VERSION_
		void setBaryLocations(SEXP RbaryLocations);
		void setLocations(SEXP Rlocations);
		void setTimeLocations(SEXP Rtime_locations);
		void setObservations(SEXP Robservations);
		void setObservationsTime(SEXP Robservations);
		void setCovariates(SEXP Rcovariates);
		void setNrealizations(SEXP Rnrealizations);
		void setDOF_matrix(SEXP RDOF_matrix);
		void setIncidenceMatrix(SEXP RincidenceMatrix);
		#endif



	public:

		RegressionData(){};

//! A basic version of the constructor.

		/*!
			It initializes the object storing the R given objects. This is the simplest of the two possible interfaces with R
			\param Rlocations an R-matrix containing the location of the observations.
			\param Robservations an R-vector containing the values of the observations.
			\param Rorder an R-integer containing the order of the approximating basis.
			\param RlambdaS an R-double containing the penalization term of the empirical evidence respect to the prior one.
			\param Rcovariates an R-matrix storing the covariates of the regression
			\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions in the model with areal data
			\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
					the other are automatically considered in Neumann Condition.
			\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in Rbindex
			\param GCV an R boolean indicating whether GCV has to be computed or not
			\param DOF an R boolean indicating whether dofs of the model have to be computed or not
	        \param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic)
	        \param Rnrealizations the number of random points used in the stochastic computation of the dofs
	        \param Rsearch an R-integer to decide the search algorithm type (tree or naive or walking search algorithm).
		*/


		#ifdef R_VERSION_
		explicit RegressionData(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS, SEXP Rcovariates,
								SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP GCV, SEXP RGCVmethod,
								SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch);
		
		explicit RegressionData(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS, SEXP RlambdaT, SEXP Rcovariates,
						SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric, SEXP GCV, SEXP RGCVmethod,
						SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch);
		#endif

		explicit RegressionData(std::vector<Point>& locations, VectorXr& observations, UInt order, std::vector<Real> lambdaS, MatrixXr& covariates, MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, bool DOF, bool GCV,  UInt search);

		explicit RegressionData(std::vector<Point>& locations, std::vector<Real>& time_locations, VectorXr& observations, UInt order,
														std::vector<Real>& lambdaS, std::vector<Real>& lambdaT, MatrixXr& covariates, MatrixXi& incidenceMatrix,
														std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, VectorXr& ic, bool flag_mass, bool flag_parabolic, bool DOF, bool GCV,  UInt search);

		explicit RegressionData(std::vector<Point>& locations, VectorXr& observations, UInt order, std::vector<Real> lambdaS, MatrixXr& covariates, VectorXr& WeightsMatrix, MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, bool DOF, bool GCV,  UInt search);

		void printObservations(std::ostream & out) const;
		void printCovariates(std::ostream & out) const;
		void printLocations(std::ostream & out) const;
		void printIncidenceMatrix(std::ostream & out) const;

		//! A method returning a reference to the observations vector
		inline VectorXr const & getObservations() const {return observations_;}
		//! A method returning a reference to the design matrix
		inline MatrixXr const & getCovariates() const {return covariates_;}
		//! A method returning a reference to the matrix of weights
		inline VectorXr const & getWeightsMatrix() const {return WeightsMatrix_;}
		//! A method returning a reference to the incidence matrix
		inline MatrixXi const & getIncidenceMatrix() const {return incidenceMatrix_;}
		//! A method returning the number of observations
		inline UInt const getNumberofObservations() const {return observations_.size();}
		//! A method returning the number of space observations
		inline UInt const getNumberofSpaceObservations() const {return observations_.size()/time_locations_.size();}
		//! A method returning the number of time observations
		inline UInt const getNumberofTimeObservations() const {return time_locations_.size();}
		//! A method returning the locations of the observations
		inline std::vector<Point> const & getLocations() const {return locations_;}
		//! A method returning the locations of the time observations
		inline std::vector<Real> const & getTimeLocations() const {return time_locations_;}
		//! A method returning the number of regions
		inline UInt const getNumberOfRegions() const {return nRegions_;}
		inline bool isLocationsByNodes() const {return locations_by_nodes_;}
		inline bool isLocationsByBarycenter() const {return locations_by_barycenter_;}
		inline bool computeDOF() const {return DOF_;}
		inline bool computeGCV() const {return GCV_;}
		inline std::vector<UInt> const & getObservationsIndices() const {return observations_indices_;}
		inline std::vector<UInt> const & getObservationsNA() const {return observations_na_;}

		//! A method returning the space penalization term
		inline std::vector<Real> const & getLambdaS() const {return lambdaS_;}
		//! A method returning the the time penalization term
		inline std::vector<Real> const & getLambdaT() const {return lambdaT_;}
		inline bool const & isSpaceTime() const {return flag_SpaceTime_;}
		inline bool const & getFlagMass() const {return flag_mass_;}
		inline bool const & getFlagParabolic() const {return flag_parabolic_;}

		//! A method returning the input order
		inline UInt const getOrder() const {return order_;}
		//! A method returning the input search
		inline UInt const getSearch() const {return search_;}
		//! A method returning the indexes of the nodes for which is needed to apply Dirichlet Conditions
		inline std::vector<UInt> const & getDirichletIndices() const {return bc_indices_;}
		//! A method returning the values to apply for Dirichlet Conditions
		inline std::vector<Real> const & getDirichletValues() const {return bc_values_;}
		//! A method returning the values to apply for Initial Conditions
		inline VectorXr const & getInitialValues() const {return ic_;}
		//! A method returning the method that should be used to compute the GCV:
		//! 1: exact calculation
		//! 2: stochastic estimation
		inline UInt const & getGCVmethod() const {return GCVmethod_;}
		//! A method returning the number of vectors to use to stochastically estimate the edf
		inline UInt const & getNrealizations() const {return nrealizations_;}
		inline MatrixXr const & getDOF_matrix() const {return dof_matrix_;}

		inline MatrixXr const & getBarycenters() const {return barycenters_;}
		inline VectorXi const & getElementIds() const {return element_ids_;}
		inline Real const & getBarycenter(int i, int j) const {return barycenters_(i,j);}
		inline UInt const & getElementId(Id i) const {return element_ids_(i);}
};


class  RegressionDataElliptic:public RegressionData
{
	private:
		Eigen::Matrix<Real,2,2> K_;
		Eigen::Matrix<Real,2,1> beta_;
		Real c_;

	public:
		//! A complete version of the constructor.
		/*!
			It initializes the object storing the R given objects. This is the simplest of the two possible interfaces with R
			\param Rlocations an R-matrix containing the location of the observations.
			\param Robservations an R-vector containing the values of the observations.
			\param Rorder an R-integer containing the order of the approximating basis.
			\param RlambdaS an R-double containing the penalization term of the empirical evidence respect to the prior one.
			\param RK an R-double 2X2 matrix containing the coefficients for a anisotropic DIFFUSION term.
			\param Rbeta an R-double 2-dim vector that contains the coefficients for the TRANSPORT coefficients.
			\param Rc an R-double that contains the coefficient of the REACTION term
			\param Rcovariates an R-matrix storing the covariates of the regression
			\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions in the model with areal data
			\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
					the other are automatically considered in Neumann Condition.
			\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in Rbindex
			\param DOF an R boolean indicating whether dofs of the model have to be computed or not
	        \param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic)
	        \param Rnrealizations the number of random points used in the stochastic computation of the dofs
	        \param Rsearch an R-integer to decide the search algorithm type (tree or naive or walking search algorithm).
		*/
		#ifdef R_VERSION_
		explicit RegressionDataElliptic(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS, SEXP RK,
				SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues,
				SEXP GCV,SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch);
		explicit RegressionDataElliptic(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS, SEXP RlambdaT, SEXP RK,
				SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric,
				SEXP GCV,SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch);
		#endif

		explicit RegressionDataElliptic(std::vector<Point>& locations, VectorXr& observations, UInt order,
										std::vector<Real> lambdaS, Eigen::Matrix<Real,2,2>& K,
										Eigen::Matrix<Real,2,1>& beta, Real c, MatrixXr& covariates,
										MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices,
										std::vector<Real>& bc_values, bool DOF, bool GCV, UInt search);

		explicit RegressionDataElliptic(std::vector<Point>& locations, std::vector<Real>& time_locations, VectorXr& observations, UInt order,
										std::vector<Real>& lambdaS, std::vector<Real>& lambdaT, Eigen::Matrix<Real,2,2>& K,	Eigen::Matrix<Real,2,1>& beta, Real c,
										MatrixXr& covariates, MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices,	std::vector<Real>& bc_values, VectorXr& ic,
										bool flag_mass, bool flag_parabolic, bool DOF, bool GCV, UInt search);

		inline Eigen::Matrix<Real,2,2> const & getK() const {return K_;}
		inline Eigen::Matrix<Real,2,1> const & getBeta() const {return beta_;}
		inline Real const getC() const {return c_;}
};

class RegressionDataEllipticSpaceVarying:public RegressionData
{
	private:
		Diffusivity K_;
		Advection beta_;
		Reaction c_;
		ForcingTerm u_;

	public:

		//! A complete version of the constructor.
		/*!
			It initializes the object storing the R given objects. This is the simplest of the two possible interfaces with R
			\param Rlocations an R-matrix containing the location of the observations.
			\param Robservations an R-vector containing the values of the observations.
			\param Rorder an R-integer containing the order of the approximating basis.
			\param RlambdaS an R-double containing the penalization term of the empirical evidence respect to the prior one.
			\param RK an R-double 2X2 matrix containing the coefficients for a anisotropic DIFFUSION term.
			\param Rbeta an R-double 2-dim vector that contains the coefficients for the TRANSPORT coefficients.
			\param Rc an R-double that contains the coefficient of the REACTION term
			\param  Ru an R-double vector of length #triangles that contaiins the forcing term integrals.
			\param Rcovariates an R-matrix storing the covariates of the regression
			\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions in the model with areal data
			\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
					the other are automatically considered in Neumann Condition.
			\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in Rbindex
			\param DOF an R boolean indicating whether dofs of the model have to be computed or not
	        \param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic)
	        \param Rnrealizations the number of random points used in the stochastic computation of the dofs
	        \param Rsearch an R-integer to decide the search algorithm type (tree or naive or walking search algorithm).

		*/
		#ifdef R_VERSION_
		explicit RegressionDataEllipticSpaceVarying(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS,
				SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices,
				SEXP RBCValues, SEXP GCV, SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix,SEXP Rsearch);
				explicit RegressionDataEllipticSpaceVarying(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS, SEXP RlambdaT, SEXP RK,
					  SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric,
					  SEXP GCV, SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch);
		#endif


		explicit RegressionDataEllipticSpaceVarying(std::vector<Point>& locations, VectorXr& observations,
													UInt order, std::vector<Real> lambdaS,
													const std::vector<Eigen::Matrix<Real,2,2>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,2> > >& K,
													const std::vector<Eigen::Matrix<Real,2,1>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,1> > >& beta,
													const std::vector<Real>& c, const std::vector<Real>& u,
													MatrixXr& covariates, MatrixXi& incidenceMatrix,
													std::vector<UInt>& bc_indices, std::vector<Real>& bc_values,
													bool DOF, bool GCV, UInt search);

		explicit RegressionDataEllipticSpaceVarying(std::vector<Point>& locations, std::vector<Real>& time_locations, VectorXr& observations, UInt order,
													std::vector<Real>& lambdaS, std::vector<Real>& lambdaT,
													const std::vector<Eigen::Matrix<Real,2,2>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,2> > >& K,
													const std::vector<Eigen::Matrix<Real,2,1>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,1> > >& beta,
													const std::vector<Real>& c, const std::vector<Real>& u,
													MatrixXr& covariates, MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices,	std::vector<Real>& bc_values, VectorXr& ic,
													bool flag_mass, bool flag_parabolic, bool DOF,bool GCV, UInt search);

		inline Diffusivity const & getK() const {return K_;}
		inline Advection const & getBeta() const {return beta_;}
		inline Reaction const & getC() const {return c_;}
		inline ForcingTerm const & getU() const {return u_;}

		void print(std::ostream & out) const;
};



/* parent template can be :

		- RegressionData
		- RegressionDataElliptic
		- RegressionDataEllipticSpaceVarying

*/

template<typename RegressionHandler>
class  RegressionDataGAM : public RegressionHandler
{
	private:
		
		VectorXr initialObservations_; // A copy of the true observations, which will not be overriden during PIRLS
	    std::vector<Real> global_lambda_; // A copy of labda vector for FPIRLS with multiple lambdas
		std::vector<UInt> initial_observations_indeces_;
	    UInt max_num_iterations_; // max n° of iterations allowed
	    Real threshold_; // limit in difference among J_k and J_k+1 for which we stop PIRLS
	    Real tune_param;
	    bool DOF_GAM; // flag for dof computation

	public:
		//! A complete version of the constructor.
		/*!
			It initializes the object storing the R given objects. This is the simplest of the two possible interfaces with R
			\param Rlocations an R-matrix containing the location of the observations.
			\param Robservations an R-vector containing the values of the observations.
			\param Rorder an R-integer containing the order of the approximating basis.
			\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
			\param RK an R-double 2X2 matrix containing the coefficients for a anisotropic DIFFUSION term.
			\param Rbeta an R-double 2-dim vector that contains the coefficients for the TRANSPORT coefficients.
			\param Rc an R-double that contains the coefficient of the REACTION term
			\param Rcovariates an R-matrix storing the covariates of the regression
			\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions in the model with areal data
			\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
					the other are automatically considered in Neumann Condition.
			\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in Rbindex
			\param DOF an R boolean indicating whether dofs of the model have to be computed or not
	        \param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic)
	        \param Rnrealizations the number of random points used in the stochastic computation of the dofs
		*/
		#ifdef R_VERSION_

		//Laplace
		explicit RegressionDataGAM(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS, 
				SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, 
				SEXP GCV, SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch, 
				SEXP Rmax_num_iteration, SEXP Rthreshold, SEXP Rtune);

		// PDE
		explicit RegressionDataGAM(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS, 
				SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, 
				SEXP GCV, SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch, 
				SEXP Rmax_num_iteration, SEXP Rthreshold, SEXP Rtune);

		// PDE SpaceVarying
		explicit RegressionDataGAM(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS, 
				SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, 
				SEXP GCV, SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch, 
				SEXP Rmax_num_iteration, SEXP Rthreshold, SEXP Rtune);

		#endif

		// Laplace
		explicit RegressionDataGAM(std::vector<Point>& locations, VectorXr& observations, UInt order, 
									std::vector<Real> lambdaS, MatrixXr& covariates, MatrixXi& incidenceMatrix, 
									std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, bool DOF, bool GCV, UInt search, 
									UInt max_num_iterations, Real threshold, Real tune);

		// PDE
		explicit RegressionDataGAM(std::vector<Point>& locations, VectorXr& observations, UInt order,
												std::vector<Real> lambdaS, Eigen::Matrix<Real,2,2>& K,
												Eigen::Matrix<Real,2,1>& beta, Real c, MatrixXr& covariates,
												MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices,
												std::vector<Real>& bc_values, bool DOF, bool GCV, UInt search, 
												UInt max_num_iterations, Real threshold, Real tune);

		// PDE SpaceVarying
		explicit RegressionDataGAM(std::vector<Point>& locations,
									VectorXr& observations, UInt order, std::vector<Real> lambdaS,
									const std::vector<Eigen::Matrix<Real,2,2>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,2> > >& K,
									const std::vector<Eigen::Matrix<Real,2,1>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,1> > >& beta,
									const std::vector<Real>& c, const std::vector<Real>& u,
									MatrixXr& covariates, MatrixXi& incidenceMatrix,
									std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, bool DOF, bool GCV, UInt search, 
									UInt max_num_iterations, Real threshold, Real tune);


		inline UInt const & get_maxiter() const {return max_num_iterations_;}
		inline Real const & get_treshold() const {return threshold_;}
		//! A method returning a reference to the observations vector
		inline VectorXr const & getInitialObservations() const {return initialObservations_;}
    	inline std::vector<Real> const & getGlobalLambda() const {return global_lambda_;}
		inline UInt const getNumberofInitialObservations() const {return initial_observations_indeces_.size();}
    	inline bool const getDOF_GAM() const {return DOF_GAM;}
    	inline Real const getTuneParam() const {return tune_param;}

		//! Update Pseudodata (observations and weights)
		void updatePseudodata(VectorXr& z_, VectorXr& P){ this-> observations_ = z_; this-> WeightsMatrix_ = P; }
    	void setCurrentLambda(UInt lambda_index){ this->lambdaS_ = std::vector<Real>(1,global_lambda_[lambda_index]);}

		void copyInitialObservations(){initialObservations_ = this->observations_;}
    	void copyLambdaVector(){global_lambda_ = this->lambdaS_;}
    	void setDOFflag(){DOF_GAM = this->DOF_; this->DOF_ = false;}

};


typedef RegressionDataGAM<RegressionData> GAMDataLaplace;

typedef RegressionDataGAM<RegressionDataElliptic> GAMDataElliptic;

typedef RegressionDataGAM<RegressionDataEllipticSpaceVarying> GAMDataEllipticSpaceVarying;


#include "regressionData_imp.h"

#endif
