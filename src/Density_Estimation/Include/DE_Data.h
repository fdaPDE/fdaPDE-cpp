#ifndef __DE_DATA_H__
#define __DE_DATA_H__

#include "../../FdaPDE.h"
#include "../../Mesh/Include/Mesh_Objects.h"
#include <array>

// This file contains the R/C++ data conversion for the Density Estimation problem

/*!  @brief An IO handler class for objects passed from R.
 * This class, given the data from R, convert them in a C++ format, offering a
 * series of method for their access, so isolating the more the possible the specific
 * code for R/C++ data conversion.
*/
template <UInt ndim>
class  DEData{
	private:
		// Data = locations.
		std::vector<Point<ndim> > data_;
		// Number of observations.
		UInt n_;
		// Finite element order.
		UInt order_;
		// Initial coefficients for the density.
		VectorXr fvec_;
		// Time step parameter for the heat diffusion process.
		Real heatStep_;
		// Number of iterations for the heat diffusion process.
		UInt heatIter_;
		// Penalization parameters. The best one is chosen with k fold cross validation.
		std::vector<Real> lambda_;
		// Number of folds for cross validatiom.
		UInt Nfolds_;
		// Number of simulations of the optimization algorithm.
		UInt nsim_;
		// Optimization parameters for fixed step methods.
		std::vector<Real> stepProposals_;
		// Tolerances for optimization algorithm termination criteria.
		Real tol1_;
		Real tol2_;
		// A boolean that is true if the user wants to see the value of the functional printed during the optimization descent.
		bool print_;

		// Integer specifying the search algorithm type (tree or naive search algorithm)
		UInt search_;

		// Auxiliary methods used in the constructor.

		void setData(SEXP Rdata);
		void setFvec(SEXP Rfvec);
		void setLambda(SEXP Rlambda);
		void setStepProposals(SEXP RstepProposals);


	public:
		// Constructors
		DEData(){};

		explicit DEData(const std::vector<Point<ndim> >& data, const UInt& order, const VectorXr& fvec, Real heatStep, UInt heatIter,
			const std::vector<Real>& lambda, const UInt& nfolds, const UInt& nsim, const std::vector<Real>& stepProposals,
			 Real tol1, Real tol2, bool print, UInt search);


		/*! Costructor useful for the R C++ interface.
			It initializes the object storing the R given objects.
			\param Rdata an R-matrix containing the data.
			\param Rorder an R-integer containing the order of the approximating basis.
			\param Rfvec an R-vector containing the the initial solution coefficients given by the user.
			\param RheatStep an R-double containing the step for the heat equation initialization.
			\para, RheatIter an R-integer containing the number of iterations to perfrom the heat equation initialization.
			\param Rlambda an R-vector containing the penalization terms.
			\param Rnfolds an R-integer specifying the number of folds for cross validation.
			\param Rnsim an R-integer specifying the number of iterations to use in the optimization algorithm.
			\param RstepProposals an R-vector containing the step parameters useful for the descent algotihm.
			\param Rtol1 an R-double specifying the tolerance to use for the termination criterion based on the percentage differences.
			\param Rtol2 an R-double specifying the tolerance to use for the termination criterion based on the norm of the gradient.
			\param Rprint and R-integer specifying if print on console.
			\param Rsearch an R-integer to decide the search algorithm type (tree or naive search algorithm).
		*/
		explicit DEData(SEXP Rdata, SEXP Rorder, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter, SEXP Rlambda, SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals,
			SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP Rsearch);


		// Setters
		//! A method to set new data (needed for projection).
		void setNewData(const std::vector<Point<ndim> >& );
		void setDatum(const Point<ndim>& , UInt );
		void updateN(UInt );

		// Getters
		//! A method returning the data.
		inline std::vector<Point<ndim> > getData() const {return data_;}
		//! A method returning a datum.
		inline Point<ndim> getDatum(UInt i) const {return data_[i];}
		//! A method returning the number of observations.
		inline UInt getNumberofData() const {return n_;}
		//! A method returning the the input order.
		inline UInt getOrder() const {return order_;}
		//! A method returning the initial coefficients for the density.
		inline VectorXr getFvec() const {return fvec_;}
		//! A method returning the heat diffusion process alpha parameter.
		inline Real getHeatStep() const {return heatStep_;}
		//! A method returning the number of iterations for the heat diffusion process.
		inline UInt getHeatIter() const {return heatIter_;}
		//! A method returning a bool which says if there is a user's initial density.
		inline bool isFvecEmpty() const {return fvec_.size()==0;}
		//! A method returning the penalization parameters.
		inline Real getLambda(UInt i) const {return lambda_[i];}
		//! A method returning the number of lambdas.
		inline UInt getNlambda()  const {return lambda_.size();}
		//! A method returning the number of folds for CV.
		inline UInt getNfolds()  const {return Nfolds_;}
		//! A method returning the number of iterations to use in the optimization algorithm.
		inline UInt getNsimulations() const {return nsim_;}
		//! A method returning the number of parameters for fixed step methods.
		inline UInt getNstepProposals() const {return stepProposals_.size();}
		//! A method returning a parameter for fixed step methods.
		inline Real getStepProposals(UInt i) const {return stepProposals_[i];}
		//! A method returning the tolerance for optimization algorithm first termination criteria.
		inline Real getTol1() const {return tol1_;}
		//! A method returning the tolerance for optimization algorithm second termination criteria.
		inline Real getTol2() const {return tol2_;}
		//! A method returning the boolean print member.
		inline bool Print() const {return print_;}
		//! A method returning the integer that specifies the search algorithm type.
		inline UInt getSearch() const {return search_;}

		// Print
		//! A method printing data.
		void printData(std::ostream & out) const;

};

#include "DE_Data_imp.h"
#endif
