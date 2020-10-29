#ifndef __REGRESSION_DATA_IMP_H__
#define __REGRESSION_DATA_IMP_H__


// Laplace
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder,
	SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch,
	SEXP Rmax_num_iteration, SEXP Rthreshold):
	RegressionData(Rlocations, RbaryLocations, Robservations, Rorder, Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch)
{
	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ =  REAL(Rthreshold)[0];
	initialObservations_ = this->observations_;
	this->isGAM = true;
}

// PDE
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder,
	SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues,
	SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch, SEXP Rmax_num_iteration, SEXP Rthreshold):
	RegressionDataElliptic(Rlocations, RbaryLocations, Robservations, Rorder, RK, Rbeta, Rc,
		Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch)
{
	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ =  REAL(Rthreshold)[0];
	initialObservations_ = this->observations_;
	this->isGAM = true;
}

// PDE SpaceVarying
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder,
	SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues,
	SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch, SEXP Rmax_num_iteration, SEXP Rthreshold):
	RegressionDataEllipticSpaceVarying(Rlocations, RbaryLocations, Robservations, Rorder, RK, Rbeta, Rc, Ru,
		Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch)
{
	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ =  REAL(Rthreshold)[0];
	initialObservations_ = this->observations_;
	this->isGAM = true;
}

#endif
