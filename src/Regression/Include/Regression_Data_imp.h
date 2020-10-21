#ifndef __REGRESSION_DATA_IMP_H__
#define __REGRESSION_DATA_IMP_H__

// -- GAM CONSTRUCTORS --
// Laplace
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(std::vector<Point> & locations, VectorXr & observations, UInt order,
	MatrixXr & covariates, std::vector<UInt> & bc_indices, std::vector<Real> & bc_values,
	MatrixXi & incidenceMatrix, bool arealDataAvg, UInt search, UInt max_num_iterations, Real threshold):
	RegressionData(locations, observations, order, covariates, bc_indices, bc_values, incidenceMatrix, arealDataAvg, search),
	max_num_iterations_(max_num_iterations), threshold_(threshold), initialObservations_(observations)
{this->isGAM = true;}

// PDE
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(std::vector<Point> & locations, VectorXr & observations, UInt order,
	Eigen::Matrix<Real,2,2> & K, Eigen::Matrix<Real,2,1> & beta, Real c, MatrixXr & covariates,
	std::vector<UInt> & bc_indices, std::vector<Real> & bc_values, MatrixXi& incidenceMatrix, bool arealDataAvg,
	UInt search, UInt max_num_iterations, Real threshold):
	RegressionDataElliptic(locations, observations, order, K, beta, c, covariates, bc_indices, bc_values, incidenceMatrix, arealDataAvg, search),
	max_num_iterations_(max_num_iterations), threshold_(threshold), initialObservations_(observations)
{this->isGAM = true;}

// PDE SpaceVarying
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(std::vector<Point> & locations, VectorXr & observations, UInt order,
	const std::vector<Eigen::Matrix<Real,2,2>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,2> > > & K,
	const std::vector<Eigen::Matrix<Real,2,1>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,1> > > & beta,
	const std::vector<Real> & c, const std::vector<Real> & u,
	MatrixXr& covariates, std::vector<UInt> & bc_indices, std::vector<Real> & bc_values,
	MatrixXi & incidenceMatrix, bool arealDataAvg, UInt search, UInt max_num_iterations, Real threshold):
	RegressionDataEllipticSpaceVarying(locations, observations, order, K, beta, c, u, covariates, bc_indices, bc_values, incidenceMatrix, arealDataAvg, search),
	max_num_iterations_(max_num_iterations), threshold_(threshold), initialObservations_(observations)
{this->isGAM = true;}

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
