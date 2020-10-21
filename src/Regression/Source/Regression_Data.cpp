#include "../Include/Regression_Data.h"

RegressionData::RegressionData(std::vector<Point> & locations, VectorXr & observations, UInt order, MatrixXr & covariates,
	std::vector<UInt> & bc_indices, std::vector<Real> & bc_values, MatrixXi & incidenceMatrix, bool arealDataAvg, UInt search):
	locations_(locations), observations_(observations), arealDataAvg_(arealDataAvg), order_(order),
	bc_values_(bc_values), bc_indices_(bc_indices), covariates_(covariates), incidenceMatrix_(incidenceMatrix),
	flag_SpaceTime_(false), search_(search)
{
	nRegions_ = incidenceMatrix_.rows();
	if(locations_.size()==0 && nRegions_==0)
	{
		locations_by_nodes_ = true;
		for(UInt i=0; i<observations_.size();++i)
			observations_indices_.push_back(i);
	}
	else
	{
		locations_by_nodes_ = false;
	}
}

RegressionData::RegressionData(std::vector<Point> & locations, std::vector<Real> & time_locations, VectorXr & observations, UInt order,
	MatrixXr & covariates, std::vector<UInt> & bc_indices, std::vector<Real> & bc_values,  MatrixXi & incidenceMatrix, bool arealDataAvg,
	bool flag_mass, bool flag_parabolic, VectorXr & ic, UInt search):
	locations_(locations), observations_(observations), arealDataAvg_(arealDataAvg), time_locations_(time_locations),
	order_(order), bc_values_(bc_values), bc_indices_(bc_indices), ic_(ic), covariates_(covariates), incidenceMatrix_(incidenceMatrix),
	flag_mass_(flag_mass), flag_parabolic_(flag_parabolic), flag_SpaceTime_(true), search_(search)
{
	nRegions_ = incidenceMatrix_.rows();
	if(locations_.size()==0 && nRegions_==0)
	{
		locations_by_nodes_ = true;
		for(UInt i=0; i<observations_.size(); ++i)
			observations_indices_.push_back(i);
	}
	else
	{
		locations_by_nodes_ = false;
	}
}

// costructor with WeightMatrix
RegressionData::RegressionData(std::vector<Point> & locations, VectorXr & observations, UInt order, MatrixXr & covariates,
	 VectorXr & WeightsMatrix, std::vector<UInt> & bc_indices, std::vector<Real> & bc_values,  MatrixXi & incidenceMatrix, bool arealDataAvg, UInt search):
	locations_(locations), observations_(observations), arealDataAvg_(arealDataAvg), WeightsMatrix_(WeightsMatrix),
	order_(order), bc_values_(bc_values), bc_indices_(bc_indices), covariates_(covariates), incidenceMatrix_(incidenceMatrix),
	flag_SpaceTime_(false), search_(search)
{
	nRegions_ = incidenceMatrix_.rows();
	if(locations_.size()==0 && nRegions_==0)
	{
		locations_by_nodes_ = true;
		for(UInt i=0; i<observations_.size();++i)
			observations_indices_.push_back(i);
	}
	else
	{
		locations_by_nodes_ = false;
	}
}

RegressionDataElliptic::RegressionDataElliptic(std::vector<Point> & locations, VectorXr & observations, UInt order,
	Eigen::Matrix<Real,2,2> & K, Eigen::Matrix<Real,2,1> & beta, Real c, MatrixXr & covariates,
	std::vector<UInt> & bc_indices, std::vector<Real> & bc_values, MatrixXi & incidenceMatrix, bool arealDataAvg, UInt search):
	RegressionData(locations, observations, order, covariates, bc_indices, bc_values, incidenceMatrix, arealDataAvg, search),
	K_(K), beta_(beta), c_(c)
{;}

RegressionDataElliptic::RegressionDataElliptic(std::vector<Point> & locations, std::vector<Real> & time_locations, VectorXr & observations,
	UInt order, Eigen::Matrix<Real,2,2> & K, Eigen::Matrix<Real,2,1> & beta, Real c, MatrixXr & covariates,
	std::vector<UInt> & bc_indices, std::vector<Real> & bc_values, MatrixXi & incidenceMatrix, bool arealDataAvg,
	bool flag_mass, bool flag_parabolic, VectorXr & ic, UInt search):
	RegressionData(locations, time_locations, observations, order, covariates, bc_indices, bc_values,
		incidenceMatrix, arealDataAvg, flag_mass, flag_parabolic, ic, search), K_(K), beta_(beta), c_(c)
{;}

RegressionDataEllipticSpaceVarying::RegressionDataEllipticSpaceVarying(std::vector<Point> & locations, VectorXr & observations, UInt order,
	const std::vector<Eigen::Matrix<Real,2,2>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,2> > > & K,
	const std::vector<Eigen::Matrix<Real,2,1>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,1> > > & beta,
	const std::vector<Real> & c, const std::vector<Real> & u, MatrixXr & covariates,
	std::vector<UInt> & bc_indices, std::vector<Real> & bc_values ,MatrixXi & incidenceMatrix, bool arealDataAvg, UInt search):
	RegressionData(locations, observations, order, covariates, bc_indices, bc_values, incidenceMatrix, arealDataAvg, search),
	K_(K), beta_(beta), c_(c), u_(u)
{;}

RegressionDataEllipticSpaceVarying::RegressionDataEllipticSpaceVarying(std::vector<Point> & locations, std::vector<Real> & time_locations, VectorXr & observations, UInt order,
	const std::vector<Eigen::Matrix<Real,2,2>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,2> > > & K,
	const std::vector<Eigen::Matrix<Real,2,1>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,1> > > & beta,
	const std::vector<Real> & c, const std::vector<Real> & u, MatrixXr& covariates,
	std::vector<UInt> & bc_indices,	std::vector<Real> & bc_values, MatrixXi & incidenceMatrix, bool arealDataAvg,
	bool flag_mass, bool flag_parabolic, VectorXr & ic, UInt search):
	RegressionData(locations, time_locations, observations, order, covariates,
		bc_indices, bc_values, incidenceMatrix, arealDataAvg, flag_mass, flag_parabolic, ic, search),
	K_(K), beta_(beta), c_(c), u_(u)
{;}


RegressionData::RegressionData(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder, SEXP Rcovariates,
	SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch)
{
	flag_SpaceTime_ = false;

	setLocations(Rlocations);
	setBaryLocations(RbaryLocations);
	setIncidenceMatrix(RincidenceMatrix);
	setObservations(Robservations);
	setCovariates(Rcovariates);

	order_ =  INTEGER(Rorder)[0];
	search_ =  INTEGER(Rsearch)[0];

	UInt length_indexes = Rf_length(RBCIndices);
	bc_indices_.assign(INTEGER(RBCIndices), INTEGER(RBCIndices) +  length_indexes);
	bc_values_.assign(REAL(RBCValues),REAL(RBCValues) + Rf_length(RBCIndices));

	arealDataAvg_ = INTEGER(RarealDataAvg)[0];
}

RegressionData::RegressionData(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rorder, SEXP Rcovariates,
	SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric, SEXP Rsearch)
{
	flag_SpaceTime_ = true;

	setLocations(Rlocations);
	setTimeLocations(Rtime_locations);
	setBaryLocations(RbaryLocations);
	setIncidenceMatrix(RincidenceMatrix);
	setObservationsTime(Robservations);
	setCovariates(Rcovariates);

	order_ =  INTEGER(Rorder)[0];
	search_ =  INTEGER(Rsearch)[0];
	flag_mass_ = INTEGER(Rflag_mass)[0];
	flag_parabolic_ = INTEGER(Rflag_parabolic)[0];

	UInt length_indexes = Rf_length(RBCIndices);
	bc_indices_.assign(INTEGER(RBCIndices), INTEGER(RBCIndices) +  length_indexes);
	bc_values_.assign(REAL(RBCValues),REAL(RBCValues) + Rf_length(RBCIndices));

	arealDataAvg_ = INTEGER(RarealDataAvg)[0];

	UInt length_ic = Rf_length(Ric);
	ic_.resize(length_ic);
	for (UInt i=0; i<length_ic; ++i)
	{
		ic_(i)=REAL(Ric)[i];
	}
}


RegressionDataElliptic::RegressionDataElliptic(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder,
	SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues,
	SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch):
	RegressionData(Rlocations, RbaryLocations, Robservations, Rorder, Rcovariates,
		RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch)
{
	K_.resize(2, 2);
	for(auto i=0; i<2; ++i)
	{
		for(auto j=0; j<2 ; ++j)
		{
			K_(i,j)=REAL(RK)[i+ 2*j];
		}
	}

	beta_.resize(2);
	for(auto i=0; i<2 ; ++i)
	{
		beta_(i)=REAL(Rbeta)[i];
	}

	c_ =  REAL(Rc)[0];
}

RegressionDataElliptic::RegressionDataElliptic(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rorder,
	SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues,
	SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric, SEXP Rsearch):
	RegressionData(Rlocations, RbaryLocations, Rtime_locations, Robservations, Rorder, Rcovariates,
		RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rflag_mass, Rflag_parabolic, Ric, Rsearch)
{
	K_.resize(2, 2);
	for(auto i=0; i<2; ++i)
	{
		for(auto j=0; j<2 ; ++j)
		{
			K_(i,j)=REAL(RK)[i+ 2*j];
		}
	}

	beta_.resize(2);
	for(auto i=0; i<2 ; ++i)
	{
		beta_(i)=REAL(Rbeta)[i];
	}

	c_ =  REAL(Rc)[0];
}

RegressionDataEllipticSpaceVarying::RegressionDataEllipticSpaceVarying(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder,
	SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues,
	SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch):
	RegressionData(Rlocations, RbaryLocations, Robservations, Rorder, Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch),
	K_(RK), beta_(Rbeta), c_(Rc), u_(Ru)
{;}

RegressionDataEllipticSpaceVarying::RegressionDataEllipticSpaceVarying(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rorder,
	SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix, SEXP RarealDataAvg,
	SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric, SEXP Rsearch):
	RegressionData(Rlocations, RbaryLocations, Rtime_locations, Robservations, Rorder, Rcovariates,
		RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rflag_mass, Rflag_parabolic, Ric, Rsearch),
	K_(RK), beta_(Rbeta), c_(Rc), u_(Ru)
{;}

void RegressionDataEllipticSpaceVarying::print(std::ostream & out) const
{
	for (auto i=0; i<18; i++)
		out << K_(i);
	for (auto i=0; i<18; i++)
		out << beta_(i);
	for (auto i=0; i<18; i++)
		out << c_(i);
}

void RegressionData::setObservations(SEXP Robservations)
{
	UInt n_obs_ = Rf_length(Robservations);
	observations_.resize(n_obs_);
	observations_indices_.reserve(n_obs_);

	UInt count = 0;
	if(locations_.size() == 0 && nRegions_ == 0)
	{
		locations_by_nodes_ = true;
		for(auto i=0;i<n_obs_;++i)
		{
			if(!ISNA(REAL(Robservations)[i]))
			{
				observations_[count] = REAL(Robservations)[i];
				count++;
				observations_indices_.push_back(i);
			}
		}
		observations_.conservativeResize(count, Eigen::NoChange);
	}
	else // locations_.size() > 0 NOR nRegions_ > 0
	{
		locations_by_nodes_ = false;
		for(auto i=0;i<n_obs_;++i)
		{
			observations_[i] = REAL(Robservations)[i];
		}
	}

	//std::cout<<"Observations #"<<observations_.size()<<std::endl<<observations_<<std::endl;
	//for(auto i=0;i<observations_indices_.size();++i)	std::cout<<observations_indices_[i]<<std::endl;
}

void RegressionData::setObservationsTime(SEXP Robservations)
{
	UInt n_obs_ = Rf_length(Robservations);
	observations_.resize(n_obs_);
	observations_indices_.reserve(n_obs_);

	UInt count = 0;
	locations_by_nodes_ = (locations_.size() == 0 && nRegions_ == 0) ? true : false;

	for(auto i=0;i<n_obs_;++i)
	{
		if(!ISNA(REAL(Robservations)[i]))
		{
			observations_(i) = REAL(Robservations)[i];
			observations_indices_.push_back(i);
		}
		else
		{
			observations_(i) = 0.0;
			observations_na_.push_back(i);
		}
	}
	//std::cout<<"Observations #"<<observations_.size()<<std::endl<<observations_<<std::endl;
	//for(auto i=0;i<observations_indices_.size();++i)	std::cout<<observations_indices_[i]<<std::endl;
}

void RegressionData::setCovariates(SEXP Rcovariates)
{
	n_ = INTEGER(Rf_getAttrib(Rcovariates, R_DimSymbol))[0];
	p_ = INTEGER(Rf_getAttrib(Rcovariates, R_DimSymbol))[1];
	UInt k=0;

	covariates_.resize(n_, p_);

	for(auto i=0; i<n_; ++i)
	{
		for(auto j=0; j<p_ ; ++j)
		{
			if(observations_na_.size()>k && i==observations_na_[k])
			{
				covariates_(i,j)=0;
				k++;
			}
			else
			{
				covariates_(i,j)=REAL(Rcovariates)[i+ n_*j];
			}
		}
	}
}


void RegressionData::setLocations(SEXP Rlocations)
{
	n_ = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[0];
	if(n_>0)
	{
		int ndim = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[1];

		if(ndim == 2)
		{
			for(auto i=0; i<n_; ++i)
			{
				locations_.emplace_back(REAL(Rlocations)[i+ n_*0],REAL(Rlocations)[i+ n_*1]);
			}
		}
		else
		{ //ndim == 3
			for(auto i=0; i<n_; ++i)
			{
				locations_.emplace_back(REAL(Rlocations)[i+ n_*0],REAL(Rlocations)[i+ n_*1],REAL(Rlocations)[i+ n_*2]);
			}
		}
	}
}

void RegressionData::setTimeLocations(SEXP Rtime_locations)
{
	UInt n_time_loc_ = Rf_length(Rtime_locations);
  	time_locations_.resize(n_time_loc_);

	for(auto i=0;i<n_time_loc_;++i)
	{
		time_locations_[i] = REAL(Rtime_locations)[i];
	}
}

void RegressionData::setBaryLocations(SEXP RbaryLocations)
{
	//RECIEVE BARYCENTER INFORMATION FROM R
	if (TYPEOF(RbaryLocations) != 0)
	{ //TYPEOF(RbaryLocations) == 0 means SEXPTYPE is NILSXP (Description is NULL)
		UInt* id_ 	= INTEGER(VECTOR_ELT(RbaryLocations, 1));
		Real* bary_ 	= REAL(VECTOR_ELT(RbaryLocations, 2));

		UInt n_ = INTEGER(Rf_getAttrib(VECTOR_ELT(RbaryLocations, 2), R_DimSymbol))[0];
		UInt p_ = INTEGER(Rf_getAttrib(VECTOR_ELT(RbaryLocations, 2), R_DimSymbol))[1]; //barycenter column dimension

		barycenters_.resize(n_, p_);
		element_ids_.resize(n_);

		if(n_>0)
		{
			for(auto i=0; i<n_; ++i)
			{
				for (auto j=0; j<p_; ++j)
				{
					barycenters_(i,j)= bary_[i+ n_*j];
				}
				element_ids_(i) = id_[i];
			}
		}
		locations_by_barycenter_ =true;
	}
	else
	{
		locations_by_barycenter_ =false;
	}
}

void RegressionData::setIncidenceMatrix(SEXP RincidenceMatrix)
{
	nRegions_ = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[0];
	UInt p = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[1];

	incidenceMatrix_.resize(nRegions_, p);

	for(auto i=0; i<nRegions_; ++i)
	{
		for(auto j=0; j<p; ++j)
		{
			incidenceMatrix_(i,j) = INTEGER(RincidenceMatrix)[i+nRegions_*j];
		}
	}
}


void RegressionData::printObservations(std::ostream & out) const
{

	for(auto i=0;i<observations_.size(); i++)
	{
		out << i << "\t" << observations_(i) << std::endl;
	}
}

void RegressionData::printCovariates(std::ostream & out) const
{

	for(auto i=0;i<covariates_.rows(); i++)
	{
		for(auto j=0; j<covariates_.cols(); j++)
		{
			out << covariates_(i,j) << "\t";
		}
		out << std::endl;
	}
}

void RegressionData::printLocations(std::ostream & out) const
{

	for(std::vector<Point>::size_type i=0;i<locations_.size(); i++)
	{
		locations_[i].print(out);
		//std::cout<<std::endl;
	}
}

void RegressionData::printIncidenceMatrix(std::ostream & out) const
{
	for (auto i=0; i<incidenceMatrix_.rows(); i++)
	{
		for (auto j=0; j<incidenceMatrix_.cols(); j++)
		{
			out << incidenceMatrix_(i,j) << "\t";
		}
		out << std::endl;
	}
}
