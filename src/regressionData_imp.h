#ifndef __REGRESSIONDATA_IMP_HPP__
#define __REGRESSIONDATA_IMP_HPP__

RegressionData::RegressionData(Real* locations, UInt n_locations, UInt ndim, VectorXr& observations, UInt order, std::vector<Real> lambdaS, MatrixXr& covariates, MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, bool DOF, bool GCV, UInt search):
					locations_(locations, n_locations, ndim), observations_(observations), covariates_(covariates), incidenceMatrix_(incidenceMatrix),
					order_(order), lambdaS_(lambdaS), bc_values_(bc_values), bc_indices_(bc_indices), DOF_(DOF), GCV_(GCV), flag_SpaceTime_(false), search_(search)
{
	nRegions_ = incidenceMatrix_.rows();
	if(locations_.nrows()==0 && nRegions_==0)
	{
		locations_by_nodes_ = true;
		for(int i = 0; i<observations_.size();++i) observations_indices_.push_back(i);
	}
	else
	{
		locations_by_nodes_ = false;
	}
}

RegressionData::RegressionData(Real* locations, UInt n_locations, UInt ndim, std::vector<Real>& time_locations, VectorXr& observations, UInt order, std::vector<Real>& lambdaS, std::vector<Real>& lambdaT, MatrixXr& covariates, MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, VectorXr& ic, bool flag_mass, bool flag_parabolic, bool DOF, bool GCV, UInt search):
					locations_(locations, n_locations, ndim), time_locations_(time_locations), observations_(observations), covariates_(covariates), incidenceMatrix_(incidenceMatrix),
					order_(order), lambdaS_(lambdaS), lambdaT_(lambdaT), bc_values_(bc_values), bc_indices_(bc_indices), ic_(ic), flag_mass_(flag_mass), flag_parabolic_(flag_parabolic), DOF_(DOF), GCV_(GCV), flag_SpaceTime_(true), search_(search)
{
	nRegions_ = incidenceMatrix_.rows();
	if(locations_.nrows()==0 && nRegions_==0)
	{
		locations_by_nodes_ = true;
		for(int i = 0; i<observations_.size();++i) observations_indices_.push_back(i);
	}
	else
	{
		locations_by_nodes_ = false;
	}
}

RegressionDataElliptic::RegressionDataElliptic(Real* locations, UInt n_locations, UInt ndim, VectorXr& observations, UInt order,
												std::vector<Real> lambdaS, Real* K,
												Real* beta, Real c, MatrixXr& covariates,
												MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices,
												std::vector<Real>& bc_values, bool DOF, bool GCV, UInt search) :
		 RegressionData(locations, n_locations, ndim, observations, order, lambdaS, covariates, incidenceMatrix, bc_indices, bc_values, DOF, GCV, search), K_(K), beta_(beta), c_(c)
{}

RegressionDataElliptic::RegressionDataElliptic(Real* locations, UInt n_locations, UInt ndim, std::vector<Real>& time_locations, VectorXr& observations, UInt order,
												std::vector<Real>& lambdaS, std::vector<Real>& lambdaT, Real* K,
												Real* beta, Real c, MatrixXr& covariates,
												MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices,
												std::vector<Real>& bc_values, VectorXr& ic, bool flag_mass, bool flag_parabolic, bool DOF, bool GCV, UInt search):
		 RegressionData(locations, n_locations, ndim, time_locations, observations, order, lambdaS, lambdaT, covariates, incidenceMatrix, bc_indices, bc_values, ic, flag_mass, flag_parabolic, DOF, GCV, search),
		 										K_(K), beta_(beta), c_(c)
{}

RegressionDataEllipticSpaceVarying::RegressionDataEllipticSpaceVarying(Real* locations, UInt n_locations, UInt ndim,
									VectorXr& observations, UInt order, std::vector<Real> lambdaS,
									Real* K, Real* beta, Real* c, Real* u,
									MatrixXr& covariates, MatrixXi& incidenceMatrix,
									std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, bool DOF, bool GCV, UInt search):
		RegressionData(locations, n_locations, ndim, observations, order, lambdaS, covariates, incidenceMatrix, bc_indices, bc_values, DOF, GCV, search), K_(K), beta_(beta), c_(c), u_(u)
{}

RegressionDataEllipticSpaceVarying::RegressionDataEllipticSpaceVarying(Real* locations, UInt n_locations, UInt ndim, std::vector<Real>& time_locations, VectorXr& observations, UInt order,
												std::vector<Real>& lambdaS, std::vector<Real>& lambdaT,
												Real* K, Real* beta, Real* c, Real* u,
												MatrixXr& covariates, MatrixXi& incidenceMatrix,
												std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, VectorXr& ic, bool flag_mass, bool flag_parabolic, bool DOF, bool GCV, UInt search):
		 RegressionData(locations, n_locations, ndim, time_locations, observations, order, lambdaS, lambdaT, covariates, incidenceMatrix, bc_indices, bc_values, ic, flag_mass, flag_parabolic, DOF, GCV, search),
											K_(K), beta_(beta), c_(c), u_(u)
{}


#ifdef R_VERSION_
RegressionData::RegressionData(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS, SEXP Rcovariates,
							SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP GCV, SEXP RGCVmethod,
							SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch) :
	locations_(Rlocations)
{
	flag_SpaceTime_=false;

	setBaryLocations(RbaryLocations);
	setIncidenceMatrix(RincidenceMatrix);
	setObservations(Robservations);
	setCovariates(Rcovariates);
	setNrealizations(Rnrealizations);
	setDOF_matrix(RDOF_matrix);

	GCVmethod_ = INTEGER(RGCVmethod)[0];

	order_ =  INTEGER(Rorder)[0];
	search_ =  INTEGER(Rsearch)[0];
	GCV_ = INTEGER(GCV)[0];
	DOF_ = INTEGER(DOF)[0];
	UInt length_indexes = Rf_length(RBCIndices);

	bc_indices_.assign(INTEGER(RBCIndices), INTEGER(RBCIndices) +  length_indexes);

	bc_values_.assign(REAL(RBCValues),REAL(RBCValues) + Rf_length(RBCIndices));

  UInt length_lambda = Rf_length(RlambdaS);
  for (UInt i = 0; i<length_lambda; ++i)  lambdaS_.push_back(REAL(RlambdaS)[i]);

	lambdaT_.push_back(0);
}

RegressionData::RegressionData(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS, SEXP RlambdaT, SEXP Rcovariates,
						SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric, SEXP GCV, SEXP RGCVmethod,
						SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch) :
	locations_(Rlocations)
{
	flag_SpaceTime_ = true;

	setTimeLocations(Rtime_locations);
	setBaryLocations(RbaryLocations);
	setIncidenceMatrix(RincidenceMatrix);
	setObservationsTime(Robservations);
	setCovariates(Rcovariates);
	setNrealizations(Rnrealizations);
	setDOF_matrix(RDOF_matrix);

	GCVmethod_ = INTEGER(RGCVmethod)[0];

	order_ =  INTEGER(Rorder)[0];
	search_ =  INTEGER(Rsearch)[0];
	flag_mass_ = INTEGER(Rflag_mass)[0];
	flag_parabolic_ = INTEGER(Rflag_parabolic)[0];

	GCV_ = INTEGER(GCV)[0];
	DOF_ = INTEGER(DOF)[0];
	UInt length_indexes = Rf_length(RBCIndices);

	bc_indices_.assign(INTEGER(RBCIndices), INTEGER(RBCIndices) +  length_indexes);

	bc_values_.assign(REAL(RBCValues),REAL(RBCValues) + Rf_length(RBCIndices));

	UInt length_lambdaS = Rf_length(RlambdaS);
	for (UInt i = 0; i<length_lambdaS; ++i)  lambdaS_.push_back(REAL(RlambdaS)[i]);

	UInt length_lambdaT = Rf_length(RlambdaT);
	for (UInt i = 0; i<length_lambdaT; ++i)  lambdaT_.push_back(REAL(RlambdaT)[i]);

	UInt length_ic = Rf_length(Ric);
	ic_.resize(length_ic);
	for (UInt i = 0; i<length_ic; ++i)  ic_(i)=REAL(Ric)[i];

}


RegressionDataElliptic::RegressionDataElliptic(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS, SEXP RK, SEXP Rbeta,
				 SEXP Rc, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP GCV, SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch) :
	RegressionData(Rlocations, RbaryLocations, Robservations, Rorder, RlambdaS, Rcovariates, RincidenceMatrix,
					 			   RBCIndices, RBCValues, GCV, RGCVmethod, Rnrealizations, DOF, RDOF_matrix, Rsearch),
									 		K_(RK), beta_(Rbeta), c_(REAL(Rc)[0]) {}


RegressionDataElliptic::RegressionDataElliptic(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS, SEXP RlambdaT, SEXP RK,
		SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric,	SEXP GCV,SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch):
	RegressionData(Rlocations, RbaryLocations, Rtime_locations, Robservations, Rorder, RlambdaS, RlambdaT,
		Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, Rflag_mass, Rflag_parabolic, Ric, GCV,
			RGCVmethod, Rnrealizations, DOF, RDOF_matrix, Rsearch), K_(RK), beta_(Rbeta), c_(REAL(Rc)[0]) {}

RegressionDataEllipticSpaceVarying::RegressionDataEllipticSpaceVarying(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS, SEXP RK, SEXP Rbeta,
				 SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP GCV, SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch):
					 RegressionData(Rlocations, RbaryLocations, Robservations, Rorder, RlambdaS, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, GCV,RGCVmethod, Rnrealizations, DOF, RDOF_matrix, Rsearch),
					 K_(RK), beta_(Rbeta), c_(Rc), u_(Ru)
{}

RegressionDataEllipticSpaceVarying::RegressionDataEllipticSpaceVarying(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS, SEXP RlambdaT, SEXP RK,
		SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric, SEXP GCV, SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch):
					 RegressionData(Rlocations, RbaryLocations, Rtime_locations, Robservations, Rorder, RlambdaS, RlambdaT, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, Rflag_mass, Rflag_parabolic, Ric, GCV, RGCVmethod, Rnrealizations, DOF, RDOF_matrix, Rsearch),
					 K_(RK), beta_(Rbeta), c_(Rc), u_(Ru)
{}

// void RegressionDataEllipticSpaceVarying::print(std::ostream & out) const
// {
// 	for (auto i=0;i<18;i++)
// 		out<<K_(i);
// 	for (auto i=0;i<18;i++)
// 		out<<beta_(i);
// 	for (auto i=0;i<18;i++)
// 		out<<c_(i);
// }

void RegressionData::setObservations(SEXP Robservations)
{
	UInt n_obs_ = Rf_length(Robservations);
	observations_.resize(n_obs_);
	observations_indices_.reserve(n_obs_);

	UInt count = 0;
	if(locations_.nrows() == 0 && nRegions_ == 0)
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
	locations_by_nodes_ = (locations_.nrows()== 0 && nRegions_ == 0) ? true : false;

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

// void RegressionData::setCovariates(SEXP Rcovariates)
// {
// 	n_ = INTEGER(Rf_getAttrib(Rcovariates, R_DimSymbol))[0];
// 	p_ = INTEGER(Rf_getAttrib(Rcovariates, R_DimSymbol))[1];
//
// 	covariates_.resize(n_, p_);
//
// 	for(auto i=0; i<n_; ++i)
// 	{
// 		for(auto j=0; j<p_ ; ++j)
// 		{
// 			covariates_(i,j)=REAL(Rcovariates)[i+ n_*j];
// 		}
// 	}
// }

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

void RegressionData::setDOF_matrix(SEXP RDOF_matrix)
{
	n_ = INTEGER(Rf_getAttrib(RDOF_matrix, R_DimSymbol))[0];
	p_ = INTEGER(Rf_getAttrib(RDOF_matrix, R_DimSymbol))[1];

	dof_matrix_.resize(n_, p_);

	for(auto i=0; i<n_; ++i)
	{
		for(auto j=0; j<p_ ; ++j)
		{
			dof_matrix_(i,j)=REAL(RDOF_matrix)[i+ n_*j];
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
	if (TYPEOF(RbaryLocations) != 0) { //TYPEOF(RbaryLocations) == 0 means SEXPTYPE is NILSXP (Description is NULL)
		UInt* id_ 	= INTEGER(VECTOR_ELT(RbaryLocations, 1));
		Real* bary_ 	= REAL(VECTOR_ELT(RbaryLocations, 2));

		UInt n_ = INTEGER(Rf_getAttrib(VECTOR_ELT(RbaryLocations, 2), R_DimSymbol))[0];
		UInt p_ = INTEGER(Rf_getAttrib(VECTOR_ELT(RbaryLocations, 2), R_DimSymbol))[1]; //barycenter column dimension

		barycenters_.resize(n_, p_);
		element_ids_.resize(n_);

		if(n_>0){
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
	} else {
		locations_by_barycenter_ =false;
	}
}

void RegressionData::setNrealizations(SEXP Rnrealizations) {
	nrealizations_ = INTEGER(Rnrealizations)[0];
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

#endif

void RegressionData::printObservations(std::ostream & out) const
{

	for(auto i=0;i<observations_.size(); i++)
	{
		out<<i<<"\t"<<observations_(i)<<std::endl;
	}
}

void RegressionData::printCovariates(std::ostream & out) const
{

	for(auto i=0;i<covariates_.rows(); i++)
	{
		for(auto j=0; j<covariates_.cols(); j++)
		{
			out<<covariates_(i,j)<<"\t";
		}
		out<<std::endl;
	}
}

void RegressionData::printLocations(std::ostream & out) const
{
	if(locations_.ncols()==2)
		for(UInt i=0; i<locations_.nrows(); i++)
			out<<getLocations<2>(i)<<std::endl;
	else
		for(UInt i=0; i<locations_.nrows(); i++)
			out<<getLocations<3>(i)<<std::endl;

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

#endif
