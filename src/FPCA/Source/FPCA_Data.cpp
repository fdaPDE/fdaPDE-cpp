#include "../Include/FPCA_Data.h"

FPCAData::FPCAData(Real* locations, UInt n_locations, UInt ndim, MatrixXr& datamatrix, UInt order, MatrixXi& incidenceMatrix,
					std::vector<Real> lambda, UInt nPC, UInt nFolds, UInt search):
					locations_(locations, n_locations, ndim), datamatrix_(datamatrix), order_(order),
					incidenceMatrix_(incidenceMatrix), lambda_(lambda),  nPC_(nPC),
					nFolds_(nFolds), search_(search)
{
	nRegions_ = incidenceMatrix.rows();
	if(locations_.nrows()==0 && nRegions_==0)
	{
		locations_by_nodes_ = true;
		for(int i = 0; i<datamatrix_.cols();++i) observations_indices_.push_back(i);
	} else
		locations_by_nodes_ = false;
}

FPCAData::FPCAData(SEXP Rlocations, SEXP RbaryLocations, SEXP Rdatamatrix, SEXP Rorder, SEXP RincidenceMatrix, SEXP Rlambda,
					SEXP RnPC, SEXP RnFolds,SEXP RGCVmethod, SEXP Rnrealizations, SEXP Rsearch) :
	locations_(Rlocations)
{


	setBaryLocations(RbaryLocations);
	setIncidenceMatrix(RincidenceMatrix);
	setDatamatrix(Rdatamatrix);
	setNrealizations(Rnrealizations);

	GCVmethod_ = INTEGER(RGCVmethod)[0];

	order_ =  INTEGER(Rorder)[0];
	search_ =  INTEGER(Rsearch)[0];

	UInt length_lambda = Rf_length(Rlambda);
	for (UInt i = 0; i<length_lambda; ++i) lambda_.push_back(REAL(Rlambda)[i]);

	nPC_ = INTEGER(RnPC)[0];

	nFolds_=INTEGER(RnFolds)[0];
}


void FPCAData::setBaryLocations(SEXP RbaryLocations)
{
	if (TYPEOF(RbaryLocations) != 0) { //TYPEOF(RbaryLocations) == 0 means SEXPTYPE is NILSXP (Description is NULL)
		Real* bary_ 	= REAL(VECTOR_ELT(RbaryLocations, 0));
		UInt* id_ 	= INTEGER(VECTOR_ELT(RbaryLocations, 1));
		UInt n_ = INTEGER(Rf_getAttrib(VECTOR_ELT(RbaryLocations, 0), R_DimSymbol))[0];
		UInt p_ = INTEGER(Rf_getAttrib(VECTOR_ELT(RbaryLocations, 0), R_DimSymbol))[1]; //barycenter column dimension

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


void FPCAData::setDatamatrix(SEXP Rdatamatrix)
{
	n_ = INTEGER(Rf_getAttrib(Rdatamatrix, R_DimSymbol))[0];
	p_ = INTEGER(Rf_getAttrib(Rdatamatrix, R_DimSymbol))[1];
	datamatrix_.resize(n_,p_);
	observations_indices_.reserve(p_);
	VectorXr auxiliary_row_;
	auxiliary_row_.resize(p_);

	nRegions_ = incidenceMatrix_.rows();

	if(locations_.nrows() == 0 && nRegions_==0)
	{
		locations_by_nodes_ = true;
		for(auto i=0; i<n_; ++i)
		{
			UInt count=0;
			for(auto j=0; j<p_ ; ++j)
			{
				if(!ISNA(REAL(Rdatamatrix)[i+n_*j]))
				{
					auxiliary_row_[count]=REAL(Rdatamatrix)[i+n_*j];
					count++;
					if(i==0) observations_indices_.push_back(j);
				}
			}
			datamatrix_.row(i)=auxiliary_row_;
		}
		datamatrix_.conservativeResize(Eigen::NoChange,observations_indices_.size());
	} else {
		locations_by_nodes_ = false;
		for(auto i=0; i<n_; ++i)
		{
			for(auto j=0; j<p_ ; ++j)
			{
				datamatrix_(i,j)=REAL(Rdatamatrix)[i+n_*j];
			}
		}
	}

}


void FPCAData::setNrealizations(SEXP Rnrealizations) {
	nrealizations_ = INTEGER(Rnrealizations)[0];
}

void FPCAData::setIncidenceMatrix(SEXP RincidenceMatrix)
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


void FPCAData::printDatamatrix(std::ostream & out) const
{

	for(auto i=0;i<datamatrix_.rows(); i++)
	{
		for(auto j=0;j<datamatrix_.cols();j++)
		{
		out<<datamatrix_(i,j)<<"\t";
		}
		out<<std::endl;
	}
}


void FPCAData::printLocations(std::ostream & out) const
{
	if(locations_.ncols()==2)
		for(UInt i=0; i<locations_.nrows(); i++)
			out<<getLocations<2>(i)<<std::endl;
	else
		for(UInt i=0; i<locations_.nrows(); i++)
			out<<getLocations<3>(i)<<std::endl;

}


void FPCAData::printIncidenceMatrix(std::ostream & out) const
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

