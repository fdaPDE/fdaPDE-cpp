#include "../Include/DE_Data.h"

DEData::DEData(SEXP Rdata, SEXP Rorder, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter, SEXP Rlambda, SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals,
  SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP Rsearch)
{
  setData(Rdata);

  order_ = INTEGER(Rorder)[0];

  setFvec(Rfvec);

  heatStep_ = REAL(RheatStep)[0];

  heatIter_ = INTEGER(RheatIter)[0];

  setLambda(Rlambda);

  Nfolds_ = INTEGER(Rnfolds)[0];

  nsim_ = INTEGER(Rnsim)[0];

  setStepProposals(RstepProposals);

  tol1_ = REAL(Rtol1)[0];

  tol2_ = REAL(Rtol2)[0];

  print_ = INTEGER(Rprint)[0];

  search_ = INTEGER(Rsearch)[0];

}


DEData::DEData(const std::vector<Point>& data, const UInt& order, const VectorXr& fvec, Real heatStep, UInt heatIter, const std::vector<Real>& lambda,
               const UInt& nfolds, const UInt& nsim, const std::vector<Real>& stepProposals, Real tol1, Real tol2,
               bool print, UInt search):
                data_(data), order_(order), fvec_(fvec), heatStep_(heatStep), heatIter_(heatIter), lambda_(lambda), Nfolds_(nfolds),
                nsim_(nsim), stepProposals_(stepProposals), tol1_(tol1), tol2_(tol2), print_(print), search_(search)
{
    n_ = data.size();
}



void DEData::setData(SEXP Rdata)
{
  n_ = INTEGER(Rf_getAttrib(Rdata, R_DimSymbol))[0];
  data_.reserve(n_);
	if(n_>0){
		UInt ndim = INTEGER(Rf_getAttrib(Rdata, R_DimSymbol))[1];

	  if (ndim == 2){
			for(UInt i=0; i<n_; ++i)
			{
				data_.emplace_back(REAL(Rdata)[i+ n_*0],REAL(Rdata)[i+ n_*1]);
      }
    }
    else {
			for(UInt i=0; i<n_; ++i)
			{
				data_.emplace_back(REAL(Rdata)[i+ n_*0],REAL(Rdata)[i+ n_*1],REAL(Rdata)[i+ n_*2]);
			}
		}
	}
}

void DEData::setFvec(SEXP Rfvec)
{
  UInt dimc = Rf_length(Rfvec);
  fvec_.resize(dimc);
  for(UInt i=0; i<dimc; i++)
  {
      fvec_[i] = REAL(Rfvec)[i];
  }
}

void DEData::setLambda(SEXP Rlambda)
{
  UInt diml = Rf_length(Rlambda);
  lambda_.reserve(diml);
  for(UInt i=0; i<diml; i++)
  {
      lambda_.push_back(REAL(Rlambda)[i]);
  }
}

void DEData::setStepProposals(SEXP RstepProposals)
{
  UInt dimPG = Rf_length(RstepProposals);
  stepProposals_.reserve(dimPG);
  for(UInt i=0; i<dimPG; i++)
  {
      stepProposals_.push_back(REAL(RstepProposals)[i]);
  }
}



void DEData::setNewData(const std::vector<Point>& p)
{
  data_.resize(p.size());
  for(UInt i = 0; i < p.size(); i++){
    data_[i] = p[i];
  }
}

void DEData::setDatum(const Point& p, UInt i)
{
  data_[i] = p;
}

void DEData::updateN(UInt n){
  n_ = n;
}

void DEData::printData(std::ostream & out) const
{
  for(std::vector<Point>::size_type i=0;i<data_.size(); i++)
	{
		data_[i].print(out);
	}
}
