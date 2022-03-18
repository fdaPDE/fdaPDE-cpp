#ifndef __DATA_PROBLEM_IMP_H__
#define __DATA_PROBLEM_IMP_H__

template<UInt ORDER, UInt mydim, UInt ndim>
DataProblem<ORDER, mydim, ndim>::DataProblem(SEXP Rdata, SEXP Rorder, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter,
  SEXP Rlambda, SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals, SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP Rsearch, SEXP Rmesh):
  deData_(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rsearch),
   mesh_(Rmesh, INTEGER(Rsearch)[0]){

    std::vector<Point<ndim> >& data = deData_.data();

    // PROJECTION
    if((mydim == 2 && ndim == 3) || (mydim == 1 && ndim == 2)){
      Rprintf("##### DATA PROJECTION #####\n");
      projection<ORDER, mydim, ndim> projection(mesh_, data);
      data = projection.computeProjection();
    }
    
    // REMOVE POINTS NOT IN THE DOMAIN
    for(auto it = data.begin(); it != data.end(); ){
      Element<EL_NNODES, mydim, ndim> tri_activated = mesh_.findLocation(data[it - data.begin()]); 
      if(tri_activated.getId() == Identifier::NVAL)
      {
        it = data.erase(it);
        Rprintf("WARNING: an observation is not in the domain. It is removed and the algorithm proceeds.\n");
      }
      else {
        it++;
      }
    }

    // FILL MATRICES
    fillFEMatrices();
    fillPsiQuad();

    std::vector<UInt> v(deData_.dataSize());
    std::iota(v.begin(),v.end(),0);
    GlobalPsi_ = computePsi(v);
}


template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem<ORDER, mydim, ndim>::fillFEMatrices(){

  //fill R0 and R1
  FiniteElement<ORDER, mydim, ndim> fe;
  typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
  typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
  Assembler::operKernel(mass, mesh_, fe, R0_);
  Assembler::operKernel(stiff, mesh_, fe, R1_);

  //fill P
  Eigen::SparseLU<SpMat> solver;
	solver.compute(R0_);
	auto X2 = solver.solve(R1_);
	P_ = R1_.transpose()* X2;
}


template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem<ORDER, mydim, ndim>::fillPsiQuad(){
	for(UInt i=0; i<Integrator::NNODES; ++i)
	   PsiQuad_.row(i)=reference_eval_point<EL_NNODES, mydim>(Integrator::NODES[i]);
}


template<UInt ORDER, UInt mydim, UInt ndim>
Real DataProblem<ORDER, mydim, ndim>::FEintegrate_exponential(const VectorXr& g) const{

  using EigenMap2WEIGHTS = Eigen::Map<const Eigen::Matrix<Real, Integrator::NNODES, 1> >;

  Real total_sum = 0.;

  for(UInt triangle=0; triangle<mesh_.num_elements(); ++triangle){

    Element<EL_NNODES, mydim, ndim> tri_activated = mesh_.getElement(triangle);

// (3) -------------------------------------------------
    Eigen::Matrix<Real,EL_NNODES,1> sub_g;
    for (UInt i=0; i<EL_NNODES; i++){
      sub_g[i]=g[tri_activated[i].getId()];
    }

// (4) -------------------------------------------------
    Eigen::Matrix<Real,Integrator::NNODES,1> expg = (PsiQuad_*sub_g).array().exp();

    total_sum+=expg.dot(EigenMap2WEIGHTS(&Integrator::WEIGHTS[0]))*tri_activated.getMeasure();

  }

  return total_sum;
}


template<UInt ORDER, UInt mydim, UInt ndim>
SpMat
DataProblem<ORDER, mydim, ndim>::computePsi(const std::vector<UInt>& indices) const{
  
  static constexpr Real eps = std::numeric_limits<Real>::epsilon(),
     tolerance = 100 * eps;

  UInt nnodes = mesh_.num_nodes();
	UInt nlocations = indices.size();
	SpMat psi(nlocations, nnodes);

	std::vector<coeff> triplets;
	triplets.reserve(EL_NNODES*nlocations);

	for(auto it = indices.cbegin(); it != indices.cend(); it++)
	{

    Element<EL_NNODES, mydim, ndim> tri_activated = mesh_.findLocation(deData_.data(*it));

		if(tri_activated.getId() == Identifier::NVAL)
		{
			Rprintf("WARNING: the following observation is not in the domain\n");
			//(deData_.getDatum(*it)).print(std::cout);
		}
    else
    {
			for(UInt node = 0; node < EL_NNODES ; ++node)
			{
				Real evaluator = tri_activated.evaluate_point(deData_.data(*it), Eigen::Matrix<Real,EL_NNODES,1>::Unit(node));
				triplets.emplace_back(it-indices.cbegin(), tri_activated[node].getId(), evaluator);
			}
		}
	}

	psi.setFromTriplets(triplets.begin(),triplets.end());

	psi.prune(tolerance);
	psi.makeCompressed();

	return psi;
}

#endif
