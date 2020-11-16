#ifndef __DATA_PROBLEM_IMP_H__
#define __DATA_PROBLEM_IMP_H__

template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
DataProblem<Integrator_noPoly, ORDER, mydim, ndim>::DataProblem(SEXP Rdata, SEXP Rorder, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter,
  SEXP Rlambda, SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals, SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP Rsearch, SEXP Rmesh):
  deData_(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rsearch),
   mesh_(Rmesh, INTEGER(Rsearch)[0]){

    // PROJECTION

    if(mydim == 2 && ndim == 3){
        Rprintf("##### DATA PROJECTION #####\n");
    }
    
    std::vector<Point<ndim> > data = deData_.getData();
    projection<ORDER, mydim, ndim> projection(mesh_, data);
    std::vector<Point<ndim> > new_data = projection.computeProjection();
    deData_.setNewData(new_data);


    // REMOVE POINTS NOT IN THE DOMAIN
    data = deData_.getData(); // for the 2.5D case
    constexpr UInt Nodes = (mydim==2) ? 3*ORDER : 6*ORDER-2;
    Element<Nodes, mydim, ndim> tri_activated;

    bool check = false;
    for(auto it = data.begin(); it != data.end(); ){
      tri_activated = mesh_.findLocationNaive(data[it - data.begin()]); // cambiare ricerca
      if(tri_activated.getId() == Identifier::NVAL)
      {
        it = data.erase(it);
        Rprintf("WARNING: an observation is not in the domain. It is removed and the algorithm proceeds.\n");

        if(!check) check = true;
      }
      else {
        it++;
      }
    }

    if(check){
      deData_.setNewData(data);
      deData_.updateN(data.size());
    }


    // FILL MATRICES
    fillFEMatrices();
    fillPsiQuad();

    std::vector<UInt> v(deData_.getNumberofData());
    std::iota(v.begin(),v.end(),0);
    GlobalPsi_ = computePsi(v);
}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
void DataProblem<Integrator_noPoly, ORDER, mydim, ndim>::fillFEMatrices(){

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


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
void DataProblem<Integrator_noPoly, ORDER, mydim, ndim>::fillPsiQuad(){

  constexpr UInt Nodes = mydim==2? 3*ORDER : 6*ORDER-2;

	PsiQuad_.resize(Integrator_noPoly::NNODES, Nodes);

	for(UInt i=0; i<Integrator_noPoly::NNODES; i++)
	{
    Eigen::Matrix<Real, Nodes, 1> evaluator=reference_eval_point<Nodes, mydim>(Integrator_noPoly::NODES[i]);
		for(UInt node = 0; node < Nodes ; ++node)
		{
			PsiQuad_(i, node) = evaluator[node];
		}
	}
}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
Real DataProblem<Integrator_noPoly, ORDER, mydim, ndim>::FEintegrate_exponential(const VectorXr& g) const{

  using EigenMap2WEIGHTS = Eigen::Map<const Eigen::Matrix<Real, Integrator_noPoly::NNODES, 1> >;
  constexpr UInt Nodes = mydim==2? 3*ORDER : 6*ORDER-2;

  Real total_sum = 0.;

  for(UInt triangle=0; triangle<mesh_.num_elements(); triangle++){

    Element<Nodes, mydim, ndim> tri_activated = mesh_.getElement(triangle);

// (3) -------------------------------------------------
    VectorXr sub_g(Nodes);
    for (UInt i=0; i<Nodes; i++){
      sub_g[i]=g[tri_activated[i].getId()];
    }

// (4) -------------------------------------------------
    VectorXr expg = (PsiQuad_*sub_g).array().exp();

    total_sum+=expg.dot(EigenMap2WEIGHTS(&Integrator_noPoly::WEIGHTS[0]))*tri_activated.getMeasure();

  }

  return total_sum;
}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
SpMat
DataProblem<Integrator_noPoly, ORDER, mydim, ndim>::computePsi(const std::vector<UInt>& indices) const{

  UInt nnodes = mesh_.num_nodes();
	UInt nlocations = indices.size();
	SpMat psi(nlocations, nnodes);

	Real eps = 2.2204e-016,
		   tolerance = 100 * eps;

	// Constexpr is used for selecting the right number of nodes to pass as a template parameter to the Element object.In case of planar domain(i.e. mydim==2), we have that the number of nodes is 3*ORDER. In case of volumetric domain (i.e. mydim==3), we have that the number of nodes is 4 nodes if ORDER==1 and 10 nodes if ORDER==2, so the expression is 6*ORDER-2. ORDER==2 if mydim==3 is not yet implemented.
	constexpr UInt Nodes = mydim ==2? 3*ORDER : 6*ORDER-2;
	Element<Nodes, mydim, ndim> tri_activated;
	Eigen::Matrix<Real,Nodes,1> coefficients;

	Real evaluator;
	std::vector<coeff> triplets;
	triplets.reserve(Nodes*nlocations);

	for(auto it = indices.cbegin(); it != indices.cend(); it++)
	{

    if (deData_.getSearch() == 1) { //use Naive search
      tri_activated = mesh_.findLocationNaive(deData_.getDatum(*it));
    } else if (deData_.getSearch() == 2) { //use Tree search (default)
      tri_activated = mesh_.findLocationTree(deData_.getDatum(*it));
    }

		if(tri_activated.getId() == Identifier::NVAL)
		{
			Rprintf("WARNING: the following observation is not in the domain\n");
			//(deData_.getDatum(*it)).print(std::cout);
		}
    else
    {
			for(UInt node = 0; node < Nodes ; ++node)
			{
				coefficients = Eigen::Matrix<Real,Nodes,1>::Zero();
				coefficients(node) = 1; //Activates only current base
				evaluator = tri_activated.evaluate_point(deData_.getDatum(*it), coefficients);
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
