#ifndef _DATA_PROBLEM_IMP_HPP_
#define _DATA_PROBLEM_IMP_HPP_

template<typename Integrator, typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
DataProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim>::DataProblem(SEXP Rdata, SEXP Rorder, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter,
  SEXP Rlambda, SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals, SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP Rsearch, SEXP Rmesh):
  deData_(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rsearch),
   mesh_(Rmesh){

    // PROJECTION
    #ifdef R_VERSION_
      if(mydim == 2 && ndim == 3){
        Rprintf("##### DATA PROJECTION #####\n");
      }
    #endif
    std::vector<Point> data = deData_.getData();
    projection<ORDER, mydim, ndim> projection(mesh_, data);
    std::vector<Point> new_data = projection.computeProjection();
    deData_.setNewData(new_data);


    // REMOVE POINTS NOT IN THE DOMAIN
    data = deData_.getData(); // for the 2.5D case
    constexpr UInt Nodes = mydim ==2? 3*ORDER : 6*ORDER-2;
    Element<Nodes, mydim, ndim> tri_activated;

    bool check = false;
    for(auto it = data.begin(); it != data.end(); ){
      tri_activated = mesh_.findLocationNaive(data[it - data.begin()]); // cambiare ricerca
      if(tri_activated.getId() == Identifier::NVAL)
      {
        it = data.erase(it);
        #ifdef R_VERSION_
        Rprintf("WARNING: an observation is not in the domain. It is removed and the algorithm proceeds.\n");
        #else
        std::cout << "WARNING: an observation is not in the domain. It is removed and the algorithm proceeds.\n";
        #endif

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


template<typename Integrator, typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
void DataProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim>::fillFEMatrices(){

  //fill R0 and R1
  FiniteElement<Integrator, ORDER, mydim, ndim> fe;
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


template<typename Integrator, typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
void DataProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim>::fillPsiQuad(){

  constexpr UInt Nodes = mydim==2? 3*ORDER : 6*ORDER-2;

	PsiQuad_.resize(Integrator_noPoly::NNODES, Nodes);

	//Set the properties of the reference element
	std::vector<Point> reference_nodes;

  if (mydim==2){
  	reference_nodes.push_back(Point(0,0));
  	reference_nodes.push_back(Point(1,0));
  	reference_nodes.push_back(Point(0,1));
  }
  else if (mydim==3){
    reference_nodes.push_back(Point(0,0,0));
  	reference_nodes.push_back(Point(1,0,0));
  	reference_nodes.push_back(Point(0,1,0));
    reference_nodes.push_back(Point(0,0,1));
  }

	Element<Nodes,mydim,ndim> referenceElem = Element<Nodes,mydim,ndim> (Id(0), reference_nodes);

	Eigen::Matrix<Real,Nodes,1> coefficients;

	Real evaluator;

	for(UInt i=0; i<Integrator_noPoly::NNODES; i++)
	{
			for(UInt node = 0; node < Nodes ; ++node)
			{
				coefficients = Eigen::Matrix<Real,Nodes,1>::Zero();
				coefficients(node) = 1; //Activates only current base
				evaluator = evaluate_point<Nodes, mydim, ndim>(referenceElem, Integrator_noPoly::NODES[i], coefficients);
				PsiQuad_(i, node) = evaluator;
			}
	}
}


template<typename Integrator, typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
Real DataProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim>::FEintegrate_exponential(const VectorXr& g) const{

  Real total_sum = 0.;

  constexpr UInt Nodes = mydim==2? 3*ORDER : 6*ORDER-2;

  for(UInt triangle=0; triangle<mesh_.num_elements(); triangle++){

    FiniteElement<Integrator_noPoly, ORDER, mydim, ndim> fe;
    Element<Nodes, mydim, ndim> tri_activated = mesh_.getElement(triangle);
    fe.updateElement(tri_activated);

// (3) -------------------------------------------------
    VectorXr sub_g(Nodes);
    for (UInt i=0; i<Nodes; i++){
      sub_g[i]=g[tri_activated[i].getId()];
    }

// (4) -------------------------------------------------
    VectorXr expg = (PsiQuad_*sub_g).array().exp();

    // mind we are using quadrature rules whom weights sum to the element measure.
    if (ndim==2){
      total_sum+=expg.dot(Integrator_noPoly::WEIGHTS)*std::abs(fe.getDet());
    }
    else if (ndim==3){
      total_sum+=expg.dot(Integrator_noPoly::WEIGHTS)*std::sqrt(std::abs(fe.getDet()));
    }
  }

  return total_sum;
}


template<typename Integrator, typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
SpMat
DataProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim>::computePsi(const std::vector<UInt>& indices) const{

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
			#ifdef R_VERSION_
			Rprintf("WARNING: the following observation is not in the domain\n");
			#else
			std::cout << "WARNING: The following observation is not in the domain\n";
      #endif
			//(deData_.getDatum(*it)).print(std::cout);
		}
    else
    {
			for(UInt node = 0; node < Nodes ; ++node)
			{
				coefficients = Eigen::Matrix<Real,Nodes,1>::Zero();
				coefficients(node) = 1; //Activates only current base
				evaluator = evaluate_point<Nodes, mydim, ndim>(tri_activated, deData_.getDatum(*it), coefficients);
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
