#ifndef __MATRIX_ASSEMBLER_IMP_H__
#define __MATRIX_ASSEMBLER_IMP_H__


template<UInt ORDER, UInt mydim, UInt ndim, typename A>
void Assembler::operKernel(EOExpr<A> oper, const MeshHandler<ORDER,mydim,ndim>& mesh,
	                     FiniteElement<ORDER,mydim,ndim>& fe, SpMat& OpMat)
{
	static constexpr Real eps = std::numeric_limits<Real>::epsilon(),
		 tolerance = 10 * eps;
	static constexpr UInt NBASES = FiniteElement<ORDER,mydim,ndim>::NBASES;
	using Integrator = typename FiniteElement<ORDER, mydim, ndim>::Integrator;

	std::vector<coeff> triplets;
	triplets.reserve(NBASES*NBASES*mesh.num_elements());

	std::vector<UInt> identifiers;
	identifiers.reserve(NBASES);

  for(int t=0; t<mesh.num_elements(); ++t){

		fe.updateElement(mesh.getElement(t));

		// Vector of vertices indices (link local to global indexing system)
		for(int i=0; i<NBASES; ++i)
			identifiers.push_back(fe[i].id());

		for(int i=0; i<NBASES; ++i)
			for(int j=0; j<NBASES; ++j)
				{
					Real s=0;
					for(int iq = 0; iq < Integrator::NNODES; ++iq)
						s += oper(fe, iq, i, j) * Integrator::WEIGHTS[iq];
					s *= fe.getMeasure();
					triplets.emplace_back(identifiers[i], identifiers[j], s);
				}

		identifiers.clear();
	}

  const UInt nnodes = mesh.num_nodes();
  OpMat.resize(nnodes, nnodes);
	OpMat.setFromTriplets(triplets.begin(),triplets.end());
	OpMat.prune(tolerance);
}



template<UInt DEGREE, UInt ORDER_DERIVATIVE>
void Assembler::operKernel(const Spline<DEGREE, ORDER_DERIVATIVE>& spline, SpMat& OpMat)
{
	using Integrator=typename Spline<DEGREE, ORDER_DERIVATIVE>::Integrator;
    const UInt M = spline.num_knots()-DEGREE-1;
  	OpMat.resize(M, M);

    for (UInt i = 0; i < M; ++i){
        for (UInt j = 0; j <= i; ++j){
            Real s = 0;
            for(UInt k = i; k <= j+DEGREE; ++k){
                Real a = spline.getKnot(k);
                Real b = spline.getKnot(k+1);
                for (UInt l = 0; l < Integrator::NNODES; ++l)
                    s += spline.time_mass_impl(i, j, (b-a)/2*Integrator::NODES[l]+(b+a)/2) * Integrator::WEIGHTS[l] * (b-a)/2;
            }

        	if(s!=0){
        		OpMat.coeffRef(i,j) = s;
        		if(i!=j) 
         			OpMat.coeffRef(j,i) = s;
        	}
       	}
    }
}

template<UInt ORDER, UInt mydim, UInt ndim>
void Assembler::forcingTerm(const MeshHandler<ORDER,mydim,ndim>& mesh,
	                     FiniteElement<ORDER,mydim,ndim>& fe, const ForcingTerm& u, VectorXr& forcingTerm)
{
	static constexpr UInt NBASES = FiniteElement<ORDER,mydim,ndim>::NBASES;

	forcingTerm = VectorXr::Zero(mesh.num_nodes());

  	for(int t=0; t<mesh.num_elements(); ++t){

		fe.updateElement(mesh.getElement(t));

		for(int i=0; i<NBASES; ++i)
			forcingTerm[fe[i].id()] += u.integrate(fe, i) * fe.getMeasure();
	}
}

#endif
