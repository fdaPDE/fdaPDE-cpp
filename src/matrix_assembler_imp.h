#ifndef MATRIX_ASSEMBLER_IMP_H_
#define MATRIX_ASSEMBLER_IMP_H_


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

	Eigen::Matrix<Real,NBASES,NBASES> loc_matr = Eigen::Matrix<Real,NBASES,NBASES>::Zero();

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
					triplets.push_back(coeff(identifiers[i],identifiers[j],s));
				}

		identifiers.clear();
	}

  const UInt nnodes = mesh.num_nodes();
  OpMat.resize(nnodes, nnodes);
	OpMat.setFromTriplets(triplets.begin(),triplets.end());
	OpMat.prune(tolerance);
}

template<UInt ORDER, UInt mydim, UInt ndim>
void Assembler::forcingTerm(const MeshHandler<ORDER,mydim,ndim>& mesh,
	                     FiniteElement<ORDER,mydim,ndim>& fe, const ForcingTerm& u, VectorXr& forcingTerm)
{
	static constexpr UInt NBASES = FiniteElement<ORDER,mydim,ndim>::NBASES;
	using Integrator = typename FiniteElement<ORDER, mydim, ndim>::Integrator;

	forcingTerm = VectorXr::Zero(mesh.num_nodes());

  for(int t=0; t<mesh.num_elements(); ++t){

		fe.updateElement(mesh.getElement(t));

		for(int i=0; i<NBASES; ++i){
			Real s=0;
			for(int iq = 0; iq < Integrator::NNODES; ++iq){
				UInt globalIndex = fe.getGlobalIndex(iq);
				s +=  fe.getPhi(i,iq) * u[globalIndex] * Integrator::WEIGHTS[iq];
			}
			forcingTerm[fe[i].id()] += s * fe.getMeasure();
		}
	}
}

#endif
