#ifndef __MATRIX_ASSEMBLER_IMP_H__
#define __MATRIX_ASSEMBLER_IMP_H__


template<UInt ORDER, typename Integrator, typename A>
void Assembler::operKernel(EOExpr<A> oper,const MeshHandler<ORDER,2,2>& mesh,
	                     FiniteElement<Integrator, ORDER,2,2>& fe, SpMat& OpMat)
{
	Real eps = 2.2204e-016,
		 tolerance = 10 * eps;
	std::vector<coeff> triplets;


  	for(auto t=0; t<mesh.num_elements(); t++)
  	{
		fe.updateElement(mesh.getElement(t));

		// Vector of vertices indices (link local to global indexing system)
		std::vector<UInt> identifiers;
		identifiers.resize(3*ORDER);
		for( auto q=0; q<3*ORDER; q++)
			identifiers[q]=mesh.getElement(t)[q].id();

		//localM=localMassMatrix(currentelem);
		for(int i = 0; i < 3*ORDER; i++)
		{
			for(int j = 0; j < 3*ORDER; j++)
			{
				Real s=0;

				for(int l = 0;l < Integrator::NNODES; l++)
				{
					s += oper(fe,i,j,l) * fe.getDet() * fe.getAreaReference() * Integrator::WEIGHTS[l];
				}
			  triplets.push_back(coeff(identifiers[i],identifiers[j],s));
			}
		}
	}

  	UInt nnodes = mesh.num_nodes();
  	OpMat.resize(nnodes, nnodes);
	OpMat.setFromTriplets(triplets.begin(),triplets.end());
	OpMat.prune(tolerance);
}

template<UInt ORDER, typename Integrator>
void Assembler::forcingTerm(const MeshHandler<ORDER,2,2>& mesh,
	                     FiniteElement<Integrator, ORDER,2,2>& fe, const ForcingTerm& u, VectorXr& forcingTerm)
{

	forcingTerm = VectorXr::Zero(mesh.num_nodes());

  	for(auto t=0; t<mesh.num_elements(); t++)
  	{
		fe.updateElement(mesh.getElement(t));

		// Vector of vertices indices (link local to global indexing system)
		std::vector<UInt> identifiers;
				identifiers.resize(3*ORDER);

		for( auto q=0; q<3*ORDER; q++)
			identifiers[q]=mesh.getElement(t)[q].id();


		//localM=localMassMatrix(currentelem);
		for(int i = 0; i < 3*ORDER; i++)
		{
			Real s=0;

			for(int iq = 0;iq < Integrator::NNODES; iq++)
			{
				UInt globalIndex = fe.getGlobalIndex(iq);
				s +=  fe.phiMaster(i,iq)* u(globalIndex) * fe.getDet() * fe.getAreaReference()* Integrator::WEIGHTS[iq];//(*)
			}
			forcingTerm[identifiers[i]] += s;
		}

	}
}


//! Surface mesh implementation

template<UInt ORDER, typename Integrator, typename A>
void Assembler::operKernel(EOExpr<A> oper,const MeshHandler<ORDER,2,3>& mesh,
	                     FiniteElement<Integrator, ORDER,2,3>& fe, SpMat& OpMat)
{
	Real eps = 2.2204e-016,
		 tolerance = 10 * eps;
	std::vector<coeff> triplets;


  	for(auto t=0; t<mesh.num_elements(); t++)
  	{
		fe.updateElement(mesh.getElement(t));

		// Vector of vertices indices (link local to global indexing system)
		std::vector<UInt> identifiers;
		identifiers.resize(3*ORDER);
		for( auto q=0; q<3*ORDER; q++)
			identifiers[q]=mesh.getElement(t)[q].id();

		//localM=localMassMatrix(currentelem);
		for(int i = 0; i < 3*ORDER; i++)
		{
			for(int j = 0; j < 3*ORDER; j++)
			{
				Real s=0;

				for(int l = 0;l < Integrator::NNODES; l++)
				{
					s += oper(fe,i,j,l) * std::sqrt(fe.getDet()) * fe.getAreaReference()* Integrator::WEIGHTS[l];
				}
			  triplets.push_back(coeff(identifiers[i],identifiers[j],s));
			}
		}

	}

  	UInt nnodes = mesh.num_nodes();
  	OpMat.resize(nnodes, nnodes);
	OpMat.setFromTriplets(triplets.begin(),triplets.end());
	OpMat.prune(tolerance);
}

template<UInt DEGREE, UInt ORDER_DERIVATIVE, typename Integrator, typename A>
void Assembler::operKernel(EOExpr<A> oper, Spline<Integrator, DEGREE, ORDER_DERIVATIVE>& spline, SpMat& OpMat)
{
    UInt M = spline.num_knots()-DEGREE-1;
  	OpMat.resize(M, M);

    for (UInt i = 0; i < M; ++i)
    {
        for (UInt j = 0; j <= i; ++j)
        {
            Real s = 0;

            for(UInt k = i; k <= j+DEGREE; ++k)
            {
                Real a = spline.getKnot(k);
                Real b = spline.getKnot(k+1);

                for (UInt l = 0; l < Integrator::NNODES; ++l)
                    s += oper(spline, i, j, (b-a)/2*Integrator::NODES[l]+(b+a)/2) * Integrator::WEIGHTS[l] * (b-a)/2;
            }

         if(s!=0) OpMat.coeffRef(i,j) = s;
         if(i!=j && s!=0) OpMat.coeffRef(j,i) = s;
        }
    }
}


template<UInt ORDER, typename Integrator>
void Assembler::forcingTerm(const MeshHandler<ORDER,2,3>& mesh,
	                     FiniteElement<Integrator, ORDER,2,3>& fe, const ForcingTerm& u, VectorXr& forcingTerm)
{

	forcingTerm = VectorXr::Zero(mesh.num_nodes());

  	for(auto t=0; t<mesh.num_elements(); t++)
  	{
		fe.updateElement(mesh.getElement(t));

		// Vector of vertices indices (link local to global indexing system)
		std::vector<UInt> identifiers;
				identifiers.resize(3*ORDER);

		for( auto q=0; q<3*ORDER; q++)
			identifiers[q]=mesh.getElement(t)[q].id();


		//localM=localMassMatrix(currentelem);
		for(int i = 0; i < 3*ORDER; i++)
		{
			Real s=0;

			for(int iq = 0;iq < Integrator::NNODES; iq++)
			{
				UInt globalIndex = fe.getGlobalIndex(iq);
				s +=  fe.phiMaster(i,iq)* u(globalIndex) * std::sqrt(fe.getDet()) * fe.getAreaReference()* Integrator::WEIGHTS[iq];//(*)
			}
			forcingTerm[identifiers[i]] += s;
		}

	}

}

//! Volume mesh implementation

template<UInt ORDER, typename Integrator, typename A>
void Assembler::operKernel(EOExpr<A> oper,const MeshHandler<ORDER,3,3>& mesh,
	                     FiniteElement<Integrator, ORDER,3,3>& fe, SpMat& OpMat)
{
	Real eps = 2.2204e-016,
		 tolerance = 10 * eps;
	std::vector<coeff> triplets;


  	for(auto t=0; t<mesh.num_elements(); t++)
  	{
		fe.updateElement(mesh.getElement(t));

		// Vector of vertices indices (link local to global indexing system)
		std::vector<UInt> identifiers;
		identifiers.resize(6*ORDER-2);
		for( auto q=0; q<6*ORDER-2; q++)
			identifiers[q]=mesh.getElement(t)[q].id();

		//localM=localMassMatrix(currentelem);
		for(int i = 0; i < 6*ORDER-2; i++)
		{
			for(int j = 0; j < 6*ORDER-2; j++)
			{
				Real s=0;

				for(int l = 0;l < Integrator::NNODES; l++)
				{
					s += oper(fe,i,j,l) * std::sqrt(fe.getDet()) * fe.getVolumeReference()* Integrator::WEIGHTS[l];
				}
			  triplets.push_back(coeff(identifiers[i],identifiers[j],s));
			}
		}

	}

  	UInt nnodes = mesh.num_nodes();
  	OpMat.resize(nnodes, nnodes);
	OpMat.setFromTriplets(triplets.begin(),triplets.end());
	OpMat.prune(tolerance);
}



template<UInt ORDER, typename Integrator>
void Assembler::forcingTerm(const MeshHandler<ORDER,3,3>& mesh,
	                     FiniteElement<Integrator, ORDER,3,3>& fe, const ForcingTerm& u, VectorXr& forcingTerm)
{

	forcingTerm = VectorXr::Zero(mesh.num_nodes());

  	for(auto t=0; t<mesh.num_elements(); t++)
  	{
		fe.updateElement(mesh.getElement(t));

		// Vector of vertices indices (link local to global indexing system)
		std::vector<UInt> identifiers;
				identifiers.resize(6*ORDER-2);

		for( auto q=0; q<6*ORDER-2; q++)
			identifiers[q]=mesh.getElement(t)[q].id();


		//localM=localMassMatrix(currentelem);
		for(int i = 0; i < 6*ORDER-2; i++)
		{
			Real s=0;

			for(int iq = 0;iq < Integrator::NNODES; iq++)
			{
				UInt globalIndex = fe.getGlobalIndex(iq);
				s +=  fe.phiMaster(i,iq)* u(globalIndex) * std::sqrt(fe.getDet()) * fe.getVolumeReference()* Integrator::WEIGHTS[iq];//(*)
			}
			forcingTerm[identifiers[i]] += s;
		}

	}

}




#endif
