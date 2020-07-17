#ifndef __MIXEDFEREGRESSION_IMP_HPP__
#define __MIXEDFEREGRESSION_IMP_HPP__

#include <iostream>
#include <chrono>
#include <random>
#include <fstream>

#include "R_ext/Print.h"

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::addDirichletBC()
{
	UInt id1,id3;

	UInt nnodes = N_*M_;

	const std::vector<UInt>& bc_indices = regressionData_.getDirichletIndices();
	const std::vector<Real>& bc_values = regressionData_.getDirichletValues();
	UInt nbc_indices = bc_indices.size();

	Real pen=10e20;

	for( auto i=0; i<nbc_indices; i++)
	 {
			id1=bc_indices[i];
			id3=id1+nnodes;

			matrixNoCov_.coeffRef(id1,id1)=pen;
			matrixNoCov_.coeffRef(id3,id3)=pen;


			_rightHandSide(id1)=bc_values[i]*pen;
			_rightHandSide(id3)=0;
	 }

	matrixNoCov_.makeCompressed();
}

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::addNA()
{

	const std::vector<UInt>& observations_na= regressionData_.getObservationsNA();

	for(UInt id:observations_na)
	{
		for(UInt j=0; j<psi_.cols(); ++j)
		{
			if(psi_.coeff(id,j)!=0)
				psi_.coeffRef(id, j) = 0;
		}
	}
	psi_.pruned();
	psi_.makeCompressed();
}

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::setPsi()
{
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.isSpaceTime() ? regressionData_.getNumberofSpaceObservations() : regressionData_.getNumberofObservations();

	psi_.resize(nlocations, nnodes);
	if (regressionData_.isLocationsByNodes() & !regressionData_.isLocationsByBarycenter()) //pointwise data
	{
		std::vector<coeff> tripletAll;
		if(!regressionData_.isSpaceTime())
		{
			auto k = regressionData_.getObservationsIndices();
			tripletAll.reserve(k.size());
			for (int i = 0; i< k.size(); ++i){
				tripletAll.push_back(coeff(i,k[i],1.0));
			}
		}
		else
		{
			tripletAll.reserve(nlocations);
			for (int i = 0; i< nlocations; ++i){
				tripletAll.push_back(coeff(i,i,1.0));
			}
		}
		psi_.setFromTriplets(tripletAll.begin(),tripletAll.end());
		psi_.makeCompressed();
	}
	else if (regressionData_.isLocationsByBarycenter() && (regressionData_.getNumberOfRegions() == 0)) //pointwise data
	{
		constexpr UInt Nodes = mydim==2 ? 3*ORDER : 6*ORDER-2;
		Element<Nodes, mydim, ndim> tri_activated;
		Real evaluator;

		for(UInt i=0; i<nlocations;i++)
		{
			tri_activated = mesh_.getElement(regressionData_.getElementId(i));

			if(tri_activated.getId() == Identifier::NVAL)
			{
				#ifdef R_VERSION_
				Rprintf("ERROR: Point %d is not in the domain, remove point and re-perform smoothing\n", i+1);
				#else
				std::cout << "ERROR: Point " << i+1 <<" is not in the domain\n";
				#endif
			}else //tri_activated.getId() found
			{
				for(UInt node = 0; node < Nodes ; ++node)
				{
					evaluator = regressionData_.getBarycenter(i,node);
					psi_.insert(i, tri_activated[node].getId()) = evaluator;
				}
			}
		} //end of for loop
		psi_.makeCompressed();
	}
	else if ((!regressionData_.isLocationsByBarycenter()) && (regressionData_.getNumberOfRegions() == 0))
	{
		constexpr UInt Nodes = mydim==2 ? 3*ORDER : 6*ORDER-2;
		Element<Nodes, mydim, ndim> tri_activated;
		Eigen::Matrix<Real,Nodes,1> coefficients;

		Real evaluator;
		this->barycenters_.resize(nlocations, Nodes);
		this->element_ids_.resize(nlocations);
		for(UInt i=0; i<nlocations;i++)
		{
			if (regressionData_.getSearch() == 1) { //use Naive search
				tri_activated = mesh_.findLocationNaive(regressionData_.getLocations()[i]);
			} else if (regressionData_.getSearch() == 2) { //use Tree search (default)
				tri_activated = mesh_.findLocationTree(regressionData_.getLocations()[i]);
			}

			if(tri_activated.getId() == Identifier::NVAL)
			{
				#ifdef R_VERSION_
				Rprintf("ERROR: Point %d is not in the domain, remove point and re-perform smoothing\n", i+1);
				#else
				std::cout << "ERROR: Point " << i+1 <<" is not in the domain\n";
				#endif
			}else //tri_activated.getId() found
			{
				element_ids_(i)=tri_activated.getId();
				for(UInt node = 0; node < Nodes ; ++node)
				{
					coefficients = Eigen::Matrix<Real,Nodes,1>::Zero();
					coefficients(node) = 1; //Activates only current base
					evaluator = evaluate_point<Nodes,mydim,ndim>(tri_activated, regressionData_.getLocations()[i], coefficients);
					barycenters_(i,node)=tri_activated.getBaryCoordinates(regressionData_.getLocations()[i])[node];
					psi_.insert(i, tri_activated[node].getId()) = evaluator;
				}
			}
		} //end of for loop
		psi_.makeCompressed();
	}
	else //areal data
	{
		constexpr UInt Nodes = mydim==2 ? 3*ORDER : 6*ORDER-2;

		Real *tab; //Psi_i
		tab = (Real*) malloc(sizeof(Real)*nnodes);
		for(UInt i=0; i<nlocations;i++) //nlocations = number of regions
		{
			for (UInt k=0; k<nnodes; k++) {tab[k]=0;}
			for (UInt j=0; j<mesh_.num_elements(); j++)
			{
				if (regressionData_.getIncidenceMatrix()(i,j) == 1) //element j is in region i
				{
					Element<Nodes, mydim, ndim> tri = mesh_.getElement(j);
					for (UInt k=0; k<Nodes; k++)
					{
						tab[tri[k].getId()] += integratePsi(tri,k); // integral over tri of psi_k
					}
				}
			}
			for (int k=0; k<nnodes; k++)
			{
				if (tab[k] != 0)
				{
					psi_.insert(i,k) = tab[k]/A_(i); //divide by |D_i|
				}
			}
		}
		free(tab);
		psi_.makeCompressed();
	}
}
template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
MatrixXr MixedFERegressionBase<InputHandler,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::LeftMultiplybyQ(const MatrixXr& u)
{
	if (regressionData_.getCovariates().rows() == 0){
		return u;
	}
	else{
		MatrixXr W(this->regressionData_.getCovariates());
		if (isWTWfactorized_ == false ){
			WTW_.compute(W.transpose()*W);
			isWTWfactorized_=true;
		}
		MatrixXr Pu= W*WTW_.solve(W.transpose()*u);
		return u-Pu;
	}

}

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::buildMatrixNoCov(const SpMat& DMat,  const SpMat& R1,  const SpMat& R0)
{

	UInt nnodes = N_*M_;

	std::vector<coeff> tripletAll;
	tripletAll.reserve(DMat.nonZeros() + 2*R1.nonZeros() + R0.nonZeros());

	for (int k=0; k<DMat.outerSize(); ++k)
		for (SpMat::InnerIterator it(DMat,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row(), it.col(),it.value()));
		}
	for (int k=0; k<R0.outerSize(); ++k)
		for (SpMat::InnerIterator it(R0,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row()+nnodes, it.col()+nnodes,it.value()));
		}
	for (int k=0; k<R1.outerSize(); ++k)
	  for (SpMat::InnerIterator it(R1,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.col(), it.row()+nnodes,it.value()));
	  }
	for (int k=0; k<R1.outerSize(); ++k)
	  for (SpMat::InnerIterator it(R1,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.row()+nnodes, it.col(), it.value()));
	  }

	matrixNoCov_.setZero();
	matrixNoCov_.resize(2*nnodes,2*nnodes);
	matrixNoCov_.setFromTriplets(tripletAll.begin(),tripletAll.end());
	matrixNoCov_.makeCompressed();
}

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::system_factorize() {

	UInt nnodes = N_*M_;

	// First phase: Factorization of matrixNoCov
	matrixNoCovdec_.compute(matrixNoCov_);

	if (regressionData_.getCovariates().rows() != 0) {
		// Second phase: factorization of matrix  G =  C + [V * matrixNoCov^-1 * U]= C + D

		// Definition of matrix U = [ psi^T * A * W | 0 ]^T and V= [ W^T*psi| 0]

		MatrixXr W(this->regressionData_.getCovariates());

		U_ = MatrixXr::Zero(2*nnodes, W.cols());

		V_=MatrixXr::Zero(W.cols(),2*nnodes);
		V_.leftCols(nnodes)=W.transpose()*psi_;

		if(regressionData_.getNumberOfRegions()==0){ // pointwise data
		  U_.topRows(nnodes) = psi_.transpose()*W;
		}
		else{                                          //areal data
		  U_.topRows(nnodes) = psi_.transpose()*A_.asDiagonal()*W;
        }

		MatrixXr D = V_*matrixNoCovdec_.solve(U_);

		// G = C + D
		MatrixXr G = -W.transpose()*W + D;
		Gdec_.compute(G);

	}
}

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
template<typename Derived>
MatrixXr MixedFERegressionBase<InputHandler,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::system_solve(const Eigen::MatrixBase<Derived> &b)
{

	// Resolution of the system matrixNoCov * x1 = b
	MatrixXr x1 = matrixNoCovdec_.solve(b);

	if (regressionData_.getCovariates().rows() != 0) {
		// Resolution of G * x2 = V * x1

		MatrixXr x2 = Gdec_.solve(V_*x1);

		// Resolution of the system matrixNoCov * x3 = U * x2
		x1 -= matrixNoCovdec_.solve(U_*x2);
	}

	return x1;
}


template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::setA()
{
	UInt nRegions = regressionData_.getNumberOfRegions();
	UInt m = regressionData_.isSpaceTime() ? regressionData_.getNumberofTimeObservations():1;
	A_=VectorXr::Zero(m*nRegions);
	for (int i=0; i<nRegions; i++)
	{
		for (int j=0; j<regressionData_.getIncidenceMatrix().cols(); j++)
		{
			if (regressionData_.getIncidenceMatrix()(i,j) == 1)
			{
				A_(i)+=mesh_.elementMeasure(j);
			}
		}
		for(int k=1; k<m; k++)
		{
			A_(i+k*nRegions)=A_(i);
		}
	}
}

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::getRightHandData(VectorXr& rightHandData)
{
	UInt nnodes = N_*M_;
	UInt nlocations = regressionData_.getNumberofObservations();
	rightHandData = VectorXr::Zero(nnodes);

	if (regressionData_.getCovariates().rows() == 0) //no covariate
	{
		if (regressionData_.isLocationsByNodes() && !regressionData_.isSpaceTime())
		{
				for (auto i=0; i<nlocations;++i)
				{
					auto index_i = regressionData_.getObservationsIndices()[i];
					rightHandData(index_i) = regressionData_.getObservations()[i];
				}
		}
		else if (regressionData_.isLocationsByNodes() && regressionData_.isSpaceTime() && regressionData_.getFlagParabolic())
		{
			for (auto i=0; i<regressionData_.getObservationsIndices().size();++i)
			{
				auto index_i = regressionData_.getObservationsIndices()[i];
				rightHandData(index_i) = regressionData_.getObservations()[index_i];
			}
		}
		else if (regressionData_.getNumberOfRegions() == 0) //pointwise data
		{
			rightHandData=psi_.transpose()*regressionData_.getObservations();
		}
		else //areal data
		{
			rightHandData=psi_.transpose()*A_.asDiagonal()*regressionData_.getObservations();
		}
	}
	else if (regressionData_.getNumberOfRegions() == 0) //with covariates, pointwise data
	{
		rightHandData=psi_.transpose()*LeftMultiplybyQ(regressionData_.getObservations());
	}
	else //with covariates, areal data
	{
		rightHandData=psi_.transpose()*A_.asDiagonal()*LeftMultiplybyQ(regressionData_.getObservations());
	}
}


template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::computeDegreesOfFreedom(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT)
{
	int GCVmethod = regressionData_.getGCVmethod();
	switch (GCVmethod) {
		case 1:
			computeDegreesOfFreedomExact(output_indexS, output_indexT, lambdaS, lambdaT);
			break;
		case 2:
			computeDegreesOfFreedomStochastic(output_indexS, output_indexT, lambdaS, lambdaT);
			break;
	}
}

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler, IntegratorSpace, ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::computeGeneralizedCrossValidation(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT)
{
	VectorXr dataHat;
	VectorXr z = regressionData_.getObservations();
	if(regressionData_.getCovariates().rows()==0) //Data estimated from the model
		dataHat = psi_*_solution(output_indexS,output_indexT).topRows(psi_.cols());
	else
		dataHat = z - LeftMultiplybyQ(z) + LeftMultiplybyQ(psi_*_solution(output_indexS,output_indexT).topRows(psi_.cols()));
	UInt n = dataHat.rows();
	if(regressionData_.isSpaceTime())
		{
			const std::vector<UInt>& observations_na= regressionData_.getObservationsNA();
			for(UInt id:observations_na)
			{
				dataHat[id]=0;
			}
			n-=observations_na.size();
		}

//! GCV computation
	_GCV(output_indexS,output_indexT) = (n / ((n-_dof(output_indexS,output_indexT)) * (n-_dof(output_indexS,output_indexT)))) * (z-dataHat).dot(z-dataHat);
	if (_GCV(output_indexS,output_indexT) < _bestGCV)
	{
		bestLambdaS_ = output_indexS;
		bestLambdaT_ = output_indexT;
		_bestGCV = _GCV(output_indexS,output_indexT);
	}
}


template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::computeDegreesOfFreedomExact(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT)
{

	UInt nnodes =N_*M_;
	UInt nlocations = regressionData_.getNumberofObservations();
	Real degrees=0;


	MatrixXr X1;
	if (regressionData_.getNumberOfRegions() == 0){ //pointwise data
		X1 = psi_.transpose() * LeftMultiplybyQ(psi_);
	}else{ //areal data
		X1 = psi_.transpose() * A_.asDiagonal() * LeftMultiplybyQ(psi_);
	}

	if (isRcomputed_ == false)
	{
		isRcomputed_ = true;
		//take R0 from the final matrix since it has already applied the dirichlet boundary conditions
		SpMat R0 = matrixNoCov_.bottomRightCorner(nnodes,nnodes)/lambdaS;
		R0dec_.compute(R0);
		if(!regressionData_.isSpaceTime() || !regressionData_.getFlagParabolic())
		{
			MatrixXr X2 = R0dec_.solve(R1_);
			R_ = R1_.transpose() * X2;
		}
	}

	MatrixXr P;
	MatrixXr X3=X1;

	//define the penalization matrix: note that for separable smoothin should be P=lambdaS*Psk+lambdaT*Ptk
	// but the second term has been added to X1 for dirichlet boundary conditions
	if (regressionData_.isSpaceTime() && regressionData_.getFlagParabolic())
	{
		SpMat X2 = R1_+lambdaT*LR0k_;
		P = lambdaS*X2.transpose()*R0dec_.solve(X2);
	}
	else
	{
		P = lambdaS*R_;
	}

	if(regressionData_.isSpaceTime() && !regressionData_.getFlagParabolic())
		X3 += lambdaT*Ptk_;

	//impose dirichlet boundary conditions if needed
	if(regressionData_.getDirichletIndices().size()!=0)
	{
		const std::vector<UInt>& bc_indices = regressionData_.getDirichletIndices();
		UInt nbc_indices = bc_indices.size();

		Real pen=10e20;
		for(UInt i=0; i<nbc_indices; i++)
		{
			UInt id = bc_indices[i];
			X3(id,id)=pen;
		}
	}

	X3 -= P;
	Eigen::PartialPivLU<MatrixXr> Dsolver(X3);

	auto k = regressionData_.getObservationsIndices();

	if(!regressionData_.isSpaceTime() && regressionData_.isLocationsByNodes() && regressionData_.getCovariates().rows() != 0) {
		degrees += regressionData_.getCovariates().cols();

		// Setup rhs B
		MatrixXr B;
		B = MatrixXr::Zero(nnodes,nlocations);
		// B = I(:,k) * Q
		for (auto i=0; i<nlocations;++i) {
			VectorXr ei = VectorXr::Zero(nlocations);
			ei(i) = 1;
			VectorXr Qi = LeftMultiplybyQ(ei);
			for (int j=0; j<nlocations; ++j) {
				B(k[i], j) = Qi(j);
			}
		}
		// Solve the system TX = B
		MatrixXr X;
		X = Dsolver.solve(B);
		// Compute trace(X(k,:))
		for (int i = 0; i < k.size(); ++i) {
			degrees += X(k[i], i);
		}
	}

	if (regressionData_.isSpaceTime() || !regressionData_.isLocationsByNodes())
	{
		MatrixXr X;
		X = Dsolver.solve(X1);

		if (regressionData_.getCovariates().rows() != 0) {
			degrees += regressionData_.getCovariates().cols();
		}
		for (int i = 0; i<nnodes; ++i) {
			degrees += X(i,i);
		}
	}

	_dof(output_indexS,output_indexT) = degrees;
}

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::computeDegreesOfFreedomStochastic(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT)
{

	UInt nnodes = N_*M_;
	UInt nlocations = regressionData_.getNumberofObservations();

	// std::random_device rd;
	auto seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	// Creation of the random matrix
	std::bernoulli_distribution distribution(0.5);
	UInt nrealizations = regressionData_.getNrealizations();
	MatrixXr u(nlocations, nrealizations);
	for (int j=0; j<nrealizations; ++j) {
		for (int i=0; i<nlocations; ++i) {
			if (distribution(generator)) {
				u(i,j) = 1.0;
			}
			else {
				u(i,j) = -1.0;
			}
		}
	}

	// Define the first right hand side : | I  0 |^T * psi^T * A * Q * u
	MatrixXr b = MatrixXr::Zero(2*nnodes,u.cols());
	if (regressionData_.getNumberOfRegions() == 0){
		b.topRows(nnodes) = psi_.transpose() * LeftMultiplybyQ(u);
	}else{
		b.topRows(nnodes) = psi_.transpose() * A_.asDiagonal() * LeftMultiplybyQ(u);
	}

	// Resolution of the system
	MatrixXr x = system_solve(b);

	MatrixXr uTpsi = u.transpose()*psi_;
	VectorXr edf_vect(nrealizations);
	Real q = 0;

	// Degrees of freedom = q + E[ u^T * psi * | I  0 |* x ]
	if (regressionData_.getCovariates().rows() != 0) {
		q = regressionData_.getCovariates().cols();
	}
	// For any realization we compute the degrees of freedom
	for (int i=0; i<nrealizations; ++i) {

		edf_vect(i) = uTpsi.row(i).dot(x.col(i).head(nnodes)) + q;
	}

	// Estimates: sample mean, sample variance
	Real mean = edf_vect.sum()/nrealizations;
	_dof(output_indexS,output_indexT) = mean;
}

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler, IntegratorSpace, ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::buildSpaceTimeMatrices()
{
	SpMat IM(M_,M_);
	SpMat phi;

	if(regressionData_.getFlagParabolic()) // Parabolic case
	{
		MixedFDRegression <InputHandler> FiniteDifference(mesh_time_,regressionData_);
		FiniteDifference.setDerOperator();
		SpMat L = FiniteDifference.getDerOpL(); // Matrix of finite differences
		IM.setIdentity();
		LR0k_ = kroneckerProduct(L,R0_);
		phi = IM;
		//! right hand side correction for the initial condition:
		rhs_ic_correction_ = (1/(mesh_time_[1]-mesh_time_[0]))*(R0_*regressionData_.getInitialValues());
	}
	else	// Separable case
	{
		MixedSplineRegression <InputHandler, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE> Spline(mesh_time_,regressionData_);
		SpMat IN(N_,N_);
		Spline.setPhi();
		Spline.setTimeMass();
		Spline.smoothSecondDerivative();
		if(regressionData_.getFlagMass()) // Mass penalization
		{
			IM = Spline.getTimeMass();
			IN = R0_;
		}
		else	// Identity penalization
		{
			IM.setIdentity();
			IN.setIdentity();
		}
		phi = Spline.getPhi();
		SpMat Pt = Spline.getPt();
		Ptk_ = kroneckerProduct(Pt,IN);
	}
	// Make the Kronecker product to tensorize the system
	SpMat psi_temp =  psi_;
	SpMat R1_temp = R1_;
	SpMat R0_temp = R0_;
	psi_.resize(N_*M_,N_*M_);
	psi_ = kroneckerProduct(phi,psi_temp);
	addNA();
	R1_.resize(N_*M_,N_*M_);
	R1_ = kroneckerProduct(IM,R1_temp);
	R1_.makeCompressed();
	R0_.resize(N_*M_,N_*M_);
	R0_ = kroneckerProduct(IM,R0_temp);
	R0_.makeCompressed();

	//! right hand side correction for the forcing term:

	if(this->isSpaceVarying)
	{
		VectorXr forcingTerm = rhs_ft_correction_;
		rhs_ft_correction_.resize(M_*N_);
		for(UInt i=0; i<N_; i++)
		{
			for(UInt j=0; j<M_; j++)
			{
				rhs_ft_correction_(i+j*N_) = forcingTerm(i);
			}
		}
	}
}

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
template<typename A>
void MixedFERegressionBase<InputHandler,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::apply(EOExpr<A> oper, const ForcingTerm & u)
{
	UInt nnodes = N_*M_;
	FiniteElement<IntegratorSpace, ORDER, mydim, ndim> fe;

	if(regressionData_.getNumberOfRegions()>0)
		setA();

	setPsi();

	typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
	Assembler::operKernel(oper, mesh_, fe, R1_);
	Assembler::operKernel(mass, mesh_, fe, R0_);

	if(this->isSpaceVarying)
	{
        //u=ForcingTerm(VectorXr::Zero(nnodes));
		Assembler::forcingTerm(mesh_, fe, u, rhs_ft_correction_);

	}

	if (regressionData_.isSpaceTime())
	{
		buildSpaceTimeMatrices();
	}

	VectorXr rightHandData;
	getRightHandData(rightHandData); //updated
	this->_rightHandSide = VectorXr::Zero(2*nnodes);
	this->_rightHandSide.topRows(nnodes)=rightHandData;

	this->_solution.resize(regressionData_.getLambdaS().size(),regressionData_.getLambdaT().size());
	this->_dof.resize(regressionData_.getLambdaS().size(),regressionData_.getLambdaT().size());
	this->_GCV.resize(regressionData_.getLambdaS().size(),regressionData_.getLambdaT().size());
	if(regressionData_.getCovariates().rows()!=0)
	{
		this->_beta.resize(regressionData_.getLambdaS().size(),regressionData_.getLambdaT().size());
	}

	VectorXr rhs= _rightHandSide;

	for(UInt s = 0; s<regressionData_.getLambdaS().size(); ++s)
	{
		for(UInt t = 0; t<regressionData_.getLambdaT().size(); ++t)
		{
			_rightHandSide=rhs;
			Real lambdaS = regressionData_.getLambdaS()[s];
			Real lambdaT = regressionData_.getLambdaT()[t];
			SpMat R0_lambda = (-lambdaS)*R0_; // build the SouthEast block of the matrix
			SpMat R1_lambda = (-lambdaS)*R1_;
			if(regressionData_.isSpaceTime() && regressionData_.getFlagParabolic())
				R1_lambda -= lambdaS*(lambdaT*LR0k_); // build the SouthWest block of the matrix (also the NorthEast block transposed)

			SpMat NWblock;
			if(regressionData_.getNumberOfRegions()==0) // pointwise data
				NWblock=psi_.transpose()*psi_;
			else                                        // areal data: need to add the diag(|D_1|,...,|D_N|)
			  NWblock=psi_.transpose()*A_.asDiagonal()*psi_;

			if(regressionData_.isSpaceTime() && !regressionData_.getFlagParabolic())
				NWblock+=lambdaT*Ptk_;

			this->buildMatrixNoCov(NWblock, R1_lambda, R0_lambda);

			//! right-hand side correction for space varying PDEs
			if(this->isSpaceVarying)
			{
			    _rightHandSide.bottomRows(nnodes)= (-lambdaS)*rhs_ft_correction_;
			}

			//! righ-hand side correction for initial condition in parabolic case
			if(regressionData_.isSpaceTime() && regressionData_.getFlagParabolic())
			{
				for(UInt i = 0; i<regressionData_.getInitialValues().rows(); i++)
				{
					_rightHandSide(nnodes+i) -= lambdaS*rhs_ic_correction_(i);
				}
			}
			//Applying boundary conditions if necessary
			if(regressionData_.getDirichletIndices().size() != 0)  // if areal data NO BOUNDARY CONDITIONS
				addDirichletBC();
			//factorization of the system for woodbury decomposition
			system_factorize();

			// system solution
			_solution(s,t) = this->template system_solve(this->_rightHandSide);

			if(regressionData_.computeGCV())
			{
				if (regressionData_.computeDOF())
					computeDegreesOfFreedom(s,t,lambdaS,lambdaT);
				computeGeneralizedCrossValidation(s,t,lambdaS,lambdaT);
			}
			else
			{
				_dof(s,t) = -1;
				_GCV(s,t) = -1;
			}

			// covariates computation
			if(regressionData_.getCovariates().rows()!=0)
			{
				MatrixXr W(this->regressionData_.getCovariates());
				VectorXr beta_rhs = W.transpose()*(regressionData_.getObservations() - psi_*_solution(s,t).topRows(psi_.cols()));
				_beta(s,t) = WTW_.solve(beta_rhs);
			}
		}
	}
}

template<typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
class MixedFERegression<RegressionData,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim> : public MixedFERegressionBase<RegressionData,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim>& mesh, const RegressionData& regressionData):MixedFERegressionBase<RegressionData,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>(mesh, regressionData){};
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim>& mesh, const std::vector<Real>& mesh_time, const RegressionData& regressionData):MixedFERegressionBase<RegressionData,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>(mesh, mesh_time, regressionData){};

	void apply()
	{
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	    MixedFERegressionBase<RegressionData,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::apply(stiff, ForcingTerm(std::vector<Real>(1)));
	}
};

template<typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
class MixedFERegression<RegressionDataElliptic,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim> : public MixedFERegressionBase<RegressionDataElliptic,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim>& mesh, const RegressionDataElliptic& regressionData):MixedFERegressionBase<RegressionDataElliptic,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>(mesh, regressionData){};
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim>& mesh, const std::vector<Real>& mesh_time, const RegressionDataElliptic& regressionData):MixedFERegressionBase<RegressionDataElliptic,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>(mesh, mesh_time, regressionData){};

	void apply()
	{
	if(mydim!=2 || ndim !=2){

	#ifdef R_VERSION_
		Rprintf("ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2");
	#else
		std::cout << "ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2\n";
	#endif

	}else{
		typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
		typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

	    const Real& c = this->regressionData_.getC();
	    const Eigen::Matrix<Real,2,2>& K = this->regressionData_.getK();
	    const Eigen::Matrix<Real,2,1>& b = this->regressionData_.getBeta();

	    MixedFERegressionBase<RegressionDataElliptic,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::apply(c*mass+stiff[K]+dot(b,grad), ForcingTerm(std::vector<Real>(1)));
	}
	}
};

template<typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
class MixedFERegression<RegressionDataEllipticSpaceVarying,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim> : public MixedFERegressionBase<RegressionDataEllipticSpaceVarying,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim>& mesh, const RegressionDataEllipticSpaceVarying& regressionData):MixedFERegressionBase<RegressionDataEllipticSpaceVarying,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>(mesh, regressionData){};
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim>& mesh, const std::vector<Real>& mesh_time, const RegressionDataEllipticSpaceVarying& regressionData):MixedFERegressionBase<RegressionDataEllipticSpaceVarying,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>(mesh, mesh_time, regressionData){};

	void apply()
	{
	if(mydim!=2 || ndim !=2){

	#ifdef R_VERSION_
		Rprintf("ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2");
	#else
		std::cout << "ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2\n";
	#endif

	}else{
		typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
		typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

		const Reaction& c = this->regressionData_.getC();
		const Diffusivity& K = this->regressionData_.getK();
		const Advection& b = this->regressionData_.getBeta();
		const ForcingTerm& u= this->regressionData_.getU();

		this->isSpaceVarying=TRUE;

		MixedFERegressionBase<RegressionDataEllipticSpaceVarying,IntegratorSpace,ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::apply(c*mass+stiff[K]+dot(b,grad), u);
	}
	}
};


// Temporal part
template<typename InputHandler, typename Integrator, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE>
void MixedSplineRegression<InputHandler, Integrator, SPLINE_DEGREE, ORDER_DERIVATIVE>::setPhi(){

		Spline<Integrator, SPLINE_DEGREE, ORDER_DERIVATIVE> spline(mesh_time_);
		UInt M = spline.num_knots()-SPLINE_DEGREE-1;
		UInt m = regressionData_.getNumberofTimeObservations();

		phi_.resize(m, M);
		Real value;

    for (UInt i = 0; i < m; ++i)
        for (UInt j = 0; j < M; ++j)
        {
					value = spline.BasisFunction(SPLINE_DEGREE, j, this->regressionData_.getTimeLocations()[i]);
					if (value!=0)
					{
						phi_.coeffRef(i,j) = value;
					}
				}
    phi_.makeCompressed();
}

template<typename InputHandler, typename Integrator, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE>
void MixedSplineRegression<InputHandler, Integrator, SPLINE_DEGREE, ORDER_DERIVATIVE>::setTimeMass(){

    using ETTimeMass = EOExpr<TimeMass>;

    Spline<Integrator, SPLINE_DEGREE, 0> spline(mesh_time_);

    TimeMass ETimeMass;
    ETTimeMass timeMass(ETimeMass);

    Assembler::operKernel(timeMass, spline, timeMass_);
}

template<typename InputHandler, typename Integrator, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE>
void MixedSplineRegression<InputHandler, Integrator, SPLINE_DEGREE, ORDER_DERIVATIVE>::smoothSecondDerivative(){

    using ETTimeMass = EOExpr<TimeMass>;

    Spline<Integrator, SPLINE_DEGREE, ORDER_DERIVATIVE> spline(mesh_time_);

    TimeMass ETimeMass;
    ETTimeMass timeMass(ETimeMass);

    Assembler::operKernel(timeMass, spline, Pt_);
}

// Parabolic
template<typename InputHandler>
void MixedFDRegression<InputHandler>::setDerOperator(){
	//Build the matrix of finite differences [1/dt1 0 .......... 0]
	//																			 [-1/dt1 1/dt1 0 ....0]
	//																										...
	//																			 [0 ....0 -1/dtM 1/dtM]
	UInt M = mesh_time_.size()-1;
	derOpL_.resize(M, M);

	// set the first and the last rows
	Real delta = mesh_time_[1] - mesh_time_[0];
	derOpL_.coeffRef(0,0) = 1/delta;

	delta = mesh_time_[M-1] - mesh_time_[M-2];
	derOpL_.coeffRef(M-1,M-1) = 1/delta;
	derOpL_.coeffRef(M-1,M-2) = -1/delta;

	for (UInt i = 1; i < M-1; ++i)
	{
		delta = mesh_time_[i] - mesh_time_[i-1];
		derOpL_.coeffRef(i,i-1) = -1/delta;
		derOpL_.coeffRef(i,i) 	= 1/delta;
	}

	derOpL_.makeCompressed();

}

#endif
