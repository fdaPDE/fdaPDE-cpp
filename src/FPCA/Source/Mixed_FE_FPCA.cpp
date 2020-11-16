#include "../Include/Mixed_FE_FPCA.h"

void MixedFEFPCABase::buildCoeffMatrix(const SpMat& DMat,  const SpMat& AMat,  const SpMat& MMat)
{
	UInt nnodes = this->nnodes_;

	std::vector<coeff> tripletAll;
	tripletAll.reserve(DMat.nonZeros() + 2*AMat.nonZeros() + MMat.nonZeros());

	for (int k=0; k<DMat.outerSize(); ++k)
		for (SpMat::InnerIterator it(DMat,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row(), it.col(),it.value()));
		}
	for (int k=0; k<MMat.outerSize(); ++k)
		for (SpMat::InnerIterator it(MMat,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row()+nnodes, it.col()+nnodes,it.value()));
		}
	for (int k=0; k<AMat.outerSize(); ++k)
		for (SpMat::InnerIterator it(AMat,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.col(), it.row()+nnodes,it.value()));
		}
	for (int k=0; k<AMat.outerSize(); ++k)
		for (SpMat::InnerIterator it(AMat,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row()+nnodes, it.col(), it.value()));
		}

	coeffmatrix_.setZero();
	coeffmatrix_.resize(2*nnodes,2*nnodes);
	coeffmatrix_.setFromTriplets(tripletAll.begin(),tripletAll.end());
	coeffmatrix_.makeCompressed();
}

//construct NW block of the system matrix when basis evaluation is necessary
//!! Depends on computeBasisEvaluations
void MixedFEFPCABase::computeDataMatrix(SpMat& DMat)
{
	UInt nnodes = this->nnodes_;
	DMat.resize(nnodes,nnodes);
	if (fpcaData_.getNumberOfRegions() > 0) //areal data
		DMat = Psi_.transpose()*Delta_.asDiagonal()*Psi_;
	else //pointwise data
		DMat = Psi_.transpose()*Psi_;
}

//construct NW block of the system matrix in Ramsay when locations of observations are
//a subset of the meshe's nodes
//!! Depends on computeBasisEvaluations
void MixedFEFPCABase::computeDataMatrixByIndices(SpMat& DMat)
{
	UInt nnodes = this->nnodes_;
	UInt nlocations = fpcaData_.getNumberofObservations();

	DMat.resize(nnodes,nnodes);

	DMat.reserve(nlocations);
	for (auto i = 0; i<nlocations; ++i)
	{
		auto index = fpcaData_.getObservationsIndices()[i];
		DMat.insert(index,index) = 1;
	}
}

void MixedFEFPCABase::computeRightHandData(VectorXr& rightHandData,FPCAObject& FPCAinput)
{
	UInt nnodes = this->nnodes_;
	UInt nlocations = fpcaData_.getNumberofObservations();
	rightHandData = VectorXr::Zero(nnodes);

	if (fpcaData_.isLocationsByNodes()) //pointwise data
	{

		for (auto i=0; i<nlocations;++i)
		{
			auto index_i = fpcaData_.getObservationsIndices()[i];
			rightHandData(index_i) = FPCAinput.getObservationData()[i];
		}
	}
	else if (fpcaData_.getNumberOfRegions()==0)
	{
		rightHandData=Psi_.transpose()*FPCAinput.getObservationData();
	}
	else //areal data
	{
		rightHandData=Psi_.transpose()*Delta_.asDiagonal()*FPCAinput.getObservationData();
	}
}

void MixedFEFPCABase::computeVarianceExplained()
{

	MatrixXr U_not_normalized(scores_mat_[0].size(),scores_mat_.size());
	for(UInt i=0;i<scores_mat_.size();i++)
		U_not_normalized.col(i)=scores_mat_[i];
	Eigen::HouseholderQR<MatrixXr> qr(U_not_normalized);
	MatrixXr R=qr.matrixQR().triangularView<Eigen::Upper>();
	variance_explained_.resize(fpcaData_.getNPC());
	for(UInt i=0;i<variance_explained_.size();i++)
	variance_explained_[i]=(R.diagonal()*R.diagonal().transpose()).diagonal()[i]/scores_mat_[0].size();
}

void MixedFEFPCABase::computeCumulativePercentageExplained()
{
	Eigen::BDCSVD<MatrixXr> svd(fpcaData_.getDatamatrix(),Eigen::ComputeThinU|Eigen::ComputeThinV);
	MatrixXr U_ALL(fpcaData_.getDatamatrix().rows(),fpcaData_.getDatamatrix().rows());
	for(UInt i=0;i<svd.singularValues().rows();i++)
		U_ALL.col(i)=svd.matrixU().col(i)*svd.singularValues().diagonal()[i]*std::sqrt((svd.matrixV().col(i).transpose()*MMat_*svd.matrixV().col(i)>0.0)?svd.matrixV().col(i).transpose()*MMat_*svd.matrixV().col(i):2.2204e-016*100);

	Real TotVar=(U_ALL.transpose()*U_ALL).trace()/fpcaData_.getDatamatrix().rows();

	cumsum_percentage_.resize(fpcaData_.getNPC());

	std::partial_sum(variance_explained_.begin(),variance_explained_.end(), cumsum_percentage_.begin());
	std::for_each(cumsum_percentage_.begin(), cumsum_percentage_.end(), [&TotVar](Real& i){i=i/TotVar;});
}

void MixedFEFPCABase::computeIterations(MatrixXr & datamatrixResiduals_, FPCAObject & FPCAinput, UInt lambda_index, UInt nnodes)
{

	Real lambda = fpcaData_.getLambda()[lambda_index];
	SpMat AMat_lambda = (-lambda)*AMat_;
	SpMat MMat_lambda = (-lambda)*MMat_;

	buildCoeffMatrix(DMat_, AMat_lambda, MMat_lambda);

	sparseSolver_.analyzePattern(coeffmatrix_);
	sparseSolver_.factorize(coeffmatrix_);
	solution_[lambda_index].resize(coeffmatrix_.rows());

	UInt niter=30;

	for(auto j=0;j<niter;j++)
	{
		FPCAinput.setObservationData(datamatrixResiduals_);
		VectorXr rightHandData;
		computeRightHandData(rightHandData,FPCAinput);
		b_ = VectorXr::Zero(2*nnodes);
		b_.topRows(nnodes)=rightHandData;

		solution_[lambda_index]=sparseSolver_.solve(b_);
		if(sparseSolver_.info()!=Eigen::Success)
		{
	        Rprintf("Solving system failed!!!\n");
		}

		if(fpcaData_.isLocationsByNodes())
			FPCAinput.setLoadings(nnodes, solution_[lambda_index], fpcaData_.getObservationsIndices());
		else
			FPCAinput.setLoadingsPsi(nnodes, solution_[lambda_index], Psi_);

		FPCAinput.setScores(datamatrixResiduals_);
	}

	if(fpcaData_.isLocationsByNodes())
	{
		UInt nlocations=nnodes;
		FPCAinput.finalizeLoadings(fpcaData_.getObservationsIndices(),nlocations);
	}

}

///CLASS MIXEDFEFPCAGCV
void MixedFEFPCAGCV::computeDegreesOfFreedom(UInt output_index, Real lambda)
{
	int GCVmethod = this->fpcaData_.getGCVmethod();
	switch (GCVmethod) {
		case 1:
			computeDegreesOfFreedomExact(output_index, lambda);
			break;
		case 2:
			computeDegreesOfFreedomStochastic(output_index, lambda);
			break;
	}
}

void MixedFEFPCAGCV::computeDegreesOfFreedomExact(UInt output_index, Real lambda)
{
	UInt nnodes = this->nnodes_;
	UInt nlocations = this->fpcaData_.getNumberofObservations();
	Real degrees=0;

	MatrixXr X1 = this->Psi_.transpose() * this->Psi_;

	if (this->isRcomputed_ == false ){
		this->isRcomputed_ = true;
		Sparse_LU solver;
		solver.compute(this->MMat_);
		auto X2 = solver.solve(this->AMat_);
		this->R_ = this->AMat_.transpose() * X2;
	}

	MatrixXr X3 = X1 + lambda * this->R_;
	Eigen::LDLT<MatrixXr> Dsolver(X3);

	auto k = this->fpcaData_.getObservationsIndices();

	if (this->fpcaData_.isLocationsByNodes()){
		MatrixXr X;		
		MatrixXr B;
		B = MatrixXr::Zero(nnodes,nlocations);
		for (auto i=0; i<nlocations; ++i)
			B.row(k[i]) = Eigen::VectorXd::Unit(nlocations,i);
		X = Dsolver.solve(B);
		for (auto i = 0; i < k.size(); ++i) {
			degrees += X(k[i], i);
		}

	} else {
		MatrixXr X;
		X = Dsolver.solve(MatrixXr(X1));
		for (int i = 0; i<nnodes; ++i) {
			degrees += X(i,i);
		}
	}

	dof_[output_index] = degrees;
	this->var_[output_index] = 0;
}

void MixedFEFPCAGCV::computeDegreesOfFreedomStochastic(UInt output_index, Real lambda)
{
	UInt nnodes = this->nnodes_;
	UInt nlocations = this->fpcaData_.getNumberofObservations();

	std::default_random_engine generator;

	// Creation of the random matrix
	std::bernoulli_distribution distribution(0.5);
	UInt nrealizations = this->fpcaData_.getNrealizations();
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

	// Define the first right hand side : | I  0 |^T * psi^T * Q * u
	MatrixXr b = MatrixXr::Zero(2*nnodes,u.cols());
	b.topRows(nnodes) = this->Psi_.transpose()* u;
	// Resolution of the system
	//MatrixXr x = system_solve(b);
	Sparse_LU solver;
	solver.compute(this->coeffmatrix_);
	auto x = solver.solve(b);
	MatrixXr uTpsi = u.transpose()*this->Psi_;
	VectorXr edf_vect(nrealizations);
	Real q = 0;
	Real var = 0;
	// For any realization we calculate the degrees of freedom
	for (int i=0; i<nrealizations; ++i)
	{
		edf_vect(i) = uTpsi.row(i).dot(x.col(i).head(nnodes)) + q;
		var += edf_vect(i)*edf_vect(i);
	}
	// Estimates: sample mean, sample variance
	Real mean = edf_vect.sum()/nrealizations;
	dof_[output_index] = mean;
	var /= nrealizations;
	var -= mean*mean;
	this->var_[output_index]=var;
}


void MixedFEFPCAGCV::computeGCV(FPCAObject& FPCAinput,UInt output_index)
{
	UInt s;
	VectorXr zhat;
	if(this->fpcaData_.isLocationsByNodes())
	{
		s= this->nnodes_;
		zhat=VectorXr::Zero(s);
		for(auto i=0;i<this->fpcaData_.getObservationsIndices().size();i++)
			zhat(this->fpcaData_.getObservationsIndices()[i])=FPCAinput.getObservationData()[i];
	} else {
		s= this->fpcaData_.getNumberofObservations();
		zhat=FPCAinput.getObservationData();
	}
	Real norm_squared=(zhat-FPCAinput.getLoadings()).transpose()*(zhat-FPCAinput.getLoadings());
	
	if(s-dof_[output_index]<0)
	{
			Rprintf("WARNING: Some values of the trace of the matrix S('lambda') are inconsistent. This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'.Value of 'lambda' that produces an error is: %d \n", this->fpcaData_.getLambda()[output_index]);
	}
	Real stderror=norm_squared/(s-dof_[output_index]);
	GCV_[output_index]=(s/(s-dof_[output_index]))*stderror;
	//GCV_[output_index]=norm_squared/(s-(1-1/s*dof_[output_index])*(1-1/s*dof_[output_index]));
}

void MixedFEFPCAGCV::computeIterationsGCV(MatrixXr & datamatrixResiduals_, UInt nnodes, UInt np)
{
	UInt niter=20;
//	FPCAObject FPCAinput(this->datamatrixResiduals_);
	UInt best_GCV=0;
	std::vector<SpMat> AMat_lambda_vec;
	std::vector<SpMat> MMat_lambda_vec;

	AMat_lambda_vec.resize(this->fpcaData_.getLambda().size());
	MMat_lambda_vec.resize(this->fpcaData_.getLambda().size());
	for (auto i=0; i<this->fpcaData_.getLambda().size(); ++i)
	{

		FPCAObject FPCAinput(this->datamatrixResiduals_);
		Real lambda = this->fpcaData_.getLambda()[i];
		AMat_lambda_vec[i] = (-lambda)*this->AMat_;
		MMat_lambda_vec[i] = (-lambda)*this->MMat_;
		for (auto j=0; j<niter; ++j)
		{
			FPCAinput.setObservationData(datamatrixResiduals_);
			this->buildCoeffMatrix(this->DMat_, AMat_lambda_vec[i], MMat_lambda_vec[i]);
			if (j==0){
				this->sparseSolver_.analyzePattern(this->coeffmatrix_);
			}

			this->sparseSolver_.factorize(this->coeffmatrix_);
			this->solution_[i].resize(this->coeffmatrix_.rows());
			VectorXr rightHandData;
			this->computeRightHandData(rightHandData,FPCAinput);
			this->b_ = VectorXr::Zero(2*nnodes);
			this->b_.topRows(nnodes) = rightHandData;

			this->solution_[i]=this->sparseSolver_.solve(this->b_);
			if (this->sparseSolver_.info()!=Eigen::Success)
			{
            	Rprintf("solving failed!!! \n");
			}
			if (this->fpcaData_.isLocationsByNodes())
				FPCAinput.setLoadings(nnodes, this->solution_[i], this->fpcaData_.getObservationsIndices());
			else
				FPCAinput.setLoadingsPsi(nnodes, this->solution_[i],this->Psi_);

			FPCAinput.setScores(this->datamatrixResiduals_);
			loadings_lambda_[i]=FPCAinput.getLoadings();
			scores_lambda_[i]=FPCAinput.getScores();
		}
		if (np==0){
			computeDegreesOfFreedom(i,lambda);
		}
		computeGCV(FPCAinput,i);
	}
	best_GCV = std::distance(GCV_.begin(),std::min_element(GCV_.begin(),GCV_.end()));

	// aggiungiamo pezzo per calcolare la soluzione una volta scelto il miglior gcv

	FPCAObject FPCAinput(this->datamatrixResiduals_);

	MixedFEFPCABase::computeIterations(this->datamatrixResiduals_,FPCAinput,best_GCV,this->nnodes_);

//	if(this->fpcaData_.isLocationsByNodes())  � gi� dentro compute_iterations
//	{
//		UInt nlocations=nnodes;
//		FPCAinput.finalizeLoadings(this->fpcaData_.getObservationsIndices(),nlocations);
//	}

	this->scores_mat_[np]=FPCAinput.getScores();
	this->loadings_mat_[np]=FPCAinput.getLoadings();
	this->lambda_PC_[np]=this->fpcaData_.getLambda()[best_GCV];

	//Devo settare la datamatrix togliendo i risultati ottenuti
	this->datamatrixResiduals_=this->datamatrixResiduals_-this->scores_mat_[np]*this->loadings_mat_[np].transpose();

	//Change for locations
	if(!this->fpcaData_.isLocationsByNodes())
	this->loadings_mat_[np]=this->solution_[best_GCV].topRows(this->nnodes_);

	//Normalize the loadings and unnormalize the scores
	Real load_norm=std::sqrt(this->loadings_mat_[np].transpose()*this->MMat_*this->loadings_mat_[np]);

	this->loadings_mat_[np]=this->loadings_mat_[np]/load_norm;

	this->scores_mat_[np]=this->scores_mat_[np]*load_norm;
}

///CLASS MIXEDFEFPCAKFOLD
void MixedFEFPCAKFold::computeKFolds(MatrixXr & datamatrixResiduals_, UInt lambda_index, UInt nnodes, UInt nFolds)
{

	Real lambda = this->fpcaData_.getLambda()[lambda_index];
	SpMat AMat_lambda = (-lambda)*this->AMat_;
	SpMat MMat_lambda = (-lambda)*this->MMat_;

	this->buildCoeffMatrix(this->DMat_, AMat_lambda, MMat_lambda);
	this->sparseSolver_.analyzePattern(this->coeffmatrix_);
	this->sparseSolver_.factorize(this->coeffmatrix_);
	this->solution_[lambda_index].resize(this->coeffmatrix_.rows());

	UInt niter=20;
	std::vector<UInt> indices_valid;

	for(auto k=0; k<nFolds; k++)
	{
		UInt length_chunk = floor(static_cast<double>(datamatrixResiduals_.rows()/nFolds));
		indices_valid.resize(datamatrixResiduals_.rows());

		std::iota(indices_valid.begin(),indices_valid.begin()+length_chunk,k*length_chunk);
		if (k==0)
			std::iota(indices_valid.begin()+length_chunk,indices_valid.end(),(k+1)*length_chunk);
		else
		{
			std::iota(indices_valid.begin()+length_chunk,indices_valid.begin()+(k+1)*length_chunk,0);
			std::iota(indices_valid.begin()+(k+1)*length_chunk,indices_valid.end(),(k+1)*length_chunk);
		}

		VectorXi indices_v=Eigen::Map<VectorXi,Eigen::Unaligned> (indices_valid.data(),indices_valid.size());

		Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> perm(indices_v);

    	        MatrixXr X_clean_train=(perm*datamatrixResiduals_).bottomRows(datamatrixResiduals_.rows()-length_chunk);
		MatrixXr X_valid=(perm*datamatrixResiduals_).topRows(length_chunk);

		VectorXr X_clean_train_mean=X_clean_train.colwise().mean();
		VectorXr X_valid_mean=X_valid.colwise().mean();

		VectorXr ones=VectorXr::Constant(X_clean_train.rows(),1,1);
		X_clean_train=X_clean_train-ones*X_clean_train_mean.transpose();
		X_valid=X_valid-ones*X_valid_mean.transpose();

		FPCAObject FPCAinputKF(X_clean_train);

		for(auto j=0; j<niter; j++)
		{

			FPCAinputKF.setObservationData(X_clean_train);

			VectorXr rightHandData;
			this->computeRightHandData(rightHandData,FPCAinputKF);
			this->b_ = VectorXr::Zero(2*nnodes);
			this->b_.topRows(nnodes)=rightHandData;

			this->solution_[lambda_index]=this->sparseSolver_.solve(this->b_);
			if(this->sparseSolver_.info()!=Eigen::Success)
			{
            	Rprintf("solving failed!!! \n");
			}

			if(this->fpcaData_.isLocationsByNodes())
				FPCAinputKF.setLoadings(nnodes, this->solution_[lambda_index],this->fpcaData_.getObservationsIndices());
			else
				FPCAinputKF.setLoadingsPsi(nnodes, this->solution_[lambda_index],this->Psi_);

			FPCAinputKF.setScores(X_clean_train);
		}
		if(this->fpcaData_.isLocationsByNodes())
		{
			UInt nlocations=nnodes;
			FPCAinputKF.finalizeLoadings(this->fpcaData_.getObservationsIndices(),nlocations);
		}

		Real U_hat_const=FPCAinputKF.getLoadings().squaredNorm() + lambda* (this->solution_[lambda_index].bottomRows(nnodes)).transpose()*this->MMat_*this->solution_[lambda_index].bottomRows(nnodes);
		VectorXr U_hat_valid=(X_valid*FPCAinputKF.getLoadings())/U_hat_const;

 		Real diffCV=(X_valid-U_hat_valid*FPCAinputKF.getLoadings().transpose()).squaredNorm()/(X_valid.rows()*X_valid.cols());

		Real sumCV=(X_valid+U_hat_valid*FPCAinputKF.getLoadings().transpose()).squaredNorm()/(X_valid.rows()*X_valid.cols());

		KFold_[lambda_index]+=std::min(diffCV,sumCV);
	}

}

///CLASS MIXEDFEFPCA
void MixedFEFPCA::apply()
{
	for(auto np=0; np<this->fpcaData_.getNPC(); np++)
	{
		UInt i=0;
		FPCAObject FPCAinput(this->datamatrixResiduals_);
		MixedFEFPCABase::computeIterations(this->datamatrixResiduals_,FPCAinput,i,this->nnodes_);

		this->scores_mat_[np]=FPCAinput.getScores();
		this->loadings_mat_[np]=FPCAinput.getLoadings();
		this->lambda_PC_[np]=this->fpcaData_.getLambda()[i];

		//Devo settare la datamatrix togliendo i risultati ottenuti
		this->datamatrixResiduals_=this->datamatrixResiduals_-this->scores_mat_[np]*this->loadings_mat_[np].transpose();

		//Change for locations
		if(!this->fpcaData_.isLocationsByNodes())
			this->loadings_mat_[np]=this->solution_[i].topRows(this->nnodes_);

		//Normalize the loadings and unnormalize the scores
		Real load_norm=std::sqrt(this->loadings_mat_[np].transpose()*this->MMat_*this->loadings_mat_[np]);

		this->loadings_mat_[np]=this->loadings_mat_[np].transpose()/load_norm;

		this->scores_mat_[np]=this->scores_mat_[np]*load_norm;
	}
	MixedFEFPCABase::computeVarianceExplained();
	MixedFEFPCABase::computeCumulativePercentageExplained();
}


///CLASS MIXEDFEFPCAGCV
void MixedFEFPCAGCV::apply()
{
	dof_.resize(this->fpcaData_.getLambda().size());
	GCV_.resize(this->fpcaData_.getLambda().size());
	this->var_.resize(this->fpcaData_.getLambda().size());
	loadings_lambda_.resize(this->fpcaData_.getLambda().size());
	scores_lambda_.resize(this->fpcaData_.getLambda().size());

	for(auto np=0; np<this->fpcaData_.getNPC(); np++)
	{
		computeIterationsGCV(this->datamatrixResiduals_,this->nnodes_,np);
	}
	MixedFEFPCABase::computeVarianceExplained();
	MixedFEFPCABase::computeCumulativePercentageExplained();

	

}


///CLASS MIXEDFEFPCAKFOLD
void MixedFEFPCAKFold::apply()
{
	nFolds=this->fpcaData_.getNFolds();
	KFold_.resize(this->fpcaData_.getLambda().size());

	for(auto np=0; np<this->fpcaData_.getNPC(); np++)
	{
		std::fill(KFold_.begin(),KFold_.end(),0);

		for(auto i = 0; i<this->fpcaData_.getLambda().size(); ++i)
		{

			FPCAObject FPCAinput(this->datamatrixResiduals_);

			computeKFolds(this->datamatrixResiduals_, i, this->nnodes_, nFolds);
		}

		UInt index_best_KF = std::distance(KFold_.begin(),std::min_element(KFold_.begin(),KFold_.end()));

		FPCAObject FPCAinput(this->datamatrixResiduals_);

		MixedFEFPCABase::computeIterations(this->datamatrixResiduals_,FPCAinput,index_best_KF,this->nnodes_);

		this->scores_mat_[np]=FPCAinput.getScores();
		this->loadings_mat_[np]=FPCAinput.getLoadings();
		this->lambda_PC_[np]=this->fpcaData_.getLambda()[index_best_KF];

		//Devo settare la datamatrix togliendo i risultati ottenuti

		this->datamatrixResiduals_=this->datamatrixResiduals_-this->scores_mat_[np]*this->loadings_mat_[np].transpose();

		//Change for locations

		if(!this->fpcaData_.isLocationsByNodes())
		this->loadings_mat_[np]=this->solution_[index_best_KF].topRows(this->nnodes_);

		//Normalize the loadings and unnormalize the scores

		Real load_norm=std::sqrt(this->loadings_mat_[np].transpose()*this->MMat_*this->loadings_mat_[np]);

		this->loadings_mat_[np]=this->loadings_mat_[np]/load_norm;

		this->scores_mat_[np]=this->scores_mat_[np]*load_norm;
	}

	MixedFEFPCABase::computeVarianceExplained();
	MixedFEFPCABase::computeCumulativePercentageExplained();

}
