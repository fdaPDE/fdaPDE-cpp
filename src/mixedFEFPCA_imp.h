#ifndef __MIXEDFEFPCA_IMP_HPP__
#define __MIXEDFEFPCA_IMP_HPP__

#include <iostream>
#include<iterator>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <random>
#include "timing.h"
#include <fstream>

#include "R_ext/Print.h"


template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<ORDER, mydim, ndim>::computeDelta()
{
	UInt nRegions = fpcaData_.getNumberOfRegions();
	Delta_.resize(nRegions,1);
	for (int i=0; i<nRegions; i++)
	{
		Delta_(i)=0;
		for (int j=0; j<fpcaData_.getIncidenceMatrix().cols(); j++)
		{
			if (fpcaData_.getIncidenceMatrix()(i,j) == 1)
			{
				Delta_(i)+=mesh_.elementMeasure(j);
			}
		}
	}
}

template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<ORDER, mydim, ndim>::computeBasisEvaluations()
{
	//std::cout<<"Data Matrix Computation by Basis Evaluation.."<<std::endl;
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = fpcaData_.getNumberofObservations();
	Real eps = 2.2204e-016,
		 tolerance = 100 * eps;

	Psi_.resize(nlocations, nnodes);
	if (fpcaData_.isLocationsByNodes()) //pointwise data
	{
		std::vector<coeff> tripletAll;
		auto k = fpcaData_.getObservationsIndices();

		tripletAll.reserve(k.size());
		for (int i = 0; i< k.size(); ++i)
		{
			tripletAll.push_back(coeff(i,k[i],1.0));
		}
		Psi_.setFromTriplets(tripletAll.begin(),tripletAll.end());
		Psi_.makeCompressed();
	}
	else if (fpcaData_.isLocationsByBarycenter() && (fpcaData_.getNumberOfRegions()==0))
	{
		//Constexpr is used for selecting the right number of nodes to pass as a template parameter to the Element object.In case of planar domain(i.e. mydim==2), we have that the number of nodes is 3*ORDER. In case of volumetric domain (i.e. mydim==3), we have that the number of nodes is 4 nodes if ORDER==1 and 10 nodes if ORDER==2, so the expression is 6*ORDER-2. ORDER==2 if mydim==3 is not yet implemented.
		constexpr UInt Nodes = mydim==2? 3*ORDER : 6*ORDER-2;
		Element<Nodes, mydim, ndim> tri_activated;
		Real evaluator;

		for(UInt i=0; i<nlocations;i++)
		{

			tri_activated = mesh_.getElement(fpcaData_.getElementId(i));

			if(tri_activated.getId() == Identifier::NVAL)
			{
				#ifdef R_VERSION_
				Rprintf("WARNING: Observation %d is not in the domain, remove point and re-perform smoothing\n", i+1);
				#else
				std::cout << "WARNING: Observation " << i+1 <<" is not in the domain\n";
				#endif
			}
			else
			{
				for(UInt node = 0; node < Nodes ; ++node)
				{
					evaluator = fpcaData_.getBarycenter(i,node);
					Psi_.insert(i, tri_activated[node].getId()) = evaluator;
				}
			}
		} //end of for loop

		Psi_.prune(tolerance);
		Psi_.makeCompressed();
	}
	else if ((!fpcaData_.isLocationsByBarycenter()) && fpcaData_.getNumberOfRegions()==0)
	{
		//Constexpr is used for selecting the right number of nodes to pass as a template parameter to the Element object.In case of planar domain(i.e. mydim==2), we have that the number of nodes is 3*ORDER. In case of volumetric domain (i.e. mydim==3), we have that the number of nodes is 4 nodes if ORDER==1 and 10 nodes if ORDER==2, so the expression is 6*ORDER-2. ORDER==2 if mydim==3 is not yet implemented.
		constexpr UInt Nodes = mydim==2? 3*ORDER : 6*ORDER-2;
		Element<Nodes, mydim, ndim> tri_activated;
		Eigen::Matrix<Real,Nodes,1> coefficients;

		Real evaluator;
		this->barycenters_.resize(nlocations, Nodes);
		this->element_ids_.resize(nlocations);
		for(UInt i=0; i<nlocations;i++)
		{

			if (fpcaData_.getSearch() == 1) { //use Naive search
				tri_activated = mesh_.findLocationNaive(fpcaData_.template getLocations<ndim>(i));
			} else if (fpcaData_.getSearch() == 2) { //use Tree search (default)
				tri_activated = mesh_.findLocationTree(fpcaData_.template getLocations<ndim>(i));
			}

			if(tri_activated.getId() == Identifier::NVAL)
			{
				#ifdef R_VERSION_
				Rprintf("WARNING: Observation %d is not in the domain, remove point and re-perform smoothing\n", i+1);
				#else
				std::cout << "WARNING: Observation " << i+1 <<" is not in the domain\n";
				#endif
			}
			else
			{
				element_ids_(i)=tri_activated.getId();
				for(UInt node = 0; node < Nodes ; ++node)
				{
					coefficients = Eigen::Matrix<Real,Nodes,1>::Zero();
					coefficients(node) = 1; //Activates only current base
					evaluator = tri_activated.evaluate_point(fpcaData_.template getLocations<ndim>(i), coefficients);
					barycenters_(i,node)=evaluator;
					Psi_.insert(i, tri_activated[node].getId()) = evaluator;
				}
			}
		} //end of for loop

		Psi_.prune(tolerance);
		Psi_.makeCompressed();
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
				if (fpcaData_.getIncidenceMatrix()(i,j) == 1) //element j is in region i
				{
					Element<Nodes, mydim, ndim> tri = mesh_.getElement(j); //can also be a tetrahedron
					for (UInt k=0; k<Nodes; k++)
					{
						tab[tri[k].getId()] += tri.integrate(Eigen::Matrix<Real,Nodes,1>::Unit(k)); // integral over tri of psi_k
					}
				}
			}
			for (int k=0; k<nnodes; k++)
			{
				if (tab[k] != 0)
				{
					Psi_.insert(i,k) = tab[k]/Delta_(i);
				}
			}
		}
		free(tab);
		Psi_.makeCompressed();
	}
}

template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<ORDER, mydim, ndim>::buildCoeffMatrix(const SpMat& DMat,  const SpMat& AMat,  const SpMat& MMat)
{
	UInt nnodes = mesh_.num_nodes();

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
template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<ORDER, mydim, ndim>::computeDataMatrix(SpMat& DMat)
{
	UInt nnodes = mesh_.num_nodes();
	DMat.resize(nnodes,nnodes);
	if (fpcaData_.getNumberOfRegions() > 0) //areal data
		DMat = Psi_.transpose()*Delta_.asDiagonal()*Psi_;
	else //pointwise data
		DMat = Psi_.transpose()*Psi_;
}

//construct NW block of the system matrix in Ramsay when locations of observations are
//a subset of the meshe's nodes
//!! Depends on computeBasisEvaluations
template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<ORDER, mydim, ndim>::computeDataMatrixByIndices(SpMat& DMat)
{
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = fpcaData_.getNumberofObservations();

	DMat.resize(nnodes,nnodes);

	DMat.reserve(nlocations);
	for (auto i = 0; i<nlocations; ++i)
	{
		auto index = fpcaData_.getObservationsIndices()[i];
		DMat.insert(index,index) = 1;
	}
}


template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<ORDER, mydim, ndim>::computeRightHandData(VectorXr& rightHandData,FPCAObject& FPCAinput)
{
	UInt nnodes = mesh_.num_nodes();
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


template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<ORDER, mydim, ndim>::computeVarianceExplained()
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

template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<ORDER, mydim, ndim>::computeCumulativePercentageExplained()
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



template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<ORDER, mydim, ndim>::computeIterations(MatrixXr & datamatrixResiduals_, FPCAObject & FPCAinput, UInt lambda_index, UInt nnodes)
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
			#ifdef R_VERSION_
	        Rprintf("Solving system failed!!!\n");
	        #endif
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



template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<ORDER, mydim, ndim>::SetAndFixParameters()
{
	FiniteElement<ORDER, mydim, ndim> fe;

	computeDelta();
	computeBasisEvaluations(); //compute Psi
	computeDataMatrix(DMat_); //NW block

	typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
	typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	Assembler::operKernel(stiff, mesh_, fe, AMat_);
	Assembler::operKernel(mass, mesh_, fe, MMat_);


	/*const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,Eigen::DontAlignCols,", ","\n");

	std::string mass_name("Mass3D.csv");
	std::string stiff_name("Stiff3D.csv");

	std::ofstream file_mass(mass_name.c_str());
	std::ofstream file_stiff(stiff_name.c_str());

	file_mass << MatrixXr(MMat_).format(CSVFormat);
	file_stiff << MatrixXr(AMat_).format(CSVFormat);*/

	scores_mat_.resize(fpcaData_.getNPC());
	loadings_mat_.resize(fpcaData_.getNPC());
	lambda_PC_.resize(fpcaData_.getNPC());


	datamatrixResiduals_ = fpcaData_.getDatamatrix();
	solution_.resize(fpcaData_.getLambda().size());
}


///CLASS MIXEDFEFPCA
template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCA<ORDER, mydim, ndim>::apply()
{
	MixedFEFPCABase<ORDER, mydim, ndim>::SetAndFixParameters();

	for(auto np=0; np<this->fpcaData_.getNPC(); np++)
	{
		UInt i=0;
		FPCAObject FPCAinput(this->datamatrixResiduals_);
		MixedFEFPCABase<ORDER, mydim, ndim>::computeIterations(this->datamatrixResiduals_,FPCAinput,i,this->mesh_.num_nodes());

		this->scores_mat_[np]=FPCAinput.getScores();
		this->loadings_mat_[np]=FPCAinput.getLoadings();
		this->lambda_PC_[np]=this->fpcaData_.getLambda()[i];

		//Devo settare la datamatrix togliendo i risultati ottenuti
		this->datamatrixResiduals_=this->datamatrixResiduals_-this->scores_mat_[np]*this->loadings_mat_[np].transpose();

		//Change for locations
		if(!this->fpcaData_.isLocationsByNodes())
			this->loadings_mat_[np]=this->solution_[i].topRows(this->mesh_.num_nodes());

		//Normalize the loadings and unnormalize the scores
		Real load_norm=std::sqrt(this->loadings_mat_[np].transpose()*this->MMat_*this->loadings_mat_[np]);

		this->loadings_mat_[np]=this->loadings_mat_[np].transpose()/load_norm;

		this->scores_mat_[np]=this->scores_mat_[np]*load_norm;
	}
	MixedFEFPCABase<ORDER, mydim, ndim>::computeVarianceExplained();
	MixedFEFPCABase<ORDER, mydim, ndim>::computeCumulativePercentageExplained();
}


///CLASS MIXEDFEFPCAGCV

template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCAGCV<ORDER, mydim, ndim>::computeDegreesOfFreedom(UInt output_index, Real lambda)
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

template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCAGCV<ORDER, mydim, ndim>::computeDegreesOfFreedomExact(UInt output_index, Real lambda)
{
	UInt nnodes = this->mesh_.num_nodes();
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

	if (!this->fpcaData_.isLocationsByNodes()){
		MatrixXr X;
		X = Dsolver.solve(MatrixXr(X1));
		for (int i = 0; i<nnodes; ++i) {
			degrees += X(i,i);
		}
	}

	dof_[output_index] = degrees;
	this->var_[output_index] = 0;
}

template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCAGCV<ORDER, mydim, ndim>::computeDegreesOfFreedomStochastic(UInt output_index, Real lambda)
{
	UInt nnodes = this->mesh_.num_nodes();
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


template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCAGCV<ORDER, mydim, ndim>::computeDegreesOfFreedom(UInt output_index)
{
	UInt nnodes = this->mesh_.num_nodes();
	UInt nlocations = this->fpcaData_.getNumberofObservations();

	SpMat I(this->coeffmatrix_.rows(),this->coeffmatrix_.cols());
	I.setIdentity();
	SpMat coeff_inv = this->sparseSolver_.solve(I);


	Real degrees=0;

	if(this->fpcaData_.isLocationsByNodes())
	{
		VectorXr d = coeff_inv.diagonal();

		for (auto i=0; i<nlocations;++i)
		{
			auto index_i = this->fpcaData_.getObservationsIndices()[i];
			degrees+=d(index_i);
		}
	}
	else
	{
		MatrixXr An(coeff_inv.topLeftCorner(nnodes, nnodes));
		MatrixXr S = this->Psi_*An*this->Psi_.transpose();
		for (auto i=0; i<nlocations;++i)
		{
			degrees+=S(i,i);
		}
	}

	//std::cout<<"TRACE "<<degrees<<std::endl;

	dof_[output_index] = degrees;
}

template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCAGCV<ORDER, mydim, ndim>::computeGCV(FPCAObject& FPCAinput,UInt output_index)
{
	UInt s;
	VectorXr zhat;
	if(this->fpcaData_.isLocationsByNodes())
	{
		s= this->mesh_.num_nodes();
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
		#ifdef R_VERSION_
			Rprintf("WARNING: Some values of the trace of the matrix S('lambda') are inconsistent. This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'.Value of 'lambda' that produces an error is: %d \n", this->fpcaData_.getLambda()[output_index]);
		#else
			std::cout << "WARNING: Some values of the trace of the matrix S('lambda') are inconsistent. This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'.Value of 'lambda' that produces an error is:" << this->fpcaData_.getLambda()[output_index] <<"\n";
		#endif
	}
	Real stderror=norm_squared/(s-dof_[output_index]);
	GCV_[output_index]=(s/(s-dof_[output_index]))*stderror;
	//GCV_[output_index]=norm_squared/(s-(1-1/s*dof_[output_index])*(1-1/s*dof_[output_index]));
}

template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCAGCV<ORDER, mydim, ndim>::computeIterationsGCV(MatrixXr & datamatrixResiduals_, UInt nnodes, UInt np)
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
				#ifdef R_VERSION_
            	Rprintf("solving failed!!! \n");
	            #endif
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

	MixedFEFPCABase<ORDER, mydim, ndim>::computeIterations(this->datamatrixResiduals_,FPCAinput,best_GCV,this->mesh_.num_nodes());

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
	this->loadings_mat_[np]=this->solution_[best_GCV].topRows(this->mesh_.num_nodes());

	//Normalize the loadings and unnormalize the scores
	Real load_norm=std::sqrt(this->loadings_mat_[np].transpose()*this->MMat_*this->loadings_mat_[np]);

	this->loadings_mat_[np]=this->loadings_mat_[np]/load_norm;

	this->scores_mat_[np]=this->scores_mat_[np]*load_norm;
}

template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCAGCV<ORDER, mydim, ndim>::apply()
{
	MixedFEFPCABase<ORDER, mydim, ndim>::SetAndFixParameters();
	this->computeBasisEvaluations();
	dof_.resize(this->fpcaData_.getLambda().size());
	GCV_.resize(this->fpcaData_.getLambda().size());
	this->var_.resize(this->fpcaData_.getLambda().size());
	loadings_lambda_.resize(this->fpcaData_.getLambda().size());
	scores_lambda_.resize(this->fpcaData_.getLambda().size());

	for(auto np=0; np<this->fpcaData_.getNPC(); np++)
	{
		computeIterationsGCV(this->datamatrixResiduals_,this->mesh_.num_nodes(),np);
	}
	MixedFEFPCABase<ORDER, mydim, ndim>::computeVarianceExplained();
	MixedFEFPCABase<ORDER, mydim,ndim>::computeCumulativePercentageExplained();
}


///CLASS MIXEDFEFPCAKFOLD
template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCAKFold<ORDER, mydim, ndim>::computeKFolds(MatrixXr & datamatrixResiduals_, UInt lambda_index, UInt nnodes, UInt nFolds)
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
		     	#ifdef R_VERSION_
            	Rprintf("solving failed!!! \n");
	            #endif
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

template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCAKFold<ORDER, mydim, ndim>::apply()
{

	MixedFEFPCABase<ORDER, mydim, ndim>::SetAndFixParameters();
	nFolds=this->fpcaData_.getNFolds();
	KFold_.resize(this->fpcaData_.getLambda().size());

	for(auto np=0; np<this->fpcaData_.getNPC(); np++)
	{
		std::fill(KFold_.begin(),KFold_.end(),0);

		for(auto i = 0; i<this->fpcaData_.getLambda().size(); ++i)
		{

			FPCAObject FPCAinput(this->datamatrixResiduals_);

			computeKFolds(this->datamatrixResiduals_, i, this->mesh_.num_nodes(), nFolds);
		}

		UInt index_best_KF = std::distance(KFold_.begin(),std::min_element(KFold_.begin(),KFold_.end()));

		FPCAObject FPCAinput(this->datamatrixResiduals_);

		MixedFEFPCABase<ORDER, mydim, ndim>::computeIterations(this->datamatrixResiduals_,FPCAinput,index_best_KF,this->mesh_.num_nodes());

		this->scores_mat_[np]=FPCAinput.getScores();
		this->loadings_mat_[np]=FPCAinput.getLoadings();
		this->lambda_PC_[np]=this->fpcaData_.getLambda()[index_best_KF];

		//Devo settare la datamatrix togliendo i risultati ottenuti

		this->datamatrixResiduals_=this->datamatrixResiduals_-this->scores_mat_[np]*this->loadings_mat_[np].transpose();

		//Change for locations

		if(!this->fpcaData_.isLocationsByNodes())
		this->loadings_mat_[np]=this->solution_[index_best_KF].topRows(this->mesh_.num_nodes());

		//Normalize the loadings and unnormalize the scores

		Real load_norm=std::sqrt(this->loadings_mat_[np].transpose()*this->MMat_*this->loadings_mat_[np]);

		this->loadings_mat_[np]=this->loadings_mat_[np]/load_norm;

		this->scores_mat_[np]=this->scores_mat_[np]*load_norm;
	}

	MixedFEFPCABase<ORDER, mydim, ndim>::computeVarianceExplained();
	MixedFEFPCABase<ORDER, mydim,ndim>::computeCumulativePercentageExplained();

}


#endif
