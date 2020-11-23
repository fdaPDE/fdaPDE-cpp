#ifndef __MIXED_FE_REGRESSION_IMP_H__
#define __MIXED_FE_REGRESSION_IMP_H__

#include <iostream>
#include <chrono>
#include <random>
#include <fstream>

#include "R_ext/Print.h"

//----------------------------------------------------------------------------//
// Dirichlet BC

template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::addDirichletBC() //adds boundary conditions to all
{
	UInt id1,id3;

	UInt nnodes = N_*M_;

	const std::vector<UInt> * bc_indices = regressionData_.getDirichletIndices();
	const std::vector<Real> * bc_values = regressionData_.getDirichletValues();
	UInt nbc_indices = bc_indices->size();

	Real pen=10e20;

	for( auto i=0; i<nbc_indices; i++)
	 {
			id1=(*bc_indices)[i];
			id3=id1+nnodes;

			matrixNoCov_.coeffRef(id1,id1)=pen;
			matrixNoCov_.coeffRef(id3,id3)=pen;


			_rightHandSide(id1)=(*bc_values)[i]*pen;
			_rightHandSide(id3)=0;
	 }

	matrixNoCov_.makeCompressed();
}


template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::addDirichletBC_matrix() //adds boundary conditions to all
{
	UInt id1,id3;

	UInt nnodes = N_*M_;

	const std::vector<UInt> * bc_indices = regressionData_.getDirichletIndices();
	const std::vector<Real> * bc_values = regressionData_.getDirichletValues();
	UInt nbc_indices = bc_indices->size();

	Real pen=10e20;

	for( auto i=0; i<nbc_indices; i++)
	 {
			id1=(*bc_indices)[i];
			id3=id1+nnodes;

			matrixNoCov_.coeffRef(id1,id1)=pen;
			matrixNoCov_.coeffRef(id3,id3)=pen;
	 }

	matrixNoCov_.makeCompressed();
}




//Add NA
template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::addNA()
{

	const std::vector<UInt> * observations_na= regressionData_.getObservationsNA();

	for(UInt id:*observations_na)
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


//----------------------------------------------------------------------------//
// Setters

template<typename InputHandler>
template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler>::setPsi(const MeshHandler<ORDER, mydim, ndim> & mesh_)
{
	// Psi is a nlocations x nnodes  matrix, first fetch the dimensions
	UInt nnodes = mesh_.num_nodes();	// Define the number of nodes
	// Set the number of locations depending on presence or not of temporal data
	UInt nlocations = regressionData_.isSpaceTime() ? regressionData_.getNumberofSpaceObservations() : regressionData_.getNumberofObservations();
	psi_.resize(nlocations, nnodes);	// Resize the matrix

	// Optimized strategies according to the presence of locations
	if(regressionData_.isLocationsByNodes() & !regressionData_.isLocationsByBarycenter()) // .Pointwise data -- no barycenters
	{
		std::vector<coeff> tripletAll;
		if(!regressionData_.isSpaceTime()) // Just spatial case
		{
			// THEORETICAL REMARK:
			// If isLocationsByNodes is active it entails that nodes are used as locations
			// However, maybe, NOT IN ALL nodes evaluations are performed: some might be utility nodes not
			// associated with values:
			// e.g. think about a triangle with 6 nodes, probably only the main 3 might be associated
			// with a value with the others not, but this is only one possible choice:
			// The efficent fact is that anyway Phi_j(p_i) will be either 1 if i is the Lagrange
			// node of j or 0 if it is not. Note that since Phi are as many as the nodes, but only
			// some of the nodes are locations some columns of Phi = [Phi_j(p_i)]_ij might be of
			// all zeros: i.e. that Phi_j is linked to a node which is not a location since it's value is NA.
			// In such a case Phi is NOT a sqare matrix and there are more culumns than rows
			// If all nodes == locations then Phi == square Identity, but this is just a particular case
			// and, even though isLocationsByNodes is true, this might not be veriefied
			const std::vector<UInt> * k = regressionData_.getObservationsIndices();
			UInt k_size = k->size();
			tripletAll.reserve(k_size);

			for (UInt i=0; i< k_size; ++i)
			{
				// Add a value 1 for each valid index in row i
				// and column k[i] (under psi_k[i], the associated node)
				tripletAll.push_back(coeff(i,(*k)[i],1.0));
			}
		}
		else
		{
			tripletAll.reserve(nlocations);
			for (UInt i=0; i<nlocations; ++i)
			{
				tripletAll.push_back(coeff(i,i,1.0));
			}
		}
		psi_.setFromTriplets(tripletAll.begin(),tripletAll.end());
	}
	else if(regressionData_.isLocationsByBarycenter() && (regressionData_.getNumberOfRegions() == 0)) // Pointwise data -- by barycenter
	{
		// Exploit isLocationsByBarycenter simplyfication
		constexpr UInt Nodes = mydim==2 ? 3*ORDER : 6*ORDER-2;
		Element<Nodes, mydim, ndim> tri_activated;	// Dummy for element search
		Real evaluator;					// Dummy for evaluation storage

		for(UInt i=0; i<nlocations;i++)
		{ // Update Psi looping on all locations
			// We already know the id of the element containing the point
			tri_activated = mesh_.getElement(regressionData_.getElementId(i));

			if(tri_activated.getId() == Identifier::NVAL)
			{ // Invald id --> error
				Rprintf("ERROR: Point %d is not in the domain, remove point and re-perform smoothing\n", i+1);
			}
			else
			{ // tri_activated.getId() found, it's action might be felt a priori by all the psi of the element, one for each node
				for(UInt node=0; node<Nodes ; ++node)
				{// Loop on all the nodes of the found element and update the related entries of Psi
					evaluator = regressionData_.getBarycenter(i,node); // We already know the value to add
					// Insert the value in the column given by the GLOBAL indexing of the evaluated NODE
					psi_.insert(i, tri_activated[node].getId()) = evaluator;
				}
			}
		} // End of for loop
	}
	else if((!regressionData_.isLocationsByBarycenter()) && (regressionData_.getNumberOfRegions() == 0))
	{
		// THEORETICAL REMARK:
		// If isLocationsByNodes && isLocationsByBarycenters are false
		// locations are unspecified points, so we have to evaluate them directly
		constexpr UInt Nodes = mydim==2 ? 3*ORDER : 6*ORDER-2;
		Element<Nodes, mydim, ndim> tri_activated;	// Dummy for element search
		Eigen::Matrix<Real,Nodes,1> coefficients;	// Dummy for point evaluation
		Real evaluator;					// Dummy for evaluation storage

		// Resize for info storage
		this->barycenters_.resize(nlocations, Nodes);
		this->element_ids_.resize(nlocations);

		for(UInt i=0; i<nlocations;i++)
		{ // Update Psi looping on all locations
			// [[GM missing a defaulted else, raising a WARNING!]]
			if(regressionData_.getSearch() == 1)
			{ // Use Naive search
				tri_activated = mesh_.findLocationNaive(regressionData_.template getLocations<ndim>(i));
			}
			else if(regressionData_.getSearch() == 2)
			{ // Use Tree search (default)
				tri_activated = mesh_.findLocationTree(regressionData_.template getLocations<ndim>(i));
			}

			// Search the element containing the point
			if(tri_activated.getId() == Identifier::NVAL)
			{ // If not found
				Rprintf("ERROR: Point %d is not in the domain, remove point and re-perform smoothing\n", i+1);
			}
			else
			{ // tri_activated.getId() found, it's action might be felt a priori by all the psi of the element, one for each node
				element_ids_(i) = tri_activated.getId(); // Save the id of the ELEMENT containing the location

				for(UInt node=0; node<Nodes ; ++node)
				{// Loop on all the nodes of the found element and update the related entries of Psi
					// Define vector of all zeros but "node" component (necessary for function evaluate_point)
					coefficients = Eigen::Matrix<Real,Nodes,1>::Zero();
					coefficients(node) = 1; //Activates only current base-node
					// Evaluate psi in the node
					evaluator = tri_activated.evaluate_point(regressionData_.template getLocations<ndim>(i), coefficients);
					// Save barycenter information
					barycenters_(i,node)=tri_activated.getBaryCoordinates(regressionData_.template getLocations<ndim>(i))[node];
					// Insert the value in the column given by the GLOBAL indexing of the evaluated NODE
					psi_.insert(i, tri_activated[node].getId()) = evaluator;
				}
			}
		} // End of for loop
	}
	else // Areal data
	{
		constexpr UInt Nodes = mydim==2 ? 3*ORDER : 6*ORDER-2;
		Real * tab; // Psi_i temporary storage
		tab = (Real*) malloc(sizeof(Real)*nnodes);

		for(UInt i=0; i<nlocations;i++) //nlocations = number of regions
		{ // Loop to update by rows Psi matrix
			for(UInt k=0; k<nnodes; k++)
				tab[k] = 0;	//Iniialize at 0
			for(UInt j=0; j<mesh_.num_elements(); j++)
			{
				if((*(regressionData_.getIncidenceMatrix()))(i,j) == 1) // Element j is in region i
				{ // Location is related to that mesh element
					Element<Nodes, mydim, ndim> tri = mesh_.getElement(j); // Identify the element
					for(UInt k=0; k<Nodes; k++)
					{ // Add contribution of the area to right location
						tab[tri[k].getId()] += tri.integrate(Eigen::Matrix<Real,Nodes,1>::Unit(k)); // integral over tri of psi_k
					}
				}
			}
			for(UInt k=0; k<nnodes; k++)
			{ // Loop update Psi matrix i-th row
				if(tab[k] != 0) // Relevant change
				{ // Update the right position
					psi_.insert(i,k) = tab[k]/A_(i); // Divide by |D_i|
				}
			}
		}
		free(tab);	// Deallocate dynamic memory
	}
	psi_.makeCompressed();	// Compress for optimization
}

template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::setH(void)
{
	// Debug print
	// Rprintf("Computing Projection Matrix\n");

	UInt nlocations = regressionData_.getNumberofObservations();
	const MatrixXr * Wp(this->regressionData_.getCovariates());
	bool ilbn = regressionData_.isLocationsByNodes();

	if(ilbn)
	{
		const std::vector<UInt> * k = regressionData_.getObservationsIndices();

		// Some rows might be discarded [[we probably have data for every node not only the points ???]]
		UInt n_Wcols = Wp->cols();
		MatrixXr * Wp_reduced  = new MatrixXr;
		Wp_reduced->resize(regressionData_.getNumberofObservations(), n_Wcols);

		for (UInt i=0; i<nlocations; ++i)
		{
		       UInt index_i = (*k)[i];
		       for (auto j=0; j<n_Wcols; ++j)
		       {
			       (*Wp_reduced)(i,j) = (*Wp)(index_i,j);
		       }
		}
		Wp = Wp_reduced;
	}

	MatrixXr Wt(Wp->transpose());		// Compute W^t
	H_ = (*Wp)*(Wt*(*Wp)).ldlt().solve(Wt);	// using cholesky LDLT decomposition for computing hat matrix

	if(ilbn)
		delete Wp;
}

template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::setQ(void)
{
	// Debug text
 	// std::cout << "Computing Orthogonal Space Projection Matrix" << std::endl;

	// Remember Q = I-H
	Q_.resize(H_.rows(), H_.cols());	// Resizing dimensions as H
	Q_ = -H_;
	for (UInt i=0; i<H_.rows(); ++i)
	{
		Q_(i,i) += 1;			// Summing the identity by rows (or columns)
	}
}

template<typename InputHandler>
template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler>::setA(const MeshHandler<ORDER, mydim, ndim> & mesh_)
{
	UInt nRegions = regressionData_.getNumberOfRegions();	// Number of regions for areal partition
	// If the problem is temporal, m stores the number of temporal nodes, else is defaulted as 1
	UInt m = regressionData_.isSpaceTime() ? regressionData_.getNumberofTimeObservations():1;
	if(!this->regressionData_.isArealDataAvg())
	{ //areal data for FPIRLS
		A_ = VectorXr::Ones(m*nRegions);	// vector of pure ones
	}
	else
	{
		A_ = VectorXr::Zero(m*nRegions);	// neutral vector to be filled
		for(UInt i=0; i<nRegions; i++)		// fill the vector
		{
			const MatrixXi * imp = regressionData_.getIncidenceMatrix();
			for(UInt j=0; j<imp->cols(); j++)
			{
				if((*imp)(i,j) == 1)	// Valid input for region
				{
					A_(i) += mesh_.elementMeasure(j); // Add area
				}
			}
			for(UInt k=1; k<m; k++) // if m=1 we avoid the step
			{
				A_(i+k*nRegions) = A_(i); // Replicate the vector m times
			}
		}
	}
}

template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::setpsi_t_(void)
{
	//Additional storage of the transpose, which changes in case of spacetime
	psi_t_ = SpMat(psi_.transpose());
	psi_t_.makeCompressed(); // Compress sparse matrix
}

template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::setDMat(void)
{
	if(regressionData_.getWeightsMatrix()->size() == 0) // no weights
		DMat_ = psi_;
	else
		DMat_ = regressionData_.getWeightsMatrix()->asDiagonal()*psi_;


	if(regressionData_.getNumberOfRegions() == 0) // pointwise data
		DMat_ = psi_t_*DMat_;
	else                                        // areal data: need to add the diag(|D_1|,...,|D_N|)
		DMat_ = psi_t_*A_.asDiagonal()*DMat_;
}
//----------------------------------------------------------------------------//
// Utilities [[GM NOT VERY OPTMIZED, SENSE??, we have Q and P...]]

template<typename InputHandler>
MatrixXr MixedFERegressionBase<InputHandler>::LeftMultiplybyQ(const MatrixXr& u)
{
	// Weight matrix is used for GAM problems, it is also automatically added to the utility
	const VectorXr * P = this->regressionData_.getWeightsMatrix();

	if(regressionData_.getCovariates()->rows() == 0)
	{
		// Q is the projecton on Col(W) if W == 0 => Q = Identity
		if(P->size()==0)
			return u;
		else
			return P->asDiagonal()*u;
	}
	else
	{
		MatrixXr W(*(this->regressionData_.getCovariates()));
		// Check factorization, if not present factorize the matrix W^t*W
		if(isWTWfactorized_ == false)
		{
			if(P->size() == 0)
				WTW_.compute(W.transpose()*W);
			else
				WTW_.compute(W.transpose()*P->asDiagonal()*W);
			isWTWfactorized_=true; // Flag to no repeat the operation next time
		}

		MatrixXr Hu;
		// Compute H (or I-Q) the projection on Col(W) and multiply it times u
		if(P->size() == 0)
			Hu = W*WTW_.solve(W.transpose()*u);
		else
			Hu = W*WTW_.solve(W.transpose()*P->asDiagonal()*u);

		// Return the result
		if(P->size()==0)
			return u-Hu;
		else
			return P->asDiagonal()*(u - Hu);
	}

}

//----------------------------------------------------------------------------//
// Builders
template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::buildSpaceTimeMatrices()
{
	SpMat IM(M_,M_); // Matrix temporl_nodes x temporal_nodes
	SpMat phi;	 // Dummy for updte, old Psi will be overwritten by Psi_tilde

	// Distinguish between two problem classes
	if(regressionData_.getFlagParabolic())
	{ // Parabolic case
		MixedFDRegression <InputHandler> FiniteDifference(mesh_time_,regressionData_);
		FiniteDifference.setDerOperator();
		SpMat L = FiniteDifference.getDerOpL(); // Matrix of finite differences
		IM.setIdentity(); // Set as identity matrix
		LR0k_ = kroneckerProduct(L,R0_); // REMARK --> HAS TO BE ADDED TO R1 (that is R1 tilde) in the system
		phi = IM;
		// Right hand side correction for the initial condition:
		rhs_ic_correction_ = (1/(mesh_time_[1]-mesh_time_[0]))*(R0_*(*(regressionData_.getInitialValues())));
	}
	else
	{ // Separable case
		MixedSplineRegression <InputHandler> Spline(mesh_time_,regressionData_);
		SpMat IN(N_,N_);
		Spline.setPhi();
		Spline.setTimeMass();
		Spline.smoothSecondDerivative();
		if(regressionData_.getFlagMass())
		{ // Mass penalization
			IM = Spline.getTimeMass();
			IN = R0_;
		}
		else
		{ // Identity penalization
			IM.setIdentity();
			IN.setIdentity();
		}
		phi = Spline.getPhi();
		SpMat Pt = Spline.getPt();
		Ptk_ = kroneckerProduct(Pt,IN);
	}

	// Make the Kronecker product to tensorize the system, overwriting old matrices
	// Update Psi_
	SpMat psi_temp = psi_;
	psi_.resize(N_*M_,N_*M_);
	psi_ = kroneckerProduct(phi,psi_temp);
	addNA();

	// Update R1_
	SpMat R1_temp  = R1_;
	R1_.resize(N_*M_,N_*M_);
	R1_ = kroneckerProduct(IM,R1_temp);
	R1_.makeCompressed();

	// Update R0
	SpMat R0_temp  = R0_;
	R0_.resize(N_*M_,N_*M_);
	R0_ = kroneckerProduct(IM,R0_temp);
	R0_.makeCompressed();

	// right hand side correction for the forcing term:
	if(this->isSpaceVarying)
	{ // otherwie no forcing term needed
		VectorXr forcingTerm = rhs_ft_correction_; // Store old data
		rhs_ft_correction_.resize(M_*N_); // New size
		for(UInt i=0; i<N_; i++) // Update forcing term (i.e. repeat identically M_ times)
		{
			for(UInt j=0; j<M_; j++)
			{
				rhs_ft_correction_(i+j*N_) = forcingTerm(i);
			}
		}
	}
}

template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::getRightHandData(VectorXr& rightHandData)
{
	UInt nnodes = N_*M_;	// Total number of spatio-temporal nodes
	UInt nlocations = regressionData_.getNumberofObservations();	// Count number of locations
	const VectorXr * obsp = regressionData_.getObservations();	// Get the value of the observations

	rightHandData = VectorXr::Zero(nnodes);	// Initialize the rhd at 0
	if(regressionData_.getCovariates()->rows() == 0)
	{ // No covariates case
		if(regressionData_.isLocationsByNodes() && !regressionData_.isSpaceTime())
		{ // Regresionbyodes + just space --> Psi^t*z [simplified Psi] (or Psi^t*Q*z in GAM)
			// does nothing since Q==I unless in GAM [where multiplication by Q also involves P (weight matrix)]
			VectorXr tmp = LeftMultiplybyQ(*obsp);
			for(UInt i=0; i<nlocations; ++i)
			{ // Simplified multiplication by Psi^t
				auto index_i = (*(regressionData_.getObservationsIndices()))[i];
				rightHandData(index_i) = tmp(i);
			}
		}
		else if(regressionData_.isLocationsByNodes() && regressionData_.isSpaceTime() && regressionData_.getFlagParabolic())
		{ // Regressionbynodes + parabolic --> Psi^t*z [simplified Psi]
			for(UInt i=0; i<regressionData_.getObservationsIndices()->size(); ++i)
			{ // Simplified multiplication by Psi^t
				auto index_i = (*(regressionData_.getObservationsIndices()))[i];
				rightHandData(index_i) = (*obsp)[index_i];
			}
		}
		else if(regressionData_.getNumberOfRegions() == 0)
		{ // Generic pointwise pata, no optimization allowed --> Psi^t*z [or Psi^t*P*z in GAM]
			// LeftMultiplybyQ does nothing since Q==I unless in GAM [where multipliction by Q also involves P (weight matrix)]
			rightHandData = psi_.transpose()*LeftMultiplybyQ(*obsp);
		}
		else
		{ // Areal data, no optimization allowed --> Psi^t*A*z [or Psi^t*A*P*z in GAM]
			// LeftMultiplybyQ does nothing since Q==I unless in GAM [where multipliction by Q also involves P (weight matrix)]
			rightHandData = psi_.transpose()*A_.asDiagonal()*LeftMultiplybyQ(*obsp);
		}
	}
	else if(regressionData_.getNumberOfRegions() == 0)
	{ // With covariates, pointwise data, no optimization --> Psi^t*Q*z [in GAM Q=Q(P)]
		rightHandData = psi_.transpose()*LeftMultiplybyQ(*obsp);
	}
	else
	{ // With covariates, areal data, no optimization --> Psi^t*A*Q*z [in GAM Q=Q(P)]
		rightHandData = psi_.transpose()*A_.asDiagonal()*LeftMultiplybyQ(*obsp);
	}
}

template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::buildMatrixNoCov(const SpMat & NWblock, const SpMat & SWblock,  const SpMat & SEblock)
{
	UInt nnodes = N_*M_; // Note that is only space M_=1
	// Vector to be filled with the triplets used to build _coeffmatrix (reserved with the right dimension)
	std::vector<coeff> tripletAll;
	tripletAll.reserve(NWblock.nonZeros() + 2*SWblock.nonZeros() + SEblock.nonZeros());

	// Parsing all matrices, reading the values to be put inside _coeffmatrix, coordinates according to the rules
	for(UInt k=0; k<NWblock.outerSize(); ++k)
		for(SpMat::InnerIterator it(NWblock,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row(), it.col(),it.value()));
		}
	for(UInt k=0; k<SEblock.outerSize(); ++k)
		for(SpMat::InnerIterator it(SEblock,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row()+nnodes, it.col()+nnodes,it.value()));
		}
	for(UInt k=0; k<SWblock.outerSize(); ++k)
		for(SpMat::InnerIterator it(SWblock,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.col(), it.row()+nnodes,it.value()));
		}
	for (UInt k=0; k<SWblock.outerSize(); ++k)
		for (SpMat::InnerIterator it(SWblock,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row()+nnodes, it.col(), it.value()));
		}

	// Define, resize, fill and compress
	matrixNoCov_.setZero();
	matrixNoCov_.resize(2*nnodes,2*nnodes);
	matrixNoCov_.setFromTriplets(tripletAll.begin(),tripletAll.end());
	matrixNoCov_.makeCompressed();
}

//----------------------------------------------------------------------------//
// Factorizer & Solver
template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::system_factorize()
{
	UInt nnodes = N_*M_;	// Note that is only space M_=1
	const VectorXr * P = regressionData_.getWeightsMatrix(); // Matrix of weights for GAM

	// First phase: Factorization of matrixNoCov
	matrixNoCovdec_.compute(matrixNoCov_);



	if(regressionData_.getCovariates()->rows() != 0)
	{ // Needed only if there are covariates, else we can stop before
		// Second phase: factorization of matrix  G =  C + [V * matrixNoCov^-1 * U]= C + D
		// Definition of matrix U = [ psi^T * A * W | 0 ]^T and V= [ W^T*psi| 0]

		MatrixXr W(*(this->regressionData_.getCovariates()));


		U_ = MatrixXr::Zero(2*nnodes,W.cols());
		V_ = MatrixXr::Zero(W.cols(),2*nnodes);

		if(P->size()==0)
		{
			V_.leftCols(nnodes) = W.transpose()*psi_;
		}
		else
		{
			V_.leftCols(nnodes) = W.transpose()*P->asDiagonal()*psi_;
		}


		if(regressionData_.getNumberOfRegions()==0)
		{ // pointwise data
			if(P->size() == 0)
				U_.topRows(nnodes) = psi_.transpose()*W;
			else
				U_.topRows(nnodes) = psi_.transpose()*P->asDiagonal()*W;
		}
		else
		{ //areal data
			if(P->size() == 0)
				U_.topRows(nnodes) = psi_.transpose()*A_.asDiagonal()*W;
			else
				U_.topRows(nnodes) = psi_.transpose()*A_.asDiagonal()*P->asDiagonal()*W;
    	}


		MatrixXr D = V_*matrixNoCovdec_.solve(U_);

		// G = C + D
		MatrixXr G;
		if(P->size()==0)
		{
			G = -W.transpose()*W + D;
		}
		else
		{
			G = -W.transpose()*P->asDiagonal()*W + D;
		}
		Gdec_.compute(G);
	}
}

template<typename InputHandler>
template<typename Derived>
MatrixXr MixedFERegressionBase<InputHandler>::system_solve(const Eigen::MatrixBase<Derived> & b)
{
	// Resolution of the system matrixNoCov * x1 = b
	MatrixXr x1 = matrixNoCovdec_.solve(b);
	if(regressionData_.getCovariates()->rows() != 0)
	{
		// Resolution of G * x2 = V * x1
		MatrixXr x2 = Gdec_.solve(V_*x1);
		// Resolution of the system matrixNoCov * x3 = U * x2
		x1 -= matrixNoCovdec_.solve(U_*x2);
	}
	return x1;
}

//----------------------------------------------------------------------------//
// GCV

template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::computeDegreesOfFreedom(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT)
{
	std::string GCVmethod = optimizationData_.get_DOF_evaluation();
	switch (GCVmethod == "exact") {
		case 1:
			computeDegreesOfFreedomExact(output_indexS, output_indexT, lambdaS, lambdaT);
			break;
		case 0:
			computeDegreesOfFreedomStochastic(output_indexS, output_indexT, lambdaS, lambdaT);
			break;
	}
}

template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::computeGeneralizedCrossValidation(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT)
{
	VectorXr dataHat;
	const VectorXr* z = regressionData_.getObservations();
	if(regressionData_.getCovariates()->rows()==0) //Data estimated from the model
		dataHat = psi_*_solution(output_indexS,output_indexT).topRows(psi_.cols());
	else
		dataHat = *z - LeftMultiplybyQ(*z) + LeftMultiplybyQ(psi_*_solution(output_indexS,output_indexT).topRows(psi_.cols()));
	UInt n = dataHat.rows();
	if(regressionData_.isSpaceTime())
		{
			const std::vector<UInt> * observations_na= regressionData_.getObservationsNA();
			for(UInt id:*observations_na)
			{
				dataHat[id]=0;
			}
			n-=observations_na->size();
		}
//! GCV computation
	_GCV(output_indexS,output_indexT) = (n / ((n - optimizationData_.get_tuning()*((this->getDOF())(output_indexS, output_indexT))) * (n - optimizationData_.get_tuning()*((this->getDOF())(output_indexS, output_indexT))))) * (*z-dataHat).dot(*z-dataHat);
	if (_GCV(output_indexS,output_indexT) < optimizationData_.get_best_value())
	{
		optimizationData_.set_best_lambda_S(output_indexS);
		optimizationData_.set_best_lambda_T(output_indexT);
		optimizationData_.set_best_value(_GCV(output_indexS,output_indexT));
	}
}

template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::computeDegreesOfFreedomExact(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT)
{
	std::string file_name;
	UInt nnodes = N_*M_;
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
	if(regressionData_.getDirichletIndices()->size()!=0)
	{
		const std::vector<UInt> * bc_indices = regressionData_.getDirichletIndices();
		UInt nbc_indices = bc_indices->size();

		Real pen=10e20;
		for(UInt i=0; i<nbc_indices; i++)
		{
			UInt id = (*bc_indices)[i];
			X3(id,id)=pen;
		}
	}

	X3 -= P;
	Eigen::PartialPivLU<MatrixXr> Dsolver(X3);

	const auto k = regressionData_.getObservationsIndices();

	if(!regressionData_.isSpaceTime() && regressionData_.isLocationsByNodes()) {
		if(regressionData_.getCovariates()->rows() != 0)
			degrees += regressionData_.getCovariates()->cols();

		// Setup rhs B
		MatrixXr B;
		B = MatrixXr::Zero(nnodes,nlocations);
		// B = I(:,k) * Q
		for (auto i=0; i<nlocations;++i) {
			VectorXr ei = VectorXr::Zero(nlocations);
			ei(i) = 1;
			VectorXr Qi = LeftMultiplybyQ(ei);
			for (int j=0; j<nlocations; ++j) {
				B((*k)[i], j) = Qi(j);
			}
		}
		// Solve the system TX = B
		MatrixXr X;
		X = Dsolver.solve(B);
		// Compute trace(X(k,:))
		for (int i = 0; i < k->size(); ++i) {
			degrees += X((*k)[i], i);
		}
	}

	if (regressionData_.isSpaceTime() || !regressionData_.isLocationsByNodes())
	{
		MatrixXr X;
		X = Dsolver.solve(X1);

		if (regressionData_.getCovariates()->rows() != 0) {
			degrees += regressionData_.getCovariates()->cols();
		}
		for (int i = 0; i<nnodes; ++i) {
			degrees += X(i,i);
		}
	}

	_dof(output_indexS,output_indexT) = degrees;
}

template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::computeDegreesOfFreedomStochastic(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT)
{

	UInt nnodes = N_*M_;
	UInt nlocations = regressionData_.getNumberofObservations();

	// std::random_device rd;
	auto seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	// Creation of the random matrix
	std::bernoulli_distribution distribution(0.5);
	UInt nrealizations = optimizationData_.get_nrealizations();
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
	if (regressionData_.getCovariates()->rows() != 0) {
		q = regressionData_.getCovariates()->cols();
	}
	// For any realization we compute the degrees of freedom
	for (int i=0; i<nrealizations; ++i) {

		edf_vect(i) = uTpsi.row(i).dot(x.col(i).head(nnodes)) + q;
	}

	// Estimates: sample mean, sample variance
	Real mean = edf_vect.sum()/nrealizations;
	_dof(output_indexS,output_indexT) = mean;
}

template<typename InputHandler>
template<UInt ORDER, UInt mydim, UInt ndim, typename A>
void MixedFERegressionBase<InputHandler>::preapply(EOExpr<A> oper, const ForcingTerm & u, const MeshHandler<ORDER, mydim, ndim> & mesh_)
{
	const MatrixXr * Wp = regressionData_.getCovariates();

	UInt nnodes = N_*M_;	// total number of spatio-temporal nodes
	FiniteElement<ORDER, mydim, ndim> fe;

	// Set Areal data if present and no already done
	if(regressionData_.getNumberOfRegions()>0 && !isAComputed)
	{
		this->template setA<ORDER, mydim, ndim>(mesh_);
		isAComputed = true;
	}

	// Set psi matrix if not already done
	if(!isPsiComputed){
		this->template setPsi<ORDER, mydim, ndim>(mesh_);
		isPsiComputed = true;
	}

	// If there are covariates in the model set H and Q
	if(Wp->rows() != 0)
	{
		setH();
		setQ();
	}

	typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
	if(!isR1Computed)
	{
		Assembler::operKernel(oper, mesh_, fe, R1_);
		isR1Computed = true;
	}
	if(!isR0Computed)
	{
		Assembler::operKernel(mass, mesh_, fe, R0_);
		isR0Computed = true;
	}

	if(this->isSpaceVarying)
	{
		Assembler::forcingTerm(mesh_, fe, u, rhs_ft_correction_);
	}

	if(regressionData_.isSpaceTime())
	{
		this->buildSpaceTimeMatrices();
	}

	// Set final transpose of Psi matrix
	setpsi_t_();
	// Set matrix DMat for all cases
	setDMat();

	// Define right hand data [rhs]
	VectorXr rightHandData;
	getRightHandData(rightHandData); //updated
	this->_rightHandSide = VectorXr::Zero(2*nnodes);
	this->_rightHandSide.topRows(nnodes)=rightHandData;

	// Debugging purpose
	//Rprintf("Preliminary problem matrices building phase completed\n");
}

//----------------------------------------------------------------------------//
// Composed operations

// To be general, it takes in input the value of lambda
template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::buildSystemMatrix(Real lambda_S)
{
        this->R1_lambda = (-lambda_S)*(R1_);
        this->R0_lambda = (-lambda_S)*(R0_);

        this->buildMatrixNoCov(this->DMat_, this->R1_lambda, this->R0_lambda);
}


template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::buildSystemMatrix(Real lambdaS, Real lambdaT)
{
	this->R0_lambda = (-lambdaS)*R0_; // build the SouthEast block of the matrix
	this->R1_lambda = (-lambdaS)*R1_;

	// Update the SouthWest block of the matrix (also the NorthEast block transposed) if parabolic
	if(regressionData_.isSpaceTime() && regressionData_.getFlagParabolic())
		this->R1_lambda -= lambdaS*(lambdaT*LR0k_);

	// Update NorthWest block of matrix if separable problem
	if(regressionData_.isSpaceTime() && !regressionData_.getFlagParabolic())
	{
		this->buildMatrixNoCov(this->DMat_+lambdaT*Ptk_, R1_lambda, R0_lambda);
	}
	else
	{
		this->buildMatrixNoCov(this->DMat_, R1_lambda, R0_lambda);
	}

}

//----------------------------------------------------------------------------//
// Public solvers
template<typename InputHandler>
MatrixXr MixedFERegressionBase<InputHandler>::apply_to_b(const MatrixXr & b)
{
	const Real last_lambda = optimizationData_.get_last_lS_used();
	const Real lambda_ = optimizationData_.get_current_lambdaS();
        if(lambda_ != last_lambda)
        {
                this->buildSystemMatrix(lambda_);

		if(regressionData_.getDirichletIndices()->size() > 0)  // if areal data NO BOUNDARY CONDITIONS
		{
			this->addDirichletBC_matrix();
		}

                this->system_factorize();
        }

	optimizationData_.set_last_lS_used(lambda_);

        return this->template system_solve(b);
}


template<typename InputHandler>
MatrixXv  MixedFERegressionBase<InputHandler>::apply(void)
{
	UInt nnodes = N_*M_; // Define nuber of nodes
	const VectorXr * obsp = regressionData_.getObservations(); // Get observations

	UInt sizeLambdaS;
	if (!regressionData_.isSpaceTime() && !isGAMData)
	   sizeLambdaS=1;
	else
	   sizeLambdaS = optimizationData_.get_size_S();
	UInt sizeLambdaT = optimizationData_.get_size_T();

	this->_solution.resize(sizeLambdaS,sizeLambdaT);
	this->_dof.resize(sizeLambdaS,sizeLambdaT);
	this->_GCV.resize(sizeLambdaS,sizeLambdaT);
	if(regressionData_.getCovariates()->rows() != 0)
	{
		this->_beta.resize(sizeLambdaS,sizeLambdaT);
	}

	VectorXr rhs = _rightHandSide; // Save rhs for modification

	for(UInt s=0; s<sizeLambdaS; ++s)
	{
		for(UInt t=0; t<sizeLambdaT; ++t)
		{
                        Real lambdaS;
			if(!regressionData_.isSpaceTime() && !isGAMData) //at the moment only space is implemented
				{
					lambdaS = optimizationData_.get_current_lambdaS();
				}
			else
			 	lambdaS = (optimizationData_.get_lambda_S())[s];

			Real lambdaT = (optimizationData_.get_lambda_T())[t];
			_rightHandSide = rhs;

			if(isGAMData || regressionData_.isSpaceTime() || optimizationData_.get_current_lambdaS()!=optimizationData_.get_last_lS_used())
			{
				if(!regressionData_.isSpaceTime())
				{
					buildSystemMatrix(lambdaS);
				}
                        	else
				{
					buildSystemMatrix(lambdaS, lambdaT);
				}
			}

			// Right-hand side correction for space varying PDEs
			if(this->isSpaceVarying)
			{
			    _rightHandSide.bottomRows(nnodes)= (-lambdaS)*rhs_ft_correction_;
			}

			// Right-hand side correction for initial condition in parabolic case
			if(regressionData_.isSpaceTime() && regressionData_.getFlagParabolic())
			{
				for(UInt i=0; i<regressionData_.getInitialValues()->rows(); i++)
				{
					_rightHandSide(nnodes+i) -= lambdaS*rhs_ic_correction_(i);
				}
			}

			// Applying boundary conditions if necessary
			if(regressionData_.getDirichletIndices()->size() != 0)  // if areal data NO BOUNDARY CONDITIONS				
				addDirichletBC();
			

			//f Factorization of the system for woodbury decomposition
			if(isGAMData || regressionData_.isSpaceTime() || optimizationData_.get_current_lambdaS()!=optimizationData_.get_last_lS_used())
			{
				system_factorize();
			}

			// system solution
			_solution(s,t) = this->template system_solve(this->_rightHandSide);


			if(optimizationData_.get_loss_function()=="GCV" && (!isGAMData&&regressionData_.isSpaceTime()))
			{
				if (optimizationData_.get_DOF_evaluation()!="not_required")
				{
					computeDegreesOfFreedom(s,t,lambdaS,lambdaT);
				}
				computeGeneralizedCrossValidation(s,t,lambdaS,lambdaT);
			}
			else
			{
				_dof(s,t) = -1;
				_GCV(s,t) = -1;
			}

			// covariates computation
			if(regressionData_.getCovariates()->rows()!=0)
			{
				MatrixXr W(*(this->regressionData_.getCovariates()));
				VectorXr P(*(this->regressionData_.getWeightsMatrix()));
				VectorXr beta_rhs;
				if( P.size() != 0)
				{
					beta_rhs = W.transpose()*P.asDiagonal()*(*obsp - psi_*_solution(s,t).topRows(psi_.cols()));
				}
				else
				{
					beta_rhs = W.transpose()*(*obsp - psi_*_solution(s,t).topRows(psi_.cols()));
				}
				_beta(s,t) = WTW_.solve(beta_rhs);
			}
		}
	}
	if(!(isGAMData||regressionData_.isSpaceTime()) && optimizationData_.get_current_lambdaS()!=optimizationData_.get_last_lS_used())
	{
		optimizationData_.set_last_lS_used(optimizationData_.get_current_lambdaS());
	}
	_rightHandSide = rhs; // Return rhs to original status for next apply call

	//std::cout<<_solution(0,0).size()<<std::endl; per il caso GCV semplice la solution è il valore in posizione 0,0
	return this->_solution;
}

//----------------------------------------------------------------------------//

template<>
class MixedFERegression<RegressionData>: public MixedFERegressionBase<RegressionData>
{
	public:
		MixedFERegression(const RegressionData & regressionData,  OptimizationData & optimizationData, UInt nnodes_):
			MixedFERegressionBase<RegressionData>(regressionData, optimizationData, nnodes_) {};
		MixedFERegression(const std::vector<Real> & mesh_time, const RegressionData & regressionData, OptimizationData & optimizationData, UInt nnodes_, UInt spline_degree):
			MixedFERegressionBase<RegressionData>(mesh_time, regressionData, optimizationData, nnodes_, spline_degree) {};

		template<UInt ORDER, UInt mydim, UInt ndim>
		void preapply(const MeshHandler<ORDER,mydim,ndim> & mesh)
		{
			typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
		    MixedFERegressionBase<RegressionData>::preapply(stiff, ForcingTerm(), mesh);
		}
};

template<>
class MixedFERegression<RegressionDataElliptic>: public MixedFERegressionBase<RegressionDataElliptic>
{
	public:
		MixedFERegression(const RegressionDataElliptic & regressionData,  OptimizationData & optimizationData, UInt nnodes_):
			MixedFERegressionBase<RegressionDataElliptic>(regressionData, optimizationData, nnodes_) {};
		MixedFERegression(const std::vector<Real> & mesh_time, const RegressionDataElliptic & regressionData,  OptimizationData & optimizationData, UInt nnodes_, UInt spline_degree):
			MixedFERegressionBase<RegressionDataElliptic>(mesh_time, regressionData, optimizationData, nnodes_, spline_degree) {};

		template<UInt ORDER, UInt mydim, UInt ndim>
		void preapply(const MeshHandler<ORDER,mydim,ndim> & mesh)
		{
			typedef EOExpr<Mass>  ETMass;  Mass EMass;   ETMass mass(EMass);
			typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
			typedef EOExpr<Grad>  ETGrad;  Grad EGrad;   ETGrad grad(EGrad);

	  		const Real& c = this->regressionData_.getC();
	  		const Diffusion<PDEParameterOptions::Constant>& K = this->regressionData_.getK();
	  		const Advection<PDEParameterOptions::Constant>& b = this->regressionData_.getBeta();

			MixedFERegressionBase<RegressionDataElliptic>::preapply(c*mass+stiff[K]+b.dot(grad), ForcingTerm(), mesh);
			
		}
};

template<>
class MixedFERegression<RegressionDataEllipticSpaceVarying> : public MixedFERegressionBase<RegressionDataEllipticSpaceVarying>
{
	public:
		MixedFERegression(const RegressionDataEllipticSpaceVarying & regressionData,  OptimizationData & optimizationData, UInt nnodes_):
			MixedFERegressionBase<RegressionDataEllipticSpaceVarying>(regressionData, optimizationData, nnodes_) {};
		MixedFERegression(const std::vector<Real> & mesh_time, const RegressionDataEllipticSpaceVarying & regressionData,  OptimizationData & optimizationData, UInt nnodes_, UInt spline_degree):
			MixedFERegressionBase<RegressionDataEllipticSpaceVarying>(mesh_time, regressionData, optimizationData, nnodes_, spline_degree) {};


		template< UInt ORDER, UInt mydim, UInt ndim>
		void preapply(const MeshHandler<ORDER,mydim,ndim> & mesh)
		{
			typedef EOExpr<Mass>  ETMass;  Mass EMass;   ETMass mass(EMass);
			typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
			typedef EOExpr<Grad>  ETGrad;  Grad EGrad;   ETGrad grad(EGrad);

			const Reaction& c = this->regressionData_.getC();
			const Diffusion<PDEParameterOptions::SpaceVarying>& K = this->regressionData_.getK();
			const Advection<PDEParameterOptions::SpaceVarying>& b = this->regressionData_.getBeta();
			const ForcingTerm& u= this->regressionData_.getU();

			this->isSpaceVarying = TRUE;

			MixedFERegressionBase<RegressionDataEllipticSpaceVarying>::preapply(c*mass+stiff[K]+b.dot(grad), u, mesh);
			
		}
};

// -- TEMPORAL PART --
template<typename InputHandler>
void MixedSplineRegression<InputHandler>::setPhi(void)
{

	Spline<SPLINE_DEGREE, ORDER_DERIVATIVE> spline(mesh_time_);
	UInt M = spline.num_knots()-SPLINE_DEGREE-1;
	UInt m = regressionData_.getNumberofTimeObservations();

	phi_.resize(m, M);
	Real value;

	for(UInt i=0; i<m; ++i)
	{
		for(UInt j=0; j<M; ++j)
		{
			value = spline.BasisFunction(j, this->regressionData_.getTimeLocations()[i]);
			if(value!=0)
			{
				phi_.coeffRef(i,j) = value;
			}
		}
	}
	phi_.makeCompressed();
}

template<typename InputHandler>
void MixedSplineRegression<InputHandler>::setTimeMass(void)
{
    Spline<SPLINE_DEGREE, 0> spline(mesh_time_);
    Assembler::operKernel(spline, timeMass_);
}

template<typename InputHandler>
void MixedSplineRegression<InputHandler>::smoothSecondDerivative(void)
{
    Spline<SPLINE_DEGREE, ORDER_DERIVATIVE> spline(mesh_time_);
    Assembler::operKernel(spline, Pt_);
}

// Parabolic
template<typename InputHandler>
void MixedFDRegression<InputHandler>::setDerOperator(void)
{
	// Build the matrix of finite differences [1/dt1 0 .......... 0]																	 [0 ....0 -1/dtM 1/dtM]
	UInt M = mesh_time_.size()-1;
	derOpL_.resize(M, M);

	// Set the first and the last rows
	Real delta = mesh_time_[1] - mesh_time_[0];
	derOpL_.coeffRef(0,0) = 1/delta;

	delta = mesh_time_[M-1] - mesh_time_[M-2];
	derOpL_.coeffRef(M-1,M-1) = 1/delta;
	derOpL_.coeffRef(M-1,M-2) = -1/delta;

	for(UInt i=1; i<M-1; ++i)
	{
		delta = mesh_time_[i] - mesh_time_[i-1];
		derOpL_.coeffRef(i,i-1) = -1/delta;
		derOpL_.coeffRef(i,i) 	= 1/delta;
	}

	derOpL_.makeCompressed();
}

// -- GAM PART --
template<>
class MixedFERegression<GAMDataLaplace>: public MixedFERegressionBase<RegressionData>
{
	public:
		MixedFERegression(const RegressionData & regressionData, OptimizationData & optimizationData, UInt nnodes_):
			MixedFERegressionBase<RegressionData>(regressionData, optimizationData, nnodes_) {};

		template<UInt ORDER, UInt mydim, UInt ndim>
		void preapply(const MeshHandler<ORDER,mydim,ndim> & mesh)
		{
			typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
			MixedFERegressionBase<RegressionData>::preapply(stiff, ForcingTerm(), mesh);
		}
};

template<>
class MixedFERegression<GAMDataElliptic>: public MixedFERegressionBase<RegressionDataElliptic>
{
	public:
		MixedFERegression(const RegressionDataElliptic & regressionData, OptimizationData & optimizationData, UInt nnodes_):
			MixedFERegressionBase<RegressionDataElliptic>(regressionData, optimizationData, nnodes_) {};

		template< UInt ORDER, UInt mydim, UInt ndim>
		void preapply(const MeshHandler<ORDER,mydim,ndim> & mesh)
		{
			typedef EOExpr<Mass>  ETMass;  Mass EMass;   ETMass mass(EMass);
			typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
			typedef EOExpr<Grad>  ETGrad;  Grad EGrad;   ETGrad grad(EGrad);

	  		const Real& c = this->regressionData_.getC();
	  		const Diffusion<PDEParameterOptions::Constant>& K = this->regressionData_.getK();
	  		const Advection<PDEParameterOptions::Constant>& b = this->regressionData_.getBeta();

			MixedFERegressionBase<RegressionDataElliptic>::preapply(c*mass+stiff[K]+b.dot(grad), ForcingTerm(), mesh);
			
		}
};

template<>
class MixedFERegression<GAMDataEllipticSpaceVarying>: public MixedFERegressionBase<RegressionDataEllipticSpaceVarying>
{
	public:
		MixedFERegression(const RegressionDataEllipticSpaceVarying & regressionData, OptimizationData & optimizationData, UInt nnodes_):
			MixedFERegressionBase<RegressionDataEllipticSpaceVarying>( regressionData, optimizationData, nnodes_) {};

		template<UInt ORDER, UInt mydim, UInt ndim>
		void preapply(const MeshHandler<ORDER,mydim,ndim> & mesh)
		{
			typedef EOExpr<Mass>  ETMass;  Mass EMass;   ETMass mass(EMass);
			typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
			typedef EOExpr<Grad>  ETGrad;  Grad EGrad;   ETGrad grad(EGrad);

			const Reaction& c = this->regressionData_.getC();
			const Diffusion<PDEParameterOptions::SpaceVarying>& K = this->regressionData_.getK();
			const Advection<PDEParameterOptions::SpaceVarying>& b = this->regressionData_.getBeta();
			const ForcingTerm& u= this->regressionData_.getU();

			this->isSpaceVarying = TRUE;

			MixedFERegressionBase<RegressionDataEllipticSpaceVarying>::preapply(c*mass+stiff[K]+b.dot(grad), u, mesh);
			
		}
};

#endif
