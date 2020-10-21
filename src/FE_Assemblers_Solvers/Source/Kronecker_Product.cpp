#include "../Include/Kronecker_Product.h"


//DENSE
/*
  SpMat kroneckerProduct(const SpMat& A, const SpMat& B)
{
	UInt Nr = A.rows();
	UInt Nc = A.cols();
	UInt Mr = B.rows();
	UInt Mc = B.cols();

	MatrixXr AB_dense(Nr*Mr, Nc*Mc);
	MatrixXr A_dense = Eigen::MatrixXd(A);

	for (UInt i = 0; i < Nr; ++i)
		for (UInt j = 0; j < Nc; ++j)
			AB_dense.block(i*Mr, j*Mc, Mr, Mc) =  A_dense.coeffRef(i,j)*B;

	SpMat AB = AB_dense.sparseView();

	AB.makeCompressed();

	return(AB);
}
*/

// SPARSE
/*
SpMat kroneckerProduct(const SpMat& A, const SpMat& B)
{
		UInt Nr = A.rows();
		UInt Nc = A.cols();
		UInt Mr = B.rows();
		UInt Mc = B.cols();

		SpMat AB(Nr*Mr, Nc*Mc);
		Real a;

		for (UInt i = 0; i < Nr; ++i) {
				for (UInt j = 0; j < Nc; ++j) {
						a = A.coeff(i,j);

						if(a != 0) {
								for (UInt k = 0; k < Mr; ++k)
								for (UInt l = 0; l < Mc; ++l) AB.insert(Mr*i+k, Mc*j+l) = a*B.coeff(k,l);
						}
				}
		}
		return(AB);
}
*/

SpMat kroneckerProduct(const SpMat& A, const SpMat& B)
{
	UInt Ar = A.rows();
	UInt Ac = A.cols();
	UInt Br = B.rows();
	UInt Bc = B.cols();

	const Real *Avalues = A.valuePtr();
	const UInt *Ainner  = A.innerIndexPtr();
	const UInt *Aouter  = A.outerIndexPtr();
	UInt Anz      = A.nonZeros();

	const Real *Bvalues = B.valuePtr();
	const UInt *Binner  = B.innerIndexPtr();
	const UInt *Bouter  = B.outerIndexPtr();
	UInt Bnz      = B.nonZeros();

	UInt ABr  = Ar * Br;
	UInt ABc  = Ac * Bc;
	UInt ABnz = Anz * Bnz;

	//Real *ABvalues = new Real[ABnz];
	//UInt *ABinner  = new UInt[ABnz];
	//UInt *ABouter  = new UInt[ABc+1];
	std::vector<Real> ABvalues(ABnz);
	std::vector<UInt> ABinner(ABnz);
	std::vector<UInt> ABouter(ABc+1);

	ABouter[0] = 0.0;
	UInt Acurrent = Aouter[0];
	UInt ij = 0;
	UInt iijj = 0;
	for (UInt i = 1; i <= Ac; i++) {
		UInt Bcurrent = Bouter[0];
		for (UInt j = 1; j <= Bc; j++) {
			ij++;
			ABouter[ij] = ABouter[ij-1] + (Aouter[i]-Aouter[i-1]) * (Bouter[j]-Bouter[j-1]);
			for (UInt ii = Acurrent; ii < Aouter[i]; ii++) {
				for (UInt jj = Bcurrent; jj < Bouter[j]; jj++) {
					ABinner[iijj]  = Ainner[ii] * Br + Binner[jj];
					ABvalues[iijj] = Avalues[ii] * Bvalues[jj];
					iijj++;
				}
			}
			Bcurrent = Bouter[j];
		}
		Acurrent = Aouter[i];
	}

	Eigen::MappedSparseMatrix<Real, Eigen::ColMajor, UInt> ABM(ABr, ABc, ABnz, ABouter.data(), ABinner.data(), ABvalues.data());
	SpMat AB(ABM);

//	delete[] ABvalues;
//	delete[] ABinner;
//	delete[] ABouter;

	return(AB);
}
