#ifndef __FUNCTIONAL_PROBLEM_IMP_H__
#define __FUNCTIONAL_PROBLEM_IMP_H__

template<UInt ORDER, UInt mydim, UInt ndim>
std::pair<Real,VectorXr>
FunctionalProblem<ORDER, mydim, ndim>::computeIntegrals(const VectorXr& g) const{

	using EigenMap2WEIGHTS = Eigen::Map<const Eigen::Matrix<Real, Integrator::NNODES, 1> >;

  // Initialization
	Real int1 = 0.;
	VectorXr int2 = VectorXr::Zero(dataProblem_.getNumNodes());

	for(UInt triangle=0; triangle<dataProblem_.getNumElements(); triangle++){

		Element<EL_NNODES, mydim, ndim> tri_activated = dataProblem_.getElement(triangle);
// (1) -------------------------------------------------

    Eigen::Matrix<Real,EL_NNODES,1> sub_g;
    for (UInt i=0; i<EL_NNODES; i++){
      sub_g[i]=g[tri_activated[i].getId()];
    }
// (2) -------------------------------------------------
		Eigen::Matrix<Real,Integrator::NNODES,1> expg = (dataProblem_.getPsiQuad()*sub_g).array().exp();

    Eigen::Matrix<Real,EL_NNODES,1> sub_int2;

    int1+=expg.dot(EigenMap2WEIGHTS(&Integrator::WEIGHTS[0]))*tri_activated.getMeasure();
  	sub_int2 = dataProblem_.getPsiQuad().transpose() * expg.cwiseProduct(EigenMap2WEIGHTS(&Integrator::WEIGHTS[0]))*tri_activated.getMeasure();

  	for (UInt i=0; i<EL_NNODES; i++){
  		int2[tri_activated[i].getId()]+= sub_int2[i];
  	}
	}

	return std::pair<Real, VectorXr> (int1, int2);
}


template<UInt ORDER, UInt mydim, UInt ndim>
std::tuple<Real, VectorXr, Real, Real>
FunctionalProblem<ORDER, mydim, ndim>::computeFunctional_g(const VectorXr& g, Real lambda, const SpMat& Psi) const{

  Real int1;
  VectorXr int2;
  std::tie(int1,int2) = computeIntegrals(g);

  const UInt n = Psi.rows();
  const Real llik = -(Psi*g).sum() + n*int1;
  const Real pen = g.dot(dataProblem_.getP()*g);

	VectorXr grad1 = - VectorXr::Constant(n,1).transpose()*Psi;
	VectorXr grad2 =  n*int2;
	VectorXr grad3 = 2*g.transpose()*dataProblem_.getP();

	VectorXr grad = grad1 + grad2 + lambda*grad3;

  return std::make_tuple(llik+lambda*pen, grad, llik, pen);

}


template<UInt ORDER, UInt mydim, UInt ndim>
std::pair<Real,Real>
FunctionalProblem<ORDER, mydim, ndim>::computeLlikPen_f(const VectorXr& f) const{

  Real llik = - (dataProblem_.getGlobalPsi()*f).array().log().sum() +
                  dataProblem_.dataSize()*dataProblem_.FEintegrate(f);
  VectorXr tmp = f.array().log();
  Real pen = tmp.dot(dataProblem_.getP()*tmp);

  return std::pair<Real, Real>(llik,pen);
}

#endif
