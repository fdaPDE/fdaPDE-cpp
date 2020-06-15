#ifndef _DENSITY_INITIALIZATION_IMP_HPP_
#define _DENSITY_INITIALIZATION_IMP_HPP_

#include <unordered_set>

template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
UserInitialization<Integrator_noPoly, ORDER, mydim, ndim>::UserInitialization(const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dp):
  DensityInitialization<Integrator_noPoly, ORDER, mydim, ndim>(dp){

    initialization = dp.getFvec();

}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
const VectorXr*
UserInitialization<Integrator_noPoly, ORDER, mydim, ndim>::chooseInitialization(Real lambda) const{

  return &(this->initialization);

}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
HeatProcess<Integrator_noPoly, ORDER, mydim, ndim>::HeatProcess(const DataProblem<Integrator_noPoly, ORDER, mydim, ndim>& dp,
  const FunctionalProblem<Integrator_noPoly, ORDER, mydim, ndim>& fp):
  DensityInitialization<Integrator_noPoly, ORDER, mydim, ndim>(dp), funcProblem_(fp){

    patch_areas_= VectorXr::Zero(this->dataProblem_.getNumNodes());
    alpha_=dp.getHeatStep();
    niter_=dp.getHeatIter();
    init_proposals_.resize(niter_);
    llik_.resize(niter_);
    penTerm_.resize(niter_);

    computePatchAreas();
    computeStartingDensities();

}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
void
HeatProcess<Integrator_noPoly, ORDER, mydim, ndim>::computePatchAreas(){

  constexpr UInt Nodes = mydim==2? 3*ORDER : 6*ORDER-2;

  for(UInt t=0; t<this->dataProblem_.getNumElements(); t++){
    Element<Nodes, mydim, ndim> current_element = this->dataProblem_.getElement(t);

    // for(UInt i=0; i<Nodes; i++){
    //   if constexpr(mydim == 2) patch_areas_[current_element[i].id()] += std::abs(current_element.getArea());
    //   if constexpr(mydim == 3) patch_areas_[current_element[i].id()] += std::abs(current_element.getVolume());
    //
    // }
    for(UInt i=0; i<Nodes; i++){
      patch_areas_[current_element[i].id()] += std::abs(this->dataProblem_.getMesh().elementMeasure(t));
    }

  }
}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
VectorXr
HeatProcess<Integrator_noPoly, ORDER, mydim, ndim>::computeDensityOnlyData(){

	constexpr UInt Nodes = mydim==2? 3*ORDER : 6*ORDER-2;

	VectorXr x = VectorXr::Zero(this->dataProblem_.getNumNodes());

	for(UInt i=0; i<this->dataProblem_.getNumberofData(); i++){

    Element<Nodes, mydim, ndim> current_element;
    if(this->dataProblem_.getSearch() == 1) { //use Naive search
      current_element = this->dataProblem_.findLocationNaive(this->dataProblem_.getDatum(i));
    } else if (this->dataProblem_.getSearch() == 2) { //use Tree search (default)
      current_element = this->dataProblem_.findLocationTree(this->dataProblem_.getDatum(i));
    }

		for(UInt j=0; j<Nodes; j++){
				x[current_element[j].id()] += 1;
		}
	}

	x.array() /= patch_areas_.array();

	return (x.array() / this->dataProblem_.FEintegrate(x));

}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
void
HeatProcess<Integrator_noPoly, ORDER, mydim, ndim>::computeStartingDensities(){

	constexpr UInt Nodes = mydim==2? 3*ORDER : 6*ORDER-2;

	VectorXr x(this->dataProblem_.getNumNodes());
	x = computeDensityOnlyData();

	std::vector<std::unordered_set<UInt>> neighbours_nodes(this->dataProblem_.getNumNodes());
	// it saves in i-th position the set of the neighboor nodes id that have node i

	for(UInt t = 0; t < this->dataProblem_.getNumElements(); t++){
		Element<Nodes, mydim, ndim> current_element = this->dataProblem_.getElement(t);
		for(UInt i=0; i<Nodes;i++){
			for(UInt j=i+1; j<Nodes; j++){
					neighbours_nodes[current_element[i].id()].insert(current_element[j].id());
					neighbours_nodes[current_element[j].id()].insert(current_element[i].id());
			}
		}
	}

	for(UInt j=0; j < niter_; j++){
		VectorXr x_new(this->dataProblem_.getNumNodes());
		for(UInt k = 0; k < this->dataProblem_.getNumNodes(); k++){
			Real mean = 0.;
			for(UInt elem : neighbours_nodes[k]){
				mean += x[elem];
			}
			mean /= neighbours_nodes[k].size();

			x_new[k] = x[k] + alpha_*(mean - x[k]);
		}

		init_proposals_[j] = x_new.array() + epsilon_;  // modify initial density

    std::tie(llik_[j], penTerm_[j]) = funcProblem_.computeLlikPen_f(init_proposals_[j]);

		x.swap(x_new);
	}
}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
const VectorXr*
HeatProcess<Integrator_noPoly, ORDER, mydim, ndim>::chooseInitialization(Real lambda) const{

  VectorXr sum = llik_ + lambda*penTerm_;

  UInt index_min;
  sum.minCoeff(&index_min);

  #ifdef R_VERSION_
  Rprintf("The initialization selected for lambda %f is the number %d\n", lambda, index_min);
  #endif

  return &(init_proposals_[index_min]);
}


#endif
