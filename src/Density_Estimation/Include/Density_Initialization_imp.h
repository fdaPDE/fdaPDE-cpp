#ifndef __DENSITY_INITIALIZATION_IMP_H__
#define __DENSITY_INITIALIZATION_IMP_H__

#include <unordered_set>

template<UInt ORDER, UInt mydim, UInt ndim>
UserInitialization<ORDER, mydim, ndim>::UserInitialization(const DataProblem<ORDER, mydim, ndim>& dp):
  DensityInitialization<ORDER, mydim, ndim>(dp){

    initialization = dp.getFvec();

}


template<UInt ORDER, UInt mydim, UInt ndim>
const VectorXr*
UserInitialization<ORDER, mydim, ndim>::chooseInitialization(Real lambda) const{

  return &(this->initialization);

}


template<UInt ORDER, UInt mydim, UInt ndim>
HeatProcess<ORDER, mydim, ndim>::HeatProcess(const DataProblem<ORDER, mydim, ndim>& dp,
  const FunctionalProblem<ORDER, mydim, ndim>& fp):
  DensityInitialization<ORDER, mydim, ndim>(dp), funcProblem_(fp){

    patch_areas_= VectorXr::Zero(this->dataProblem_.getNumNodes());
    alpha_=dp.getHeatStep();
    niter_=dp.getHeatIter();
    init_proposals_.resize(niter_);
    llik_.resize(niter_);
    penTerm_.resize(niter_);

    data_index_.resize(this->dataProblem_.dataSize());
    std::iota(data_index_.begin(),data_index_.end(),0);

    computePatchAreas();
    computeStartingDensities();

}


template<UInt ORDER, UInt mydim, UInt ndim>
void
HeatProcess<ORDER, mydim, ndim>::computePatchAreas(){

  for(UInt t=0; t<this->dataProblem_.getNumElements(); ++t){
    Element<EL_NNODES, mydim, ndim> current_element = this->dataProblem_.getElement(t);
    for(const auto& node : current_element)
      patch_areas_[node.id()] += current_element.getMeasure();
  }
}


template<UInt ORDER, UInt mydim, UInt ndim>
VectorXr
HeatProcess<ORDER, mydim, ndim>::computeDensityOnlyData(){

	VectorXr x = VectorXr::Zero(this->dataProblem_.getNumNodes());

  for(UInt i : data_index_){
    Element<EL_NNODES, mydim, ndim> current_element = this->dataProblem_.findLocation(this->dataProblem_.data(i));
    for(const auto& node : current_element)
				x[node.id()] += 1;
	}

	x.array() /= patch_areas_.array();

	return (x.array() / this->dataProblem_.FEintegrate(x));

}


template<UInt ORDER, UInt mydim, UInt ndim>
void
HeatProcess<ORDER, mydim, ndim>::computeStartingDensities(){

	VectorXr x = computeDensityOnlyData();

	std::vector<std::unordered_set<UInt>> neighbours_nodes(this->dataProblem_.getNumNodes());
	// it saves in i-th position the set of the neighboor nodes id that have node i

	for(UInt t = 0; t < this->dataProblem_.getNumElements(); ++t){
		Element<EL_NNODES, mydim, ndim> current_element = this->dataProblem_.getElement(t);
		for(UInt i=0; i<EL_NNODES;i++){
			for(UInt j=i+1; j<EL_NNODES; j++){
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


template<UInt ORDER, UInt mydim, UInt ndim>
const VectorXr*
HeatProcess<ORDER, mydim, ndim>::chooseInitialization(Real lambda) const{

  VectorXr sum = llik_ + lambda*penTerm_;

  UInt index_min;
  sum.minCoeff(&index_min);

  Rprintf("The initialization selected for lambda %f is the number %d\n", lambda, index_min);

  return &(init_proposals_[index_min]);
}

template<UInt ORDER, UInt mydim, UInt ndim>
Heat_CV<ORDER, mydim, ndim>::Heat_CV(const DataProblem<ORDER, mydim, ndim>& dp,
  const FunctionalProblem<ORDER, mydim, ndim>& fp, UInt K):
  HeatProcess<ORDER, mydim, ndim>(dp, fp), error_(dp), nFolds_(K){

    cv_errors_.resize(this->niter_, 0);

    K_folds_.resize(dp.dataSize());

    perform_init_cv();

}

template<UInt ORDER, UInt mydim, UInt ndim>
void
Heat_CV<ORDER, mydim, ndim>::perform_init_cv(){

    UInt N = this->dataProblem_.dataSize();
    UInt K = nFolds_;

    for (UInt i = 0; i< N; i++){
      UInt length = ((i % K) <= (N % K))? (i % K)*(N/K +1) : (N % K) + (i % K)*(N/K);
      K_folds_[length + i/K] = i;
    }

    // cycle on the folds
    for (UInt i = 0; i < K; i++){

      std::vector<UInt> x_valid, x_train;

      if (i < N % K){ // fold grossi
        std::set_union(K_folds_.cbegin(), K_folds_.cbegin()+ i*(N/K +1), K_folds_.cbegin()+ (i + 1)*(N/K +1), K_folds_.cend(), std::back_inserter(x_train));
        std::copy(K_folds_.cbegin()+ i*(N/K +1), K_folds_.cbegin()+ (i + 1)*(N/K +1), std::back_inserter(x_valid));
      }
      else{ //fold piccoli
        std::set_union(K_folds_.cbegin(), K_folds_.cbegin()+ (N % K) + i*(N/K), K_folds_.cbegin()+ (N % K) + (i+1)*(N/K) , K_folds_.cend(), std::back_inserter(x_train));
        std::copy(K_folds_.cbegin()+ (N % K) + i*(N/K), K_folds_.cbegin()+ (N % K) + (i+1)*(N/K), std::back_inserter(x_valid));
      }

      // train
      this->data_index_ = x_train;
      this-> computeStartingDensities();

      // error
      SpMat Psi_valid = this->dataProblem_.computePsi(x_valid);

      for(UInt j=0; j<this->niter_; j++){
        cv_errors_[j] += this->error_(Psi_valid, this->init_proposals_[j]);
      }

    }

    init_best_ = std::distance(cv_errors_.cbegin(), std::min_element(cv_errors_.cbegin(), cv_errors_.cend()));

    Rprintf("The initialization selected is the number %d\n", init_best_);

    // totale
    this->data_index_.resize(this->dataProblem_.dataSize());
    std::iota(this->data_index_.begin(),this->data_index_.end(),0);

    this-> computeStartingDensities();

}


template<UInt ORDER, UInt mydim, UInt ndim>
const VectorXr*
Heat_CV<ORDER, mydim, ndim>::chooseInitialization(Real lambda) const{

  return &(this->init_proposals_[init_best_]);

}

#endif
