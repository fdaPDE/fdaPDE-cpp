#ifndef __DESCENT_DIRECTION_IMP_H__
#define __DESCENT_DIRECTION_IMP_H__


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
std::unique_ptr<DirectionBase<Integrator_noPoly, ORDER, mydim, ndim>>
DirectionGradient<Integrator_noPoly, ORDER, mydim, ndim>::clone() const {

  return make_unique<DirectionGradient<Integrator_noPoly, ORDER, mydim, ndim>>(*this);

}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
VectorXr
DirectionGradient<Integrator_noPoly, ORDER, mydim, ndim>::computeDirection(const VectorXr& g, const VectorXr& grad){

  return (- grad);
}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
DirectionBFGS<Integrator_noPoly, ORDER, mydim, ndim>::DirectionBFGS(const DirectionBFGS<Integrator_noPoly, ORDER, mydim, ndim>& rhs):
DirectionBase<Integrator_noPoly, ORDER, mydim, ndim>(rhs) {

  updateH_ = false;
  HInit_ = rhs.HInit_;
  HOld_ = rhs.HInit_;

}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
std::unique_ptr<DirectionBase<Integrator_noPoly, ORDER, mydim, ndim>>
DirectionBFGS<Integrator_noPoly, ORDER, mydim, ndim>::clone() const {

  return make_unique<DirectionBFGS<Integrator_noPoly, ORDER, mydim, ndim>>(*this);

}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
VectorXr
DirectionBFGS<Integrator_noPoly, ORDER, mydim, ndim>::computeDirection(const VectorXr& g, const VectorXr& grad){

  if(updateH_){
    const VectorXr delta = g - gOld_;
    const VectorXr gamma = grad - gradOld_;

    const Real dg = delta.dot(gamma);
    const VectorXr Hg = HOld_*gamma;

    HOld_ = HOld_ + (1 + (gamma.dot(Hg))/dg)*(delta*delta.transpose())/dg -
          (Hg*delta.transpose() + delta*Hg.transpose())/dg;
  }

  gOld_ = g;
  gradOld_ = grad;

  if(!updateH_) updateH_ = true;

  return (-HOld_*grad);
}


template<typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
void
DirectionBFGS<Integrator_noPoly, ORDER, mydim, ndim>::resetParameters(){
  updateH_ = false;
  HOld_ = HInit_;
}

#endif
