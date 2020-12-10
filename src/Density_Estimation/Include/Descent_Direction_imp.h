#ifndef __DESCENT_DIRECTION_IMP_H__
#define __DESCENT_DIRECTION_IMP_H__


template<UInt ORDER, UInt mydim, UInt ndim>
std::unique_ptr<DirectionBase<ORDER, mydim, ndim>>
DirectionGradient<ORDER, mydim, ndim>::clone() const {

  return make_unique<DirectionGradient<ORDER, mydim, ndim>>(*this);

}


template<UInt ORDER, UInt mydim, UInt ndim>
VectorXr
DirectionGradient<ORDER, mydim, ndim>::computeDirection(const VectorXr& g, const VectorXr& grad){

  return (- grad);
}


template<UInt ORDER, UInt mydim, UInt ndim>
DirectionBFGS<ORDER, mydim, ndim>::DirectionBFGS(const DirectionBFGS<ORDER, mydim, ndim>& rhs):
DirectionBase<ORDER, mydim, ndim>(rhs) {

  updateH_ = false;
  HInit_ = rhs.HInit_;
  HOld_ = rhs.HInit_;

}


template<UInt ORDER, UInt mydim, UInt ndim>
std::unique_ptr<DirectionBase<ORDER, mydim, ndim>>
DirectionBFGS<ORDER, mydim, ndim>::clone() const {

  return make_unique<DirectionBFGS<ORDER, mydim, ndim>>(*this);

}


template<UInt ORDER, UInt mydim, UInt ndim>
VectorXr
DirectionBFGS<ORDER, mydim, ndim>::computeDirection(const VectorXr& g, const VectorXr& grad){

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


template<UInt ORDER, UInt mydim, UInt ndim>
void
DirectionBFGS<ORDER, mydim, ndim>::resetParameters(){
  updateH_ = false;
  HOld_ = HInit_;
}

#endif
