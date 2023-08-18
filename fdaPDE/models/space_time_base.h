// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __SPACE_TIME_BASE_H__
#define __SPACE_TIME_BASE_H__

#include <fdaPDE/utils.h>
#include "model_base.h"

namespace fdapde {
namespace models {

// abstract base interface for any *space-time* fdaPDE statistical model. This class is not directly usable from
// models since the type of time regularization is not yet defined at this point
template <typename Model> class SpaceTimeBase : public ModelBase<Model> {
    static_assert(is_space_time<Model>::value);
   protected:
    typedef typename model_traits<Model>::PDE PDE;   // PDE used for regularization in space
    typedef ModelBase<Model> Base;
    using Base::lambda_;   // vector of smoothing parameters
    using Base::model;     // underlying model object
    using Base::pde_;      // regularizing PDE

    DVector<double> time_;   // time domain [0, T]
   public:
    // constructor
    SpaceTimeBase() = default;
    SpaceTimeBase(const PDE& pde, const DVector<double>& time) : ModelBase<Model>(pde), time_(time) {};

    // setters
    void set_lambda_D(double lambda_D) { lambda_[0] = lambda_D; }
    void set_lambda_T(double lambda_T) { lambda_[1] = lambda_T; }
    void set_time_domain(const DVector<double>& time) { time_ = time; }
    // getters
    inline double lambda_D() const { return lambda_[0]; }
    inline double lambda_T() const { return lambda_[1]; }
    const DVector<double>& time_domain() const { return time_; }           // number of nodes in time
    const DVector<double>& time_locs() const { return time_; }             // time locations where we have observations
    inline std::size_t n_temporal_locs() const { return time_.rows(); }    // number of time instants
    std::size_t n_spatial_basis() const { return pde_->domain().dof(); }   // number of basis functions in space

    // destructor
    virtual ~SpaceTimeBase() = default;
};

}   // namespace models
}   // namespace fdapde

#endif   // __SPACE_TIME_BASE_H__
