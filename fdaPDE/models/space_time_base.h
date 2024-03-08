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
template <typename Model, typename RegularizationType>
class SpaceTimeBase : public ModelBase<Model> {
   public:
    using Base = ModelBase<Model>;
    static constexpr int n_lambda = n_smoothing_parameters<RegularizationType>::value;
    using Base::lambda;       // dynamic sized smoothing parameter vector
    using Base::model;        // underlying model object
    using Base::set_lambda;   // dynamic sized setter for \lambda
    // constructor
    SpaceTimeBase() = default;
    SpaceTimeBase(const DVector<double>& time) : Base(), time_(time) {};
    // setters
    void set_lambda(const SVector<n_lambda>& lambda) {
        if(lambda_ == lambda) return;
        model().runtime().set(runtime_status::is_lambda_changed);
        lambda_ = lambda;
    }
    void set_lambda_D(double lambda_D) { set_lambda(SVector<n_lambda>(lambda_D, lambda_[1])); }
    void set_lambda_T(double lambda_T) { set_lambda(SVector<n_lambda>(lambda_[0], lambda_T)); }
    void set_time_domain(const DVector<double>& time) { time_ = time; }
    // getters
    SVector<n_lambda> lambda() const { return lambda_; }
    inline double lambda_D() const { return lambda_[0]; }
    inline double lambda_T() const { return lambda_[1]; }
    const DVector<double>& time_domain() const { return time_; }           // number of nodes in time
    const DVector<double>& time_locs() const { return time_; }             // time locations where we have observations
    inline int n_temporal_locs() const { return time_.rows(); }            // number of time instants
    // destructor
    virtual ~SpaceTimeBase() = default;
   protected:
    DVector<double> time_;   // time domain [0, T]
    SVector<n_lambda> lambda_ = SVector<n_lambda>::Zero();
};

}   // namespace models
}   // namespace fdapde

#endif   // __SPACE_TIME_BASE_H__
