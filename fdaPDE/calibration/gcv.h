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

#ifndef __CALIBRATION_GCV_H__
#define __CALIBRATION_GCV_H__

#include "../core/fdaPDE/optimization/optimizer.h"
#include "../models/regression/gcv.h"
#include "calibration_base.h"

namespace fdapde {
namespace calibration {

// GCV calibrator (fulfills the calibration strategy concept)
template <typename RegularizationType> class GCV : public CalibratorBase<GCV<RegularizationType>> {
   private:
    models::GCV gcv_ {};
    core::Optimizer<models::GCV> opt_ {};
    template <typename T> struct is_void {
        static constexpr bool value = std::is_same_v<T, void>;
    };
   public:
    // constructor
    GCV() {
        if constexpr (!is_void<RegularizationType>::value) {
            constexpr int n_lambda = std::is_same_v<RegularizationType, models::SpaceOnly> ? 1 : 2;
            gcv_.resize(n_lambda);
        }
    }
    template <typename Optimizer_, typename EDFStrategy_> GCV(Optimizer_&& opt, EDFStrategy_&& edf) : GCV() {
        opt_ = opt;
        gcv_.set_edf_strategy(edf);
    }
    // selects best smoothing parameter of regression model by minimization of GCV index
    template <typename ModelType_, typename... Args> DVector<double> fit(ModelType_& model, Args&&... args) {
        fdapde_assert(gcv_.inner_size() != 0 && bool(opt_) == true);
        gcv_.set_model(model);
        return opt_.optimize(gcv_, std::forward<Args>(args)...);
    }
    DVector<double> optimum() { return opt_.optimum(); }              // optimal \lambda found
    const std::vector<double>& edfs() const { return gcv_.edfs(); }   // equivalent degrees of freedom q + Tr[S]
    const std::vector<double>& gcvs() const { return gcv_.gcvs(); }   // computed values of GCV index
    // setters
    void set_step(double step) { gcv_.set_step(step); }
    template <typename Optimizer_> void set_optimizer(Optimizer_&& opt) { opt_ = opt; }
    template <typename RegularizationType_> void set() {   // set GCV's domain dimension
        fdapde_static_assert(is_void<RegularizationType>::value, THIS_METHOD_IS_FOR_VOID_REGULARIZATION_ONLY);
        gcv_.resize(models::n_smoothing_parameters<RegularizationType_>::value);
    }
    template <typename EDFStrategy_> void set_edf_strategy(EDFStrategy_&& edf) { gcv_.set_edf_strategy(edf); }
};
// template argument deduction rule
template <typename Optimizer_, typename EDFStrategy_> GCV(Optimizer_&& opt, EDFStrategy_&& edf) -> GCV<void>;

}   // namespace calibration
}   // namespace fdapde

#endif   // __CALIBRATION_GCV_H__
