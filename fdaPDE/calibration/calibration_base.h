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

#ifndef __CALIBRATION_BASE_H__
#define __CALIBRATION_BASE_H__

#include <fdaPDE/utils.h>

namespace fdapde {
namespace calibration {

template <typename Calibrator, typename... Args_> class ConfiguredCalibrator {
   private:
    using CalibratorType = std::decay_t<Calibrator>;
    using ArgsPackType = std::tuple<std::decay_t<Args_>...>;
    ArgsPackType args_;   // parameter pack forwarded to calibration strategy at fit time
    CalibratorType c_;
    template <typename ModelType_, std::size_t... I> auto call_(ModelType_& model, std::index_sequence<I...> seq) {
        return c_.fit(model, std::get<I>(args_)...);   // forward stored parameters to calibrator
    }

    // templated member detection trait for set_model()
    template <typename T, typename M, typename = void> struct callback_require_model : std::false_type { };
    template <typename T, typename M>
    struct callback_require_model<
      T, M, std::void_t<decltype(std::declval<T>().template set_model<M>(std::declval<M>()))>> : std::true_type { };
   public:
    ConfiguredCalibrator(const CalibratorType& c, Args_&&... args) :
        c_(c), args_(std::make_tuple(std::forward<Args_>(args)...)) {};
    template <typename ModelType_> DVector<double> fit(ModelType_& model) {
        // forward model instance to callbacks which require it
        auto set_model = [&model](auto&& arg) {
            if constexpr (callback_require_model<std::decay_t<decltype(arg)>, ModelType_>::value) {
                arg.set_model(model);
            }
        };
        std::apply([&](auto&&... arg) { (set_model(arg), ...); }, args_);
        return call_(model, std::make_index_sequence<sizeof...(Args_)> {});
    }
};

template <typename T> struct CalibratorBase {
    template <typename... Args> ConfiguredCalibrator<T, Args...> operator()(Args&&... args) const {
        return ConfiguredCalibrator<T, Args...>(static_cast<const T&>(*this), std::forward<Args>(args)...);
    }
};

// a type-erasure wrapper for a (configured) calibration algorithm for models of type ModelType
template <typename ModelType_> struct Calibrator__ {
    template <typename T> using fn_ptrs = fdapde::mem_fn_ptrs<&T::template fit<ModelType_>>;
    DVector<double> fit(ModelType_& model) { return fdapde::invoke<DVector<double>, 0>(*this, model); }
};
template <typename ModelType_> using Calibrator = fdapde::erase<fdapde::heap_storage, Calibrator__<ModelType_>>;

}   // namespace calibration
}   // namespace fdapde

#endif   // __CALIBRATION_BASE_H__
