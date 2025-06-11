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

#ifndef __FUNCTIONAL_DTI_H__
#define __FUNCTIONAL_DTI_H__

#include "header_check.h"

namespace fdapde {

template <typename VariationalSolver>
    requires(std::is_same_v<typename VariationalSolver::solver_category, it_solver>)
class FDTI {
   private:
    using solver_t = std::decay_t<VariationalSolver>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using binary_t = BinaryMatrix<Dynamic, Dynamic>;
    static constexpr int n_lambda = solver_t::n_lambda;
   public:
    FDTI() noexcept = default;
    template <typename GeoFrame, typename Penalty>
    FDTI(const vector_t& b, const matrix_t& g, const GeoFrame& gf, Penalty&& penalty) noexcept :
        solver_(), geo_category_(gf[0].category().begin(), gf[0].category().end()) {
        fdapde_assert(gf.n_layers() == 1);
        solver_ = solver_t(b, g, gf, penalty.get());
    }
    template <typename... Args> auto fit(Args&&... args) {
        solver_.fit(std::forward<Args>(args)...);
        return L();
    }
    // observers
    const matrix_t& L() const { return solver_.L(); }
    matrix_t Ln() const {
        matrix_t Ln = solver_.Ln();
        return Ln;
    }
    const matrix_t& D() const {
        // implement cache
        D_ = solver_.L();   // compute D starting from D
        return D_;
    }
    int n_obs() const { return solver_->n_obs_; }
    int n_locs() const { return solver_->n_locs_; }
    // double edf(int r = 100, int seed = random_seed) { return solver_.edf(r, seed); }
    matrix_t fitted() const { return solver_.Psi() * solver_.D(); }
    // matrix_t residuals() const { return Y_ - fitted().replicate(1, Y_.rows()).transpose(); }

    // Generalized Cross Validation index
    /*
    struct gcv_t : public ScalarFieldBase<n_lambda, gcv_t> {
        using Base = ScalarFieldBase<1, gcv_t>;
        static constexpr int StaticInputSize = n_lambda;
        static constexpr int NestAsRef = 0;
        static constexpr int XprBits = 0;
        using Scalar = double;
        using InputType = Vector<Scalar, StaticInputSize>;
        using edf_cache_t = std::unordered_map<
          std::array<double, StaticInputSize>, double, internals::std_array_hash<double, StaticInputSize>>;

        gcv_t() noexcept = default;
        gcv_t(FDTI* model, const edf_cache_t& edf_cache) :
            model_(model),
            n_(model->n_obs()),
            k_(model->n_stat_units()),
            edf_cache_(edf_cache),
            r_(100),
            seed_(random_seed) { }
        gcv_t(FDTI* model, const edf_cache_t& edf_cache, int r, int seed) :
            model_(model), n_(model->n_obs()), k_(model->n_stat_units()), edf_cache_(edf_cache), r_(r), seed_(seed) { }
        gcv_t(FDTI* model) : gcv_t(model, edf_cache_t()) { }
        gcv_t(FDTI* model, int r, int seed) : gcv_t(model, edf_cache_t(), r, seed) { }

        template <typename InputType_>
            requires(internals::is_subscriptable<InputType_, int>)
        constexpr double operator()(const InputType_& lambda) {
            return internals::apply_index_pack<n_lambda>([&]<int... Ns_>() { return operator()(lambda[Ns_]...); });
        }
        template <typename... LambdaT>
            requires(std::is_convertible_v<LambdaT, double> && ...) && (sizeof...(LambdaT) == StaticInputSize)
        constexpr double operator()(LambdaT... lambda) {
            model_->fit(static_cast<double>(lambda)...);
            std::array<double, StaticInputSize> lambda_vec {lambda...};
            if (edf_cache_.find(lambda_vec) == edf_cache_.end()) {   // cache Tr[S]
                edf_cache_[lambda_vec] = model_->edf(r_, seed_);
            }
            double dor = n_ * k_ - edf_cache_.at(lambda_vec);   // residual degrees of freedom
            matrix_t residuals = (~model_->nan_pattern_).select(model_->residuals(), 0);
            double gcv = ((n_ * k_) / std::pow(dor, 2)) * residuals.squaredNorm();
            // al posto di (n_) forse ci vuole il numero di non nan in nan_pattern_?
            return gcv;
        }
        // observers
        const edf_cache_t& edf_cache() const { return edf_cache_; }
        edf_cache_t& edf_cache() { return edf_cache_; }
       private:
        FDTI* model_;
        int n_ = 0, k_ = 0;
        edf_cache_t edf_cache_;
        // stochastic edf approximation parameter
        int r_, seed_;
    };
    gcv_t gcv() { return gcv_t(this); }
    gcv_t gcv(const typename gcv_t::edf_cache_t& edf_cache) { return gcv_t(this, edf_cache); }
    gcv_t gcv(int r, int seed) { return gcv_t(this, r, seed); }
    gcv_t gcv(const typename gcv_t::edf_cache_t& edf_cache, int r, int seed) { return gcv_t(this, edf_cache, r, seed); }
    */
   private:
    solver_t solver_;
    std::vector<ltype> geo_category_;
    matrix_t D_;
    binary_t nan_pattern_;
};

// deduction guide
template <typename GeoFrame, typename Penalty>
FDTI(
  const Eigen::Matrix<double, Dynamic, 1>& b, const Eigen::Matrix<double, Dynamic, Dynamic>& g, const GeoFrame& gf,
  Penalty&& solver) -> FDTI<typename Penalty::solver_t>;

}   // namespace fdapde

#endif   //  __FUNCTIONAL_DTI_H__
