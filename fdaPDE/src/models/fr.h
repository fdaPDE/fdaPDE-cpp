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

#ifndef __FUNCTIONAL_REGRESSION_H__
#define __FUNCTIONAL_REGRESSION_H__

#include "header_check.h"

namespace fdapde {

template <typename VariationalSolver>
    requires(std::is_same_v<typename VariationalSolver::solver_category, ls_solver>)
class FRPDE {
   private:
    using solver_t = std::decay_t<VariationalSolver>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using binary_t = BinaryMatrix<Dynamic, Dynamic>;
    static constexpr int n_lambda = solver_t::n_lambda;
   public:
    FRPDE() noexcept = default;
    template <typename GeoFrame, typename Penalty>
    FRPDE(const std::string& formula, const GeoFrame& gf, Penalty&& penalty) noexcept :
        solver_(), geo_category_(gf[0].category().begin(), gf[0].category().end()) {
        fdapde_assert(gf.n_layers() == 1);
        Formula formula_(formula);
        n_obs_ = gf[0].rows();
        n_stat_units_ = gf[0].template col<double>("Y").as_matrix().cols();
        n_covs_ = 0;
        for (const std::string& token : formula_.rhs()) {
            if (gf.contains(token)) { n_covs_++; }
        }
        if (n_covs_) { std::cerr << "COVRIATES ARE NOT ALLOWED WITH THIS SCHEME!" << std::endl; }
        // copy data
        Y_.resize(n_stat_units_, n_obs_);
        for (int i = 0; i < n_stat_units_; ++i)
            Y_.row(i) = gf[0].template col<double>("Y").as_matrix().col(i);   // trasposing
        nan_pattern_ = na_matrix(Y_);
        // create local geoframe
        GeoFrame gf_ = gf;
        gf_[0].data().append_blk("z", vector_t::Zero(n_obs_));
        // initialize solver
        solver_ = solver_t("z ~ f", gf_, penalty.get());
    }
    template <typename... Args> auto fit(Args&&... args) {
        vector_t z;
        vector_t w {vector_t::Ones(n_stat_units_)};
        z = (~nan_pattern_).select(Y_, 0).transpose() * w;
        for (std::size_t i = 0; i < z.size(); ++i) {
            z(i, 0) /= (~nan_pattern_).col(i).select(w, 0).sum();
            // da capire bene cosa ci va, prima c'era squaredSum, perchè?
        }
        solver_.update_response(z);
        solver_.fit(std::forward<Args>(args)...);
        return solver_.f();
    }
    // observers
    const matrix_t& F() const { return solver_.f(); }
    int n_covs() const { return n_covs_; }
    int n_obs() const { return n_obs_; }
    int n_stat_units() const { return n_stat_units_; }
    double edf(int r = 100, int seed = random_seed) { return solver_.edf(r, seed); }
    const matrix_t& response() const { return Y_; }
    const matrix_t& design_matrix() const { return solver_.design_matrix(); }
    const sparse_matrix_t& weights() const { return solver_.weights(); }
    matrix_t fitted() const { return solver_.Psi() * solver_.f(); }
    matrix_t residuals() const { return Y_ - fitted().replicate(1, Y_.rows()).transpose(); }

    // Generalized Cross Validation index
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
        gcv_t(FRPDE* model, const edf_cache_t& edf_cache) :
            model_(model),
            n_(model->n_obs()),
            k_(model->n_stat_units()),
            edf_cache_(edf_cache),
            r_(100),
            seed_(random_seed) { }
        gcv_t(FRPDE* model, const edf_cache_t& edf_cache, int r, int seed) :
            model_(model), n_(model->n_obs()), k_(model->n_stat_units()), edf_cache_(edf_cache), r_(r), seed_(seed) { }
        gcv_t(FRPDE* model) : gcv_t(model, edf_cache_t()) { }
        gcv_t(FRPDE* model, int r, int seed) : gcv_t(model, edf_cache_t(), r, seed) { }

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
        FRPDE* model_;
        int n_ = 0, k_ = 0;
        edf_cache_t edf_cache_;
        // stochastic edf approximation parameter
        int r_, seed_;
    };
    gcv_t gcv() { return gcv_t(this); }
    gcv_t gcv(const typename gcv_t::edf_cache_t& edf_cache) { return gcv_t(this, edf_cache); }
    gcv_t gcv(int r, int seed) { return gcv_t(this, r, seed); }
    gcv_t gcv(const typename gcv_t::edf_cache_t& edf_cache, int r, int seed) { return gcv_t(this, edf_cache, r, seed); }
   private:
    solver_t solver_;
    int n_obs_ = 0, n_covs_ = 0, n_stat_units_ = 0;
    std::vector<ltype> geo_category_;
    vector_t f_;
    matrix_t Y_;
    binary_t nan_pattern_;
};

// deduction guide
template <typename GeoFrame, typename Penalty>
FRPDE(const std::string& formula, const GeoFrame& gf, Penalty&& solver) -> FRPDE<typename Penalty::solver_t>;

}   // namespace fdapde

#endif   //  __SPATIAL_REGRESSION_H__
