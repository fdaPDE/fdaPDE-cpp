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

#ifndef __VECTOR_SPATIAL_REGRESSION_H__
#define __VECTOR_SPATIAL_REGRESSION_H__

#include "header_check.h"

namespace fdapde {

template <typename VariationalSolver>
    requires(std::is_same_v<typename VariationalSolver::solver_category, ls_solver>)
class VSRPDE {
   private:
    using solver_t = std::decay_t<VariationalSolver>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using binary_t = BinaryMatrix<Dynamic, Dynamic>;
    static constexpr int n_lambda = solver_t::n_lambda;
   public:
    VSRPDE() noexcept = default;
    template <typename GeoFrame, typename Penalty>
    VSRPDE(const std::string& formula, const GeoFrame& gf, Penalty&& penalty) noexcept :
        solver_(), geo_category_(gf[0].category().begin(), gf[0].category().end()) {
        fdapde_assert(gf.n_layers() == 1);
        Formula formula_(formula);
        n_obs_ = gf[0].rows();
        n_comp_ = gf[0].cols();
        n_covs_ = 0;
        for (const std::string& token : formula_.rhs()) {
            if (gf.contains(token)) { n_covs_++; }
        }
        if (n_covs_) { std::cerr << "COVRIATES ARE NOT ALLOWED WITH THIS SCHEME!" << std::endl; }
        // room for results
        F_.resize(std::get<0>(penalty.get().penalty).n_dofs(), n_comp_);
        // copy data
        Y_.resize(n_obs_, n_comp_);
        for (int i = 0; i < n_comp_; ++i) gf[0].template col<double>(i).assign_to(Y_.col(i));
        nan_pattern_ = na_matrix(Y_);   // compute missingness pattern
        // create local geoframe
        vector_t zero {vector_t::Zero(n_obs_)};
        GeoFrame gf_ = gf;
        geo_cast<POINT>(gf_[0]).load_vec("yi", zero);
        // initialize solver
        solver_ = solver_t("yi ~ fi", gf, penalty.get());
    }
    template <typename... Args> auto fit(Args&&... args) {
        for (int i = 0; i < n_comp_; ++i) {
            solver_.update_response(Y_.col(i));
            F_.col(i) = solver_.fit(std::forward<Args>(args)...).first;
        }
        return F_;
    }
    // observers
    const matrix_t& F() const { return F_; }
    int n_covs() const { return n_covs_; }
    int n_obs() const { return n_obs_; }
    int n_comp() const { return n_comp_; }
    double edf(int r = 100, int seed = random_seed) { return solver_.edf(r, seed); }
    const matrix_t& response() const { return Y_; }
    const matrix_t& design_matrix() const { return solver_.design_matrix(); }
    const sparse_matrix_t& weights() const { return solver_.weights(); }
    matrix_t fitted() const {
        matrix_t fitted_ = solver_.Psi() * F_;
        return fitted_;
    }

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
        gcv_t(VSRPDE* model, const edf_cache_t& edf_cache) :
            model_(model),
            n_(model->n_obs()),
            k_(model->n_comp()),
            edf_cache_(edf_cache),
            r_(100),
            seed_(random_seed) { }
        gcv_t(VSRPDE* model, const edf_cache_t& edf_cache, int r, int seed) :
            model_(model), n_(model->n_obs()), k_(model->n_comp()), edf_cache_(edf_cache), r_(r), seed_(seed) { }
        gcv_t(VSRPDE* model) : gcv_t(model, edf_cache_t()) { }
        gcv_t(VSRPDE* model, int r, int seed) : gcv_t(model, edf_cache_t(), r, seed) { }

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
            double dor = n_ - edf_cache_.at(lambda_vec);   // residual degrees of freedom
            matrix_t residuals = (~model_->nan_pattern_).select(model_->fitted() - model_->response(), 0);
            double gcv = ((n_ * k_) / std::pow(dor, 2)) * (residuals * residuals.transpose()).trace();
            // al posto di (n_ * k_) forse ci vuole il numero di non nan in nan_pattern_?
            return gcv;
        }
        // observers
        const edf_cache_t& edf_cache() const { return edf_cache_; }
        edf_cache_t& edf_cache() { return edf_cache_; }
       private:
        VSRPDE* model_;
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
    int n_obs_ = 0, n_covs_ = 0, n_comp_ = 0;
    std::vector<ltype> geo_category_;
    matrix_t F_;
    matrix_t Y_;
    binary_t nan_pattern_;
};

// deduction guide
template <typename GeoFrame, typename Penalty>
VSRPDE(const std::string& formula, const GeoFrame& gf, Penalty&& solver) -> VSRPDE<typename Penalty::solver_t>;

}   // namespace fdapde

#endif   //  __SPATIAL_REGRESSION_H__
