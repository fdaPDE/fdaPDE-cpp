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

#ifndef __SPATIAL_REGRESSION_H__
#define __SPATIAL_REGRESSION_H__

#include "header_check.h"

namespace fdapde {

template <typename VariationalSolver> class SRPDE {
   private:
    using solver_t = std::decay_t<VariationalSolver>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    static constexpr int n_lambda = solver_t::n_lambda;
   public:
    SRPDE() noexcept = default;
    template <typename GeoFrame, typename Penalty>
    SRPDE(const std::string& formula, const GeoFrame& gf, Penalty&& penalty) noexcept : solver_() {
        fdapde_assert(gf.n_layers() == 1);	
        Formula formula_(formula);
	n_obs_  = gf[0].rows();
	n_covs_ = 0;
        for (const std::string& token : formula_.rhs()) {
            if (gf.contains(token)) { n_covs_++; }
        }
        if constexpr (requires(Penalty p) { p.get(); }) {
            solver_ = solver_t(formula, gf, penalty.get());
        } else {
            solver_ = solver_t(formula, gf, penalty(gf.template triangulation<0>()).get());
        }
    }
    template <typename... LambdaT>
        requires(std::is_convertible_v<LambdaT, double> && ...) ||
                (sizeof...(LambdaT) == 1 && (internals::is_vector_like_v<LambdaT> && ...))
    void fit(LambdaT... lambda) {
        solver_.fit(lambda...);
    }
    // observers
    const vector_t& f() const { return solver_.f(); }
    const vector_t& beta() const { return solver_.beta(); }
    int n_covs() const { return n_covs_; }
    int n_obs() const { return n_obs_; }
    double edf(int r = 100, int seed = random_seed) { return solver_.edf(r, seed); }
    const vector_t& response() const { return solver_.response(); }
    vector_t fitted() const {
        vector_t fitted_ = solver_.fn();
        if constexpr (requires(solver_t s) { s.design_matrix(); }) {
            if (n_covs_ != 0) { fitted_ += solver_.design_matrix() * beta(); }
        }
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
        gcv_t(SRPDE* model, const edf_cache_t& edf_cache) :
            model_(model),
            n_(model->n_obs()),
            q_(model->n_covs()),
            edf_cache_(edf_cache),
            r_(100),
            seed_(random_seed) { }
        gcv_t(SRPDE* model, const edf_cache_t& edf_cache, int r, int seed) :
            model_(model), n_(model->n_obs()), q_(model->n_covs()), edf_cache_(edf_cache), r_(r), seed_(seed) { }
        gcv_t(SRPDE* model) : gcv_t(model, edf_cache_t()) { }
        gcv_t(SRPDE* model, int r, int seed) : gcv_t(model, edf_cache_t(), r, seed) { }

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
            double dor = n_ - (q_ + edf_cache_.at(lambda_vec));   // residual degrees of freedom	    
            return (n_ / std::pow(dor, 2)) * (model_->fitted() - model_->response()).squaredNorm();
        }
        // observers
        const edf_cache_t& edf_cache() const { return edf_cache_; }
        edf_cache_t& edf_cache() { return edf_cache_; }
       private:
        SRPDE* model_;
        int n_ = 0, q_ = 0;
        edf_cache_t edf_cache_;
        // stochastic edf approximation parameter
        int r_, seed_;
    };
    gcv_t gcv() { return gcv_t(this); }
    gcv_t gcv(const typename gcv_t::edf_cache_t& edf_cache) { return gcv_t(this, edf_cache); }
    gcv_t gcv(int r, int seed) { return gcv_t(this, r, seed); }
    gcv_t gcv(const typename gcv_t::edf_cache_t& edf_cache, int r, int seed) { return gcv_t(this, edf_cache, r, seed); }

    // inference
    class wald_t {
        using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
        using vector_t = Eigen::Matrix<double, Dynamic, 1>;

        template <typename Distribution>
        std::pair<vector_t, vector_t> confint_beta_(double alpha, const matrix_t& C, Distribution&& distr) const {
            fdapde_assert(C.cols() == q_);
            int p = C.rows();
            double q = distr.quantile(alpha);
            vector_t a = C * m_->beta();
            vector_t b(p);
            for (int i = 0; i < p; ++i) { b[i] = std::sqrt(q * (C.row(i) * (*V_) * C.row(i).transpose()).value()); }
            // build confidence interval
            return std::make_pair(a - b, a + b);
        }
       public:
        wald_t() noexcept = default;
        wald_t(const SRPDE* m, bool approx) : m_(m), q_(m->n_covs()) {
            // compute model's variance-covariance matrix
            const auto& s = m_->solver_;
            const matrix_t& X = s.design_matrix();
            const Eigen::SparseMatrix<double>& W = s.weights();

            matrix_t XtWX = X.transpose() * W * X;
            matrix_t invT;
            if (approx) {
                Eigen::SparseMatrix<double> E = s.PsiNA().transpose() * W * s.PsiNA() + s.P(s.lambda(), FSPAI(s.mass()));
                FSPAI invE(E);   // compute approximate inverse
                int n_dofs = s.n_dofs();

                invT = woodbury_system_solve(
                  invE, s.U().topRows(n_dofs), -XtWX, s.V().leftCols(n_dofs), matrix_t::Identity(n_dofs, n_dofs));
            } else {
                matrix_t E = s.PsiNA().transpose() * W * s.PsiNA() + s.P(s.lambda());
                Eigen::PartialPivLU<matrix_t> invE(E);
                int n_dofs = s.n_dofs();

                invT = woodbury_system_solve(
                  invE, s.U().topRows(n_dofs), -XtWX, s.V().leftCols(n_dofs), matrix_t::Identity(n_dofs, n_dofs));
            }
            // request matrix Q = W(I - H)
            matrix_t Q = s.Q();
            matrix_t S = s.PsiNA() * invT * s.PsiNA().transpose() * Q;

            // compute variance estimator \sigma^2
            vector_t eps = s.response() - m_->fitted();
            matrix_t invSigma = XtWX.inverse();
            matrix_t H = X * invSigma * X.transpose() * W;

            double sigma_squared = (eps.transpose() * W * eps).value() / (m_->n_obs() - q_ - S.trace());
            matrix_t e = (invSigma * X.transpose() * W * S).transpose();
            Eigen::SparseLU<Eigen::SparseMatrix<double>> invW(W);
            V_ = sigma_squared * (invSigma + e.transpose() * invW.solve(e));
        }

        // parametric confidence intervals
        std::pair<vector_t, vector_t> confint_sim_beta(double alpha, const matrix_t& C) const {
            return confint_beta_(1 - alpha, C, chi_squared_distribution(C.rows()));
        }
        auto confint_sim_beta(double alpha) const { return confint_sim_beta(alpha, matrix_t::Identity(q_, q_)); }
        std::pair<vector_t, vector_t> confint_bon_beta(double alpha, const matrix_t& C) const {
            return confint_beta_(1 - alpha, C, normal_distribution(1 - alpha / (2 * C.rows())));
        }
        auto confint_bon_beta(double alpha) const { return confint_bon_beta(alpha, matrix_t::Identity(q_, q_)); }
        std::pair<vector_t, vector_t> confint_oat_beta(double alpha, const matrix_t& C) const {
            return confint_beta_(1 - alpha, C, normal_distribution(1 - alpha / 2));
        }
        auto confint_oat_beta(double alpha) const { return confint_oat_beta(alpha, matrix_t::Identity(q_, q_)); }
        // parametric testing
        template <typename BetaT>
            requires(internals::is_vector_like_v<BetaT>)
        double test_sim_beta(const BetaT& beta0, const matrix_t& C) const {
            fdapde_assert(beta0.size() == q_);
            vector_t beta0_(q_);
            for (int i = 0; i < q_; ++i) { beta0_[i] = beta0[i]; }
            matrix_t Sigma = C * (*V_) * C.transpose();
            Eigen::PartialPivLU<matrix_t> invSigma(Sigma);
            double stat = ((C * m_->beta() - beta0_).transpose() * invSigma.solve(C * m_->beta() - beta0_)).value();
            return 1.0 - chi_squared_distribution(q_).cdf(stat);   // return p-value
        }
        template <typename BetaT> double test_sim_beta(const BetaT& beta0) const {
            return test_sim_beta(beta0, matrix_t::Identity(q_, q_));
        }
        double test_sim_beta(const std::initializer_list<double>& beta0, const matrix_t& C) const {
            return test_sim_beta(std::vector<double> {beta0.begin(), beta0.end()}, C);
        }
        double test_sim_beta(const std::initializer_list<double>& beta0) const {
            return test_sim_beta(std::vector<double> {beta0.begin(), beta0.end()}, matrix_t::Identity(q_, q_));
        }
        template <typename BetaT>
            requires(internals::is_vector_like_v<BetaT>)
        vector_t test_oat_beta(const BetaT& beta0, const matrix_t& C) const {
            fdapde_assert(beta0.size() == q_);
            vector_t pvalue(q_);
            for (int i = 0; i < q_; ++i) {
                double sigma = (C.row(i) * (*V_) * C.col(i)).value();
                double stat = (C.row(i).dot(m_->beta()) - beta0[i]) / std::sqrt(sigma);
                pvalue[i] = 2 * normal_distribution(0, 1).cdf(-std::abs(stat));   // compute p-value
            }
            return pvalue;
        }
        template <typename BetaT> vector_t test_oat_beta(const BetaT& beta0) const {
            return test_oat_beta(beta0, matrix_t::Identity(q_, q_));
        }
        vector_t test_oat_beta(const std::initializer_list<double>& beta0, const matrix_t& C) const {
            return test_oat_beta(std::vector<double> {beta0.begin(), beta0.end()}, C);
        }
        vector_t test_oat_beta(const std::initializer_list<double>& beta0) const {
            return test_oat_beta(std::vector<double> {beta0.begin(), beta0.end()}, matrix_t::Identity(q_, q_));
        }
       private:
        mutable std::optional<matrix_t> V_;
        const SRPDE* m_;
        int q_;
    };
    wald_t wald(bool approx = true) const { return wald_t(this, approx); }
   private:
    solver_t solver_;
    int n_obs_ = 0, n_covs_ = 0;
};

// deduction guide
template <typename GeoFrame, typename Penalty>
SRPDE(const std::string& formula, const GeoFrame& gf, Penalty&& solver) -> SRPDE<typename Penalty::solver_t>;

}   // namespace fdapde

#endif //  __SPATIAL_REGRESSION_H__
