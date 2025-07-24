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

template <typename VariationalSolver>
    requires(std::is_same_v<typename VariationalSolver::solver_category, ls_solver>)
class SRPDE {
   private:
    using solver_t = std::decay_t<VariationalSolver>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    static constexpr int n_lambda = solver_t::n_lambda;
   public:
    SRPDE() noexcept = default;
    template <typename GeoFrame, typename Penalty>
    SRPDE(const std::string& formula, const GeoFrame& gf, Penalty&& penalty) noexcept :
        solver_(), geo_category_(gf[0].category().begin(), gf[0].category().end()) {
        discretize(penalty.get());
        analyze_data(formula, gf);
    }
    // modifiers
    template <typename... Args> void discretize(Args&&... args) { solver_.discretize(std::forward<Args>(args)...); }
    template <typename GeoFrame, typename WeightMatrix>
    void analyze_data(const std::string& formula, const GeoFrame& gf, const WeightMatrix& W) {
        fdapde_assert(gf.n_layers() == 1);
        Formula formula_(formula);
	n_obs_  = gf[0].rows();
	n_covs_ = 0;
        for (const std::string& token : formula_.rhs()) {
            if (gf.contains(token)) { n_covs_++; }
        }
        solver_.analyze_data(formula, gf, W);
    }
    template <typename GeoFrame> void analyze_data(const std::string& formula, const GeoFrame& gf) {
        analyze_data(formula, gf, vector_t::Ones(gf[0].rows()).asDiagonal());
    }
    // fitting
    template <typename... Args> auto fit(Args&&... args) { return solver_.fit(std::forward<Args>(args)...); }
    // observers
    const vector_t& f() const { return solver_.f(); }
    const vector_t& beta() const { return solver_.beta(); }
    const vector_t& misfit() const { return solver_.misfit(); }
    int n_covs() const { return n_covs_; }
    int n_obs() const { return n_obs_; }
    double edf(int r = 100, int seed = random_seed) { return solver_.edf(r, seed); }
    const vector_t& response() const { return solver_.response(); }
    const matrix_t& design_matrix() const { return solver_.design_matrix(); }
    const sparse_matrix_t& weights() const { return solver_.weights(); }
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
        using sparse_matrix_t = Eigen::SparseMatrix<double>;

        // compute model's variance-covariance matrix for the parametric component
        void compute_Vbeta_() const {
            const matrix_t& X = m_->design_matrix();
            const sparse_matrix_t& W = m_->weights();

            matrix_t invSigma = XtWX_.inverse();
            matrix_t H = X * invSigma * X.transpose() * W;

            matrix_t e = (invSigma * X.transpose() * W * S_).transpose();
            Eigen::SparseLU<sparse_matrix_t> invW(W);
            Vbeta_ = sigma_squared_ * (invSigma + e.transpose() * invW.solve(e));
	    return;
        }
        template <typename Distribution>
        std::pair<vector_t, vector_t> confint_beta_(double alpha, const matrix_t& C, Distribution&& distr) const {
            fdapde_assert(C.cols() == q_);
            if (!Vbeta_.has_value()) { compute_Vbeta_(); }
            int p = C.rows();
            double q = distr.quantile(alpha);
            vector_t a = C * m_->beta();
            vector_t b(p);
            for (int i = 0; i < p; ++i) { b[i] = std::sqrt(q * (C.row(i) * (*Vbeta_) * C.row(i).transpose()).value()); }
            // build confidence interval
            return std::make_pair(a - b, a + b);
        }
        // compute model's variance-covariance matrix for the non-parametric component      
        void compute_Vf_() const {
            const auto& s = m_->solver_;
            matrix_t e = invT_ * s.PsiNA().transpose();
            Vf_ = sigma_squared_ * e * Q_ * e.transpose();
	    return;
        }
        template <typename Distribution>
        std::pair<vector_t, vector_t> confint_f_(double alpha, const sparse_matrix_t& Psi, Distribution&& distr) const {
            const auto& s = m_->solver_;
            fdapde_assert(Psi.rows() > 0 && Psi.cols() == s.n_dofs());
            if (!Vf_.has_value()) { compute_Vf_(); }
	    sparse_matrix_t Vf__ = Psi * (*Vf_) * Psi.transpose();
	    double q = distr.quantile(alpha);
	    vector_t a = Psi * s.f();
	    vector_t b = q * (Vf__.diagonal().array()).sqrt();
            // build confidence interval
            return std::make_pair(a - b, a + b);    
        }
       public:
        wald_t() noexcept = default;
        wald_t(const SRPDE* m, bool approx) : m_(m), q_(m->n_covs()) {
            const auto& s = m_->solver_;
            const matrix_t& X = m_->design_matrix();
            const sparse_matrix_t& W = m_->weights();
            XtWX_ = X.transpose() * W * X;

            if (approx) {
                sparse_matrix_t E =
                  s.PsiNA().transpose() * W * s.PsiNA() + s.P(s.lambda(), FSPAI(s.mass()));
                FSPAI invE(E);   // compute approximate inverse
                int n_dofs = s.n_dofs();

                invT_ = woodbury_system_solve(
                  invE, s.U().topRows(n_dofs), -XtWX_, s.V().leftCols(n_dofs), matrix_t::Identity(n_dofs, n_dofs));
            } else {
                matrix_t E = s.PsiNA().transpose() * W * s.PsiNA() + s.P(s.lambda());
                Eigen::PartialPivLU<matrix_t> invE(E);
                int n_dofs = s.n_dofs();

                invT_ = woodbury_system_solve(
                  invE, s.U().topRows(n_dofs), -XtWX_, s.V().leftCols(n_dofs), matrix_t::Identity(n_dofs, n_dofs));
            }
            Q_ = s.Q();   // request matrix Q = W(I - H)
            S_ = s.PsiNA() * invT_ * s.PsiNA().transpose() * Q_;

            // compute variance estimator \sigma^2
            vector_t eps = s.response() - m_->fitted();
            sigma_squared_ = (eps.transpose() * W * eps).value() / (m_->n_obs() - q_ - S_.trace());
        }

        // parametric confidence intervals
        // simultaneous
        std::pair<vector_t, vector_t> confint_sim_beta(double alpha, const matrix_t& C) const {
            return confint_beta_(1 - alpha, C, chi_squared_distribution(C.rows()));
        }
        auto confint_sim_beta(double alpha) const { return confint_sim_beta(alpha, matrix_t::Identity(q_, q_)); }
        // bonferroni corrected
        std::pair<vector_t, vector_t> confint_bon_beta(double alpha, const matrix_t& C) const {
            return confint_beta_(1 - alpha, C, normal_distribution(1 - alpha / (2 * C.rows())));
        }
        auto confint_bon_beta(double alpha) const { return confint_bon_beta(alpha, matrix_t::Identity(q_, q_)); }
        // one-at-a-time
        std::pair<vector_t, vector_t> confint_oat_beta(double alpha, const matrix_t& C) const {
            return confint_beta_(1 - alpha, C, normal_distribution(1 - alpha / 2));
        }
        auto confint_oat_beta(double alpha) const { return confint_oat_beta(alpha, matrix_t::Identity(q_, q_)); }
        // parametric testing
        // simultaneous
        template <typename BetaT>
            requires(internals::is_vector_like_v<BetaT>)
        double test_sim_beta(const BetaT& beta0, const matrix_t& C) const {
            fdapde_assert(beta0.size() == q_);
            if (!Vbeta_.has_value()) { compute_Vbeta_(); }
            vector_t beta0_(q_);
            for (int i = 0; i < q_; ++i) { beta0_[i] = beta0[i]; }
            matrix_t Sigma = C * (*Vbeta_) * C.transpose();
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
        // one-at-a-time
        template <typename BetaT>
            requires(internals::is_vector_like_v<BetaT>)
        vector_t test_oat_beta(const BetaT& beta0, const matrix_t& C) const {
            fdapde_assert(beta0.size() == q_);
            if (!Vbeta_.has_value()) { compute_Vbeta_(); }
            vector_t pvalue(q_);
            for (int i = 0; i < q_; ++i) {
                double sigma = (C.row(i) * (*Vbeta_) * C.col(i)).value();
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

        // non-parametric confidence interval
        template <typename... DataLocs>
        std::pair<vector_t, vector_t> confint_oat_f(double alpha, const DataLocs&... locs) const {
            fdapde_assert(
		std::all_of(m_->geo_category_.begin() FDAPDE_COMMA m_->geo_category_.end() FDAPDE_COMMA
		    [&](auto t) { return t == ltype::point; })
	    );

	    // eval_basis_at must be made public
	    // what if more than one layer and just one matrix of joint space-time points?
	    // what in case of parabolic penalty?
	    
            sparse_matrix_t Psi_p = m_->solver_.eval_basis_at(locs...);

	    
	    
	    return confint_f_(1 - alpha, Psi_p, normal_distribution(1 - alpha / 2));
        }
        auto confint_oat_f(double alpha) const { return confint_f_(alpha, m_->solver_.PsiNA()); }

        // non-parametric testing
        double test_sim_f(const vector_t& f0, const sparse_matrix_t& Psi, double tol = 1e-4) const {
            const auto& s = m_->solver_;
            fdapde_assert(f0.rows() == Psi.rows() && Psi.cols() == s.n_dofs());
            // compute pseudoinverse of matrix Vf_
            if (!Vf_.has_value()) { compute_Vf_(); }
            if (!invVf_.has_value()) {   // compute Vf pseudoinverse
                Eigen::SelfAdjointEigenSolver<matrix_t> eigenVf(Vf_);
                // discard eigenvalues smaller than tol (search first eigenvalues greater than tol)
                vector_t eigval = eigenVf.eigenvalues();
                int i = 0, n = eigval.size();
                for (; i < n && eigval[i] > tol; ++i);
                // build (rank r) pseudoinverse as V_r * D_r^{-1} * V_r^\top (V_r: eigenvectors, D_r: eigenvalues)
                r_ = n - i + 1;
                vector_t inv_eigval = eigval.tail(r_).array().inverse();
                auto eigvec = eigenVf.eigenvectors().rightCols(r_);
                invVf_ = eigvec * inv_eigval.asDiagonal() * eigvec.transpose();
            }
            vector_t fn = Psi * s.f();
            double stat = (fn - f0).transpose() * (*invVf_) * (fn - f0);
            return 1.0 - chi_squared_distribution(r_).cdf(stat);
        }
        double test_sim_f(const vector_t& f0, double tol = 1e-4) const {
            return test_oat_f(f0, m_->solver_.PsiNA(), tol);
        }
       private:
        matrix_t invT_;   // n_dofs x n_dofs matrix (\Psi^\top * Q * \Psi + P_{\lambda})^{-1}
        mutable std::optional<matrix_t> Vbeta_;
        // non-parametric testing
        mutable std::optional<matrix_t> Vf_;
        mutable std::optional<matrix_t> invVf_;
        mutable int r_;   // rank of pseudoinverse

        const SRPDE* m_;
        int q_;

      
        matrix_t XtWX_;   // q x q matrix X^\top * W * X
        double sigma_squared_;
        matrix_t Q_;
        matrix_t S_;

        // if we estimate the trace of S stochastically, and use lmbQ, we never need to assemble and store S_ nor Q_
    };
    wald_t wald(bool approx = true) const { return wald_t(this, approx); }

    class speckman_t { };
   private:
    solver_t solver_;
    int n_obs_ = 0, n_covs_ = 0;
    std::vector<ltype> geo_category_;
};

// deduction guide
template <typename GeoFrame, typename Penalty>
SRPDE(const std::string& formula, const GeoFrame& gf, Penalty&& solver)
  -> SRPDE<typename std::decay_t<Penalty>::solver_t>;

}   // namespace fdapde

#endif //  __SPATIAL_REGRESSION_H__
