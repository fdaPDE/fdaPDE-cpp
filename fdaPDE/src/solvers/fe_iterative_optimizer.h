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

#ifndef __FE_LS_ELLIPTIC_SOLVER_IT_H__
#define __FE_LS_ELLIPTIC_SOLVER_IT_H__

#include "header_check.h"

namespace fdapde {
namespace internals {

// solves \min_{f} L(f | y, W) + \lambda * P(f)
template <typename Derived> struct fe_iterative_optimizer {
   protected:
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using binary_t = BinaryMatrix<Dynamic, Dynamic>;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using diag_matrix_t = Eigen::DiagonalMatrix<double, Dynamic, Dynamic>;
    using sparse_solver_t = eigen_sparse_solver_movable_wrap<Eigen::SparseLU<sparse_matrix_t>>;
    using dense_solver_t = Eigen::PartialPivLU<matrix_t>;
    template <typename DataLocs>
    static constexpr bool is_valid_data_locs_descriptor_v =
      std::is_same_v<DataLocs, matrix_t> || std::is_same_v<DataLocs, binary_t>;
    template <typename InfoT> struct is_valid_info_t {
        static constexpr bool value = requires(InfoT info) { info.penalty; };
    };

    // access to derived class
    Derived& derived() { return static_cast<Derived&>(*this); }
    const Derived& derived() const { return static_cast<const Derived&>(*this); }

    // evaluation of basis system at spatial locations
    template <typename DataLocs>
        requires(is_valid_data_locs_descriptor_v<DataLocs>)
    void eval_basis_at_(const DataLocs& locs) {
        fdapde_assert(n_locs_ == locs.rows());
        if constexpr (std::is_same_v<DataLocs, matrix_t>) {   // pointwise sampling
            Psi_ = point_eval_(locs);
            D_ = vector_t::Ones(n_locs_).asDiagonal() * (1. / n_obs_);
        } else {   // areal sampling
            const auto& [psi, measure_vect] = areal_eval_(locs);
            Psi_ = psi;
            D_ = measure_vect.asDiagonal();
        }
        return;
    }
    // optimized basis evaluation at geoframe
    template <typename GeoFrame> void eval_basis_at_(const GeoFrame& gf) {
        switch (gf.category(0)[0]) {
        case ltype::point: {
            const auto& spatial_index = geo_index_cast<0, POINT>(gf[0]);
            if (spatial_index.points_at_dofs()) {
                Psi_.resize(n_locs_, n_dofs_);
                Psi_.setIdentity();
            } else {
                Psi_ = point_eval_(spatial_index.coordinates());
            }
            D_ = vector_t::Ones(n_locs_).asDiagonal() * (1. / n_obs_);
            break;
        }
        case ltype::areal: {
            const auto& spatial_index = geo_index_cast<0, POLYGON>(gf[0]);
            const auto& [psi, measure_vect] = areal_eval_(spatial_index.incidence_matrix());
            Psi_ = psi;
            D_ = measure_vect.asDiagonal();
            break;
        }
        }
        return;
    }
   public:
    static constexpr int n_lambda = 1;
    using solver_category = ls_solver;

    // default constructor
    fe_iterative_optimizer() noexcept = default;
    // construct from formula + geoframe
    template <typename GeoFrame, typename InfoT, typename WeightMatrix>
        requires(is_valid_info_t<InfoT>::value)
    fe_iterative_optimizer(const std::string& formula, const GeoFrame& gf, InfoT&& info, const WeightMatrix& W) :
        W_(W) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        using BilinearForm = std::tuple_element_t<0, std::decay_t<decltype(info.penalty)>>;
        using FeSpace = typename BilinearForm::TrialSpace;
        using DofHandler = typename FeSpace::DofHandlerType;
        using Triangulation = typename FeSpace::Triangulation;
        constexpr int embed_dim = Triangulation::embed_dim;
        fdapde_assert(gf.n_layers() == 1);
        n_obs_ = gf[0].rows();
        n_locs_ = n_obs_;

        const Triangulation& triangulation = gf.template triangulation<0>();
        const FeSpace& fe_space = std::get<0>(info.penalty).trial_space();
        const DofHandler& dof_handler = fe_space.dof_handler();

        discretize(info.penalty);
        analyze_data(formula, gf, W);
    }
    template <typename GeoFrame, typename InfoT>
        requires(is_valid_info_t<InfoT>::value)
    fe_iterative_optimizer(const std::string& formula, const GeoFrame& gf, InfoT&& info) :
        fe_iterative_optimizer(formula, gf, info, vector_t::Ones(gf[0].rows()).asDiagonal()) { }
    // construct with no data
    template <typename GeoFrame, typename InfoT, typename WeightMatrix>
        requires(is_valid_info_t<InfoT>::value)
    fe_iterative_optimizer(const GeoFrame& gf, InfoT&& info, const WeightMatrix& W) : W_(W) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        using BilinearForm = std::tuple_element_t<0, std::decay_t<decltype(info.penalty)>>;
        using FeSpace = typename BilinearForm::TrialSpace;
        using DofHandler = typename FeSpace::DofHandlerType;
        using Triangulation = typename FeSpace::Triangulation;
        constexpr int embed_dim = Triangulation::embed_dim;
        fdapde_assert(gf.n_layers() == 1 && gf[0].category()[0] == ltype::point);
        n_obs_ = gf[0].rows();
        n_locs_ = n_obs_;

        const Triangulation& triangulation = gf.template triangulation<0>();
        const FeSpace& fe_space = std::get<0>(info.penalty).trial_space();
        const DofHandler& dof_handler = fe_space.dof_handler();

        discretize(info.penalty);
        eval_basis_at_(gf);
    }
    template <typename GeoFrame, typename InfoT>
        requires(is_valid_info_t<InfoT>::value)
    fe_iterative_optimizer(const GeoFrame& gf, InfoT&& info) :
        fe_iterative_optimizer(gf, info, vector_t::Ones(gf[0].rows()).asDiagonal()) { }

    // perform finite element based numerical discretization
    template <typename Penalty> void discretize(Penalty&& penalty) {
        fdapde_static_assert(internals::is_valid_penalty_pair_v<Penalty>, INVALID_PENALTY_DESCRIPTION);
        using BilinearForm = std::tuple_element_t<0, std::decay_t<Penalty>>;
        using LinearForm = std::tuple_element_t<1, std::decay_t<Penalty>>;
        using FeSpace = typename BilinearForm::TrialSpace;
        // discretization
        const BilinearForm& bilinear_form = std::get<0>(penalty);
        const LinearForm& linear_form = std::get<1>(penalty);
        n_dofs_ = bilinear_form.n_dofs();   // number of basis functions over physical domain
        internals::fe_mass_assembly_loop<FeSpace> mass_assembler(bilinear_form.trial_space());
        R0_ = mass_assembler.assemble();
        R1_ = bilinear_form.assemble();
        u_ = linear_form.assemble();
        // store handles for basis system evaluation at locations
        point_eval_ = [fe_space = bilinear_form.trial_space()](const matrix_t& locs) -> decltype(auto) {
            return internals::point_basis_eval(fe_space, locs);
        };
        areal_eval_ = [fe_space = bilinear_form.trial_space()](const binary_t& locs) -> decltype(auto) {
            return internals::areal_basis_eval(fe_space, locs);
        };
        return;
    }

    // analyze_data
    template <typename DataLocs, typename WeightMatrix>
        requires(std::is_same_v<DataLocs, matrix_t> || std::is_same_v<DataLocs, binary_t>)
    void analyze_data(const DataLocs& locs, const matrix_t& y, const WeightMatrix& W) {
        fdapde_assert(
          locs.rows() > 0 && y.rows() == locs.rows() && y.cols() == 1 && W.rows() == locs.rows() &&
          W.rows() == W.cols());
        n_obs_ = locs.rows();
        n_locs_ = n_obs_;
        n_covs_ = 0;
        eval_basis_at_(locs);   // update \Psi matrix
        update_response_and_weights(y, W);
        return;
    }
    template <typename GeoFrame, typename WeightMatrix>
    void analyze_data(const std::string& formula, const GeoFrame& gf, const WeightMatrix& W) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        fdapde_assert(gf.n_layers() == 1);
        n_obs_ = gf[0].rows();
        n_locs_ = n_obs_;
        eval_basis_at_(gf);   // update \Psi matrix

        // parse formula, extract response vector and design matrix
        Formula formula_(formula);
        std::vector<std::string> covs;
        for (const std::string& token : formula_.rhs()) {
            if (gf.contains(token)) { covs.push_back(token); }
        }
        n_covs_ = covs.size();
        if (n_covs_) { std::cerr << "COVRIATES ARE NOT ALLOWED WITH THIS SCHEME!" << std::endl; }
        const auto& y_data = gf[0].data().template col<double>(formula_.lhs());
        y_.resize(n_locs_, y_data.blk_sz());
        y_data.assign_to(y_);

        update_response_and_weights(y_, W);   // this updates also design_matrix releated matrices
        return;
    }

    // modifiers
    void update_response(const vector_t& y) {
        fdapde_assert(Psi_.rows() > 0 && y.rows() == n_locs_ && y.cols() == 1);
        y_ = y;
        // correct \Psi for missing observations
        nan_pattern_ = na_matrix(y);
        int old_n_obs = n_obs_;
        if (nan_pattern_.any()) {
            n_obs_ = n_locs_ - nan_pattern_.count();
            B_ = nan_pattern_.repeat(1, n_dofs_).select(Psi_, 0);
            y_ = nan_pattern_.select(y_, 0);
        }
        if (old_n_obs != n_obs_) { W_ *= (double)old_n_obs / n_obs_; }
        return;
    }
    template <typename WeightMatrix> void update_weights(const WeightMatrix& W) {
        fdapde_assert(Psi_.rows() > 0 && W.rows() == n_locs_ && W.rows() == W.cols());
        W_ = W;
        W_changed_ = true;
        return;
    }
    template <typename WeightMatrix> void update_response_and_weights(const vector_t& y, const WeightMatrix& W) {
        fdapde_assert(
          Psi_.rows() > 0 && y.rows() == n_locs_ && y.cols() == 1 && W.rows() == W.cols() && W.rows() == n_locs_);
        y_ = y;
        // correct \Psi for missing observations
        nan_pattern_ = na_matrix(y);
        if (nan_pattern_.any()) {
            n_obs_ = n_locs_ - nan_pattern_.count();
            B_ = nan_pattern_.repeat(1, n_dofs_).select(Psi_, 0);
            y_ = nan_pattern_.select(y_, 0);
        }
        update_weights(W);
        return;
    }

    // main fit entry point
    auto fit(double lambda, double tol = 1e-15) {
        fdapde_assert(lambda > 0 && n_dofs_ > 0 && n_obs_ > 0);
        // check if P has already been built
        if (!P_built_) derived().build_P();
        // update tolerance
        tol_ = tol;
        // optimize
        BFGS<Dynamic, BacktrackingLineSearch> opt {50000, tol_, 1e-2};
        // GradientDescent<Dynamic, BacktrackingLineSearch> opt {50000, tol_, 1e-2};
        f_ = opt.optimize(
          typename Derived::ls_t(derived(), lambda), vector_t::Random(n_dofs_)
          // , [](auto value) { std::cout << "obj value: " << value << std ::endl;  }
        );

        lambda_saved_ = lambda;
        return std::make_pair(f_, beta_);
    }
    template <typename LambdaT>
        requires(internals::is_vector_like_v<LambdaT>)
    auto fit(LambdaT&& lambda, double tol = 1e-15) {
        fdapde_assert(lambda.size() == n_lambda);
        return fit(lambda[0], tol);
    }
   private:
    template <typename ResponseT> auto fit_(ResponseT&& response, double lambda) {
        vector_t y = y_;
        update_response(response);
        auto [f, beta] = fit(lambda);
        update_response(y);
        return std::make_pair(f, beta);
    }
   public:
    // hutchinson approximation for Tr[S]
    double edf(int r = 100, int seed = random_seed) {
        fdapde_assert(lambda_saved_.has_value());
        if (!Us_.has_value()) {
            int seed_ = (seed == random_seed) ? std::random_device()() : seed;
            std::mt19937 rng(seed_);
            rademacher_distribution rademacher;
            Us_->resize(n_locs_, r);
            for (int i = 0; i < n_locs_; ++i) {
                for (int j = 0; j < r; ++j) { Us_->operator()(i, j) = rademacher(rng); }
            }
        }
        // Tr[S] \approx \sum_{i=1}^r (u_i^\top * S * u_i)
        double trS = 0;
        // std::cout << std::endl;
        // std::cout << lambda_saved_.value() << " ";
        for (int i = 0; i < r; ++i) {
            auto [f, beta] = fit_(Us_->col(i), lambda_saved_.value());
            vector_t fn(n_locs_);
            fn = Psi_ * f;
            trS += Us_->col(i).dot(fn);
        }
        // std::cout << trS << " " << r << " " << trS / r << std::endl;
        fit(lambda_saved_.value());
        return trS / r;
    }

    // observers
    int n_dofs() const { return n_dofs_; }
    const sparse_matrix_t& mass() const { return R0_; }
    const sparse_matrix_t& stiff() const { return R1_; }
    const sparse_matrix_t& Psi() const { return Psi_; }
    const sparse_matrix_t& PsiNA() const { return B_.has_value() ? *B_ : Psi_; }
    const matrix_t& P(double lambda) const { return lambda * P_; }
    double ftPf(double lambda) {
        if (lambda_saved_.value() != lambda || W_changed_) { fit(lambda); }
        return lambda * f_.dot(P_ * f_);
    }
    const vector_t& force() const { return u_; }
    const vector_t& f() const { return f_; }
    const vector_t fn() const { return Psi_ * f_; }
    const vector_t& beta() const { return beta_; }
    const matrix_t& design_matrix() const { return X_; }
    const vector_t& response() const { return y_; }
    const binary_t& nan_pattern() const { return nan_pattern_; }
    const sparse_matrix_t& weights() const { return W_; }
    double lambda() const { return *lambda_saved_; }
   protected:
    // matrices for hutchinson stochastic estimation of Tr[S]
    std::optional<matrix_t> Us_;
    std::optional<double> lambda_saved_ = -1;
    int n_dofs_ = 0, n_locs_ = 0, n_obs_ = 0, n_covs_ = 0;
    sparse_matrix_t R0_;    // n_dofs x n_dofs matrix [R0]_{ij} = \int_D \psi_i * \psi_j
    sparse_matrix_t R1_;    // n_dofs x n_dofs matrix [R1]_{ij} = \int_D a(\psi_i, \psi_j)
    sparse_matrix_t Psi_;   // n_obs x n_dofs matrix [Psi]_{ij} = \psi_j(p_i)
    vector_t u_;            // n_dofs x 1 vector u_i = \int_D u * \psi_i
    diag_matrix_t D_;       // vector of regions' measures (areal sampling)
    mutable sparse_solver_t invR0_;
    std::optional<sparse_matrix_t> B_;   // \Psi matrix corrected for missing observations

    matrix_t P_, R0invP_;   // n_dofs x n_dofs penalty matrix P_ = R1^\top * (R0)^{-1} * R1
    vector_t f_, beta_;
    // basis system evaluation handle
    std::function<sparse_matrix_t(const matrix_t& locs)> point_eval_;
    std::function<std::pair<sparse_matrix_t, vector_t>(const binary_t& locs)> areal_eval_;

    matrix_t X_;   // n_obs x n_covs design matrix
    vector_t y_;   // n_obs x 1 observation vector
    binary_t nan_pattern_;
    sparse_matrix_t W_;   // n_obs x n_obs matrix of observation weights
    bool W_changed_;
    bool P_built_ = false;

    double tol_ = 1e-15;
};

// Derived classes

struct fe_it_ls_elliptic : fe_iterative_optimizer<fe_it_ls_elliptic> {
    struct ls_t {
        // constructor
        ls_t(fe_it_ls_elliptic& m, double lambda) : m_(std::addressof(m)), lambda_(lambda) { }

        // penalized negative log-likelihood at point
        double operator()(const vector_t& f) {
            vector_t res = (m_->y_ - m_->Psi_ * f).array();
            return res.dot(m_->D_ * m_->W_ * res) + lambda_ * f.dot(m_->P_ * f);
        }
        // gradient functor
        std::function<vector_t(const vector_t&)> derive() {
            return [this](const vector_t& f) {
                vector_t res = (m_->y_ - m_->Psi_ * f).array();
                return vector_t(-2 * m_->Psi_.transpose() * m_->D_ * m_->W_ * res + 2 * lambda_ * m_->P_ * f);
            };
        }
        // injected optimization stopping criterion
        template <typename Optimizer> bool stop_if(Optimizer& opt) {
            double loss_old = operator()(opt.x_old);
            double loss_new = operator()(opt.x_new);
            return std::abs((loss_new - loss_old) / loss_old) < m_->tol_;
        }
       private:
        fe_it_ls_elliptic* m_;
        double lambda_;
    };

    // penalty matrix builder
    void build_P() {
        // sparse_matrix_t invR0 = lump(R0_);
        // for (int k = 0; k < invR0.outerSize(); ++k)
        //     for (sparse_matrix_t::InnerIterator it(invR0, k); it; ++it) { it.valueRef() = 1. / it.value(); }
        sparse_solver_t invR0;
        invR0.compute(R0_);
        P_ = R1_.transpose() * invR0.solve(R1_);
        P_built_ = true;
    }

    fe_it_ls_elliptic() noexcept = default;
    // construct from formula + geoframe
    template <typename GeoFrame, typename InfoT, typename WeightMatrix>
        requires(is_valid_info_t<InfoT>::value)
    fe_it_ls_elliptic(const std::string& formula, const GeoFrame& gf, InfoT&& info, const WeightMatrix& W) :
        fe_iterative_optimizer(formula, gf, info, W) { }
    template <typename GeoFrame, typename InfoT>
        requires(is_valid_info_t<InfoT>::value)
    fe_it_ls_elliptic(const std::string& formula, const GeoFrame& gf, InfoT&& info) :
        fe_iterative_optimizer(formula, gf, info) { }
    // construct with no data
    template <typename GeoFrame, typename InfoT, typename WeightMatrix>
        requires(is_valid_info_t<InfoT>::value)
    fe_it_ls_elliptic(const GeoFrame& gf, InfoT&& info, const WeightMatrix& W) : fe_iterative_optimizer(gf, info, W) { }
    template <typename GeoFrame, typename InfoT>
        requires(is_valid_info_t<InfoT>::value)
    fe_it_ls_elliptic(const GeoFrame& gf, InfoT&& info) : fe_iterative_optimizer(gf, info) { }
};

struct fe_it_ls_dirichlet : fe_iterative_optimizer<fe_it_ls_dirichlet> {
    struct ls_t {
        // constructor
        ls_t(fe_it_ls_dirichlet& m, double lambda) : m_(std::addressof(m)), lambda_(lambda) { }

        // penalized negative log-likelihood at point
        double operator()(const vector_t& f) {
            vector_t res = (m_->y_ - m_->Psi_ * f).array();
            return res.dot(m_->D_ * m_->W_ * res) + lambda_ * f.dot(m_->P_ * f);
        }
        // gradient functor
        std::function<vector_t(const vector_t&)> derive() {
            return [this](const vector_t& f) {
                vector_t res = (m_->y_ - m_->Psi_ * f).array();
                return vector_t(-2 * m_->Psi_.transpose() * m_->D_ * m_->W_ * res + 2 * lambda_ * m_->P_ * f);
            };
        }
        // injected optimization stopping criterion
        template <typename Optimizer> bool stop_if(Optimizer& opt) {
            double loss_old = operator()(opt.x_old);
            double loss_new = operator()(opt.x_new);
            return std::abs((loss_new - loss_old) / loss_old) < m_->tol_;
        }
       private:
        fe_it_ls_dirichlet* m_;
        double lambda_;
    };

    // penalty matrix builder
    void build_P() {
        P_ = R1_;
        P_built_ = true;
    }

    fe_it_ls_dirichlet() noexcept = default;
    // construct from formula + geoframe
    template <typename GeoFrame, typename InfoT, typename WeightMatrix>
        requires(is_valid_info_t<InfoT>::value)
    fe_it_ls_dirichlet(const std::string& formula, const GeoFrame& gf, InfoT&& info, const WeightMatrix& W) :
        fe_iterative_optimizer(formula, gf, info, W) { }
    template <typename GeoFrame, typename InfoT>
        requires(is_valid_info_t<InfoT>::value)
    fe_it_ls_dirichlet(const std::string& formula, const GeoFrame& gf, InfoT&& info) :
        fe_iterative_optimizer(formula, gf, info) { }
    // construct with no data
    template <typename GeoFrame, typename InfoT, typename WeightMatrix>
        requires(is_valid_info_t<InfoT>::value)
    fe_it_ls_dirichlet(const GeoFrame& gf, InfoT&& info, const WeightMatrix& W) :
        fe_iterative_optimizer(gf, info, W) { }
    template <typename GeoFrame, typename InfoT>
        requires(is_valid_info_t<InfoT>::value)
    fe_it_ls_dirichlet(const GeoFrame& gf, InfoT&& info) : fe_iterative_optimizer(gf, info) { }
};

}   // namespace internals

// solver factory
template <typename BilinearForm, typename LinearForm> struct fe_it_ls_elliptic {
    using solver_t = internals::fe_it_ls_elliptic;
   private:
    struct info_t {
        std::tuple<BilinearForm, LinearForm> penalty;
    };
   public:
    fe_it_ls_elliptic(const BilinearForm& bilinear_form, const LinearForm& linear_form) :
        info_(std::make_tuple(bilinear_form, linear_form)) { }
    const info_t& get() const { return info_; }
   private:
    info_t info_;
};

template <typename BilinearForm, typename LinearForm> struct fe_it_ls_dirichlet {
    using solver_t = internals::fe_it_ls_dirichlet;
   private:
    struct info_t {
        std::tuple<BilinearForm, LinearForm> penalty;
    };
   public:
    fe_it_ls_dirichlet(const BilinearForm& bilinear_form, const LinearForm& linear_form) :
        info_(std::make_tuple(bilinear_form, linear_form)) { }
    const info_t& get() const { return info_; }
   private:
    info_t info_;
};

}   // namespace fdapde

#endif   // __FE_LS_ELLIPTIC_SOLVER_IT_H__
