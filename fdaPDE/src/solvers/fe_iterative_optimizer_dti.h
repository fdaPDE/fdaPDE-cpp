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

#ifndef __FE_ITERATIVE_OPTIMIZER_DTI__
#define __FE_ITERATIVE_OPTIMIZER_DTI__

#include "header_check.h"

using namespace Eigen::indexing;

namespace fdapde {
namespace internals {

struct dwi_data {
   protected:
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
   public:
    vector_t b;
    matrix_t g;
    vector_t S0;
    matrix_t S;
};

template <typename Derived> struct fe_it_opt_dti {
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
    using solver_category = it_solver;

    // default constructor
    fe_it_opt_dti() noexcept = default;
    // construct from geoframe
    template <typename GeoFrame, typename InfoT>
        requires(is_valid_info_t<InfoT>::value)
    fe_it_opt_dti(const vector_t& b, const matrix_t& g, const GeoFrame& gf, InfoT&& info) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        using BilinearForm = std::tuple_element_t<0, std::decay_t<decltype(info.penalty)>>;
        using FeSpace = typename BilinearForm::TrialSpace;
        using DofHandler = typename FeSpace::DofHandlerType;
        using Triangulation = typename FeSpace::Triangulation;
        constexpr int embed_dim = Triangulation::embed_dim;
        fdapde_assert(gf.n_layers() == 1);
        n_locs_ = gf[0].rows();
        n_obs_ = gf[0].template col<double>("S").as_matrix().cols();
        n_dim_ = embed_dim;

        data_ = dwi_data {
          b,
          g,
          gf[0].template col<double>("S0").as_matrix(),
          gf[0].template col<double>("S").as_matrix(),
        };

        const Triangulation& triangulation = gf.template triangulation<0>();
        const FeSpace& fe_space = std::get<0>(info.penalty).trial_space();
        const DofHandler& dof_handler = fe_space.dof_handler();

        discretize(info.penalty);
        eval_basis_at_(gf);
    }

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
    template <typename GeoFrame> void analyze_data(const GeoFrame& gf) {
        fdapde_static_assert(GeoFrame::Order == 1, THIS_CLASS_IS_FOR_ORDER_ONE_GEOFRAMES_ONLY);
        fdapde_assert(gf.n_layers() == 1);
        gf = gf;
        n_locs_ = gf[0].rows();
        n_obs_ = gf[0].template col<double>("S").as_matrix().cols();
        // n_dim_ = ; ?? serve
        eval_basis_at_(gf);   // update \Psi matrix
        return;
    }

    vector_t get_row(const vector_t& V, int row, int size) const {
        int n = V.size() / size;
        vector_t v(n);
        for (int i = 0; i < n; ++i) { v[i] = V[i * size + row]; }
        return v;
    }

    void write_row(const vector_t& v, vector_t& V, int row, int size) const {
        int n = v.size();
        for (int i = 0; i < n; ++i) { V[i * size + row] = v[i]; }
    }

    vector_t get_component(const vector_t& V, int k, int size) const {
        vector_t v(size);
        v = V(seqN(k * size, size));
        return v;
    }

    matrix_t to_matrix(const vector_t& V, int size) const {
        int n = V.size() / size;
        matrix_t M(size, n);
        for (int k = 0; k < n; ++k) M.col(k) = get_component(V, k, n_dofs_);
        return M;
    }

    vector_t to_vector(const matrix_t& V) const {
        int n = V.cols();
        int size = V.rows();
        vector_t v(size * n);
        for (int k = 0; k < n; ++k) v(seqN(k * size, size)) = V.col(k);
        return v;
    }

    // matrix view of a symmetrix matrix in vector format
    matrix_t matrix_view(const vector_t& s) const {
        // assert to check that the lenght of s il ok
        matrix_t S(n_dim_, n_dim_);
        for (int i = 0; i < n_dim_; ++i) S(i, i) = s[i];
        int index = 0;
        for (int i = 0; i < n_dim_; ++i) {
            for (int j = i + 1; j < n_dim_; ++j) {
                S(j, i) = S(i, j) = s[n_dim_ + index];
                index++;
            }
        }
        return S;
    }

    vector_t vector_view(const matrix_t& S) const {
        // assert to check that the lenght of s il ok
        vector_t s(n_dim_ * (n_dim_ + 1) / 2);
        for (int i = 0; i < n_dim_; ++i) s[i] = S(i, i);
        int index = 0;
        for (int i = 0; i < n_dim_; ++i) {
            for (int j = i + 1; j < n_dim_; ++j) {
                s[n_dim_ + index] = S(i, j);
                index++;
            }
        }
        return s;
    }

    // Attenzione: la formula fornita è R^T M R, ma la decomposizione spettrale di L è R S R^T.
    // Se la derivata direzionale è definita come d_G exp(L), e la formula è R^T M R,
    // significa che R è la matrice degli autovettori di L.
    // In base alla formula d_G exp.L = R^T M R
    // l'immagine mostra M = d_RGR^T exp.S
    // e poi d_G exp.L = R^T M R
    // L'uso di RGR^T per M e poi R^TMR non è immediato.
    // Se L = R S R^T, allora exp(L) = R exp(S) R^T.
    // La derivata d_G exp(L) dovrebbe essere R d_{R^T G R} exp(S) R^T.
    // Quindi M = d_{R^T G R} exp(S), e la formula finale sarebbe R M R^T.
    // Se fosse R^T M R, allora M dovrebbe essere d_G exp(L') dove L'=RSR^T.
    // Assumiamo che la formula d_G exp.L = R^T M R sia corretta e che R sia
    // la matrice degli autovettori di L.

    matrix_t computeM(const matrix_t& G, const matrix_t& R, const vector_t& S_diag) const {
        int d = S_diag.size();
        matrix_t M = R.transpose() * G * R;

        for (int l = 0; l < d; ++l) {
            for (int m = 0; m < d; ++m) {
                double s_l = S_diag(l);
                double s_m = S_diag(m);
                if (std::abs(s_l - s_m) > 1e-9) {
                    M(l, m) *= (std::exp(s_m) - std::exp(s_l)) / (s_m - s_l);
                } else {
                    M(l, m) *= std::exp(s_l);
                }
            }
        }
        return M;
    }

    vector_t dG_exp(double b, const vector_t& g, const matrix_t& L) {
        matrix_t G = b * g * g.transpose();   // controllare che b sia necessario
        Eigen::SelfAdjointEigenSolver<matrix_t> es_L(L);
        matrix_t R = es_L.eigenvectors();
        vector_t S_diag = es_L.eigenvalues();

        matrix_t M = computeM(G, R, S_diag);

        matrix_t dG_exp_L = R * M * R.transpose();

        return vector_view(dG_exp_L);
    }

    // main fit entry point
    auto fit(double lambda, double tol = 1e-15) {
        std::cout << "fit" << std::endl;
        fdapde_assert(lambda > 0 && n_dofs_ > 0 && n_obs_ > 0);
        // check if P has already been built
        if (!built_) derived().build_P();
        // update tolerance
        tol_ = tol;
        // optimize
        BFGS<Dynamic, BacktrackingLineSearch> opt {50000, tol_, 1e-2};
        // GradientDescent<Dynamic, BacktrackingLineSearch> opt {50000, tol_, 1e-2};
        std::cout << "optimizer" << std::endl;
        vector_t vec_L = opt.optimize(
          typename Derived::opt_functor_t(derived(), lambda), vector_t::Random(n_dofs_ * (n_dim_ + 1) * n_dim_ / 2),
          [](auto value) { std::cout << "obj value: " << value << std ::endl; });
        L_ = to_matrix(vec_L, n_dofs_);
        std::cout << "optimizer" << std::endl;

        lambda_saved_ = lambda;
        std::cout << "fit" << std::endl;
        return L_;
    }
    template <typename LambdaT>
        requires(internals::is_vector_like_v<LambdaT>)
    auto fit(LambdaT&& lambda, double tol = 1e-15) {
        fdapde_assert(lambda.size() == n_lambda);
        return fit(lambda[0], tol);
    }
   private:
    template <typename ResponseT> auto fit_(ResponseT&& response, double lambda) {
        // vector_t y = y_; ??
        update_response(response);
        matrix_t L_ = fit(lambda);
        // update_response(y); ??
        return L_;
    }
   public:
    // hutchinson approximation for Tr[S]
    /*
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
    */

    // observers
    int n_dofs() const { return n_dofs_; }
    int n_obs() const { return n_obs_; }
    int n_locs() const { return n_locs_; }
    int n_dim() const { return n_dim_; }
    const sparse_matrix_t& mass() const { return R0_; }
    const sparse_matrix_t& stiff() const { return R1_; }
    const sparse_matrix_t& Psi() const { return Psi_; }
    const sparse_matrix_t& PsiNA() const { return B_.has_value() ? *B_ : Psi_; }
    const matrix_t& P(double lambda) const { return lambda * P_; }
    double LtPL(double lambda) {
        if (lambda_saved_.value() != lambda) { fit(lambda); }
        double pen = 0;
        for (int k = 0; k < L_.cols(); ++k) { pen += L_.col(k).dot(P_ * L_.col(k)); }
        return lambda * pen;
    }
    const matrix_t& L() const { return L_; }
    matrix_t Ln() const {
        matrix_t Ln = Psi_ * L_;
        return Ln;
    }
    // const binary_t& nan_pattern() const { return nan_pattern_; }
    double lambda() const { return *lambda_saved_; }
   protected:
    // matrices for hutchinson stochastic estimation of Tr[S]
    std::optional<matrix_t> Us_;
    std::optional<double> lambda_saved_ = -1;
    int n_dofs_ = 0, n_locs_ = 0, n_obs_ = 0, n_dim_ = 0;
    sparse_matrix_t R0_;    // n_dofs x n_dofs matrix [R0]_{ij} = \int_D \psi_i * \psi_j
    sparse_matrix_t R1_;    // n_dofs x n_dofs matrix [R1]_{ij} = \int_D a(\psi_i, \psi_j)
    sparse_matrix_t Psi_;   // n_obs x n_dofs matrix [Psi]_{ij} = \psi_j(p_i)
    vector_t u_;            // n_dofs x 1 vector u_i = \int_D u * \psi_i
    diag_matrix_t D_;       // vector of regions' measures (areal sampling)
    mutable sparse_solver_t invR0_;
    std::optional<sparse_matrix_t> B_;
    matrix_t P_, R0invP_;

    // basis system evaluation handle
    std::function<sparse_matrix_t(const matrix_t& locs)> point_eval_;
    std::function<std::pair<sparse_matrix_t, vector_t>(const binary_t& locs)> areal_eval_;

    dwi_data data_;
    matrix_t L_;
    // binary_t nan_pattern_;
    bool built_ = false;

    double tol_ = 1e-15;
};

// Derived classes
struct fe_it_opt_dti_linearized_gaussian_dirichlet : fe_it_opt_dti<fe_it_opt_dti_linearized_gaussian_dirichlet> {
    struct opt_functor_t {
        // constructor
        opt_functor_t(fe_it_opt_dti_linearized_gaussian_dirichlet& m, double lambda) :
            m_(std::addressof(m)), lambda_(lambda) { }

        // penalized negative log-likelihood at point
        double operator()(const vector_t& L) {
            double obj = 0;
            vector_t S0 = m_->data_.S0;
            // loss
            for (int i = 0; i < m_->n_obs(); ++i) {
                vector_t Si = m_->data_.S.col(i);
                double bi = m_->data_.b[i];
                vector_t gi = m_->data_.g.col(i);
                for (int j = 0; j < m_->n_locs_; ++j) {
                    matrix_t Lj = m_->matrix_view(m_->get_row(L, j, m_->n_dofs()));
                    obj += std::pow(log(S0[j] / Si[j]) - bi * gi.dot(Lj * gi), 2);
                }
            }
            obj /= m_->n_obs() * m_->n_locs();
            // penalty
            for (int k = 0; k < m_->n_dim() * (m_->n_dim() + 1) / 2; ++k) {
                vector_t lk = m_->get_component(L, k, m_->n_dofs());
                obj += lambda_ * lk.dot(m_->P_ * lk);
            };
            return obj;
        }
        // gradient functor
        std::function<vector_t(const vector_t&)> derive() {
            return [this](const vector_t& L) {
                matrix_t gradient_locs(m_->n_locs(), m_->n_dim() * (m_->n_dim() + 1) / 2);
                vector_t S0 = m_->data_.S0;
                // loss
                for (int i = 0; i < m_->n_obs_; ++i) {
                    vector_t Si = m_->data_.S.col(i);
                    double bi = m_->data_.b[i];
                    vector_t gi = m_->data_.g.col(i);
                    for (int j = 0; j < m_->n_locs_; ++j) {
                        matrix_t Lj = m_->matrix_view(m_->get_row(L, j, m_->n_dofs()));
                        vector_t dG_exp_L = m_->dG_exp(bi, gi, Lj);
                        vector_t gradient_loc_j = -(log(S0[j] / Si[j]) - bi * gi.dot(Lj * gi)) * dG_exp_L;
                        gradient_locs.row(j) = gradient_loc_j;
                    }
                }
                matrix_t gradient = m_->Psi().transpose() * gradient_locs;
                // penalty
                gradient += lambda_ * m_->P_ * m_->to_matrix(L, m_->n_dofs());
                return m_->to_vector(gradient);
            };
        }
        // injected optimization stopping criterion
        template <typename Optimizer> bool stop_if(Optimizer& opt) {
            double loss_old = operator()(opt.x_old);
            double loss_new = operator()(opt.x_new);
            return std::abs((loss_new - loss_old) / loss_old) < m_->tol_;
        }
       private:
        fe_it_opt_dti_linearized_gaussian_dirichlet* m_;
        double lambda_;
    };

    // penalty matrix builder
    void build_P() {
        P_ = R1_;
        built_ = true;
    }

    fe_it_opt_dti_linearized_gaussian_dirichlet() noexcept = default;
    // construct from formula + geoframe
    template <typename GeoFrame, typename InfoT>
        requires(is_valid_info_t<InfoT>::value)
    fe_it_opt_dti_linearized_gaussian_dirichlet(
      const vector_t& b, const matrix_t& g, const GeoFrame& gf, InfoT&& info) :
        fe_it_opt_dti(b, g, gf, info) { }
};

}   // namespace internals

template <typename BilinearForm, typename LinearForm> struct fe_it_opt_dti_linearized_gaussian_dirichlet {
    using solver_t = internals::fe_it_opt_dti_linearized_gaussian_dirichlet;
   private:
    struct info_t {
        std::tuple<BilinearForm, LinearForm> penalty;
    };
   public:
    fe_it_opt_dti_linearized_gaussian_dirichlet(const BilinearForm& bilinear_form, const LinearForm& linear_form) :
        info_(std::make_tuple(bilinear_form, linear_form)) { }
    const info_t& get() const { return info_; }
   private:
    info_t info_;
};

}   // namespace fdapde

#endif   // __FE_ITERATIVE_OPTIMIZER_DTI__