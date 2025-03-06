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

#ifndef __FE_SEPARABLE_DRIVER_H__
#define __FE_SEPARABLE_DRIVER_H__

#include "header_check.h"

namespace fdapde {
namespace internals {

template <typename Strategy> class fe_separable_driver;
  
// solves \min_{f, \beta} \| W^{1/2} * (y_i - x_i^\top * \beta - f(p_i, t_j)) \|_2^2 + \int_D \int_T (L_D(f) - u_D)^2 +
// \int_T \int_D (L_T(f) - u_T)^2
template <> class fe_separable_driver<monolithic> {
   private:
    template <typename Tuple> struct function_space_tuple {
        using type = decltype([]<size_t... Is_>(std::index_sequence<Is_...>) {
            return std::make_tuple(typename std::tuple_element_t<Is_, Tuple>::TrialSpace {}...);
        }(std::make_index_sequence<std::tuple_size_v<Tuple>>()));
    };
    // select one between arg1 and arg2 based on condition f
    template <typename F, typename Arg1, typename Arg2>
    const auto& select_arg_(F&& f, const Arg1& arg1, const Arg2& arg2) const {
        if constexpr (f(arg1, arg2)) {
            return arg1;
        } else {
            return arg2;
        }
    }
    template <typename FuncSpace>
    static constexpr bool is_fe_space_v =
      std::is_same_v<typename FuncSpace::discretization_category, finite_element_tag>;
    template <typename Penalty1, typename Penalty2>
    const auto& fe_penalty_(const Penalty1& penalty1, const Penalty2& penalty2) const {
        return select_arg_([]<typename Pen1_, typename Pen2_>(const Pen1_&, const Pen2_&) {
            return  is_fe_space_v<typename std::tuple_element_t<0, Pen1_>::TrialSpace>;
	  }, penalty1, penalty2);
    }
    template <typename Penalty1, typename Penalty2>
    const auto& bs_penalty_(const Penalty1& penalty1, const Penalty2& penalty2) const {
        return select_arg_([]<typename Pen1_, typename Pen2_>(const Pen1_&, const Pen2_&) {
            return !is_fe_space_v<typename std::tuple_element_t<0, Pen1_>::TrialSpace>;
	  }, penalty1, penalty2);
    }
    using vector_t        = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t        = Eigen::Matrix<double, Dynamic, Dynamic>;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using diag_matrix_t   = Eigen::DiagonalMatrix<double, Dynamic, Dynamic>;
    using sparse_solver_t = Eigen::SparseLU<sparse_matrix_t>;
    using dense_solver_t  = Eigen::PartialPivLU<matrix_t>;

    template <typename GeoFrame, typename Penalty1, typename Penalty2, typename WeightMatrix>
    void init_(
      const std::string& formula, const GeoFrame& gf, Penalty1&& penalty1, Penalty2&& penalty2, const WeightMatrix& W) {
        fdapde_static_assert(
          internals::is_valid_penalty_pair_v<Penalty1> && internals::is_valid_penalty_pair_v<Penalty2>,
          INVALID_PENALTY_DESCRIPTION);
        using BilinearForms =
          std::tuple<std::tuple_element_t<0, std::decay_t<Penalty1>>, std::tuple_element_t<0, std::decay_t<Penalty2>>>;
        using FunctionSpaces = typename function_space_tuple<BilinearForms>::type;
        using FS1 = std::tuple_element_t<0, FunctionSpaces>;
        using FS2 = std::tuple_element_t<1, FunctionSpaces>;
        // one penalty must be on a FeSpace
        fdapde_static_assert(is_fe_space_v<FS1> || is_fe_space_v<FS2>, NO_FINITE_ELEMENT_SPACE_DETECTED);
	constexpr int fe_space_index = is_fe_space_v<FS1> ? 0 : 1;
	using FeSpace = std::tuple_element_t<fe_space_index, FunctionSpaces>;
	constexpr int bs_space_index = is_fe_space_v<FS1> ? 1 : 0;
	using BsSpace = std::tuple_element_t<bs_space_index, FunctionSpaces>;
        // we enforce a space-time (or SpaceMajor) expansion of the field by reordering the forms so that, index 0
        // always refer to the spatial finite element discretization
	const auto& fe_penalty = fe_penalty_(penalty1, penalty2);
	const auto& bs_penalty = bs_penalty_(penalty1, penalty2);
        // get references to bilinear and linear forms
        auto bilinear_form = std::tie(std::get<0>(fe_penalty), std::get<0>(bs_penalty));
        auto linear_form   = std::tie(std::get<1>(fe_penalty), std::get<1>(bs_penalty));
        {
            const BsSpace& bs_space = std::get<bs_space_index>(bilinear_form).trial_space();
            fdapde_assert(bs_space.sobolev_regularity() > 1);
        }	
        // discretization
        auto assemble_ = [&, this]<int Index>() {
            auto& space = std::get<Index>(bilinear_form).trial_space();
            // assemble mass matrix
            TrialFunction u(space);
            TestFunction  v(space);
            R0__[Index] = integral(space.triangulation())(u * v).assemble();	    
            R1__[Index] = std::get<Index>(bilinear_form).assemble();
        };
        assemble_.template operator()<0>();
        assemble_.template operator()<1>();
        // tensorization
        R0_ = kronecker(R0__[1], R0__[0]);   // R0_T \kron R0_D
        R1_ = kronecker(R0__[1], R1__[0]);   // R0_T \kron R1_D
        K_  = kronecker(R1__[1], R0__[0]);   // R1_T \kron R0_D
        // number of basis functions on physical domain
        n_dofs__[0] = std::get<0>(bilinear_form).trial_space().n_dofs();
        n_dofs__[1] = std::get<1>(bilinear_form).trial_space().n_dofs();
        n_dofs_ = n_dofs__[0] * n_dofs__[1];
	// forcing discretization
        u_.resize(n_dofs_);
        {
            vector_t u = std::get<0>(linear_form).assemble();
            for (int i = 0; i < n_dofs__[1]; ++i) { u_.segment(i * n_dofs__[0], n_dofs__[0]) = u; }
        }
        // basis system evaluation
        std::array<sparse_matrix_t, 2> Psi__;
        internals::for_each_index_in_pack<2>([&]<int Ns>() {
            switch (gf.category(0)[Ns]) {
            case ltype::point: {
                const auto& spatial_index = geo_index_cast<Ns, POINT>(gf[0]);
                // evaluate basis at locations
                Psi__[Ns] = internals::point_basis_eval(std::get<Ns>(bilinear_form).trial_space(), spatial_index);
                break;
            }
            case ltype::areal: {
                const auto& spatial_index = geo_index_cast<Ns, POLYGON>(gf[0]);
                const auto& [psi, measure_vec] =
                  internals::areal_basis_eval(std::get<Ns>(bilinear_form).trial_space(), spatial_index);
                Psi__[Ns] = psi;
                vector_t D(n_obs_);
                // for (int i = 0; i < m_; ++i) { D.segment(i * measure_vec.rows(), measure_vec.rows()) = measure_vec; }
                D_ = D.asDiagonal();
                break;
            }
            }
        });
        Psi_ = kronecker(Psi__[1], Psi__[0]);
        D_ = vector_t::Ones(n_obs_).asDiagonal();
	
        // data extraction
        Formula formula_(formula);
        std::vector<std::string> covs;
        for (const std::string& token : formula_.rhs()) {
            if (gf.contains(token)) { covs.push_back(token); }
        }
        n_covs_ = covs.size();
        y_.resize(n_obs_);
        {
            const auto& y_data = gf[0].data().template col<double>(formula_.lhs());
            y_data.assign_to(y_);
            if (y_data.has_nan()) {   // correct \Psi for missing observations
                Psi_ = y_data.nan().as_matrix().repeat(1, n_dofs_).select(Psi_);
            }
        }
        if (n_covs_ != 0) {
            // assemble design matrix
            X_.resize(n_obs_, n_covs_);
            for (int i = 0; i < n_covs_; ++i) { gf[0].template col<double>(covs[i]).assign_to(X_.col(i)); }
            XtWX_ = X_.transpose() * W * X_;
            invXtWX_ = XtWX_.partialPivLu();
            invXtWXXtW_ = invXtWX_.solve(X_.transpose() * W);   // (X^\top * W * X)^{-1} * X^\top * W
            // woodbury decomposition matrices
            U_ = matrix_t::Zero(2 * n_dofs_, n_covs_);
            U_.block(0, 0, n_dofs_, n_covs_) = Psi_.transpose() * D_ * W * X_;
            V_ = matrix_t::Zero(n_covs_, 2 * n_dofs_);
            V_.block(0, 0, n_covs_, n_dofs_) = X_.transpose() * W * Psi_;
        }
	return;
    }
   public:
    using solution_policy = monolithic;

    fe_separable_driver() noexcept = default;
    template <typename GeoFrame, typename Penalty1, typename Penalty2>
        requires(internals::is_pair_v<Penalty1> && internals::is_pair_v<Penalty2>)
    fe_separable_driver(const std::string& formula, const GeoFrame& gf, Penalty1&& penalty1, Penalty2&& penalty2) {
        fdapde_static_assert(GeoFrame::Order == 2, THIS_CLASS_IS_FOR_ORDER_TWO_GEOFRAMES_ONLY);
        n_obs_ = gf[0].rows();   // number of data locations on physical domain
        init_(formula, gf, penalty1, penalty2, Eigen::Matrix<double, Dynamic, 1>::Ones(n_obs_).asDiagonal());
    }
    template <typename GeoFrame, typename Penalty1, typename Penalty2, typename WeightMatrix>
        requires(internals::is_pair_v<Penalty1> && internals::is_pair_v<Penalty2>)
    fe_separable_driver(
      const std::string& formula, const GeoFrame& gf, Penalty1&& penalty1, Penalty2&& penalty2, const WeightMatrix& W) {
        fdapde_static_assert(GeoFrame::Order == 2, THIS_CLASS_IS_FOR_ORDER_TWO_GEOFRAMES_ONLY);
        n_obs_ = gf[0].rows();   // number of data locations on physical domain
        init_(formula, gf, penalty1, penalty2, W);
    }

    void operator()(double lambda_D, double lambda_T) {
        // assemble system matrix for the nonparameteric part
        SparseBlockMatrix<double, 2, 2> A_(
          -Psi_.transpose() * D_ * Psi_ - lambda_T * K_, lambda_D * R1_.transpose(), lambda_D * R1_, lambda_D * R0_);
        invA_.compute(A_);
        // linear system rhs
        vector_t b_(2 * n_dofs_);
        b_.block(n_dofs_, 0, n_dofs_, 1) = lambda_D * u_;
	
        vector_t x;
        if (n_covs_ == 0) {   // nonparametric case
            b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * y_;
            x = invA_.solve(b_);
            f_ = x.head(n_dofs_);
        } else {   // parametric case
            b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * internals::lmbQ(X_, invXtWX_, y_);
            x = woodbury_system_solve(invA_, U_, XtWX_, V_, b_);
            f_ = x.head(n_dofs_);
            beta_ = invXtWXXtW_ * (y_ - Psi_ * f_);
        } 
        g_ = x.tail(n_dofs_);   // PDE misfit
        return;
    }
    template <typename WeightMatrix> void operator()(double lambda_D, double lambda_T, WeightMatrix&& W) {
        // assemble system matrix for the nonparameteric part
        SparseBlockMatrix<double, 2, 2> A_(
          -Psi_.transpose() * D_ * W * Psi_ - lambda_T * K_, lambda_D * R1_.transpose(), lambda_D * R1_,
          lambda_D * R0_);
        invA_.compute(A_);
        // linear system rhs
        vector_t b_(2 * n_dofs_);
        b_.block(n_dofs_, 0, n_dofs_, 1) = lambda_D * u_;

        vector_t x;
        if (n_covs_ == 0) {   // nonparametric case
            b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * W * y_;
            x = invA_.solve(b_);
            f_ = x.head(n_dofs_);
        } else {   // parametric case
            XtWX_ = X_.transpose() * W * X_;
            invXtWX_ = XtWX_.partialPivLu();
            b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * internals::lmbQ(W, X_, invXtWX_, y_);
            // woodbury matrices
            U_.block(0, 0, n_dofs_, n_covs_) = Psi_.transpose() * D_ * W * X_;
            V_.block(0, 0, n_covs_, n_dofs_) = X_.transpose() * W * Psi_;
            // solve A * x = (A_ + U_ * (X^\top*W*X) * V_) * x = b
            x = woodbury_system_solve(invA_, U_, XtWX_, V_, b_);
            f_ = x.head(n_dofs_);
            beta_ = invXtWX_.solve(X_.transpose() * W) * (y_ - Psi_ * f_);
        }
        g_ = x.tail(n_dofs_);   // PDE misfit
        return;
    }
  
    // observers
    int n_dofs() const { return n_dofs_; }
    const sparse_matrix_t& R0() const { return R0_; }
    const sparse_matrix_t& R1() const { return R1_; }
    const sparse_matrix_t& Psi() const { return Psi_; }
    const vector_t& u() const { return u_; }
    const vector_t& f() const { return f_; }
    const vector_t& beta() const { return beta_; }
    const vector_t& g() const { return g_; }
    // penalty matrix: \lambda_D * R0_T \kron (R1_D^\top * R0_D^{-1} * R1_D) + \lambda_T * R1_T \kron R0_D
    matrix_t P(double lambda_D, double lambda_T) const {
        if (!invR0_.has_value()) { invR0_->compute(R0__[0]); }
        if (!PD_.has_value()) { PD_ = kronecker(R0__[1], R1__[0].transpose() * invR0_->solve(R1__[0])); }
        if (!PT_.has_value()) { PT_ = kronecker(R1__[1], R0__[0]); }
        return lambda_D * (*PD_) + lambda_T * (*PT_);
    }
   protected:
    int n_dofs_ = 0, n_obs_ = 0, n_covs_ = 0;
    // not tensorized quantities
    std::array<int, 2> n_dofs__;           // number of spatial and temporal degrees of freedom {n_dofs_D, n_dofs_T}
    std::array<sparse_matrix_t, 2> R0__;   // {R0_D, R0_T} = { \int_D \psi_i * \psi_j, \int_T \phi_i * \phi_j }
    std::array<sparse_matrix_t, 2> R1__;   // {R1_D, R1_T} = { \int_D a_D(\psi_i, \psi_j), \int_T a_T(\phi_T, \phi_D) }

    sparse_matrix_t R0_;    // n_dofs x n_dofs matrix R0 = R0_T \kron R0_D
    sparse_matrix_t R1_;    // n_dofs x n_dofs matrix R1 = R0_T \kron R1_D
    sparse_matrix_t K_;     // n_dofs x n_dofs matrix  K = R1_T \kron R0_D
    sparse_matrix_t Psi_;   // n_obs x n_dofs matrix Psi = Psi_T \kron Psi_D
    vector_t u_;            // (n_dofs_D * n_dofs_T) x 1 vector u = [u_1 \ldots u_n, \ldots, u_1 \ldots u_n]
    diag_matrix_t D_;       // vector of regions' measures (areal sampling)
    mutable std::optional<sparse_solver_t> invR0_;
    mutable std::optional<sparse_matrix_t> PD_;   // matrix PD = R0_T \kron (R1_D^\top * R0_D^{-1} * R1_D)
    mutable std::optional<sparse_matrix_t> PT_;   // matrix PT = R1_T \kron R0_D
    vector_t f_, beta_, g_;
    sparse_solver_t invA_;   // factorization of (2 * n_dofs) x (2 * n_dofs) nonparametric matrix

    matrix_t X_;               // n_obs x n_covs design matrix
    vector_t y_;               // n_obs x 1 observation vector
    matrix_t U_, V_;           // (2 * n_dofs) x n_covs matrices [\Psi^\top * D * W * y, 0] and [X^\top * W * \Psi, 0]
    matrix_t XtWX_;            // n_covs x n_covs matrix X^\top * W * X
    matrix_t invXtWXXtW_;      // n_covs x n_obs matrix (X^\top * X)^{-1} * (X^\top W)
    dense_solver_t invXtWX_;   // factorization of n_covs x n_covs matrix X^\top * W * X
};

template <> struct fe_separable_driver<iterative> {
   private:
    using vector_t        = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t        = Eigen::Matrix<double, Dynamic, Dynamic>;
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using diag_matrix_t   = Eigen::DiagonalMatrix<double, Dynamic, Dynamic>;
    using sparse_solver_t = Eigen::SparseLU<sparse_matrix_t>;
    using dense_solver_t  = Eigen::PartialPivLU<matrix_t>;

    class block_map_t {
        static constexpr int Order = 3;
        using Scalar = double;
        using storage_t = MdMap<Scalar, full_dynamic_extent_t<Order>, internals::layout_left>;

        storage_t data_;
        int rows_ = 0, cols_ = 0;
        int blk_rows_ = 0, blk_cols_ = 0;   // single block size
       public:
        block_map_t() noexcept = default;
        template <typename DataT>
            requires(internals::is_eigen_dense_xpr_v<DataT> &&
                     std::is_same_v<typename std::decay_t<DataT>::Scalar, double>)
        block_map_t(DataT&& data, int rows, int blk_rows, int blk_cols) :
            data_(data.data(), rows, data.cols(), (data.rows() / rows)),
            rows_(rows),
            cols_(data.cols()),
            blk_rows_(blk_rows),
            blk_cols_(blk_cols) {
            fdapde_assert(data.rows() % rows == 0 && rows % blk_rows == 0 && data.cols() % blk_cols == 0);
        }
        template <typename DataT>
            requires(internals::is_eigen_dense_xpr_v<DataT> &&
                     std::is_same_v<typename std::decay_t<DataT>::Scalar, double>)
        block_map_t(DataT&& data, int rows) :   // divide data in (data.rows() / rows) blocks of size rows x data.cols()
            data_(data.data(), rows, data.cols(), (data.rows() / rows)),
            rows_(rows),
	    cols_(data.cols()),
            blk_rows_(rows),
            blk_cols_(data.cols()) {
            fdapde_assert(data.rows() % rows == 0);
        }
        // observers
        auto operator()(int i, int k) const {   // get i-th row-block of k-th time instant
            auto slice_ = data_.template slice<2>(k);
            return slice_.as_eigen_map().block(i * blk_rows_, 0, blk_rows_, cols_);
        }
        auto topRows(int j, int k) const {   // get the first j row-blocks of the k-th time instant
            return data_.template slice<2>(k).as_eigen_map().block(0, 0, j * blk_rows_, cols_);
        }
        auto operator()(int k) const { return data_.template slice<2>(k).as_eigen_map(); }
        int size() const { return data_.size(); }
        // modifiers
        auto operator()(int i, int k) {
            auto slice_ = data_.template slice<2>(k);
            return slice_.as_eigen_map().block(i * blk_rows_, 0, blk_rows_, cols_);
        }
        auto topRows(int j, int k) {
            return data_.template slice<2>(k).as_eigen_map().block(0, 0, j * blk_rows_, cols_);
        }
        auto operator()(int k) { return data_.template slice<2>(k).as_eigen_map(); }
    };

    // J(f,g) = \sum_{k=1}^m (y^k - \Psi*f^k)^T*(y^k - \Psi*f^k) + \lambda_S*(g^k)^T*(g^k) + \lambda_T*(l^k)^T*(l^k)
    double J_(const block_map_t& y, const block_map_t& x, double lambda_D, double lambda_T) const {
        double sse = 0;
        for (int t = 0; t < m_; ++t) {
            sse += (y(t) - Psi_ * x(0, t)).squaredNorm() + lambda_D * x(1, t).squaredNorm() +
                   lambda_T * x(2, t).squaredNorm();
        }
        return sse;
    }
    template <typename GeoFrame, typename Penalty, typename WeightMatrix>
    void init_(const std::string& formula, const GeoFrame& gf, Penalty&& penalty, const WeightMatrix& W) {     
        fdapde_static_assert(internals::is_valid_penalty_pair_v<Penalty>, INVALID_PENALTY_DESCRIPTION);
	// enforce POINT indexing in time
	fdapde_assert(gf.category(0)[1] == ltype::point);
        using BilinearForm = std::tuple_element_t<0, std::decay_t<Penalty>>;
        using LinearForm = std::tuple_element_t<1, std::decay_t<Penalty>>;
        using FeSpace = typename BilinearForm::TrialSpace;
        fdapde_static_assert(
          std::is_same_v<typename FeSpace::discretization_category FDAPDE_COMMA finite_element_tag>,
          NO_FINITE_ELEMENT_SPACE_DETECTED);
        // discretization
        const BilinearForm& bilinear_form = std::get<0>(penalty);
        const LinearForm& linear_form = std::get<1>(penalty);
        n_dofs_ = bilinear_form.n_dofs();
        internals::fe_mass_assembly_loop<FeSpace> mass_assembler(bilinear_form.trial_space());
        R0_ = mass_assembler.assemble();
        R1_ = bilinear_form.assemble();
	// extract time step and number of time instants
        const auto& time_index = geo_index_cast<1, POINT>(gf[0]);
        const auto& time_coords = time_index.coordinates();
        m_ = time_coords.rows();
        fdapde_assert(m_ > 0 && time_coords.cols() == 1);
        DeltaT_ = time_coords(1, 0) - time_coords(0, 0);
        for (int i = 1; i < m_ - 1; ++i) {
            double lag_i = time_coords(i + 1, 0) - time_coords(i, 0);
            fdapde_assert(DeltaT_ > 0 && lag_i > 0 && almost_equal(DeltaT_ FDAPDE_COMMA lag_i));
        }
        u_.resize(n_dofs_ * m_);
        {
            vector_t u = linear_form.assemble();
            for (int i = 0; i < m_; ++i) { u_.segment(i * n_dofs_, n_dofs_) = u; }
        }
        // basis system evaluation
        switch (gf.category(0)[0]) {
        case ltype::point: {
            const auto& spatial_index = geo_index_cast<0, POINT>(gf[0]);
	    n_ = spatial_index.rows();
            // evaluate basis at locations
            Psi_ = internals::point_basis_eval(bilinear_form.trial_space(), spatial_index);
            D_ = vector_t::Ones(n_).asDiagonal();
            break;
        }
        case ltype::areal: {
            const auto& spatial_index = geo_index_cast<0, POLYGON>(gf[0]);
            n_ = spatial_index.rows();
            const auto& [psi, measure_vec] =
              internals::areal_basis_eval(bilinear_form.trial_space(), spatial_index);
            Psi_ = psi;
            vector_t D(n_obs_);
            for (int i = 0; i < m_; ++i) { D.segment(i * measure_vec.rows(), measure_vec.rows()) = measure_vec; }
            D_ = D.asDiagonal();
            break;
        }
        }
        // data extraction
        Formula formula_(formula);
        std::vector<std::string> covs;
        for (const std::string& token : formula_.rhs()) {
            if (gf.contains(token)) { covs.push_back(token); }
        }
        n_covs_ = covs.size();
        y_.resize(n_obs_);
        {
            const auto& y_data = gf[0].data().template col<double>(formula_.lhs());
            y_data.assign_to(y_);
            if (y_data.has_nan()) {   // correct \Psi for missing observations
                Psi_ = y_data.nan().as_matrix().repeat(1, n_dofs_).select(Psi_);
            }
        }
        return;
    }
   public:
    using solution_policy = iterative;

    fe_separable_driver() noexcept = default;
    template <typename GeoFrame, typename Penalty, typename WeightMatrix>
        requires(internals::is_pair_v<Penalty>)
    fe_separable_driver(
      const std::string& formula, const GeoFrame& gf, Penalty&& penalty, const WeightMatrix& W, double tol,
      int max_iter) :
        tol_(tol), max_iter_(max_iter) {
        fdapde_static_assert(GeoFrame::Order == 2, THIS_CLASS_IS_FOR_ORDER_TWO_GEOFRAMES_ONLY);
        n_obs_ = gf[0].rows();   // number of data locations on physical domain
        init_(formula, gf, penalty, W);
    }
    template <typename GeoFrame, typename Penalty, typename WeightMatrix>
        requires(internals::is_pair_v<Penalty>)
    fe_separable_driver(const std::string& formula, const GeoFrame& gf, Penalty&& penalty, const WeightMatrix& W) :
        fe_separable_driver(formula, gf, penalty, W, 1e-4, 50) { }
    template <typename GeoFrame, typename Penalty>
        requires(internals::is_pair_v<Penalty>)
    fe_separable_driver(const std::string& formula, const GeoFrame& gf, Penalty&& penalty) :
        fe_separable_driver(formula, gf, penalty, Eigen::Matrix<double, Dynamic, 1>::Ones(gf[0].rows()).asDiagonal()) {
    }

    void operator()(double lambda_D, double lambda_T) {
        // define auxiliary structures
        block_map_t y(y_, n_);
        block_map_t u(u_, n_dofs_);
        matrix_t x_old_buff(3 * n_dofs_ * m_, y_.cols()), x_new_buff(3 * n_dofs_ * m_, y_.cols());
        block_map_t x_old(x_old_buff, 3 * n_dofs_, n_dofs_, y_.cols());
        block_map_t x_new(x_new_buff, 3 * n_dofs_, n_dofs_, y_.cols());
        double alpha = lambda_T / std::pow(DeltaT_, 2);

        {   // compute starting point (f^(k,0), g^(k,0), l^(k,0)) k = 1 ... m
            SparseBlockMatrix<double, 2, 2> A_(
              Psi_.transpose() * D_ * Psi_, lambda_D * R1_.transpose(), lambda_D * R1_, -lambda_D * R0_);
            invA_.compute(A_);
            vector_t b_(2 * n_dofs_);
            for (int t = 0; t < m_; ++t) {
                b_ << Psi_.transpose() * D_ * y(t), lambda_D * u(t);
                x_old.topRows(2, t) = invA_.solve(b_);
            }
            x_old(2, 0).setZero();
            x_old(2, m_ - 1).setZero();
            for (int t = 1; t < m_ - 1; ++t) {
                x_old(2, t) = (x_old(0, t + 1) - 2 * x_old(0, t) + x_old(0, t - 1)) / std::pow(DeltaT_, 2);
            }
        }
        // iterative scheme initialization
        double Jold = std::numeric_limits<double>::max();
        double Jnew = J_(y, x_old, lambda_D, lambda_T);
        int i = 1;
        sparse_matrix_t Zero(n_dofs_, n_dofs_);
        SparseBlockMatrix<double, 3, 3> A_(
          Psi_.transpose() * D_ * Psi_, lambda_D * R1_.transpose(), -2 * alpha * R0_, lambda_D * R1_, -lambda_D * R0_,
          Zero, -2 * alpha * R0_, Zero, -lambda_T * R0_);
        invA_.compute(A_);
	vector_t b_(3 * n_dofs_);
        // iterative loop
        while (i < max_iter_ && std::abs((Jnew - Jold) / Jnew) > tol_) {
            // at step 0: f^(-1,i-1) = l^(-1,i-1) = 0
            b_ << Psi_.transpose() * D_ * y(0) - alpha * R0_ * x_old(2, 1), lambda_D * u(0),
              -alpha * R0_ * x_old(0, 1);
            x_new.topRows(2, 0) = invA_.solve(b_).topRows(2 * n_dofs_);   // l^(0) = 0
            // general step
            for (int t = 1; t < m_ - 1; ++t) {
                b_ << Psi_.transpose() * D_ * y(t) - alpha * R0_ * (x_old(2, t + 1) + x_old(2, t - 1)),
                  lambda_D * u(t), -alpha * R0_ * (x_old(0, t + 1) + x_old(0, t - 1));
                x_new(t) = invA_.solve(b_);
            }
            // at step m_ - 1: f^(m+1,i-1) = l^(m+1, i-1) = 0
            b_ << Psi_.transpose() * D_ * y(m_ - 1) - alpha * R0_ * x_old(2, m_ - 1), lambda_D * u(m_ - 1),
              -alpha * R0_ * x_old(2, m_ - 1);
            x_new.topRows(2, m_ - 1) = invA_.solve(b_).topRows(2 * n_dofs_);   // l^(m_ - 1) = 0
            // prepare for next iteration
            Jold = Jnew;
            Jnew = J_(y, x_new, lambda_D, lambda_T);
            x_old = x_new;
            i++;
	}
        // store result
        f_.resize(n_dofs_ * m_, y_.cols());
        g_.resize(n_dofs_ * m_, y_.cols());
        for (int i = 0; i < m_; ++i) {
            f_.middleRows(i * n_dofs_, n_dofs_) = x_new(0, i);
            g_.middleRows(i * n_dofs_, n_dofs_) = x_new(1, i);
        }
        return;
    }
    // observers
    int n_dofs() const { return n_dofs_; }
    const sparse_matrix_t& R0() const { return R0_; }
    const sparse_matrix_t& R1() const { return R1_; }
    const sparse_matrix_t& Psi() const { return Psi_; }
    const vector_t& u() const { return u_; }
    const vector_t& f() const { return f_; }
    const vector_t& beta() const { return beta_; }
    const vector_t& g() const { return g_; }
   protected:
    int max_iter_ = 50;   // maximum number of iterations
    double tol_ = 1e-4;   // convergence tolerance
    int n_dofs_ = 0, n_obs_ = 0, n_covs_ = 0;
    double DeltaT_ = 0;
    sparse_matrix_t R0_;    // n_dofs x n_dofs matrix [R0]_{ij} = \int_D \psi_i * \psi_j
    sparse_matrix_t R1_;    // n_dofs x n_dofs matrix [R1]_{ij} = \int_D a(\psi_i, \psi_j)
    sparse_matrix_t Psi_;   // n_obs x n_dofs matrix [Psi]_{ij} = \psi_j(p_i)
    vector_t u_;            // (n_dofs * m) x 1 vector u = [u_1 + + R0_*s / DeltaT, u_2, \ldots, u_n]
    diag_matrix_t D_;       // vector of regions' measures (areal sampling)
    mutable std::optional<sparse_solver_t> invR0_;
    vector_t f_, beta_, g_;
    sparse_solver_t invA_;   // factorization of nonparametric matrix

    int n_, m_;    // number of spatial and temporal locations
    vector_t y_;   // n_obs x 1 observation vector
};

}   // namespace internals
}   // namespace fdapde

#endif // __FE_SEPARABLE_DRIVER_H__
