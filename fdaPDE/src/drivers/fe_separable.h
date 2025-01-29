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

struct fe_separable_driver_base {
   private:
    template <typename Tuple> struct function_space_tuple {
        using type = decltype([]<size_t... Is_>(std::index_sequence<Is_...>) {
            return std::make_tuple(typename std::tuple_element_t<Is_, Tuple>::TrialSpace {}...);
        }(std::make_index_sequence<std::tuple_size_v<Tuple>>()));
    };
   public:
    using VectorType = Eigen::Matrix<double, Dynamic, 1>;
    using MatrixType = Eigen::Matrix<double, Dynamic, Dynamic>;
    using SparseMatrixType = Eigen::SparseMatrix<double>;
    using DiagonalMatrixType = Eigen::DiagonalMatrix<double, Dynamic, Dynamic>;
    using SparseSolverType = Eigen::SparseLU<SparseMatrixType>;
    using DenseSolverType  = Eigen::PartialPivLU<MatrixType>;
  
    fe_separable_driver_base() noexcept = default;
    template <typename BilinearForm_, typename LinearForm_, typename GeoFrame>
        requires(
          (internals::is_tuple_v<BilinearForm_> || internals::is_pair_v<BilinearForm_>) &&
          (internals::is_tuple_v<LinearForm_> || internals::is_pair_v<LinearForm_>))
    fe_separable_driver_base(const GeoFrame& gf, BilinearForm_&& bilinear_form, LinearForm_&& linear_form) {
        using BilinearForm = std::decay_t<BilinearForm_>;
        using LinearForm   = std::decay_t<LinearForm_>;
        fdapde_static_assert(
          std::tuple_size_v<BilinearForm> == 2 && std::tuple_size_v<LinearForm> == 2,
          THIS_CLASS_IS_FOR_EXACTLY_TWO_PENALTY_TERMS_ONLY);
        // check that function spaces are appropriate for the discretization
        using FunctionSpaces = typename function_space_tuple<BilinearForm>::type;
        using FeSpace = std::conditional_t<
          std::is_same_v<typename std::tuple_element_t<0, FunctionSpaces>::discretization_category, finite_element_tag>,
          std::tuple_element_t<0, FunctionSpaces>, std::tuple_element_t<1, FunctionSpaces>>;
        using Triangulation = typename FeSpace::Triangulation;
        fdapde_static_assert(
          Triangulation::embed_dim == 1 ||
            std::is_same_v<typename FeSpace::discretization_category FDAPDE_COMMA finite_element_tag>,
          SPACE_EMBED_DIMENSION_LARGER_THAN_ONE_AND_NO_FINITE_ELEMENT_SPACE_DETECTED);
        using OtherSpace = std::conditional_t<
          std::is_same_v<std::tuple_element_t<0, FunctionSpaces>, FeSpace>, std::tuple_element_t<1, FunctionSpaces>,
          std::tuple_element_t<0, FunctionSpaces>>;
        auto& fe_bilinear_form = std::get<internals::index_of<FeSpace, FunctionSpaces>::value>(bilinear_form);
        const FeSpace& fe_space = fe_bilinear_form.trial_space();
        // check not finite element space has a sufficiently high regularity
        auto& other_bilinear_form = std::get<internals::index_of<OtherSpace, FunctionSpaces>::value>(bilinear_form);
        const OtherSpace& other_space = other_bilinear_form.trial_space();
        fdapde_assert(other_space.sobolev_regularity() > 1);

        // assemble (take care of spaces position in function call, as this will influence the tensor basis expansion)
        std::array<SparseMatrixType, 2> R0__;
        std::array<SparseMatrixType, 2> R1__;
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
	// compute tensor products
        R0_ = kronecker(R0__[1], R0__[0]);
        R1_ = kronecker(R0__[1], R1__[0]);
        K_  = kronecker(R1__[1], R0__[0]);
        // number of basis functions on physical domain
        n_dofs_ = fe_bilinear_form.n_dofs() * other_bilinear_form.n_dofs();
	u_ = Eigen::Matrix<double, Dynamic, 1>::Zero(n_dofs_); // ----------------------------------- TODO

        // evaluate basis system
        std::array<SparseMatrixType, 2> Psi__;
        switch (gf.layer_category(0).value()) {
        case ltype::point: {
            const auto& layer = gf.get_as(layer_t::point, 0);
	    // O(n) unique coordinates extraction
            auto extract_unique_coords = [](int start_col, int dim, const MatrixType& coords) {
                auto cols = coords.middleCols(start_col, dim);
                // find unique coordinates
                std::unordered_set<MatrixType, eigen_matrix_hash> coords_set;
                std::vector<int> coords_idx;
                for (int i = 0; i < cols.rows(); ++i) {
                    MatrixType p(cols.row(i));
                    if (!coords_set.contains(p)) {
                        coords_set.insert(p);
                        coords_idx.push_back(i);
                    }
                }
                MatrixType coords_unique(coords_set.size(), dim);
                int i = 0;
                for (int idx : coords_idx) { coords_unique.row(i++) = cols.row(idx); }
                return coords_unique;
            };
            // derive coordinates from geoframe
            int lhs_embed_dim = std::tuple_element_t<0, FunctionSpaces>::embed_dim;
            int rhs_embed_dim = std::tuple_element_t<1, FunctionSpaces>::embed_dim;

            Eigen::Matrix<double, Dynamic, Dynamic> lhs_coords =
              extract_unique_coords(0, lhs_embed_dim, layer.coordinates());
            Eigen::Matrix<double, Dynamic, Dynamic> rhs_coords =
              extract_unique_coords(lhs_embed_dim, rhs_embed_dim, layer.coordinates());
	    // evaluate basis at locations
	    Psi__[0] = internals::point_basis_eval(std::get<0>(bilinear_form).trial_space(), lhs_coords);
	    Psi__[1] = internals::point_basis_eval(std::get<1>(bilinear_form).trial_space(), rhs_coords);		    
	    // tensorize
	    Psi_ = kronecker(Psi__[1], Psi__[0]);
            D_ = VectorType::Ones(Psi_.rows()).asDiagonal();
	    
	    // can we avoid to evaluate the basis system if finite element and nodes coincide with coordinates?
	    // we need some idea on how to represent the data at geoframe layer
            break;
        }
        // case ltype::areal: {
        //     const auto& layer = gf.get_as(layer_t::areal, 0);
        //     const auto& [psi, measure_vect] =
        //       internals::fe_areal_basis_eval(bilinear_form.trial_space(), layer.incidence_matrix());
        //     Psi_ = psi;
        //     D_ = measure_vect.asDiagonal();   // regions' measure
        //     break;
        // }
        }
    }
    // observers
    int n_dofs() const { return n_dofs_; }
    const SparseMatrixType& R0() const { return R0_; }
    const SparseMatrixType& R1() const { return R1_; }
    const SparseMatrixType& Psi() const { return Psi_; }
    const VectorType& u() const { return u_; }
    // penalty matrix
    // MatrixType P(double lambda) const {
    //     if (invR0_.info() != Eigen::Success) { invR0_.compute(R0_); }
    //     return lambda * R1_.transpose() * invR0_.solve(R1_);
    // }
   protected:
    int n_dofs_ = 0;         // number of basis on physical domain
    SparseMatrixType R0_;    // mass matrix: [R0]_{ij} = \int_D \psi_i * \psi_j
    SparseMatrixType R1_;    // discretization of bilinear form a: [R1]_{ij} = \int_D a(\psi_i, \psi_j)
    SparseMatrixType Psi_;   // evaluation of basis system \psi_1, ..., \psi_{n_dofs_} on physical domain
    SparseMatrixType K_;
    VectorType u_;           // discretized forcing term u_i = \int_D u * \psi_i
    DiagonalMatrixType D_;   // vector of regions' measures (areal sampling)
    mutable SparseSolverType invR0_;
};

template <typename SolutionStrategy> struct fe_separable_driver_impl;

// solves \min_{f, \beta} \| W^{1/2} * (y_i - x_i^\top * \beta - f(p_i, t_j)) \|_2^2 + \int_D \int_T (L_D(f) - u_D)^2 +
// \int_T \int_D (L_T(f) - u_T)^2
template <> struct fe_separable_driver_impl<monolithic_tag> : fe_separable_driver_base {
   private:
    template <typename GeoFrame, typename WeightMatrix>
    void init_(const std::string& formula, const GeoFrame& gf, const WeightMatrix& W) {
        // parse formula
        Formula formula_(formula);
        std::vector<std::string> covs;
        for (const std::string& token : formula_.rhs()) {
            if (gf.contains(token)) { covs.push_back(token); }
        }
        q_ = covs.size();
        // extract data from geoframe
	n_obs_ = gf[0].rows();
        y_.resize(n_obs_);
        gf[0].template col<double>(formula_.lhs()).data().assign_to(y_);
        if (q_ != 0) {
            // assemble design matrix
            X_.resize(gf[0].rows(), q_);
            for (int i = 0; i < q_; ++i) { gf[0].template col<double>(covs[i]).data().assign_to(X_.col(i)); }
            XtX_ = X_.transpose() * W * X_;
            invXtX_ = XtX_.partialPivLu();
            invXtXXt_ = invXtX_.solve(X_.transpose() * W);   // (X^\top * X)^{-1} * X^\top * W
            // woodbury decomposition matrices
            U_ = MatrixType::Zero(2 * n_dofs_, q_);
            U_.block(0, 0, n_dofs_, q_) = Psi_.transpose() * D_ * W * X_;
            V_ = MatrixType::Zero(q_, 2 * n_dofs_);
            V_.block(0, 0, q_, n_dofs_) = X_.transpose() * W * Psi_;
        }
	return;
    }
   public:
    fe_separable_driver_impl() noexcept = default;
    template <typename BilinearForm, typename LinearForm, typename GeoFrame, typename WeightMatrix>
    fe_separable_driver_impl(
      const std::string& formula, const GeoFrame& gf, BilinearForm&& bilinear_form,
      LinearForm&& linear_form, const WeightMatrix& W) :
        fe_separable_driver_base(gf, bilinear_form, linear_form) {
        init_(formula, gf, W);
    }
    template <typename BilinearForm, typename LinearForm, typename GeoFrame>
    fe_separable_driver_impl(
      const std::string& formula, const GeoFrame& gf, BilinearForm&& bilinear_form, LinearForm&& linear_form) :
        fe_separable_driver_base(gf, bilinear_form, linear_form) {
        init_(formula, gf, Eigen::Matrix<double, Dynamic, 1>::Ones(y_.rows()).asDiagonal());
    }

    void operator()(double lambda_D, double lambda_T) {
        // assemble system matrix for the nonparameteric part
        SparseBlockMatrix<double, 2, 2> A_(
          -Psi_.transpose() * D_ * Psi_ - lambda_T * K_, lambda_D * R1_.transpose(), lambda_D * R1_, lambda_D * R0_);
        invA_.compute(A_);
        // linear system rhs
        VectorType b_(2 * n_dofs_);
        b_.block(n_dofs_, 0, n_dofs_, 1) = lambda_D * u_;
	
        VectorType x;
        if (q_ == 0) {   // nonparametric case
            b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * y_;
            x = invA_.solve(b_);
            f_ = x.head(n_dofs_);
        } else {   // parametric case
            b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * internals::lmbQ(X_, invXtX_, y_);
            x = woodbury_system_solve(invA_, U_, XtX_, V_, b_);
            f_ = x.head(n_dofs_);
            beta_ = invXtXXt_ * (y_ - Psi_ * f_);
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
        VectorType b_(2 * n_dofs_);
        b_.block(n_dofs_, 0, n_dofs_, 1) = lambda_D * u_;

        VectorType x;
        if (q_ == 0) {   // nonparametric case
            b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * W * y_;
            x = invA_.solve(b_);
            f_ = x.head(n_dofs_);
        } else {   // parametric case
            XtX_ = X_.transpose() * W * X_;
            invXtX_ = XtX_.partialPivLu();
            b_.block(0, 0, n_dofs_, 1) = -Psi_.transpose() * D_ * internals::lmbQ(W, X_, invXtX_, y_);
            // woodbury matrices
            U_.block(0, 0, n_dofs_, q_) = Psi_.transpose() * D_ * W * X_;
            V_.block(0, 0, q_, n_dofs_) = X_.transpose() * W * Psi_;
            // solve A * x = (A_ + U_ * (X^\top*W*X) * V_) * x = b
            x = woodbury_system_solve(invA_, U_, XtX_, V_, b_);
            f_ = x.head(n_dofs_);
            beta_ = invXtX_.solve(X_.transpose() * W) * (y_ - Psi_ * f_);
        }
        g_ = x.tail(n_dofs_);   // PDE misfit
        return;
    }
    // observers
    const VectorType& f() const { return f_; }
    const VectorType& beta() const { return beta_; }
    const VectorType& g() const { return g_; }
   private:
    int q_;                       // number of covariates
    int n_obs_;
    MatrixType X_;                // n_obs x q design matrix
    VectorType y_;                // n_obs x 1 observation vector
    MatrixType U_, V_;            // (2 * n_dofs) x q matrices [\Psi^\top * D * y, 0] and [X^\top * \Psi, 0]
    MatrixType XtX_, invXtXXt_;   // q x q matrix X^\top * X and q x n_obs matrix (X^\top * X)^{-1} * X^\top
    DenseSolverType invXtX_;      // factorization of q x q matrix X^\top * X
    SparseSolverType invA_;
    VectorType f_, beta_, g_;
};

}   // namespace internals
}   // namespace fdapde

#endif // __FE_SEPARABLE_DRIVER_H__
