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

#ifndef __FDAPDE_DRIVERS_UTILITY_H__
#define __FDAPDE_DRIVERS_UTILITY_H__

namespace fdapde {
namespace internals {

// checks if BilinearForm - LinearForm pair is valid
template <typename BilinearForm, typename LinearForm> struct is_valid_penalty_pair {
   private:
    using BilinearForm_ = std::decay_t<BilinearForm>;
    using LinearForm_ = std::decay_t<LinearForm>;
   public:
    static constexpr bool value =
      requires(BilinearForm_ bilinear_form) {   // first pair element: bilinear form
          { bilinear_form.assemble() } -> std::same_as<Eigen::SparseMatrix<double>>;
      } &&
      requires(LinearForm_ linear_form) {   // second pair element: linear form
          { linear_form.assemble() } -> std::same_as<Eigen::Matrix<double, Dynamic, 1>>;
      } && 
      std::is_same_v<
	typename BilinearForm_::discretization_category FDAPDE_COMMA typename LinearForm_::discretization_category>;
};
template <typename BilinearForm, typename LinearForm>
constexpr bool is_valid_penalty_pair_v = is_valid_penalty_pair<BilinearForm, LinearForm>::value;

// efficient left multiplication Q*x, with Q = W * (I - X * (X^\top * W * X)^{-1} * X^\top * W)
template <typename WeightMatrix, typename DesignMatrix, typename InvDesignMatrix>
Eigen::Matrix<double, Dynamic, Dynamic> lmbQ(
  const WeightMatrix& W, const DesignMatrix& X, const InvDesignMatrix& invXtWX,
  const Eigen::Matrix<double, Dynamic, Dynamic>& x) {
    using MatrixType = Eigen::Matrix<double, Dynamic, Dynamic>;
    if (X.cols() == 0) return W * x;
    MatrixType v = X.transpose() * W * x;   // X^\top * W * x
    MatrixType z = invXtWX.solve(v);        // (X^\top * W * X)^{-1} * X^\top * W * x
    // compute W*x - W*X*z = W*x - (W*X*(X^\top*W*X)^{-1}*X^\top*W)*x = W(I - H)*x = Q*x
    return W * (x - X * z);
}
template <typename DesignMatrix, typename InvDesignMatrix>
Eigen::Matrix<double, Dynamic, Dynamic>
lmbQ(const DesignMatrix& X, const InvDesignMatrix& invXtX, const Eigen::Matrix<double, Dynamic, Dynamic>& x) {
    using MatrixType = Eigen::Matrix<double, Dynamic, Dynamic>;
    if (X.cols() == 0) return x;
    MatrixType v = X.transpose() * x;   // X^\top * x
    MatrixType z = invXtX.solve(v);     // (X^\top * X)^{-1} * X^\top * x
    // compute x - X*z = x - (X*(X^\top*X)^{-1}*X^\top)*x = (I - H)*x = Q*x
    return x - X * z;
}

// pointwise basis evaluation for finite element basis system
template <typename Triangulation_, typename FeType_, typename CoordsMatrix_>
    requires(internals::is_eigen_dense_xpr_v<CoordsMatrix_>)
Eigen::SparseMatrix<double> point_basis_eval(const FeSpace<Triangulation_, FeType_>& fe_space, CoordsMatrix_&& coords) {
    static constexpr int embed_dim = Triangulation_::embed_dim;
    fdapde_assert(coords.rows() > 0 && coords.cols() == embed_dim);

    int n_shape_functions = fe_space.n_shape_functions();
    int n_dofs = fe_space.n_dofs();
    int n_locs = coords.rows();
    Eigen::SparseMatrix<double> psi_(n_locs, n_dofs);
    // evaluate basis system at locations
    std::vector<fdapde::Triplet<double>> triplet_list;
    triplet_list.reserve(n_locs * n_shape_functions);

    Eigen::Matrix<int, Dynamic, 1> cell_id = fe_space.triangulation().locate(coords);
    const auto& dof_handler = fe_space.dof_handler();
    // build basis evaluation matrix
    for (int i = 0; i < n_locs; ++i) {
        if (cell_id[i] != -1) {   // point falls inside domain
            Eigen::Matrix<double, embed_dim, 1> p_i(coords.row(i));
            auto cell = dof_handler.cell(cell_id[i]);
            // update matrix
            for (int h = 0; h < n_shape_functions; ++h) {
                triplet_list.emplace_back(
                  i, cell.dofs()[h], fe_space.eval_shape_value(h, cell.invJ() * (p_i - cell.node(0))));   // \psi_j(p_i)
            }
        }
    }
    // finalize construction
    psi_.setFromTriplets(triplet_list.begin(), triplet_list.end());
    psi_.makeCompressed();
    return psi_;
}
template <typename Triangulation_, typename FeType_, typename GeoIndex_>
    requires(!internals::is_eigen_dense_xpr_v<GeoIndex_>)
Eigen::SparseMatrix<double>
point_basis_eval(const FeSpace<Triangulation_, FeType_>& fe_space, const GeoIndex_& geo_index) {
    if (geo_index.points_at_dofs()) {
        int n_dofs = fe_space.n_dofs();
        int n_locs = geo_index.rows();
        Eigen::SparseMatrix<double> psi_(n_locs, n_dofs);
        psi_.setIdentity();   // \psi_i(p_j) = 1 \iff i == j, otherwise \psi_i(p_j) = 0
        return psi_;
    }
    return point_basis_eval(fe_space, geo_index.coordinates());
}

// pointwise basis evaluation for spline basis system
template <typename Triangulation_, typename CoordsMatrix_>
    requires(internals::is_eigen_dense_xpr_v<CoordsMatrix_>)
Eigen::SparseMatrix<double> point_basis_eval(const BsSpace<Triangulation_>& bs_space, CoordsMatrix_&& coords) {
    static constexpr int embed_dim = Triangulation_::embed_dim;
    fdapde_assert(coords.rows() > 0 && coords.cols() == embed_dim);

    int n_shape_functions = bs_space.n_shape_functions();
    int n_dofs = bs_space.n_dofs();
    int n_locs = coords.rows();
    Eigen::SparseMatrix<double> psi_(n_locs, n_dofs);    
    std::vector<Triplet<double>> triplet_list;
    triplet_list.reserve(n_locs * n_shape_functions);

    Eigen::Matrix<int, Dynamic, 1> cell_id = bs_space.triangulation().locate(coords);
    const auto& dof_handler = bs_space.dof_handler();
    // build basis evaluation matrix
    for (int i = 0; i < n_locs; ++i) {
        if (cell_id[i] != -1) {   // point falls inside domain
            Eigen::Matrix<double, embed_dim, 1> p_i(coords.row(i));
            auto cell = dof_handler.cell(cell_id[i]);
            // update matrix
            for (int h = 0; h < cell.dofs().size(); ++h) {
                int active_dof = cell.dofs()[h];
                triplet_list.emplace_back(i, active_dof, bs_space.eval_cell_value(active_dof, p_i));   // \psi_j(p_i)
            }
        }
    }
    // finalize construction
    psi_.setFromTriplets(triplet_list.begin(), triplet_list.end());
    psi_.makeCompressed();
    return psi_;
}
template <typename Triangulation_, typename GeoIndex_>
    requires(!internals::is_eigen_dense_xpr_v<GeoIndex_>)
Eigen::SparseMatrix<double> point_basis_eval(const BsSpace<Triangulation_>& bs_space, const GeoIndex_& geo_index) {
    return point_basis_eval(bs_space, geo_index.coordinates());
}

// areal basis evaluation for finite element basis system
template <typename Triangulation_, typename FeType_>
std::pair<Eigen::SparseMatrix<double>, Eigen::Matrix<double, Dynamic, 1>> areal_basis_eval(
  const FeSpace<Triangulation_, FeType_>& fe_space, const BinaryMatrix<Dynamic, Dynamic>& incidence_mat) {
    using FeSpace_ = FeSpace<Triangulation_, FeType_>;
    fdapde_assert(incidence_mat.rows() > 0 && incidence_mat.cols() == fe_space.triangulation().n_cells());
    static constexpr int local_dim = Triangulation_::local_dim;
    using FeType = typename FeSpace_::FeType;
    using cell_dof_descriptor = typename FeSpace_::cell_dof_descriptor;
    using BasisType = typename cell_dof_descriptor::BasisType;
    using Quadrature = typename FeType::template cell_quadrature_t<local_dim>;
    static constexpr int n_quadrature_nodes = Quadrature::order;
    static constexpr int n_shape_functions = fe_space.n_shape_functions();
    // compile time evaluation of \int_{\hat K} \psi_i on reference element \hat K
    static constexpr Matrix<double, n_shape_functions, 1> int_table_ {[]() {
        std::array<double, n_shape_functions> int_table_ {};
        BasisType basis {cell_dof_descriptor().dofs_phys_coords()};
        for (int i = 0; i < n_shape_functions; ++i) {
            for (int k = 0; k < n_quadrature_nodes; ++k) {
                int_table_[i] += Quadrature::weights[k] * basis[i](Quadrature::nodes.row(k).transpose());
            }
        }
        return int_table_;
    }};

    int n_dofs = fe_space.n_dofs();
    int n_regions = incidence_mat.rows();
    Eigen::SparseMatrix<double> psi_(n_regions, n_dofs);
    Eigen::Matrix<double, Dynamic, 1> D(n_regions);
    std::vector<fdapde::Triplet<double>> triplet_list;
    triplet_list.reserve(n_regions * n_shape_functions);

    const auto& dof_handler = fe_space.dof_handler();
    int tail = 0;
    for (int k = 0; k < n_regions; ++k) {
        int head = 0;
        double Di = 0;   // measure of region D_i
        for (int l = 0, n_cells = incidence_mat.cols(); l < n_cells; ++l) {
            if (incidence_mat(k, l)) {   // element with ID l belongs to k-th region
                auto cell = dof_handler.cell(l);
                for (int h = 0; h < n_shape_functions; ++h) {
                    // compute \int_e \psi_h on physical element e
                    triplet_list.emplace_back(k, cell.dofs()[h], int_table_[h] * cell.measure());
                    head++;
                }
                Di += cell.measure();
            }
        }
        // divide each \int_{D_i} \psi_j by the measure of region D_i
        for (int j = 0; j < head; ++j) { triplet_list[tail + j].value() /= Di; }
        D[k] = Di;
        tail += head;
    }
    // finalize construction
    psi_.setFromTriplets(triplet_list.begin(), triplet_list.end());
    psi_.makeCompressed();
    return std::make_pair(std::move(psi_), std::move(D));
}
template <typename Triangulation_, typename FeType_, typename GeoIndex_>
std::pair<Eigen::SparseMatrix<double>, Eigen::Matrix<double, Dynamic, 1>>
areal_basis_eval(const FeSpace<Triangulation_, FeType_>& fe_space, const GeoIndex_& geo_index) {
    return areal_basis_eval(fe_space, geo_index.incidence_matrix());
}

// areal basis evaluation for spline basis system
template <typename Triangulation_, typename GeoIndex_>
std::pair<Eigen::SparseMatrix<double>, Eigen::Matrix<double, Dynamic, 1>>
areal_basis_eval(const BsSpace<Triangulation_>&, const GeoIndex_&) {
    return {};   // TODO
}

}   // namespace internals
}   // namespace fdapde

#endif // __FE_ELLIPTIC_DRIVER_H__
