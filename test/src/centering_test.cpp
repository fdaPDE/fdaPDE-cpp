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

#include <fdaPDE/core.h>
#include <gtest/gtest.h>   // testing framework

#include <cstddef>
using fdapde::core::fem_order;
using fdapde::core::FEM;
using fdapde::core::Grid;
using fdapde::core::laplacian;
using fdapde::core::PDE;
using fdapde::core::Triangulation;

#include "../../fdaPDE/models/regression/srpde.h"
#include "../../fdaPDE/models/regression/gcv.h"
#include "../../fdaPDE/models/sampling_design.h"
using fdapde::models::SRPDE;
using fdapde::models::SpaceOnly;
using fdapde::models::ExactEDF;
using fdapde::models::StochasticEDF;
using fdapde::models::Sampling;
#include "../../fdaPDE/models/functional/center.h"
using fdapde::models::center;
#include "../../fdaPDE/calibration/gcv.h"

#include "utils/constants.h"
#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MeshLoader;
using fdapde::testing::read_mtx;
using fdapde::testing::read_csv;

// test 1
//    domain:       unit square [1,1] x [1,1]
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid stochastic
TEST(centering_test, srpde_gcv_stochastic_grid) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("unit_square");
    // import data from files
    DMatrix<double> X = read_csv<double>("../data/models/centering/2D_test1/X.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> pde(domain.mesh, L, u);
    // perform centering
    DMatrix<double> lambda_grid(12, 1);
    for (int i = 0; i < 12; ++i) lambda_grid(i, 0) = std::pow(10, -6 + 0.5 * i);
    auto centered_data = center(
      X, SRPDE {pde, Sampling::mesh_nodes},
      fdapde::calibration::GCV<SpaceOnly> {Grid<fdapde::Dynamic> {}, StochasticEDF(100)}(lambda_grid));
    // test correctness
    EXPECT_TRUE(almost_equal(centered_data.fitted, "../data/models/centering/2D_test1/fitted.mtx"));
}
