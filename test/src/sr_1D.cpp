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

#include "fdaPDE/io.h"
using namespace fdapde;
using fdapde::test::almost_equal;

// test 1
//    mesh:         unit_interval
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
TEST(sr_1D, test_01) {

    // geometry
    int n_nodes = 101;
    Triangulation<1,1> T = Triangulation<1,1>::UnitInterval(n_nodes);

    // data
    int n_obs = n_nodes;
    GeoFrame data(T);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);
    l1.load_csv<double>("../data/sr/01_1D/y.csv");

    // physics
    BsSpace Bh(T, 3);   // cubic B-splines
    TrialFunction f(Bh);
    TestFunction  v(Bh);
    auto a = integral(T)(dxx(f) * dxx(v));  // curvature penalty
    ZeroField<1> u;
    auto F = integral(T)(u * v);

    // modeling
    SRPDE m("x ~ f", data, bs_ls_elliptic(a, F));
    m.fit(1e-4/n_obs);
    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/01_1D/sol.mtx"));
}

// test 1
//    mesh:         unit_interval
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
TEST(sr_1D, test_02) {

    // geometry
    int n_nodes = 101;
    Triangulation<1,1> T = Triangulation<1,1>::UnitInterval(n_nodes);

    // data
    int n_obs = n_nodes;
    GeoFrame data(T);
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);
    l1.load_csv<double>("../data/sr/01_1D/y.csv");

    // physics
    BsSpace Bh(T, 3); // cubic B-splines
    TrialFunction f(Bh);
    TestFunction  v(Bh);
    auto a = integral(T)(dxx(f) * dxx(v));
    ZeroField<1> u;
    auto F = integral(T)(u * v);

    // modeling
    SRPDE m("x ~ f", data, bs_ls_elliptic(a, F));
    m.set_trace_mode(TraceMode::Exact);

    // calibration
    std::vector<double> lambda_grid;
    for (double x = -9.0; x <= .0; x += 0.5) lambda_grid.push_back(std::pow(10, x) / n_obs);
    GridSearch<1> optimizer;
    optimizer.optimize(m.gcv(100, 476813), lambda_grid);
    EXPECT_TRUE(almost_equal<double>(optimizer.values(), "../data/sr/01_1D/gcv_scores.mtx"));

    // final fit
    m.fit(optimizer.optimum());
    EXPECT_TRUE(almost_equal<double>(m.f(), "../data/sr/01_1D/sol_gcv.mtx"));
}