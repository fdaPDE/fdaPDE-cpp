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
using fdapde::core::advection;
using fdapde::core::diffusion;
using fdapde::core::fem_order;
using fdapde::core::FEM;
using fdapde::core::Grid;
using fdapde::core::Mesh; 
using fdapde::core::laplacian;
using fdapde::core::bilaplacian;
using fdapde::core::SPLINE;
using fdapde::core::spline_order;
using fdapde::core::PDE;

#include "../../fdaPDE/models/sampling_design.h"
#include "../../fdaPDE/models/regression/qsrpde.h"
using fdapde::models::QSRPDE;
using fdapde::models::SpaceTimeSeparable;
using fdapde::models::SpaceTimeParabolic;

#include "../../fdaPDE/models/regression/gcv.h"
using fdapde::models::ExactEDF;
using fdapde::models::GCV;
using fdapde::models::StochasticEDF;
using fdapde::models::Sampling;

#include "utils/constants.h"
#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MeshLoader;
using fdapde::testing::read_mtx;
using fdapde::testing::read_csv;


// for time and memory performances
#include <chrono>
#include <iomanip>
using namespace std::chrono;
#include <unistd.h>
#include <fstream>



// test 1
//  domain:       c-shaped
//  space sampling: locations != nodes
//  time sampling:  locations != nodes
//  penalization: simple laplacian
//  missing:      yes
//  covariates:   no
//  BC:           no
//  order FE:     1
//  GCV optimization: grid exact
//  time penalization: separable (mass penalization)
TEST(gcv_qstrpde_test, laplacian_nonparametric_samplingatlocations_gridexact) {
 
    // define domain
    MeshLoader<Mesh2D> domain("c_shaped_adj");
    unsigned int M = 7;  
    Mesh<1, 1> time_mesh(0, fdapde::testing::pi, M-1);     

    // import data from files
    DMatrix<double> space_locs = read_csv<double>("../data/models/gcv/2D_test9/locs.csv");
    DMatrix<double> time_locs = read_csv<double>("../data/models/gcv/2D_test9/time_locations.csv"); 
    DMatrix<double> y = read_csv<double>("../data/models/gcv/2D_test9/y.csv");

    // define regularizing PDE in space 
    auto Ld = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
    PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

    // define regularizing PDE in time
    auto Lt = -bilaplacian<SPLINE>();
    PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

    // define model
    double alpha = 0.5;
    QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);
    model.set_spatial_locations(space_locs);
    model.set_temporal_locations(time_locs);
    // set model's data
    BlockFrame<double, int> df;
    df.stack(OBSERVATIONS_BLK, y);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    auto GCV = model.gcv<ExactEDF>();
    std::vector<DVector<double>> lambdas_d_t;
    for(double xs = -4.0; xs <= -2.0; xs +=0.5)
      for(double xt = -7.0; xt <= -5.0; xt +=1.0) 
        lambdas_d_t.push_back(SVector<2>(std::pow(10,xs), std::pow(10,xt)));

    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambdas_d_t);
    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/models/gcv/2D_test9/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/models/gcv/2D_test9/gcvs.mtx"));
  
}

// test 2
//  domain:       c-shaped
//  space sampling: locations != nodes
//  time sampling:  locations != nodes
//  penalization: simple laplacian
//  missing:      yes
//  covariates:   no
//  BC:           no
//  order FE:     1
//  GCV optimization: grid stochastic
//  time penalization: separable (mass penalization)
TEST(gcv_qstrpde_test, laplacian_nonparametric_samplingatlocations_gridstochastic) {

    // define domain
    MeshLoader<Mesh2D> domain("c_shaped_adj");
    unsigned int M = 7;  
    Mesh<1, 1> time_mesh(0, fdapde::testing::pi, M-1);
    // import data from files
    DMatrix<double> space_locs = read_csv<double>("../data/models/gcv/2D_test10/locs.csv");
    DMatrix<double> time_locs = read_csv<double>("../data/models/gcv/2D_test10/time_locations.csv"); 
    DMatrix<double> y = read_csv<double>("../data/models/gcv/2D_test10/y.csv");

    // define regularizing PDE in space 
    auto Ld = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
    PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

    // define regularizing PDE in time
    auto Lt = -bilaplacian<SPLINE>();
    PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

    // define model
    double alpha = 0.5;
    QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);
    model.set_spatial_locations(space_locs);
    model.set_temporal_locations(time_locs);
    // set model's data
    BlockFrame<double, int> df;
    df.stack(OBSERVATIONS_BLK, y);
    model.set_data(df);
    model.init();
    // define GCV function and grid of \lambda_D values
    std::size_t seed = 66546513;
    auto GCV = model.gcv<StochasticEDF>(1000, seed);
    std::vector<DVector<double>> lambdas_d_t;
    for(double xs = -4.0; xs <= -2.0; xs +=0.5)
      for(double xt = -7.0; xt <= -5.0; xt +=1.0)
        lambdas_d_t.push_back(SVector<2>(std::pow(10,xs), std::pow(10,xt)));

    // optimize GCV
    Grid<fdapde::Dynamic> opt;
    opt.optimize(GCV, lambdas_d_t);

    // test correctness
    EXPECT_TRUE(almost_equal(GCV.edfs(), "../data/models/gcv/2D_test10/edfs.mtx"));
    EXPECT_TRUE(almost_equal(GCV.gcvs(), "../data/models/gcv/2D_test10/gcvs.mtx"));

}