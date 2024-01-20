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
using fdapde::core::dt;
using fdapde::core::FEM;
using fdapde::core::fem_order;
using fdapde::core::bilaplacian;
using fdapde::core::laplacian;
using fdapde::core::PDE;
using fdapde::core::Mesh;
using fdapde::core::SPLINE;
using fdapde::core::spline_order;

#include "../../fdaPDE/models/sampling_design.h"
#include "../../fdaPDE/models/regression/qsrpde.h"

using fdapde::models::QSRPDE;
using fdapde::models::SpaceTimeSeparable;
using fdapde::models::SpaceTimeParabolic;
using fdapde::models::SpaceOnly;
using fdapde::models::Sampling;

#include "utils/constants.h"
#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MeshLoader;
using fdapde::testing::read_mtx;
using fdapde::testing::read_csv;

// test 1 
//    domain:         c-shaped
//    space sampling: locations != nodes
//    time sampling:  locations != nodes
//    missing data:   no
//    penalization:   simple laplacian
//    covariates:     no
//    BC:             no
//    order FE:       1
//    time penalization: separable (mass penalization)
TEST(qstrpde_test, laplacian_nonparametric_samplingatlocations_separable_monolithic) {
  
  // define temporal domain
  unsigned int M = 7; 
  double tf = fdapde::testing::pi;   
  Mesh<1, 1> time_mesh(0, tf, M-1);     // t0, tf, #subintervals 

  // define spatial domain and regularizing PDE
  MeshLoader<Mesh2D> domain("c_shaped_adj");

  // define regularizing PDE in space 
  auto Ld = -laplacian<FEM>();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.n_nodes(), 1);
  PDE<Mesh<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);

  // define regularizing PDE in time
  auto Lt = -bilaplacian<SPLINE>();
  PDE<Mesh<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);

  // import locs from files
  DMatrix<double> space_locs = read_csv<double>("../data/models/qstrpde/2D_test1/locs.csv");
  DMatrix<double> time_locs = read_csv<double>("../data/models/qstrpde/2D_test1/time_locations.csv");

  // load data from .csv files
  DMatrix<double> y; 
  y = read_csv<double>("../data/models/qstrpde/2D_test1/y.csv");

  // alpha
  double alpha = 0.5;

  // lambdas 
  double lambda_D = 1e-3; 
  double lambda_T = 1e-6; 

  // define model
  QSRPDE<SpaceTimeSeparable> model(space_penalty, time_penalty, Sampling::pointwise, alpha);
  model.set_lambda_D(lambda_D);
  model.set_lambda_T(lambda_T);
  model.set_spatial_locations(space_locs);
  model.set_temporal_locations(time_locs);

  // set model's data
  BlockFrame<double, int> df;
  df.stack(OBSERVATIONS_BLK, y);
  model.set_data(df);

  // solve smoothing problem
  model.init();
  model.solve();

  // test correctness
  EXPECT_TRUE(almost_equal(model.f(), "../data/models/qstrpde/2D_test1/sol.mtx"));

}

