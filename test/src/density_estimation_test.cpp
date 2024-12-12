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

#include <cstddef>
#include <gtest/gtest.h>   // testing framework

#include <fdaPDE/core.h>
using fdapde::core::FEM;
using fdapde::core::fem_order;
using fdapde::core::laplacian;
using fdapde::core::PDE;
using fdapde::core::Triangulation;
using fdapde::core::GradientDescent;
using fdapde::core::BacktrackingLineSearch;
using fdapde::core::WolfeLineSearch;
using fdapde::core::BFGS;

#include "../../fdaPDE/models/density_estimation/depde.h"
#include "../../fdaPDE/models/sampling_design.h"
using fdapde::models::DEPDE;
using fdapde::models::Sampling;
#include "../../fdaPDE/calibration/kfold_cv.h"
using fdapde::calibration::KCV;

#include "utils/constants.h"
#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MeshLoader;
using fdapde::testing::read_csv;

// test 1
//    domain:       unit square [1,1] x [1,1]
//    optimizer:    BGFS, fixed step
//    lambda:       fixed
//    order FE:     1
TEST(depde_test, fixed_smoothing_bfgs_fixed_step) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("square_density");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<Triangulation<2, 2>, decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double lambda = 0.1;
    DEPDE<SpaceOnly> model(problem);
    model.set_lambda_D(lambda);
    model.set_tolerance(1e-15);   // to let the optimization process stop for maximum number of iterations
    // set model's data
    BlockFrame<double, int> df;
    df.insert(SPACE_LOCS, read_csv<double>("../data/models/depde/2D_test1/data.csv"));
    model.set_data(df);
    model.set_optimizer(BFGS<fdapde::Dynamic> {500, 1e-5, 1e-2});
    model.set_g_init(read_csv<double>("../data/models/depde/2D_test1/f_init.csv").array().log());
    // solve density estimation problem
    model.init();
    model.solve();   // this stops at maximum number of iterations (expected)
    // test correctness
    EXPECT_TRUE(almost_equal(model.g(), "../data/models/depde/2D_test1/g.mtx"));
}

// test 2
//    domain:       unit square [1,1] x [1,1]
//    optimizer:    gradient descent, backtracking adaptive step
//    lambda:       fixed
//    order FE:     1
TEST(depde_test, fixed_smoothing_gradient_descent_backtracking_step) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("square_density");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<Triangulation<2, 2>, decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double lambda = 0.1;
    DEPDE<SpaceOnly> model(problem);
    model.set_lambda_D(lambda);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(SPACE_LOCS, read_csv<double>("../data/models/depde/2D_test2/data.csv"));
    model.set_data(df);
    model.set_optimizer(GradientDescent<fdapde::Dynamic, BacktrackingLineSearch> {1000, 1e-5, 1e-2});
    model.set_g_init(read_csv<double>("../data/models/depde/2D_test2/f_init.csv").array().log());
    // solve density estimation problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.g(), "../data/models/depde/2D_test2/g.mtx"));
}

// test 3
//    domain:       unit square [1,1] x [1,1]
//    optimizer:    BFGS, wolfe adaptive step
//    lambda:       fixed
//    order FE:     1
TEST(depde_test, fixed_smoothing_bfgs_wolfe_step) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("square_density");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<Triangulation<2, 2>, decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    double lambda = 0.1;
    DEPDE<SpaceOnly> model(problem);
    model.set_lambda_D(lambda);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(SPACE_LOCS, read_csv<double>("../data/models/depde/2D_test3/data.csv"));
    model.set_data(df);
    model.set_optimizer(BFGS<fdapde::Dynamic, WolfeLineSearch> {1000, 1e-5, 1e-2});
    model.set_g_init(read_csv<double>("../data/models/depde/2D_test3/f_init.csv").array().log());
    // solve density estimation problem
    model.init();
    model.solve();
    // test correctness
    EXPECT_TRUE(almost_equal(model.g(), "../data/models/depde/2D_test3/g.mtx"));
}

// test 4
//    domain:       unit square [1,1] x [1,1]
//    optimizer:    BFGS, backtracking adaptive step
//    lambda:       K-fold cross validation, 5 folds
//    order FE:     1
TEST(depde_test, kcv_smoothing_bfgs_backtracking_step) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("square_density");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<Triangulation<2, 2>, decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define model
    DVector<double> lambda_vect;
    lambda_vect.resize(4);
    lambda_vect << 0.001, 0.01, 0.1, 1;
    DEPDE<SpaceOnly> model(problem);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(SPACE_LOCS, read_csv<double>("../data/models/depde/2D_test4/data.csv"));
    model.set_data(df);
    model.set_optimizer(BFGS<fdapde::Dynamic, WolfeLineSearch> {1000, 1e-5, 1e-2});
    DMatrix<double> f_init = read_csv<double>("../data/models/depde/2D_test4/f_init.csv");
    // declare cross validation engine
    KCV kcv(5, 10);
    DVector<double> cv_error;
    cv_error.resize(lambda_vect.size());

    for (int i = 0; i < lambda_vect.rows(); ++i) {
        model.set_g_init(f_init.col(i).array().log());
        kcv.fit(model, SVector<1>(lambda_vect[i]), DEPDE<SpaceOnly>::CVScore(model));
        cv_error[i] = kcv.avg_scores()[0];
    }
    // test correctness
    EXPECT_TRUE(almost_equal(cv_error, "../data/models/depde/2D_test4/cv_error.mtx"));
}
