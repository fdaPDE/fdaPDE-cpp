#include <cstddef>
#include <gtest/gtest.h>   // testing framework

#include <fdaPDE/core.h>
using fdapde::core::advection;
using fdapde::core::diffusion;
using fdapde::core::dt;
using fdapde::core::FEM;
using fdapde::core::SPLINE;
using fdapde::core::bilaplacian;
using fdapde::core::laplacian;
using fdapde::core::PDE;
using fdapde::core::Triangulation;
using fdapde::core::spline_order;

#include "../../fdaPDE/models/regression/strpde.h"
#include "../../fdaPDE/models/sampling_design.h"
using fdapde::models::STRPDE;
using fdapde::models::SpaceTimeSeparable;
using fdapde::models::SpaceTimeParabolic;
using fdapde::models::Sampling;

#include "utils/constants.h"
#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MeshLoader;
using fdapde::testing::read_mtx;
using fdapde::testing::read_csv;

#include "../../fdaPDE/models/regression/wald.h"
#include "../../fdaPDE/models/regression/speckman.h"
#include "../../fdaPDE/models/regression/esf.h"

using fdapde::core::fem_order;
using fdapde::core::Newton;
using fdapde::core::laplacian;
using fdapde::core::PDE;
using fdapde::core::advection;
using fdapde::core::diffusion;
using fdapde::core::dt;
using fdapde::core::SPLINE;
using fdapde::core::bilaplacian;
using fdapde::core::spline_order;

#include "../../fdaPDE/models/regression/srpde.h"
#include "../../fdaPDE/models/regression/gcv.h"
#include "../../fdaPDE/models/sampling_design.h"
#include "../../fdaPDE/models/regression/regression_type_erasure.h"
using fdapde::models::SRPDE;
using fdapde::models::ExactEDF;
using fdapde::models::GCV;
using fdapde::models::StochasticEDF;
using fdapde::models::Sampling;
using fdapde::models::RegressionView;

#include "../../fdaPDE/calibration/gcv.h"


// test 2
//    domain:       c-shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
//    time penalization: separable (mass penalization)


TEST(inferencetime_test, Exact24) {
    // define temporal and spatial domain
    Triangulation<1, 1> time_mesh(0, fdapde::testing::pi, 4);
    MeshLoader<Triangulation<2, 2>> domain("c_shaped");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/models/strpde/2D_test2/locs.csv");
    DMatrix<double> y    = read_csv<double>("../data/models/strpde/2D_test2/y.csv");
    DMatrix<double> X    = read_csv<double>("../data/models/strpde/2D_test2/X.csv");
    // define regularizing PDE in space
    auto Ld = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<Triangulation<2, 2>, decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);
    // define regularizing PDE in time
    auto Lt = -bilaplacian<SPLINE>();
    PDE<Triangulation<1, 1>, decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);
    // define model
    double lambda_D = 0.01;
    double lambda_T = 0.01;
    STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::pointwise);
    model.set_lambda_D(lambda_D);
    model.set_lambda_T(lambda_T);
    model.set_spatial_locations(locs);
    // set model's data
    BlockFrame<double, int> df;
    df.stack(OBSERVATIONS_BLK, y);
    df.stack(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();

    // test correctness WALD
    fdapde::models::Wald<STRPDE<SpaceTimeSeparable, fdapde::monolithic>, fdapde::models::exact> inferenceW(model);
    fdapde::models::Speckman<STRPDE<SpaceTimeSeparable, fdapde::monolithic>, fdapde::models::exact> inferenceS(model);
    fdapde::models::ESF<STRPDE<SpaceTimeSeparable, fdapde::monolithic>, fdapde::models::exact> inferenceESF(model);
    int cols = model.beta().size();
    DMatrix<double> C=DMatrix<double>::Identity(cols, cols);    
    inferenceW.setC(C);
    inferenceS.setC(C);
    inferenceESF.setC(C);
    DVector<double> beta0(1);
    beta0(0)=2;
    inferenceW.setBeta0(beta0);
    inferenceS.setBeta0(beta0);
    inferenceESF.setBeta0(beta0);
    inferenceESF.setNflip(10000);

    inferenceESF.setseed(46);
    
    EXPECT_TRUE(almost_equal(inferenceW.p_value(fdapde::models::one_at_the_time)(0), 0.7660934 , 1e-7));
    EXPECT_TRUE(almost_equal(inferenceS.p_value(fdapde::models::one_at_the_time)(0), 0.715712 , 1e-7));
    //EXPECT_TRUE(almost_equal(inferenceW.f_p_value(), 0.715712 , 1e-7));
    EXPECT_TRUE(almost_equal(inferenceESF.p_value(fdapde::models::one_at_the_time)(0), 0.8168 , 1e-7));

}



TEST(inferencetime_test, spacetime25D) {
    // define temporal and spatial domain
    Triangulation<1, 1> time_mesh(0, 4, 4);
    MeshLoader<Triangulation<2, 3>> domain("hub2.5D");
    // import data from files
    DMatrix<double> y    = read_csv<double>("../data/models/strpde/25D_test1/y.csv");
    DMatrix<double> X    = read_csv<double>("../data/models/strpde/25D_test1/X.csv");
    // define regularizing PDE in space
    auto Ld = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(Ld), DMatrix<double>, FEM, fem_order<1>> space_penalty(domain.mesh, Ld, u);
    // define regularizing PDE in time
    auto Lt = -bilaplacian<SPLINE>();
    PDE<decltype(time_mesh), decltype(Lt), DMatrix<double>, SPLINE, spline_order<3>> time_penalty(time_mesh, Lt);
    // define model
    double lambda_D = 0.00001;
    double lambda_T = 0.00001;
    STRPDE<SpaceTimeSeparable, fdapde::monolithic> model(space_penalty, time_penalty, Sampling::mesh_nodes);
    model.set_lambda_D(lambda_D);
    model.set_lambda_T(lambda_T);
    // set model's data
    BlockFrame<double, int> df;
    df.stack(OBSERVATIONS_BLK, y);
    df.stack(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();

    // test correctness 
    fdapde::models::Wald<STRPDE<SpaceTimeSeparable, fdapde::monolithic>, fdapde::models::exact> inferenceW(model);
    fdapde::models::Speckman<STRPDE<SpaceTimeSeparable, fdapde::monolithic>, fdapde::models::exact> inferenceS(model);
    fdapde::models::ESF<STRPDE<SpaceTimeSeparable, fdapde::monolithic>, fdapde::models::exact> inferenceESF(model);

    // set H0
    DVector<double> beta0(1);
    beta0 << 0.45;
    inferenceW.setBeta0(beta0);
    inferenceS.setBeta0(beta0);
    inferenceESF.setBeta0(beta0);

    // set C
    DMatrix<double> C = DMatrix<double>::Identity(1, 1);
    inferenceW.setC(C);
    inferenceS.setC(C);
    inferenceESF.setC(C);

    // set N flips
    inferenceESF.setNflip(10000);
    inferenceESF.setseed(46);
    
    DVector<double> waldpval = inferenceW.p_value(fdapde::models::one_at_the_time);
    DVector<double> speckpval = inferenceS.p_value(fdapde::models::one_at_the_time);
    DVector<double> esfpval = inferenceESF.p_value(fdapde::models::one_at_the_time);

    EXPECT_TRUE(almost_equal(waldpval(0), 9.184956e-08 , 1e-7));
    EXPECT_TRUE(almost_equal(speckpval(0), 0.9566908 , 1e-7));
    //EXPECT_TRUE(almost_equal(esfpval(0), 0.5652 , 1e-3));


}


