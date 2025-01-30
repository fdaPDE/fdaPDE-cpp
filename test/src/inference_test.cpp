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


// QUESTO FILE VA SPOSTATO IN FDAPDE-CORE/TEST/SRC
// NEL MAIN.CPP DELLA REPOSITORY FDAPDE-CORE/TEST VA AGGIUNTO 
// #include inference_test.cpp
// e la parte che chiama la classe inference_test 


// questi sono da controllare 
#include <fdaPDE/core.h>
#include <gtest/gtest.h>   // testing framework
#include <cstddef>
#include <chrono>

#include <cstddef>
#include <gtest/gtest.h>   // testing framework

#include <fdaPDE/core.h>
using fdapde::core::advection;
using fdapde::core::diffusion;
using fdapde::core::FEM;
using fdapde::core::fem_order;
using fdapde::core::laplacian;
using fdapde::core::DiscretizedMatrixField;
using fdapde::core::PDE;
using fdapde::core::DiscretizedVectorField;
using fdapde::core::Triangulation;

#include "../../fdaPDE/models/regression/srpde.h"
#include "../../fdaPDE/models/sampling_design.h"
using fdapde::models::SRPDE;
using fdapde::models::Sampling;

#include "utils/constants.h"
#include "utils/mesh_loader.h"
#include "utils/utils.h"
using fdapde::testing::almost_equal;
using fdapde::testing::MeshLoader;
using fdapde::testing::read_csv;
using fdapde::testing::read_mtx;


#include "../../fdaPDE/models/regression/wald.h"
#include "../../fdaPDE/models/regression/speckman.h"
#include "../../fdaPDE/models/regression/esf.h"
#include "../../fdaPDE/models/regression/pesf.h"

//#include <../../../fdaPDE-core/fdaPDE/core.h>
using fdapde::core::DiscretizedMatrixField;
using fdapde::core::DiscretizedVectorField;


// test 
//    domain:       c-shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1



TEST(inference_test, exact27) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("c_shaped");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/models/srpde/2D_test2/locs.csv");
    DMatrix<double> y    = read_csv<double>("../data/models/srpde/2D_test2/y.csv");
    DMatrix<double> X    = read_csv<double>("../data/models/srpde/2D_test2/X.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define statistical model
    double lambda = 0.2201047;
    SRPDE model(problem, Sampling::pointwise);
    model.set_lambda_D(lambda);
    model.set_spatial_locations(locs);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    df.insert(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();

    fdapde::models::Wald<SRPDE, fdapde::models::exact> inferenceWald(model);
    fdapde::models::Speckman<SRPDE, fdapde::models::exact> inferenceSpeck(model);
    fdapde::models::ESF<SRPDE,fdapde::models::exact> inferenceESF(model);
    fdapde::models::PESF<SRPDE,fdapde::models::exact> inferencePESF(model);

    int cols = model.beta().size();
    DMatrix<double> C=DMatrix<double>::Identity(cols, cols);
    
    inferenceWald.setC(C);
    inferenceSpeck.setC(C);
    inferenceESF.setC(C);
    inferencePESF.setC(C);

    DVector<double> beta0(2);
    beta0(0)=2;
    beta0(1)=-1;
    inferenceWald.setBeta0(beta0);
    inferenceSpeck.setBeta0(beta0);
    inferenceESF.setBeta0(beta0);
    inferencePESF.setBeta0(beta0);

    int n = 100000;
    inferenceESF.setNflip(n);
    inferenceESF.setseed(46);
    inferencePESF.setNflip(n);
    inferencePESF.setseed(46);

    inferenceESF.setNflip(10000);
    DVector<int> loc_indexes(6);
    loc_indexes << 1, 5, 7, 8, 9, 10;
    inferenceWald.setLocationsF(loc_indexes);
    inferenceESF.setLocationsF(loc_indexes);

    DVector<double> pvalueswald = inferenceWald.p_value(fdapde::models::simultaneous);
    DVector<double> pvaluesspeck = inferenceSpeck.p_value(fdapde::models::one_at_the_time);
    DVector<double> pvaluesesf = inferenceESF.p_value(fdapde::models::one_at_the_time);
    //double waldstatisticf = inferenceWald.f_p_value();
    double pvalueesf_f = inferenceESF.f_p_value();

    // test correctness Wald
    EXPECT_TRUE(almost_equal(pvalueswald(0), 0.411991314607044 , 1e-7));
    
    // test correctness Speckman
    EXPECT_TRUE(almost_equal(pvaluesspeck(0), 0.0868023617435293, 1e-7));
    EXPECT_TRUE(almost_equal(pvaluesspeck(1), 0.4810795610695496, 1e-7));

    // test correctness ESF
    EXPECT_TRUE(almost_equal(pvaluesesf(0), 0.159 , 1e-7));
    EXPECT_TRUE(almost_equal(pvaluesesf(1), 0.9022 , 1e-7));

    EXPECT_TRUE(almost_equal(pvalueesf_f, 0.6249 , 1e-7));

}



TEST(inference_test, nonexact27) {
    // define domain
    MeshLoader<Triangulation<2, 2>> domain("c_shaped");
    // import data from files
    DMatrix<double> locs = read_csv<double>("../data/models/srpde/2D_test2/locs.csv");
    DMatrix<double> y    = read_csv<double>("../data/models/srpde/2D_test2/y.csv");
    DMatrix<double> X    = read_csv<double>("../data/models/srpde/2D_test2/X.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    // define statistical model
    double lambda = 0.2201047;
    SRPDE model(problem, Sampling::pointwise);
    model.set_lambda_D(lambda);
    model.set_spatial_locations(locs);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    df.insert(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();

    fdapde::models::Wald<SRPDE, fdapde::models::nonexact> inferenceWald(model);
    fdapde::models::Speckman<SRPDE, fdapde::models::nonexact> inferenceSpeck(model);
 
    int cols = model.beta().size();
    DMatrix<double> C=DMatrix<double>::Identity(cols, cols);
    
    inferenceWald.setC(C);
    inferenceSpeck.setC(C);

    DVector<double> beta0(2);
    beta0(0)=2;
    beta0(1)=-1;
    inferenceWald.setBeta0(beta0);
    inferenceSpeck.setBeta0(beta0);

    DVector<double> pvalueswald = inferenceWald.p_value(fdapde::models::simultaneous);
    DVector<double> pvaluesspeck = inferenceSpeck.p_value(fdapde::models::one_at_the_time);

    // test correctness Wald
    EXPECT_TRUE(almost_equal(pvalueswald(0), 0.2368866 , 1e-6));
    
    // test correctness Speckman
    EXPECT_TRUE(almost_equal(pvaluesspeck(0), 0.0903161, 1e-6));
    EXPECT_TRUE(almost_equal(pvaluesspeck(1), 0.6466615, 1e-7));
}





TEST(inference_test, inference25D){
    MeshLoader<Triangulation<2, 3>> domain("horsehoe2.5D");
    // import data from files
    DMatrix<double> y    = read_csv<double>("../data/models/srpde/25D_test1/y.csv");
    DMatrix<double> X    = read_csv<double>("../data/models/srpde/25D_test1/X.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 3, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    double lambda = 0.1;
    SRPDE model(problem, Sampling::mesh_nodes);
    model.set_lambda_D(lambda);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    df.insert(DESIGN_MATRIX_BLK, X);
    model.set_data(df);
    // solve smoothing problem
    model.init();
    model.solve();    
        
    fdapde::models::Wald<SRPDE, fdapde::models::exact> inferenceWald(model);
    fdapde::models::Speckman<SRPDE, fdapde::models::exact> inferenceSpeck(model);
    fdapde::models::ESF<SRPDE, fdapde::models::exact> inferenceESF(model);

    DVector<double> beta0(1);
    int cols = model.beta().size();
    DMatrix<double> C = DMatrix<double>::Identity(cols, cols);
    beta0(0) = 1;
    inferenceWald.setBeta0(beta0);
    inferenceSpeck.setBeta0(beta0);
    inferenceESF.setBeta0(beta0);

    inferenceWald.setC(C);
    inferenceSpeck.setC(C);
    inferenceESF.setC(C);

    inferenceESF.setNflip(10000);
  
    DVector<double> Wald_beta_p = inferenceWald.p_value(fdapde::models::one_at_the_time);
    DMatrix<double> Wald_beta_CI = inferenceWald.computeCI(fdapde::models::one_at_the_time);
    DVector<double> Speck_beta_p = inferenceSpeck.p_value(fdapde::models::one_at_the_time);
    DMatrix<double> Speck_beta_CI = inferenceSpeck.computeCI(fdapde::models::one_at_the_time);

    EXPECT_TRUE(almost_equal(Wald_beta_p(0), 0.01282658 , 1e-7));
    EXPECT_TRUE(almost_equal(Speck_beta_p(0), 0.02520083 , 1e-7));

}


TEST(inference_test, inference3D){

    MeshLoader<Triangulation<3,3>> domain("unit_sphere3D");
    DMatrix<double> y    = read_csv<double>("../data/models/srpde/3D_test1/y.csv");
    DMatrix<double> X    = read_csv<double>("../data/models/srpde/3D_test1/X.csv");
    // define regularizing PDE
    auto L = -laplacian<FEM>();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_cells() * 4, 1);
    PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
    double lambda = 0.01;
    SRPDE model(problem, Sampling::mesh_nodes);
    model.set_lambda_D(lambda);
    //model.set_spatial_locations(locs);
    // set model's data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, y);
    df.insert(DESIGN_MATRIX_BLK, X);
    model.set_data(df);

    // solve smoothing problem
    model.init();
    model.init_psi_esf(problem);
    model.solve();   
        
    fdapde::models::Wald<SRPDE, fdapde::models::exact> inferenceWald(model);
    fdapde::models::Speckman<SRPDE, fdapde::models::exact> inferenceSpeck(model);
    fdapde::models::ESF<SRPDE, fdapde::models::exact> inferenceESF(model);

    DVector<double> beta0(2);
    int cols = model.beta().size();
    DMatrix<double> C = DMatrix<double>::Identity(cols, cols);
    beta0(0) = 2;
    beta0(1) = -1;
    inferenceWald.setBeta0(beta0);
    inferenceSpeck.setBeta0(beta0);
    inferenceESF.setBeta0(beta0);

    inferenceWald.setC(C);
    inferenceSpeck.setC(C);
    inferenceESF.setC(C);

    inferenceESF.setNflip(10000);
    inferenceESF.setseed(46);

    DVector<int> locs_ind(10);
    locs_ind << 1, 6, 8, 11, 16, 18, 20, 21, 23, 24;
    inferenceWald.setLocationsF(locs_ind);
    inferenceESF.setMesh_loc(locs_ind);

    DVector<double> Wald_beta_p = inferenceWald.p_value(fdapde::models::simultaneous);
    //double wald_f_p_val = inferenceWald.f_p_value();
    DVector<double> Speck_beta_p = inferenceSpeck.p_value(fdapde::models::one_at_the_time); 

    double pvalueesf_f = inferenceESF.f_p_value();
    double pvalueesf_sf = inferenceESF.sign_flip_p_value();
    
    EXPECT_TRUE(almost_equal(Wald_beta_p(0), 0.9684002 , 1e-7));
    EXPECT_TRUE(almost_equal(Speck_beta_p(0), 0.6479218 , 1e-7));
    EXPECT_TRUE(almost_equal(Speck_beta_p(1), 0.4182482 , 1e-7));
    //EXPECT_TRUE(almost_equal(wald_f_p_val, 0.00003363157 , 1e-7));
    EXPECT_TRUE(almost_equal(pvalueesf_f, 0.3525 , 1e-7));
    EXPECT_TRUE(almost_equal(pvalueesf_sf, 0.5488 , 1e-7));

}


