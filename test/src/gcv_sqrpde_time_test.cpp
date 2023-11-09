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
using fdapde::core::laplacian;
using fdapde::core::MatrixDataWrapper;
using fdapde::core::PDE;
using fdapde::core::VectorDataWrapper;

#include "../../fdaPDE/models/regression/sqrpde.h"
#include "../../fdaPDE/models/sampling_design.h"
using fdapde::models::SQRPDE;
using fdapde::models::SpaceTimeSeparable;
using fdapde::models::SpaceTimeParabolic;
using fdapde::models::GeoStatMeshNodes;
using fdapde::models::GeoStatLocations;
using fdapde::models::Areal;
using fdapde::models::MonolithicSolver;
using fdapde::models::IterativeSolver;

#include "../../fdaPDE/calibration/gcv.h"
using fdapde::calibration::ExactEDF;
using fdapde::calibration::ExactGCV;
using fdapde::calibration::GCV;
using fdapde::calibration::StochasticEDF;

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




// /* test 1 SQRPDE - Time
//    domain:       unit square [0,1] x [0,1] (coarse)
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
//  */
// TEST(gcv_sqrpde_time_test, laplacian_nonparametric_samplingatnodes_spacetimeseparable_gridexact) {
  
  
//   // Parameters 
//   const std::string TestNumber = "1"; 
  
//   std::vector<double> alphas = {0.1};   // , 0.5, 0.9}; 

//   // number of simulations
//   unsigned int n_sim = 11; 

//   // Marco
//   // std::string path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/PACS_project_shared"; 
//   // Ilenia 
//   std::string path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared"; 

//   std::string path_test = path + "/space_time/Test_" + TestNumber ;


//   // define time domain
//   DVector<double> time_mesh;
//   time_mesh.resize(11);
//   std::size_t i = 0;
//   for(double x = 0; x <= 2; x+=0.2, ++i) time_mesh[i] = x;
  
//   // define spatial domain 
//   MeshLoader<Mesh2D> domain("unit_square_coarse");
//   // define regularizing PDE    
//   auto L = -laplacian<FEM>();
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.rows(), 1);
//   PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
//   for(double alpha : alphas){

//     unsigned int alpha_int = alpha*100; 

//     std::cout << "------------------------------------------alpha=" << std::to_string(alpha_int) << "%-------------------------------------------------" << std::endl; 

  
//     // define statistical model
//     SQRPDE<decltype(problem), SpaceTimeSeparable, GeoStatMeshNodes, MonolithicSolver> model(problem, time_mesh, alpha);


//     for(unsigned int sim = 11; sim <= n_sim; ++sim){

//       // import data from files
//       DMatrix<double> y = read_csv<double>(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/y.csv");

//       // set model data
//       BlockFrame<double, int> df;
//       df.stack(OBSERVATIONS_BLK, y);
//       model.set_data(df);
//       model.init(); 

//       // define GCV function and grid of \lambda_D values
//       GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
//       ScalarField<2, decltype(GCV)> obj(GCV);
//       std::vector<SVector<2>> lambdas_d_t;
//       // for (double x = -7.0; x <= -3.6; x +=0.20) lambdas_d_t.push_back(SVector<2>(std::pow(10, x), std::pow(10,x)));
//       for (double x = -7.0; x <= -6.9; x +=0.20) lambdas_d_t.push_back(SVector<2>(std::pow(10, x), std::pow(10,x)));
//       // optimize GCV
//       Grid<2> opt;
//       opt.optimize(obj, lambdas_d_t);
//       SVector<2> best_lambda = opt.optimum();
      
//       std::cout << "Best lambda = " << best_lambda << std::endl ; 

//       // Save Lambda opt
//       std::ofstream fileLambdaoptS(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/GCV/Exact/LambdaS.csv");
//       if (fileLambdaoptS.is_open()){
//       fileLambdaoptS << std::setprecision(16) << best_lambda[0];
//       fileLambdaoptS.close();
//       }

//       std::ofstream fileLambdaoptT(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/GCV/Exact/LambdaT.csv");
//       if (fileLambdaoptT.is_open()){
//       fileLambdaoptT << std::setprecision(16) << best_lambda[1];
//       fileLambdaoptT.close();
//       }
           
      
//       // Save GCV scores
//       std::ofstream fileGCV_scores(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/GCV/Exact/GCV_scores.csv");
//       for(std::size_t i = 0; i < GCV.values().size(); ++i) 
//         fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.values()[i]) << "\n" ; 

//       fileGCV_scores.close(); 

//     }

//   }


// }



// /* test 1 SQRPDE - Time
//    domain:       unit square [0,1] x [0,1] (coarse)
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid stoch
//  */
// TEST(gcv_sqrpde_time_test, laplacian_nonparametric_samplingatnodes_spacetimeseparable_gridstochastic) {
  
  
//   // Parameters 
//   const std::string TestNumber = "1"; 

//   double alpha = 0.90;
//   unsigned int alpha_int = alpha*100; 
//   const std::string alpha_string = std::to_string(alpha_int);
  
//   // Marco
//   // std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/PACS_project_shared"; 
//   // Ilenia 
//   std::string path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared"; 

//   std::string path_test = path + "/space_time/Test_" + TestNumber ;


//   // define time domain
//   DVector<double> time_mesh;
//   time_mesh.resize(11);
//   std::size_t i = 0;
//   for(double x = 0; x <= 2; x+=0.2, ++i) time_mesh[i] = x;
  
//   // define spatial domain and regularizing PDE
//   MeshLoader<Mesh2D> domain("unit_square_coarse");
//   // define regularizing PDE    
//   auto L = -laplacian<FEM>();
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3 * time_mesh.rows(), 1);
//   PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);
  
//   // define statistical model
//   SQRPDE<decltype(problem), SpaceTimeSeparable, GeoStatMeshNodes, MonolithicSolver> model(problem, time_mesh, alpha);


//   // import data from files
//   DMatrix<double> y = read_csv<double>(path_test + "/alpha_" + std::to_string(alpha_int) + "/y.csv");

//   // set model data
//   BlockFrame<double, int> df;
//   df.stack(OBSERVATIONS_BLK, y);
//   model.set_data(df);
//   model.init(); 
  

//   // define GCV function and grid of \lambda_D values
//   std::size_t seed = 438172;
//   GCV<decltype(model), StochasticEDF<decltype(model)>> GCV(model, 1000, seed);
    
//   ScalarField<2, decltype(GCV)> obj(GCV);
//   std::vector<SVector<2>> lambdas_d_t;
//   for (double x = -7.0; x <= -3.6; x +=0.20) lambdas_d_t.push_back(SVector<2>(std::pow(10, x), std::pow(10,x)));
//   // optimize GCV
//   Grid<2> opt;
//   opt.optimize(obj, lambdas_d_t);
//   SVector<2> best_lambda = opt.optimum();
  
//   std::cout << "Best lambda = " << best_lambda << std::endl ; 

//   // Save Lambda opt
//   std::ofstream fileLambdaopt(path_test + "/alpha_" + alpha_string + "/GCV/Stoch/LambdaCpp.csv");
//   for(std::size_t i = 0; i < best_lambda.size(); ++i) 
//     fileLambdaopt << std::setprecision(16) << best_lambda[i] << "\n" ; 
    
//   fileLambdaopt.close();
  
  
//   // Save GCV scores
//   std::ofstream fileGCV_scores(path_test + "/alpha_" + alpha_string + "/GCV/Stoch/GCV_scores.csv");
//   for(std::size_t i = 0; i < GCV.values().size(); ++i) 
//     fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.values()[i]) << "\n" ; 

//   fileGCV_scores.close(); 


// }


/* test 2 
   domain:       c-shaped
   sampling:     locations != nodes
   penalization: simple laplacian
   covariates:   yes
   BC:           no
   order FE:     1
   time penalization: separable (mass penalization)
 */
TEST(gcv_sqrpde_time_test, laplacian_semiparametric_samplingatlocations_separablemonolithic_gridexact) {

  // Parameters 
  const std::string TestNumber = "2"; 
  
  // Marco
  // std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/PACS_project_shared"; 
  // Ilenia 
  std::string path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared"; 

  std::string path_test = path + "/space_time/Test_" + TestNumber ;

  std::vector<double> alphas = {0.1};   //  0.5,0.9}; 

  // number of simulations
  unsigned int n_sim = 11; 

  // define time domain
  DVector<double> time_mesh;
  time_mesh.resize(5);
  for(std::size_t i = 0; i < 5; ++i)
    time_mesh[i] = (fdapde::testing::pi/4)*i;

  // define domain and regularizing PDE
  MeshLoader<Mesh2D> domain("c_shaped");
  // define regularizing PDE
  auto L = -laplacian<FEM>();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
  PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

  // import data from files
  DMatrix<double> locs = read_csv<double>(path_test + "/locs.csv");
  DMatrix<double> X = read_csv<double>(path_test + "/X.csv");
  

  for(double alpha : alphas){

    unsigned int alpha_int = alpha*100; 

    std::cout << "------------------------------------------alpha=" << std::to_string(alpha_int) << "%-------------------------------------------------" << std::endl; 

    // Define model 
    SQRPDE<decltype(problem), SpaceTimeSeparable, GeoStatLocations, MonolithicSolver> model(problem, time_mesh, alpha);

    for(unsigned int sim = 11; sim <= n_sim; ++sim){

      // load data from .csv files
      DMatrix<double> y = read_csv<double>(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/y.csv");
    
      // set model's data
      model.set_spatial_locations(locs);
      
      BlockFrame<double, int> df;
      df.stack(OBSERVATIONS_BLK, y);
      df.stack(DESIGN_MATRIX_BLK, X);
      model.set_data(df);
      model.init();

      // define GCV function and grid of \lambda_D values
      GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
      ScalarField<2, decltype(GCV)> obj(GCV);
      std::vector<SVector<2>> lambdas_d_t;
      // for(double x = -3.8; x <= -2.6; x +=0.1) lambdas_d_t.push_back(SVector<2>(std::pow(10,x), std::pow(10,x)));   // 50%
      // for(double x = -4.8; x <= -3.0; x +=0.1) lambdas_d_t.push_back(SVector<2>(std::pow(10,x), std::pow(10,x)));  // 10% and 90%
      for(double x = -4.8; x <= -4.75; x +=0.1) lambdas_d_t.push_back(SVector<2>(std::pow(10,x), std::pow(10,x)));  // 10% and 90%
      // optimize GCV
      Grid<2> opt;
      opt.optimize(obj, lambdas_d_t);
      SVector<2> best_lambda = opt.optimum();
      
      std::cout << "Best lambda = " << best_lambda << std::endl ; 

      // Save Lambda opt
      std::ofstream fileLambdaoptS(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/GCV/Exact/LambdaS.csv");
      if (fileLambdaoptS.is_open()){
      fileLambdaoptS << std::setprecision(16) << best_lambda[0];
      fileLambdaoptS.close();
      }

      std::ofstream fileLambdaoptT(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/GCV/Exact/LambdaT.csv");
      if (fileLambdaoptT.is_open()){
      fileLambdaoptT << std::setprecision(16) << best_lambda[1];
      fileLambdaoptT.close();
      }


      // Save GCV scores
      std::ofstream fileGCV_scores(path_test + "/alpha_" + std::to_string(alpha_int) + "/sim_" + std::to_string(sim) + "/GCV/Exact/GCV_scores.csv");
      for(std::size_t i = 0; i < GCV.values().size(); ++i) 
        fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.values()[i]) << "\n" ; 

      fileGCV_scores.close(); 

    }
  }
  
}




// /* test 2 
//    domain:       c-shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
//    time penalization: separable (mass penalization)
//  */
// TEST(gcv_sqrpde_time_test, laplacian_semiparametric_samplingatlocations_gridstochastic) {

//   // Parameters 
//   const std::string TestNumber = "2"; 

//   double alpha = 0.50;
//   unsigned int alpha_int = alpha*100; 
//   const std::string alpha_string = std::to_string(alpha_int);
  
//   // Marco
//   // std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/PACS_project_shared"; 
//   // Ilenia 
//   std::string path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/Thesis_shared"; 

//   std::string path_test = path + "/space_time/Test_" + TestNumber ;


//   // define time domain
//   DVector<double> time_mesh;
//   time_mesh.resize(5);
//   for(std::size_t i = 0; i < 5; ++i)
//     time_mesh[i] = (fdapde::testing::pi/4)*i;

//   // define domain and regularizing PDE
//   MeshLoader<Mesh2D> domain("c_shaped");
//   // define regularizing PDE
//   auto L = -laplacian<FEM>();
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.n_elements() * 3, 1);
//   PDE<decltype(domain.mesh), decltype(L), DMatrix<double>, FEM, fem_order<1>> problem(domain.mesh, L, u);

//   // import data from files
//   DMatrix<double> locs = read_csv<double>(path_test + "/locs.csv");
//   DMatrix<double> X = read_csv<double>(path_test + "/X.csv");
//   DMatrix<double> y = read_csv<double>(path_test + "/y.csv");

//   // Define model 
//   SQRPDE<decltype(problem), SpaceTimeSeparable, GeoStatLocations, MonolithicSolver> model(problem, time_mesh, alpha);

//   // set model's data
//   model.set_spatial_locations(locs);
//   BlockFrame<double, int> df;
//   df.insert(OBSERVATIONS_BLK, y);
//   df.insert(DESIGN_MATRIX_BLK, X);
//   model.set_data(df);
//   model.init();

//   // define GCV function and grid of \lambda_D values
//   std::size_t seed = 66546513;
//   GCV<decltype(model), StochasticEDF<decltype(model)>> GCV(model, 1000, seed);  
//   ScalarField<2, decltype(GCV)> obj(GCV);
//   std::vector<SVector<2>> lambdas_d_t;
//   // for(double x = -3.8; x <= -2.6; x +=0.1) lambdas_d_t.push_back(SVector<2>(std::pow(10,x), std::pow(10,x)));   // 50%
//   for(double x = -4.8; x <= -3.0; x +=0.1) lambdas_d_t.push_back(SVector<2>(std::pow(10,x), std::pow(10,x)));  // 10% and 90%
//   // optimize GCV
//   Grid<2> opt;
//   opt.optimize(obj, lambdas_d_t);
//   SVector<2> best_lambda = opt.optimum();
  
//   std::cout << "Best lambda = " << best_lambda << std::endl ; 

//   // Save Lambda opt
//   std::ofstream fileLambdaopt(path_test + "/alpha_" + alpha_string + "/GCV/Stoch/LambdaCpp.csv");
//   for(std::size_t i = 0; i < best_lambda.size(); ++i) 
//     fileLambdaopt << std::setprecision(16) << best_lambda[i] << "\n" ; 
    
//   fileLambdaopt.close();
  
  
//   // Save GCV scores
//   std::ofstream fileGCV_scores(path_test + "/alpha_" + alpha_string + "/GCV/Stoch/GCV_scores.csv");
//   for(std::size_t i = 0; i < GCV.values().size(); ++i) 
//     fileGCV_scores << std::setprecision(16) << std::sqrt(GCV.values()[i]) << "\n" ; 

//   fileGCV_scores.close(); 

// }
  
