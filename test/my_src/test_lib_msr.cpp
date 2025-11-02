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


#include <fdaPDE/fdapde.h>

#include <iostream>
#include <string>
#include <vector>

using namespace fdapde;

int test_02();
int test_06_fake();
int test_06_fake_missing();
int test_norm_gcv(); 

int main(){
    // test_02();
    // test_06_fake(); 
    // test_06_fake_missing();
    // test_norm_gcv();
    return 0; 
}

// test 2
//    mesh:         c-shaped
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
int test_02() {

    // geometry 
    using PointT = Eigen::Matrix<double, 2, 1>;

    Eigen::Matrix<double, Dynamic, Dynamic> points = read_csv<double>("my_data/mesh/c_shaped_242/points.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> elements = read_csv<int>("my_data/mesh/c_shaped_242/elements.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> boundary = read_csv<int>("my_data/mesh/c_shaped_242/boundary.csv").as_matrix();

    elements.array() -= 1; // non necessario

    std::string datadir = "my_data/msr/02/";
    Triangulation<2, 2> D(points, elements, boundary);

    GeoFrame data(D);

    // data 
    auto& l = data.insert_scalar_layer<POINT>("layer", MESH_NODES);
    l.load_csv<double>(datadir + "response.csv");
    l.load_csv<double>(datadir + "design_matrix.csv");
    // l.load_csv<double>(datadir + "ids_groups.csv");   // compile error se metto unsigned int...

    //std::cout << data[0] << std::endl;

    // physics 
    FeSpace Vh(D, P1<1>);   // functional space definition

    // trial and test function definition
    TrialFunction f(Vh);
    TestFunction v(Vh);

    auto a = integral(D)(dot(grad(f), grad(v)));

    // homogeneous forcing linear form
    ZeroField<2> u;
    auto F = integral(D)(u * v);

    // modeling
    // MSRPDE m("y ~ x1 + x2 + 1|g + f", data, fe_ls_elliptic(a, F)); 
    QSRPDE m("y ~ x1 + x2 + f", data, 0.25, fe_ls_elliptic(a, F)); 

    // calibration
    std::cout << "setting lambdas..." << std::endl;
    Eigen::Matrix<double, Dynamic, 1> lambda_grid(10); 
    unsigned int count = 0; 
    for(double x = -3.0; x <= +2.0; x += 0.55555555){
        lambda_grid(count) = std::pow(10, x);
        count++;
    } 
    std::cout << "end set lambdas..." << std::endl;
    GridSearch<1> opt;
    std::cout << "Lambda grid: " << lambda_grid << std::endl;
    opt.optimize(m.gcv(1000, 1234), lambda_grid);
    std::cout << "Optimal lambda: " << std::setprecision(16) << opt.optimum() << std::endl;
    double sum_abs_of_elems = 0.;
    std::vector<double> scores = opt.values(); 
    for(std::vector<double>::iterator it = scores.begin(); it != scores.end(); ++it)
        sum_abs_of_elems += std::abs(*it);
    std::cout << "sum abs scores: " << std::setprecision(16) << sum_abs_of_elems << std::endl;

    // // fit at optimal smoothing level
    // m.fit(opt.optimum());

    // double lambda = 2.154434524672918; // read_csv<double>(datadir + "lambda.csv", false, false).as_matrix()(0,0);
    // std::cout << "lambda = " << lambda << std::endl;
    // m.fit(lambda);
    
    write_csv(datadir + "f.csv", m.f());
    write_csv(datadir + "beta.csv", m.beta());

    return 0;
}


// test 6 FAKE ---> SOLO PER TESTARE LIBRERIA NUOVA!!
//    mesh:         unit square
//    sampling:     locations != nodes
//    penalization: anisotropic diffusion
//    time penalization: separable
//    covariates:   yes
//    BC:           no
//    order FE:     1
int test_06_fake() {

    // geometry 
    using PointT = Eigen::Matrix<double, 2, 1>;

    Eigen::Matrix<double, Dynamic, Dynamic> points = read_csv<double>("my_data/mesh/unit_square_reduced_censoring_476/points.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> elements = read_csv<int>("my_data/mesh/unit_square_reduced_censoring_476/elements.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> boundary = read_csv<int>("my_data/mesh/unit_square_reduced_censoring_476/boundary.csv").as_matrix();

    elements.array() -= 1; // non necessario

    std::string datadir = "my_data/msr/06/";
    const unsigned int M = 8; 
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 1, M);  // ATT qui non bisogna fare più M-1 come nella vecchia lib!! Vuole direttamente il numero di nodi, cioè M!!! 
    std::cout << "ATT check M" << std::endl;
    Triangulation<2, 2> D(points, elements, boundary);

    GeoFrame data(D, T);

    const unsigned int max_fpirls_iter = 15;

    // data 
    auto& l = data.insert_scalar_layer<POINT, POINT>("layer", std::pair{"my_data/msr/06/space_locs.csv", "my_data/msr/06/time_locs.csv"});
    l.load_csv<double>(datadir + "response.csv");
    l.load_csv<double>(datadir + "design_matrix.csv");
    l.load_csv<double>(datadir + "ids_groups.csv");

    std::cout << "geoframe:" << data[0] << std::endl;

    // physics 
    FeSpace Vh(D, P1<1>);   // functional space definition

    // trial and test function definition
    Eigen::Matrix<double, 2, 2> K;
    K << 1, 0, 0, 1;   // sto guardando sim 1, che aveva stimato isotropia
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a_D = integral(D)(dot(K * grad(f), grad(v)));
    // homogeneous forcing linear form
    ZeroField<2> u_D;
    auto F_D = integral(D)(u_D * v);

    BsSpace Bh(T, 3);   // cubic B-splines in time
    TrialFunction g(Bh);
    TestFunction  w(Bh);
    auto a_T = integral(T)(dxx(g) * dxx(w));
    ZeroField<1> u_T;
    auto F_T = integral(T)(u_T * w);    

    // modeling
    MSRPDE m("y ~ x1 + x2 + 1|g + f", data, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  

    m.set_fpirls_max_iter(max_fpirls_iter);

    // // calibration
    // std::cout << "setting lambdas..." << std::endl;
    // Eigen::Matrix<double, Dynamic, 1> lambda_grid(10); 
    // unsigned int count = 0; 
    // for(double x = -3.0; x <= +2.0; x += 0.55555555){
    //     lambda_grid(count) = std::pow(10, x);
    //     count++;
    // } 
    // std::cout << "end set lambdas..." << std::endl;
    // GridSearch<2> opt;
    // std::cout << "Lambda grid: " << lambda_grid << std::endl;
    // opt.optimize(m.gcv(1000, 1234), lambda_grid);
    // std::cout << "Optimal lambda: " << std::setprecision(16) << opt.optimum() << std::endl;
    // double sum_abs_of_elems = 0.;
    // std::vector<double> scores = opt.values(); 
    // for(std::vector<double>::iterator it = scores.begin(); it != scores.end(); ++it)
    //     sum_abs_of_elems += std::abs(*it);
    // std::cout << "sum abs scores: " << std::setprecision(16) << sum_abs_of_elems << std::endl;

    // // // fit at optimal smoothing level
    // // m.fit(opt.optimum());

    double lambdaD = 0.000562341325190349; 
    double lambdaT = 0.0001; // read_csv<double>(datadir + "lambda.csv", false, false).as_matrix()(0,0);
    std::cout << "lambdaD = " << lambdaD << std::endl;
    std::cout << "lambdaT = " << lambdaT << std::endl;
    m.fit(lambdaD, lambdaT);
    
    // write_csv(datadir + "f.csv", m.f());
    // write_csv(datadir + "beta.csv", m.beta());

    return 0;
}



// test 6 FAKE con missing ---> SOLO PER TESTARE LIBRERIA NUOVA!!
//    mesh:         unit square
//    sampling:     locations != nodes
//    penalization: anisotropic diffusion
//    time penalization: separable
//    covariates:   yes
//    BC:           no
//    order FE:     1
int test_06_fake_missing() {

    // geometry 
    using PointT = Eigen::Matrix<double, 2, 1>;

    Eigen::Matrix<double, Dynamic, Dynamic> points = read_csv<double>("my_data/mesh/unit_square_reduced_censoring_476/points.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> elements = read_csv<int>("my_data/mesh/unit_square_reduced_censoring_476/elements.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> boundary = read_csv<int>("my_data/mesh/unit_square_reduced_censoring_476/boundary.csv").as_matrix();

    elements.array() -= 1; // non necessario

    std::string datadir = "my_data/msr/06-missing/";
    const unsigned int M = 8; 
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 1, M);  // ATT qui non bisogna fare più M-1 come nella vecchia lib!! Vuole direttamente il numero di nodi, cioè M!!! 
    std::cout << "ATT check M" << std::endl;
    Triangulation<2, 2> D(points, elements, boundary);

    GeoFrame data(D, T);

    const unsigned int max_fpirls_iter = 15;

    // data 
    auto& l = data.insert_scalar_layer<POINT, POINT>("layer", std::pair{"my_data/msr/06/space_locs.csv", "my_data/msr/06/time_locs.csv"});
    l.load_csv<double>(datadir + "response.csv");
    // std::cout << "------ATT: NO missing loading..." << std::endl;
    // l.load_csv<double>(datadir + "response_full.csv");

    l.load_csv<double>(datadir + "design_matrix.csv");
    l.load_csv<double>(datadir + "ids_groups.csv");

    std::cout << "geoframe:" << data[0] << std::endl;

    // physics 
    FeSpace Vh(D, P1<1>);   // functional space definition

    // trial and test function definition
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a_D = integral(D)(dot(grad(f), grad(v)));
    // homogeneous forcing linear form
    ZeroField<2> u_D;
    auto F_D = integral(D)(u_D * v);

    BsSpace Bh(T, 3);   // cubic B-splines in time
    TrialFunction g(Bh);
    TestFunction  w(Bh);
    auto a_T = integral(T)(dxx(g) * dxx(w));
    ZeroField<1> u_T;
    auto F_T = integral(T)(u_T * w);    

    // modeling
    MSRPDE m("y ~ x1 + x2 + 1|g + f", data, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  

    m.set_fpirls_max_iter(max_fpirls_iter);

    // calibration
    std::vector<double> lambdas_d;
    std::vector<double> lambdas_t;
    std::vector<Eigen::Matrix<double, Dynamic, 1>> lambdas_d_t;
    for(double xs = -4.5; xs <= -3.0; xs += 0.5)
    lambdas_d.push_back(std::pow(10,xs));

    for(double xt = -4.0; xt <= -2.0; xt += 2.0)
        lambdas_t.push_back(std::pow(10,xt));

    for(auto i = 0; i < lambdas_d.size(); ++i)
        for(auto j = 0; j < lambdas_t.size(); ++j) 
            lambdas_d_t.push_back(Eigen::Matrix<double, 2, 1>(lambdas_d[i], lambdas_t[j]));
    std::cout << "end set lambdas..." << std::endl;
    
    Eigen::Matrix<double, Dynamic, 2> lambda_grid(lambdas_d.size()*lambdas_t.size(), 2);
    for(int i = 0; i < lambdas_d.size(); ++i) { 
        for (int j = 0; j < lambdas_t.size(); ++j) {
            lambda_grid(i * lambdas_t.size() + j, 0) = lambdas_d[i];
            lambda_grid(i * lambdas_t.size() + j, 1) = lambdas_t[j];
        }
    }
    std::cout << "lambda_grid.rows() = " << lambda_grid.rows() << std::endl;
    std::cout << "lambda_grid.cols() = " << lambda_grid.cols() << std::endl; 

    GridSearch<2> opt;
    std::cout << "Lambda grid: " << lambda_grid << std::endl;
    opt.optimize(m.gcv(1000, 1234), lambda_grid);
    std::cout << "Optimal lambda: " << std::setprecision(16) << opt.optimum() << std::endl;
    double sum_abs_of_elems = 0.;
    std::vector<double> scores = opt.values(); 
    for(std::vector<double>::iterator it = scores.begin(); it != scores.end(); ++it)
        sum_abs_of_elems += std::abs(*it);
    std::cout << "sum abs scores: " << std::setprecision(16) << sum_abs_of_elems << std::endl;

    // // fit at optimal smoothing level
    // m.fit(opt.optimum());

    // double lambdaD = 0.000562341325190349; 
    // double lambdaT = 0.0001; // read_csv<double>(datadir + "lambda.csv", false, false).as_matrix()(0,0);
    // std::cout << "lambdaD = " << lambdaD << std::endl;
    // std::cout << "lambdaT = " << lambdaT << std::endl;
    // m.fit(lambdaD, lambdaT);
    
    // write_csv(datadir + "f.csv", m.f());
    // write_csv(datadir + "beta.csv", m.beta());

    return 0;
}


// test norm gcv other models --> PER TESTARE SE RUNNA IL CALCOLO DI norm in GCV con missing 
//    mesh:         unit square
//    sampling:     locations != nodes
//    penalization: anisotropic diffusion
//    time penalization: separable
//    covariates:   yes
//    BC:           no
//    order FE:     1
int test_norm_gcv() {

    // geometry 
    using PointT = Eigen::Matrix<double, 2, 1>;

    Eigen::Matrix<double, Dynamic, Dynamic> points = read_csv<double>("my_data/mesh/unit_square_reduced_censoring_476/points.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> elements = read_csv<int>("my_data/mesh/unit_square_reduced_censoring_476/elements.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> boundary = read_csv<int>("my_data/mesh/unit_square_reduced_censoring_476/boundary.csv").as_matrix();

    elements.array() -= 1; // non necessario

    std::string datadir = "my_data/msr/06-missing/";
    const unsigned int M = 8; 
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 1, M);  // ATT qui non bisogna fare più M-1 come nella vecchia lib!! Vuole direttamente il numero di nodi, cioè M!!! 
    std::cout << "ATT check M" << std::endl;
    Triangulation<2, 2> D(points, elements, boundary);

    GeoFrame data(D, T);

    // data 
    auto& l = data.insert_scalar_layer<POINT, POINT>("layer", std::pair{"my_data/msr/06/space_locs.csv", "my_data/msr/06/time_locs.csv"});
    l.load_csv<double>(datadir + "response.csv");
    l.load_csv<double>(datadir + "design_matrix.csv");

    std::cout << "geoframe:" << data[0] << std::endl;

    // physics 
    FeSpace Vh(D, P1<1>);   // functional space definition

    // trial and test function definition
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a_D = integral(D)(dot(grad(f), grad(v)));
    // homogeneous forcing linear form
    ZeroField<2> u_D;
    auto F_D = integral(D)(u_D * v);

    BsSpace Bh(T, 3);   // cubic B-splines in time
    TrialFunction g(Bh);
    TestFunction  w(Bh);
    auto a_T = integral(T)(dxx(g) * dxx(w));
    ZeroField<1> u_T;
    auto F_T = integral(T)(u_T * w);    


    std::vector<double> lambdas_d;
    std::vector<double> lambdas_t;
    std::vector<Eigen::Matrix<double, Dynamic, 1>> lambdas_d_t;
    for(double xs = -4.5; xs <= -3.0; xs += 0.5)
        lambdas_d.push_back(std::pow(10,xs));

    for(double xt = -4.0; xt <= -2.0; xt += 2.0)
        lambdas_t.push_back(std::pow(10,xt));

    for(auto i = 0; i < lambdas_d.size(); ++i)
        for(auto j = 0; j < lambdas_t.size(); ++j) 
            lambdas_d_t.push_back(Eigen::Matrix<double, 2, 1>(lambdas_d[i], lambdas_t[j]));
    
    Eigen::Matrix<double, Dynamic, 2> lambda_grid(lambdas_d.size()*lambdas_t.size(), 2);
    for(int i = 0; i < lambdas_d.size(); ++i) { 
        for (int j = 0; j < lambdas_t.size(); ++j) {
            lambda_grid(i * lambdas_t.size() + j, 0) = lambdas_d[i];
            lambda_grid(i * lambdas_t.size() + j, 1) = lambdas_t[j];
        }
    }
    std::cout << "Lambda grid: " << lambda_grid << std::endl;

    GridSearch<2> opt;

    // std::cout << "-----------------SRPDE------------------" << std::endl;
    // // SRPDE
    // SRPDE m_sr("y ~ x1 + x2 + f", data, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));  
    // opt.optimize(m_sr.gcv(200, 1234), lambda_grid);

    // std::cout << "-----------------QSRPDE------------------" << std::endl;
    // // QSRPDE
    // QSRPDE m_qsr("y ~ x1 + x2 + f", data, 0.50, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));
    // opt.optimize(m_qsr.gcv(200, 1234), lambda_grid);

    std::cout << "-----------------GSRPDE------------------" << std::endl;
    // GSRPDE
    GSRPDE m_gsr("y ~ x1 + x2 + f", data, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));
    opt.optimize(m_gsr.gcv(200, 1234), lambda_grid);

    return 0;
}


