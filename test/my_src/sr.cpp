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

int test_01();
int test_02();
int test_03();
int test_06(); 
int test_12();

int main(){
    // test_01();  
    // test_02();
    // test_03();
    test_06();
    //test_12();
    return 0; 
}

// test 1
//    mesh:         unit_square_60
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
int test_01() {

    // geometry 
    using PointT = Eigen::Matrix<double, 2, 1>;

    Eigen::Matrix<double, Dynamic, Dynamic> points = read_csv<double>("my_data/mesh/unit_square_60/points.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> elements = read_csv<int>("my_data/mesh/unit_square_60/elements.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> boundary = read_csv<int>("my_data/mesh/unit_square_60/boundary.csv").as_matrix();
    
    elements.array() -= 1; // non necessario

    std::string datadir = "my_data/sr/01/";
    Triangulation<2, 2> D(points, elements, boundary);

    GeoFrame data(D);

    // data 
    auto& l = data.insert_scalar_layer<POINT>("layer", MESH_NODES);
    l.load_csv<double>(datadir + "response.csv");

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
    SRPDE m("y ~ f", data, fe_ls_elliptic(a, F));

    // calibration
    std::vector<double> lambda_grid = {1e-4, 1e-3, 1e-2, 1e-1};
    GridSearch<1> opt;
    opt.optimize(m.gcv(), lambda_grid);

    // fit at optimal smoothing level
    m.fit(opt.optimum());

    // m.fit(1.56206e-08);
    write_csv(datadir + "f.csv", m.f());

    return 0;
}

// test 2
//    mesh:         c_shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
int test_02() {
    
    // geometry 
    using PointT = Eigen::Matrix<double, 2, 1>;

    Eigen::Matrix<double, Dynamic, Dynamic> points = read_csv<double>("my_data/mesh/c_shaped/points.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> elements = read_csv<int>("my_data/mesh/c_shaped/elements.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> boundary = read_csv<int>("my_data/mesh/c_shaped/boundary.csv").as_matrix();
    
    elements.array() -= 1; // non necessario

    std::string datadir = "my_data/sr/02/";
    Triangulation<2, 2> D(points, elements, boundary);

    GeoFrame data(D);

    // data
    auto& l1 = data.insert_scalar_layer<POINT>("l1", "my_data/sr/02/locs.csv");
    l1.load_csv<double>("my_data/sr/02/response.csv");
    l1.load_csv<double>("my_data/sr/02/design_matrix.csv");

    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);

    // modeling
    SRPDE m("y ~ x1 + x2 + f", data, fe_ls_elliptic(a, F));

    // calibration
    std::vector<double> lambda_grid = {1e-4, 1e-3, 1e-2, 1e-1};
    GridSearch<1> opt;
    opt.optimize(m.gcv(), lambda_grid);

    // fit at optimal smoothing level
    m.fit(opt.optimum());

    write_csv(datadir + "f.csv", m.f());
    write_csv(datadir + "beta.csv", m.beta());

    return 0;
}


// test 3
//    mesh:         unit_square_60
//    sampling:     locations = nodes
//    penalization: anisotropic diffusion
//    covariates:   no
//    BC:           no
//    order FE:     1
int test_03() {
    
    // geometry 
    using PointT = Eigen::Matrix<double, 2, 1>;

    Eigen::Matrix<double, Dynamic, Dynamic> points = read_csv<double>("my_data/mesh/unit_square_60/points.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> elements = read_csv<int>("my_data/mesh/unit_square_60/elements.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> boundary = read_csv<int>("my_data/mesh/unit_square_60/boundary.csv").as_matrix();
    
    elements.array() -= 1; // non necessario

    std::string datadir = "my_data/sr/03/";
    Triangulation<2, 2> D(points, elements, boundary);

    GeoFrame data(D);

    // data
    auto& l1 = data.insert_scalar_layer<POINT>("l1", MESH_NODES);
    l1.load_csv<double>("my_data/sr/03/response.csv");

    // physics: anisotropic diffussion
    Eigen::Matrix<double, 2, 2> K;
    K << 1, 0, 0, 4;
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a = integral(D)(dot(K * grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);

    // modeling
    SRPDE m("y ~ f", data, fe_ls_elliptic(a, F));

    // calibration
    std::vector<double> lambda_grid = {1e-4, 1e-3, 1e-2, 1e-1};
    GridSearch<1> opt;
    opt.optimize(m.gcv(), lambda_grid);

    // fit at optimal smoothing level
    m.fit(opt.optimum());

    write_csv(datadir + "f.csv", m.f());

    return 0;
}

// test 6
//    mesh:         unit_square_21
//    sampling:     locations = nodes
//    space penalization: laplacian
//    time penalization: separable 
//    covariates:   no
//    BC:           no
//    order FE:     1
int test_06() {

    // geometry
    Triangulation<1, 1> T = Triangulation<1, 1>::Interval(0, 2, 11);
    Eigen::Matrix<double, Dynamic, Dynamic> points = read_csv<double>("my_data/mesh/unit_square_21/points.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> elements = read_csv<int>("my_data/mesh/unit_square_21/elements.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> boundary = read_csv<int>("my_data/mesh/unit_square_21/boundary.csv").as_matrix();
    
    elements.array() -= 1; // non necessario

    std::string datadir = "my_data/sr/06/";
    Triangulation<2, 2> D(points, elements, boundary);

    // data
    GeoFrame data(D, T);
    auto& l1 = data.insert_scalar_layer<POINT, POINT>("l1", std::pair {MESH_NODES, MESH_NODES});
    l1.load_csv<double>("my_data/sr/06/response.csv");  
    
    // physics
    FeSpace Vh(D, P1<1>);   // linear finite element in space
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    auto a_D = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u_D;
    auto F_D = integral(D)(u_D * v);

    BsSpace Bh(T, 3);   // cubic B-splines in time
    TrialFunction g(Bh);
    TestFunction  w(Bh);
    auto a_T = integral(T)(dxx(g) * dxx(w));
    ZeroField<1> u_T;
    auto F_T = integral(T)(u_T * w);

    // modeling
    SRPDE m("y ~ f", data, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}));

    // calibration
    std::vector<double> lambda_grid_d = {1e-6, 1e-5};
    std::vector<double> lambda_grid_t = {1e-4, 1e-3}; 
    Eigen::MatrixXd lambda_grid(lambda_grid_d.size()*lambda_grid_t.size(), 2);  
    for(auto i = 0; i < lambda_grid_d.size(); ++i){
        for(auto j = 0; j < lambda_grid_t.size(); ++j){
            lambda_grid(i * lambda_grid_t.size() + j, 0) = lambda_grid_d[i];
            lambda_grid(i * lambda_grid_t.size() + j, 1) = lambda_grid_t[j];
        }
    }
    GridSearch<2> opt;
    opt.optimize(m.gcv(), lambda_grid);

    // fit at optimal smoothing level
    m.fit(opt.optimum());

    write_csv(datadir + "f.csv", m.f());

    return 0;
}


// test 12
//    mesh:         unit_square_60
//    sampling:     areal
//    penalization: non-constant PDE coefficients
//    covariates:   no
//    BC:           no
//    order FE:     1
int test_12() {

    // geometry 
    using PointT = Eigen::Matrix<double, 2, 1>;

    Eigen::Matrix<double, Dynamic, Dynamic> points = read_csv<double>("my_data/mesh/quasi_circle/points.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> elements = read_csv<int>("my_data/mesh/quasi_circle/elements.csv").as_matrix();
    Eigen::Matrix<int, Dynamic, Dynamic> boundary = read_csv<int>("my_data/mesh/quasi_circle/boundary.csv").as_matrix();
    
    elements.array() -= 1; // non necessario

    std::string datadir = "my_data/sr/12/";
    Triangulation<2, 2> D(points, elements, boundary);

    GeoFrame data(D);

    // data
    auto& l1 = data.insert_scalar_layer<POLYGON>("l1", "my_data/sr/12/incidence_mat.csv");
    l1.load_csv<double>("my_data/sr/12/response.csv");

    // physics
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;  
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction  v(Vh);
    FeCoeff<2, 2, 2, matrix_t> K(read_csv<double>("my_data/sr/12/diffusion.csv").as_matrix());
    FeCoeff<2, 2, 1, matrix_t> b(read_csv<double>("my_data/sr/12/transport.csv").as_matrix());
    auto a = integral(D)(dot(K * grad(f), grad(v)) + dot(b, grad(f)) * v);
    FeCoeff<2, 1, 1, vector_t> u(read_csv<double>("my_data/sr/12/force.csv").as_matrix());
    auto F = integral(D)(u * v);

    // modeling
    SRPDE m("y ~ f", data, fe_ls_elliptic(a, F));

    // calibration
    std::vector<double> lambda_grid = {1e-4, 1e-3, 1e-2, 1e-1};
    GridSearch<1> opt;
    opt.optimize(m.gcv(), lambda_grid);

    // fit at optimal smoothing level
    m.fit(opt.optimum());

    write_csv(datadir + "f.csv", m.f());

    return 0;
}