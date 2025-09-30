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

using namespace fdapde;
using fdapde::test::almost_equal;

TEST(rgcca, test_00) {
    std::string path = "../../../../projects/cca/";

    // data
    Eigen::Matrix<double, Dynamic, Dynamic> X1 = read_csv<double>(path + "X1.csv").as_matrix();
    internals::MultivariateBlock block_1("X1", X1, 0.1);
    std::cout << block_1 << std::endl;

    Eigen::Matrix<double, Dynamic, Dynamic> X2 = read_csv<double>(path + "X2.csv").as_matrix();
    internals::MultivariateBlock block_2("X2", X2, 0.1);
    std::cout << block_2 << std::endl;

    Eigen::Matrix<double, Dynamic, Dynamic> X3 = read_csv<double>(path + "X3.csv").as_matrix();
    internals::MultivariateBlock block_3("X3", X3, 0.1);
    std::cout << block_3 << std::endl;

    Eigen::Matrix<double, Dynamic, Dynamic> X4 = read_csv<double>(path + "X4.csv").as_matrix();
    internals::MultivariateBlock block_4("X4", X4, 0.1);
    std::cout << block_4 << std::endl;

    // block_4.l_compute(Eigen::Matrix<double, Eigen::Dynamic, 1>::Ones(201));
    // std::cout << block_4.loadings_m().transpose() << "\n" << std::endl;
    // std::cout << block_4.components().transpose().leftCols(10) << "\n" << std::endl;

}

TEST(rgcca, test_01) {
    std::string path = "../../../../projects/cca/";

    // geometries
    Triangulation<1, 1> I(0, 1, 21);

    // define physics
    FeSpace Bh(I, P1<1>);
    TrialFunction f(Bh);
    TestFunction  v(Bh);
    auto a = integral(I)(dx(f) * dx(v));
    ZeroField<1> u;
    auto F = integral(I)(u * v);

    // data
    GeoFrame gf_1(I);
    {
        Eigen::Matrix<double, Dynamic, Dynamic> X1 = read_csv<double>(path + "X1.csv").as_matrix();
        auto& level = gf_1.insert_scalar_layer<POINT>("data", path + "locs_1.csv");
        level.load_blk("X1", X1.transpose());
    }
    internals::FunctionalBlock block_1("X1", gf_1, fe_ls_elliptic(a, F), 0.1);
    std::cout << block_1 << std::endl;

    GeoFrame gf_2(I);
    {
        Eigen::Matrix<double, Dynamic, Dynamic> X2 = read_csv<double>(path + "X2.csv").as_matrix();
        auto& level = gf_2.insert_scalar_layer<POINT>("data", path + "locs_2.csv");
        level.load_blk("X2", X2.transpose());
    }
    internals::FunctionalBlock block_2("X2", gf_2, fe_ls_elliptic(a, F), 0.1);
    std::cout << block_2 << std::endl;

    GeoFrame gf_3(I);
    {
        Eigen::Matrix<double, Dynamic, Dynamic> X3 = read_csv<double>(path + "X3.csv").as_matrix();
        auto& level = gf_3.insert_scalar_layer<POINT>("data", path + "locs_3.csv");
        level.load_blk("X3", X3.transpose());
    }
    internals::FunctionalBlock block_3("X3", gf_3, fe_ls_elliptic(a, F), 0.1);
    std::cout << block_3 << std::endl;

    GeoFrame gf_4(I);
    {
        Eigen::Matrix<double, Dynamic, Dynamic> X4 = read_csv<double>(path + "X4.csv").as_matrix();
        auto& level = gf_4.insert_scalar_layer<POINT>("data", path + "locs_4.csv");
        level.load_blk("X4", X4.transpose());
    }
    internals::FunctionalBlock block_4("X4", gf_4, fe_ls_elliptic(a, F), 0.1);
    block_4.set_lambda(1e-15);
    std::cout << block_4 << std::endl;
    // block_4.l_compute(Eigen::Matrix<double, Eigen::Dynamic, 1>::Ones(201));
    // std::cout << block_4.loadings_m().transpose() << "\n" << std::endl;
    // std::cout << block_4.components().transpose().leftCols(10) << "\n" << std::endl;

}


TEST(rgcca, test_02) {
    std::string path = "../../../../projects/cca/";

    // chose options
    RGCCA::Options options;
    options.tau_selection = TauSelection::Automatic;

    // model initialization
    int n_comp = 3;
    RGCCA rgcca(201, Scheme::Factorial(), options, n_comp);
    rgcca.set_noise_sigma_sqr(0.2);

    // add blocks
    Eigen::Matrix<double, Dynamic, Dynamic> X1 = read_csv<double>(path + "X1.csv").as_matrix();
    rgcca.add_multivariate_block("X1", X1);
    Eigen::Matrix<double, Dynamic, Dynamic> X2 = read_csv<double>(path + "X2.csv").as_matrix();
    rgcca.add_multivariate_block("X2", X2);
    Eigen::Matrix<double, Dynamic, Dynamic> X3 = read_csv<double>(path + "X3.csv").as_matrix();
    rgcca.add_multivariate_block("X3", X3);
    Eigen::Matrix<double, Dynamic, Dynamic> X4 = read_csv<double>(path + "X4.csv").as_matrix();
    rgcca.add_multivariate_block("X4", X4);

    // add connections
    rgcca.connect(0,1);
    rgcca.connect(0,2);
    rgcca.connect(1,3);

    /*
    // check
    for (int j = 0; j < rgcca.n_blocks(); ++j) {
        const auto& X = rgcca.blocks()[j]->data();
        Eigen::BDCSVD<RGCCA::Matrix> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
        const double s1 = svd.singularValues()(0);
        const double n  = double(rgcca.n_obs());
        const double tau = rgcca.blocks()[j]->tau();
        const double predicted = 1.0 / std::sqrt(((1.0 - tau)/double(n)) * s1 * s1 + tau);
        std::cout << "block " << j
                  << "  s1=" << s1
                  << "  predicted ||a||=" << predicted
                  << std::endl;
    }
    */

    // fit
    auto results = rgcca.fit();

    /*
    for (int j = 0; j < rgcca.n_blocks(); ++j) {
        const auto& a = rgcca.blocks()[j]->loadings().col(rgcca.h());
        std::cout << "block " << j << "  ||a||=" << a.norm()
                  << "  a^T Σ a=" << std::sqrt( (a.transpose() * rgcca.blocks()[j]->Sigma() * a)(0) )
                  << std::endl;
    }
    std::cout << std::endl;
    */

    std::cout << results << std::endl;
    std::cout << std::endl;

    for (auto& block : rgcca.blocks()) {
        write_csv(path + "loadings_"+ block -> name() + ".csv", block -> loadings_m());
        write_csv(path + "components_"+ block -> name() + ".csv", block -> components());
    }
}


TEST(rgcca, test_03) {
    std::string path = "../../../../projects/cca/";

    // geometries
    Triangulation<1, 1> I(0, 1, 21);

    // define physics (same for all the blocks)
    FeSpace Bh(I, P1<1>);
    TrialFunction f(Bh);
    TestFunction  v(Bh);
    auto a = integral(I)(dx(f) * dx(v));
    ZeroField<1> u;
    auto F = integral(I)(u * v);

    // chose options
    RGCCA::Options options;
    options.tau_selection = TauSelection::Automatic;
    options.lambda_selection = LambdaSelection::Automatic;

    // model initialization
    int n_comp = 3;
    RGCCA rgcca(201, Scheme::Factorial(), options, n_comp);
    rgcca.set_noise_sigma_sqr(0.2);

    // add blocks
    for (int i = 1; i <=4; ++i) {
        GeoFrame gf(I);
        Eigen::Matrix<double, Dynamic, Dynamic> X = read_csv<double>(path + "X"+std::to_string(i)+".csv").as_matrix();
        auto& level = gf.insert_scalar_layer<POINT>("data", path + "locs_"+std::to_string(i)+".csv");
        level.load_blk("X"+std::to_string(i), X.transpose());
        rgcca.add_functional_block("X"+std::to_string(i), gf, fe_ls_elliptic(a, F));
    }

    // add connections
    rgcca.connect(0,1);
    rgcca.connect(0,2);
    rgcca.connect(1,3);

    // set lambda parameter
    rgcca.set_lambda_all(-1);

    /*
    // check
    for (int j = 0; j < rgcca.n_blocks(); ++j) {
        const auto& X = rgcca.blocks()[j]->data();
        Eigen::BDCSVD<RGCCA::Matrix> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
        const double s1 = svd.singularValues()(0);
        const double n  = double(rgcca.n_obs());
        const double tau = rgcca.blocks()[j]->tau();
        const double predicted = 1.0 / std::sqrt(((1.0 - tau)/double(n)) * s1 * s1 + tau);
        std::cout << "block " << j
                  << "  s1=" << s1
                  << "  predicted ||a||=" << predicted
                  << std::endl;
    }
    */

    // fit
    auto results = rgcca.fit();

    /*
    for (int j = 0; j < rgcca.n_blocks(); ++j) {
        const auto& a = rgcca.blocks()[j]->loadings().col(rgcca.h());
        std::cout << "block " << j << "  ||a||=" << a.norm()
                  << "  a^T Σ a=" << std::sqrt( (a.transpose() * rgcca.blocks()[j]->Sigma() * a)(0) )
                  << std::endl;
    }
    std::cout << std::endl;
    */

    std::cout << results << std::endl;
    std::cout << std::endl;

    for (auto& block : rgcca.blocks()) {
        write_csv(path + "f_loadings_"+ block -> name() + ".csv", block -> loadings_m());
        write_csv(path + "f_components_"+ block -> name() + ".csv", block -> components());
    }
}
