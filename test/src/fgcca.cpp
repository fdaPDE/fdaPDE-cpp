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

/*
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

    // block_4.l_compute(Eigen::Matrix<double, Eigen::Dynamic, 1>::Ones(401));
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
    // block_4.l_compute(Eigen::Matrix<double, Eigen::Dynamic, 1>::Ones(401));
    // std::cout << block_4.loadings_m().transpose() << "\n" << std::endl;
    // std::cout << block_4.components().transpose().leftCols(10) << "\n" << std::endl;

}
*/

TEST(rgcca, test_02) {
    std::string path = "../../../../projects/cca/";

    // chose options
    RGCCA<IndependentSampling>::Options options;
    options.scheme = Scheme::Factorial();
    // change defaults if needed ...

    // model initialization
    int n_comp = 3;
    int n_obs = 401;
    RGCCA<IndependentSampling> rgcca(n_obs, options, n_comp);

    // set empirical noise variance
    rgcca.set_noise_variance(0.2*0.2);

    // add blocks

    for (int i = 1; i <=4; ++i) {
        Eigen::Matrix<double, Dynamic, Dynamic> X = read_csv<double>(path + "X" + std::to_string(i) + ".csv").as_matrix();
        rgcca.add_multivariate_block("X" + std::to_string(i), X);
    }

    // add connections
    rgcca.connect(0,1);
    rgcca.connect(0,2);
    rgcca.connect(1,3);

    // fit
    const auto results = rgcca.fit();
    std::cout << results << std::endl;

    // save results
    for (const auto& block : rgcca.blocks()) {
        write_csv(path + "loadings_"+ block -> name() + ".csv", block -> loadings_m());
        write_csv(path + "components_"+ block -> name() + ".csv", block -> components_m());
    }
}

TEST(rgcca, test_03) {
    std::string path = "../../../../projects/cca/";

    // geometries
    Triangulation<1, 1> I(0, 1, 31);

    // define physics (same for all the blocks)
    FeSpace Bh(I, P1<1>);
    TrialFunction f(Bh);
    TestFunction  v(Bh);
    auto a = integral(I)(dx(f) * dx(v));
    ZeroField<1> u;
    auto F = integral(I)(u * v);

    // chose options
    RGCCA<IndependentSampling>::Options options;
    options.scheme = Scheme::Factorial();
    // change defaults if needed ...

    // model initialization
    int n_comp = 3;
    int n_obs = 401;
    RGCCA<IndependentSampling> rgcca(n_obs, options, n_comp);

    // set empirical noise variance
    rgcca.set_noise_variance(0.2*0.2);

    // add blocks
    for (int i = 1; i <= 4; ++i) {
        GeoFrame gf(I);
        Eigen::Matrix<double, Dynamic, Dynamic> X = read_csv<double>(path + "X" + std::to_string(i) + ".csv").as_matrix();
        auto& level = gf.insert_scalar_layer<POINT>("data", path + "locs_" + std::to_string(i) + ".csv");
        level.load_blk("X" + std::to_string(i), X.transpose());
        rgcca.add_functional_block("X" + std::to_string(i), gf, fe_ls_elliptic(a, F));
    }

    // add connections
    rgcca.connect(0,1);
    rgcca.connect(0,2);
    rgcca.connect(1,3);

    // fit
    const auto results = rgcca.fit();
    std::cout << results << std::endl;

    // save results
    for (const auto& block : rgcca.blocks()) {
        write_csv(path + "f_loadings_"+ block -> name() + ".csv", block -> loadings_m());
        write_csv(path + "f_components_"+ block -> name() + ".csv", block -> components_m());
    }
}

TEST(rgcca, test_04) {
    std::string path = "../../../../projects/cca/";

    // geometries
    Triangulation<1, 1> T(0, 1, 151);
    Triangulation<1, 1> I(0, 1, 31);

    // define physic in space (same for all the blocks)
    FeSpace Vh(I, P1<1>);
    TrialFunction f_D(Vh);
    TestFunction v_D(Vh);
    auto a_D = integral(I)(dx(f_D) * dx(v_D));
    ZeroField<1> u;
    auto F_D = integral(I)(u * v_D);
    auto penalty = fe_ls_elliptic(a_D, F_D);

    // chose options
    RGCCA<TimeDependentSampling>::Options options;
    options.scheme = Scheme::Factorial();
    // change defaults if needed ...

    // model initialization
    int n_comp = 3;
    int n_obs = 401;
    RGCCA<TimeDependentSampling> rgcca(n_obs, T, options, n_comp);

    // set empirical noise variance
    rgcca.set_noise_variance(0.2*0.2);

    // add blocks
    for (int i = 1; i <= 3; ++i) {
        GeoFrame gf(I);
        Eigen::Matrix<double, Dynamic, Dynamic> X = read_csv<double>(path + "X" + std::to_string(i) + ".csv").as_matrix();
        Eigen::Matrix<double, Dynamic, Dynamic> times = read_csv<double>(path + "times_" + std::to_string(i)+".csv").as_matrix();
        auto& level = gf.insert_scalar_layer<POINT>("data", path + "locs_" + std::to_string(i)+".csv");
        level.load_blk("X" + std::to_string(i), X.transpose());
        rgcca.add_functional_block("X" + std::to_string(i), times, gf, fe_ls_elliptic(a_D, F_D));
    }

    GeoFrame gf(I);
    Eigen::Matrix<double, Dynamic, Dynamic> X = read_csv<double>(path + "X4_short.csv").as_matrix();
    Eigen::Matrix<double, Dynamic, Dynamic> times = read_csv<double>(path + "times_4_short.csv").as_matrix();
    auto& level = gf.insert_scalar_layer<POINT>("data", path + "locs_4.csv");
    level.load_blk("X4", X.transpose());
    rgcca.add_functional_block("X4", times, gf, fe_ls_elliptic(a_D, F_D));

    // add connections
    rgcca.connect(0,1);
    rgcca.connect(0,2);
    rgcca.connect(1,3);

    // fit
    const auto results = rgcca.fit();
    std::cout << results << std::endl;

    // save results
    for (const auto& block : rgcca.blocks()) {
        write_csv(path + "tf_loadings_"+ block -> name() + ".csv", block -> loadings_m());
        write_csv(path + "tf_components_"+ block -> name() + ".csv", block -> components_m());
    }
}

TEST(rgcca, test_secanti) {

    // ---- a, b, c (quadratic form of μ inside ρ) ----
    const double a = 1.;
    const double b = 1.;
    const double c = 2.;

    // ---- Covariances ----
    const double C_DD = 71./31.;
    const double C_TT = 1.;
    const double C_DT = 1.00007802919;
    const double C_ND = 0.874;
    const double C_NT = 0.816;
    const double C_NN = 247./310.;

    // ---- Noise ----
    const double s = 0.33;
    const double d  = s - C_NN;

    auto f = [&](double mu) -> double {
        const double rho = std::sqrt(a*mu*mu + 2.0*b*mu + c);
        if (!std::isfinite(rho)) return std::numeric_limits<double>::infinity();

        const double P = C_TT*mu*mu + 2.0*C_DT*mu + C_DD;
        const double L = C_NT*mu + C_ND;
        return (P - d*(a*mu*mu + 2.0*b*mu + c)) - 2.0*L*rho;  // target = 0
    };


    const double mu_star = internals::find_root_secant(f, 0.0, .5, .5);

    std::cout << "mu_star = " << std::setprecision(8) << mu_star << std::endl;

}