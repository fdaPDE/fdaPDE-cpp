#include <fdaPDE/fdapde.h>

#include <filesystem>

using namespace fdapde;
namespace fs = std::filesystem;

int main() {
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    using binary_t = BinaryMatrix<Dynamic, Dynamic>;

    // geometry
    Triangulation<2, 2> D = Triangulation<2, 2>::UnitSquare(10);

    // data
    matrix_t data = read_csv<double>("../data/vsr/tensors/Y.csv").as_matrix();
    // processing of tensors in the log-euclidean is simply euclidean in the logarithm domain
    // move to logarithm domain

    int n_comp = 3;   // number of components in the euclidean domain
    matrix_t log_data(data.rows(), n_comp);
    for (int i = 0; i < log_data.rows(); ++i) {
        Eigen::Matrix<double, 2, 2> m;
        m << data(i, 0), data(i, 1), data(i, 2), data(i, 3);
        Eigen::Matrix<double, 2, 2> log_m = logm(m);
        // std::cout << log_m << std::endl;
        log_data(i, 0) = log_m(0, 0);
        log_data(i, 1) = log_m(1, 1);
        log_data(i, 2) = std::sqrt(2) * log_m(0, 1);
    }
    GeoFrame gf(D);
    auto& l = gf.insert_scalar_layer<POINT>("layer", "../data/vsr/tensors/locations.csv");
    for (int i = 0; i < log_data.cols(); ++i) { l.load_vec("K" + std::to_string(i), log_data.col(i)); }
    int n_obs = gf[0].rows();   // number of observations
    // std::cout << l << std::endl;

    // physics
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    int n_dofs = Vh.n_dofs();   // number of degrees of freedom

    // modeling
    auto vsr = [&](double lambda) -> std::tuple<matrix_t, matrix_t, binary_t, double> {   // smoothing of vector field
        matrix_t coeff(n_dofs, n_comp);
        matrix_t fitted(n_obs, n_comp);
        binary_t nan_pattern(n_obs, n_comp);
        double edf = 0;
        // fit for each component
        for (int i = 0; i < n_comp; ++i) {
            SRPDE m("K" + std::to_string(i) + " ~ f", gf, fe_ls_elliptic(a, F));
            m.fit(lambda);
            fitted.col(i) = m.fitted();
            coeff.col(i) = m.f();
            nan_pattern.col(i) = m.nan_pattern();
            if (i == 0) { edf = m.edf(100, 1); }   // compute edf once
        }
        return std::tuple {coeff, fitted, nan_pattern, edf};
    };

    // calibration
    auto gcv = [&](auto lambda) -> double {
        const auto& [f, fitted, nan_pattern, edf] = vsr(lambda[0]);
        int n = n_obs, k = n_comp;
        double dor = n - edf;
        matrix_t residuals = (~nan_pattern).select(fitted - log_data, 0);
        return ((n * k) / std::pow(dor, 2)) * (residuals * residuals.transpose()).trace();
    };
    int n_lambda = 200;
    double lambda_min = -9.0, lambda_max = -2.0;
    double grid_step = std::abs(lambda_max - lambda_min) / n_lambda;
    std::vector<double> lambda_grid(n_lambda);
    for (int i = 0; i < n_lambda; ++i) { lambda_grid[i] = std::pow(10, lambda_min + grid_step * i) / n_obs; }
    GridOptimizer<1> opt;
    opt.optimize(gcv, lambda_grid);
    write_csv("gcv_curve.csv", opt.values());

    // fit
    const auto& [coeff, fitted, nan_pattern, edf] = vsr(opt.optimum()[0]);

    // export
    // evaluate estimated tensor field at quadrature nodes
    // matrix_t locs = quadrature_nodes(D, QS2DP2);
    matrix_t locs = D.nodes();
    Eigen::SparseMatrix<double> Psi = internals::point_basis_eval(Vh, locs);
    matrix_t fitted_ = Psi * coeff;

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    VSRPDE model("Y ~ F", gf, fe_ls_elliptic(a, F));
    GridOptimizer<1> opt_1;
    opt_1.optimize(model.gcv(100, 1), lambda_grid);

    model.fit(opt.optimum()[0]);

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "confronto" << std::endl;
    for (int i = 0; i < n_lambda; ++i) std::cout << opt_1.values()[i] << "\t" << opt.values()[i] << std::endl;

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "coeff modello wrappato" << std::endl;
    std::cout << model.F().topRows(10) << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "coeff codice Ale" << std::endl;
    std::cout << coeff.topRows(10) << std::endl;

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "fitted modello wrappato" << std::endl;
    std::cout << model.fitted().topRows(10) << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "fitted codice Ale" << std::endl;
    std::cout << fitted.topRows(10) << std::endl;

    /*

    // back to exponential domain
    matrix_t exp_data(locs.rows(), 4);
    for (int i = 0; i < locs.rows(); ++i) {
        Eigen::Matrix<double, 2, 2> m;
        m(0, 0) = fitted_(i, 0);
        m(1, 1) = fitted_(i, 1);
        m(0, 1) = fitted_(i, 2) / std::sqrt(2);
        m(1, 0) = m(0, 1);

        Eigen::Matrix<double, 2, 2> exp_m = expm(m);
        for (int j = 0; j < exp_m.rows(); ++j) {
            for (int h = 0; h < exp_m.cols(); ++h) { exp_data(i, j * exp_m.rows() + h) = exp_m(j, h); }
        }
    }
    write_csv("f_hat.csv", exp_data);
    */

    return 0;
}
