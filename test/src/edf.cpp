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

/// @brief checks probe reuse and resizing against independently initialized smoothers
template <typename GeoFrame, typename Penalty, typename... Lambda>
void check_edf_probes(
  const GeoFrame& data,     // observations and their spatial or temporal locations
  const Penalty& penalty,   // discretized model specification
  Lambda... lambda          // fixed smoothing parameters
) {
    SRPDE model("y ~ f", data, penalty);
    model.fit(lambda...);
    const int seed = 476813;
    const double first = model.edf(3, seed);
    // the first EDF call must initialize its probe storage and yield a finite estimate
    ASSERT_TRUE(std::isfinite(first));
    // an unchanged probe count must reuse the existing sample, even with a different seed argument
    EXPECT_DOUBLE_EQ(model.edf(3, 123), first);
    for (int r : {model.n_obs(), 2}) {
        SRPDE reference("y ~ f", data, penalty);
        reference.fit(lambda...);
        const double expected = reference.edf(r, seed);
        const double actual = model.edf(r, seed);
        // increasing to n_locs or shrinking must reproduce a fresh smoother with the requested sample
        EXPECT_NEAR(actual, expected, 1e-12);
        // the resized sample must also survive repeated EDF calls without being regenerated
        EXPECT_DOUBLE_EQ(model.edf(r, 123), actual);
    }
}

// check the elliptic EDF cache through the public regression model
TEST(edf, elliptic_probe_cache) {
    auto D = Triangulation<2, 2>::Rectangle(0, 1, 0, 1, 4, 4);
    GeoFrame data(D);
    auto& layer = data.insert_scalar_layer<POINT>("l1", MESH_NODES);
    const Eigen::MatrixXd y = Eigen::VectorXd::LinSpaced(D.n_nodes(), 0, 1);
    layer.load_blk("y", y);
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u;
    auto F = integral(D)(u * v);
    check_edf_probes(data, fe_ls_elliptic(a, F), 0.01);
}

// check both monolithic space-time solvers and the iterative separable solver against fresh EDF samples
TEST(edf, space_time_probe_cache) {
    auto D = Triangulation<2, 2>::Rectangle(0, 1, 0, 1, 3, 3);
    auto T = Triangulation<1, 1>::Interval(0, 1, 4);
    GeoFrame data(D, T);
    auto& layer = data.insert_scalar_layer<POINT, POINT>("l1", std::pair {MESH_NODES, MESH_NODES});
    const Eigen::MatrixXd y = Eigen::VectorXd::LinSpaced(D.n_nodes() * T.n_nodes(), 0, 1);
    layer.load_blk("y", y);
    FeSpace Vh(D, P1<1>);
    TrialFunction f(Vh);
    TestFunction v(Vh);
    auto a_D = integral(D)(dot(grad(f), grad(v)));
    ZeroField<2> u_D;
    auto F_D = integral(D)(u_D * v);
    BsSpace Bh(T, 3);
    TrialFunction g(Bh);
    TestFunction w(Bh);
    auto a_T = integral(T)(dxx(g) * dxx(w));
    ZeroField<1> u_T;
    auto F_T = integral(T)(u_T * w);
    const Eigen::VectorXd ic = Eigen::VectorXd::Zero(D.n_nodes());
    {
        SCOPED_TRACE("parabolic monolithic");
        check_edf_probes(data, fe_ls_parabolic_mono(std::pair {a_D, F_D}, ic), 0.01, 0.01);
    }
    {
        SCOPED_TRACE("separable monolithic");
        check_edf_probes(data, fe_ls_separable_mono(std::pair {a_D, F_D}, std::pair {a_T, F_T}), 0.01, 0.01);
    }
    {
        SCOPED_TRACE("separable iterative");
        check_edf_probes(data, fe_ls_separable_cdti(std::pair {a_D, F_D}), 0.01, 0.01);
    }
}
