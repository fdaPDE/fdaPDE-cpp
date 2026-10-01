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

#include <fdaPDE/models.h>

#include <cassert>

// check areal finite-element evaluation against exact integrals on the unit square
int main() {
    using namespace fdapde;
    auto domain = Triangulation<2, 2>::Rectangle(0, 1, 0, 1, 2, 2);
    FeSpace space(domain, FeP<1, 1>());
    BinaryMatrix<Dynamic, Dynamic> incidence(1, domain.n_cells());
    incidence.set();
    auto [basis, measure] = internals::areal_basis_eval(space, incidence);
    // the region containing every cell must have unit area
    assert(measure.size() == 1 && std::abs(measure[0] - 1.0) < 1e-12);
    Eigen::VectorXd constant = Eigen::VectorXd::Ones(space.n_dofs());
    // partition of unity must reproduce the exact average of the constant function
    assert(std::abs((basis * constant)[0] - 1.0) < 1e-12);
    Eigen::VectorXd x = domain.nodes().col(0);
    // linear finite elements must reproduce the exact average of x on the unit square
    assert(std::abs((basis * x)[0] - 0.5) < 1e-12);
}
