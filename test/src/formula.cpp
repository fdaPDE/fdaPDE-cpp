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

// preserve the legacy covariate accessor while retaining separate pipe effects
TEST(formula, legacy_rhs) {
    const Formula basic("y ~ x1 + x2 + f");
    const std::vector<std::string> expected {"x1", "x2", "f"};
    // ordinary formulas must retain every right-hand-side token in its original order
    EXPECT_EQ(basic.rhs(), expected);
    // the legacy and current accessors must refer to the same covariate collection
    EXPECT_EQ(&basic.rhs(), &basic.covs());
    const Formula grouped("y ~ x1 + x2|group");
    // a pipe effect must remain separate from the ordinary covariate collection
    EXPECT_EQ(grouped.rhs(), std::vector<std::string> {"x1"});
    // parsing a single pipe token must produce exactly one effect
    ASSERT_EQ(grouped.efxs().size(), 1);
    // the effect covariate must preserve the token to the left of the pipe
    EXPECT_EQ(grouped.efxs()[0].cov(), "x2");
    // the effect grouping term must preserve the token to the right of the pipe
    EXPECT_EQ(grouped.efxs()[0].efx(), "group");
}
