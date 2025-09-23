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

// #include <fdaPDE/models.h>   // fdaPDE
#include <gtest/gtest.h>     // testing framework

// namespace fdapde {
// namespace test {

// [[maybe_unused]] constexpr double testing_double_tolerance = 1e-7;

// // floating point comparision utilities
// template <typename Scalar>
// bool almost_equal(
//   const Eigen::Matrix<Scalar, Dynamic, Dynamic>& op1, const Eigen::Matrix<Scalar, Dynamic, Dynamic>& op2,
//   double epsilon) {
//     return (op1 - op2).template lpNorm<Eigen::Infinity>() < epsilon ||
//            (op1 - op2).template lpNorm<Eigen::Infinity>() <
//       (std::max(op1.template lpNorm<Eigen::Infinity>(), op2.template lpNorm<Eigen::Infinity>()) * epsilon);
// }
// template <typename Scalar>
// bool almost_equal(
//   const Eigen::Matrix<Scalar, Dynamic, Dynamic>& op1, const Eigen::Matrix<Scalar, Dynamic, Dynamic>& op2) {
//     return almost_equal(op1, op2, testing_double_tolerance);
// }
// template <typename Scalar> bool almost_equal(const Eigen::SparseMatrix<Scalar>& op1, std::string op2) {
//     Eigen::SparseMatrix<Scalar> mem_buff;
//     Eigen::loadMarket(mem_buff, op2);
//     return almost_equal(op1, mem_buff);
// }
// template <typename Scalar> bool almost_equal(const Eigen::Matrix<Scalar, Dynamic, Dynamic>& op1, std::string op2) {
//     Eigen::SparseMatrix<Scalar> mem_buff;
//     Eigen::loadMarket(mem_buff, op2);
//     return almost_equal(op1, Eigen::Matrix<Scalar, Dynamic, Dynamic>(mem_buff));
// }
// template <typename Scalar> bool almost_equal(const std::vector<Scalar>& op1, std::string op2) {
//     Eigen::SparseMatrix<Scalar> mem_buff;
//     Eigen::Matrix<double, Dynamic, Dynamic> values(op1.size(), 1);

//     Eigen::loadMarket(mem_buff, op2);
//     for (int i = 0; i < op1.size(); ++i) { values(i, 0) = op1[i]; }
//     return almost_equal(values, Eigen::Matrix<Scalar, Dynamic, Dynamic>(mem_buff));
// }

// }   // namespace test
// }   // namespace fdapde

// #include "src/sr.cpp"
// #include "src/gsr.cpp"
// #include "src/qsr.cpp"
// #include "src/de.cpp"

int main(/*int argc, char **argv*/){
  // start testing
  // testing::InitGoogleTest(&argc, argv);
  // return RUN_ALL_TESTS();

  return 0;
}
