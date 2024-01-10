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

#ifndef __RSVD_H__
#define __RSVD_H__

#include <fdaPDE/utils.h>

#include <Eigen/SVD>

#include "../model_traits.h"

namespace fdapde {
namespace models {

// finds a rank r matrix U minimizing \norm{X - U*\Psi^\top}_F^2 + Tr[U*P_{\lambda}(f)*U^\top]
// with P_{\lambda}(f) the penalty term and \norm{}_F the Frobenius norm
template <typename Model_> class RSVD {
   private:
    using Model = typename std::remove_reference<Model_>::type;
    static constexpr int n_lambda = Model::n_lambda;
    const Model* m_;

    // let E*\Sigma*F^\top the reduced (rank r) SVD of X*\Psi*(D^{1})^\top, with D^{-1} the inverse of the cholesky
    // factor of \Psi^\top * \Psi + P(\lambda), then
    DMatrix<double> H_;        // matrix E in the reduced SVD of X*\Psi*(D^{-1})^\top (scores)
    DMatrix<double> W_;        // \Sigma*F^\top*D^{-1} (field expansion coefficients)
    DVector<double> w_norm_;   // L^2 norm of estimated fields
   public:
    // constructors
    RSVD() = default;
    RSVD(const Model* m) : m_(m) { }

    // executes the Regularized Singular Value Decomposition on data X for a given rank r and smoothing level \lambda
    void compute(const DMatrix<double>& X, const SVector<n_lambda>& lambda, std::size_t r) {
        // compute matrix C = \Psi^\top*\Psi + P(\lambda)
        DMatrix<double> C = m_->Psi().transpose() * m_->Psi() + m_->P(lambda);
        // compute the inverse of the cholesky factor of C, D^{-1}
        DMatrix<double> D = C.llt().matrixL();
        DMatrix<double> invD = D.inverse();

        // compute SVD of X*\Psi*(D^{-1})^\top
        Eigen::JacobiSVD<DMatrix<double>> svd(
          X * m_->Psi() * invD.transpose(), Eigen::ComputeThinU | Eigen::ComputeThinV);
        // store results
        H_ = svd.matrixU().leftCols(r);
        W_ = (svd.singularValues().head(r).asDiagonal() * svd.matrixV().leftCols(r).transpose() * invD).transpose();
        w_norm_ = (W_.transpose() * m_->R0() * W_).array().sqrt();   // L2 norm of estimated fields
        return;
    }

    // getters
    const DMatrix<double>& H() const { return H_; }
    const DMatrix<double>& W() const { return W_; }
    const DVector<double>& w_norm() const { return w_norm_; }
};

}   // namespace models
}   // namespace fdapde

#endif   // __RSVD_H__
