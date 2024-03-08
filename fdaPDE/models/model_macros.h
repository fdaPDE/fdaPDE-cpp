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

#ifndef __MODEL_MACROS_H__
#define __MODEL_MACROS_H__

// imports all basic symbols a model type can expect independently on its specific module membership
#define IMPORT_MODEL_SYMBOLS                                                                                           \
    using Base::n_spatial_basis; /* number of basis function for discretization in space N */                          \
    using Base::n_locs;          /* number of locations p_1 ... p_n where data are observed */                         \
    using Base::Psi;             /* n x N matrix of spatial basis evaluations at p_1 ... p_n */                        \
    using Base::PsiTD;           /* block P^T*D, being D the matrix of subdomains' measure */                          \
                                 /* returns P^T if sampling is not areal */                                            \
    using Base::R1;              /* discretization of differential operator L (tensorized for */                       \
                                 /* space-time problems) */                                                            \
    using Base::R0;              /* mass matrix in space (tensorized for space-time problems) */                       \
    using Base::u;               /* discretization of forcing term */                                                  \
    using Base::pde;             /* differential operator Lf- u in space */                                            \
    using Base::data;            /* BlockFrame object containing data */

// imports all basic symbols a model can expect from a RegressionBase
// symbols specific for the regularization type used need to be imported via dedicated using declaration
#define IMPORT_REGRESSION_SYMBOLS                                                                                      \
    IMPORT_MODEL_SYMBOLS                                                                                               \
    using Base::y;              /* vector of observations y = [y_1 ... y_n] */                                         \
    using Base::n_obs;          /* number of observations n */                                                         \
    using Base::W;              /* matrix of observation weights W_ = diag[W_1 ... W_n] */                             \
    using Base::q;              /* number of covariates */                                                             \
    using Base::X;              /* n x q design matrix X = [X_1 ... X_q] */                                            \
    using Base::XtWX;           /* q x q dense matrix X^T*W*X */                                                       \
    using Base::invXtWX;        /* partialPivLU factorization of X^T*W*X */                                            \
    using Base::has_covariates; /* true if the model is semi-parametric */                                             \
    using Base::has_weights;    /* true if heteroscedastic observations are assumed */                                 \
    using Base::lmbQ;           /* efficient left multiplication by Q */                                               \
    using Base::f_;             /* estimate of the nonparametric part of the model */                                  \
    using Base::g_;             /* PDE misfit */                                                                       \
    using Base::beta_;          /* estimate of coefficient vector for parametric part */                               \
    using Base::U_;             /* woodbury matrix [\Psi^T*D*W*X, 0] */                                                \
    using Base::V_;             /* woodbury matrix [X^T*W*\Psi,   0] */

#define FDAPDE_DELETE_INIT                                                                                             \
    void init() { return; }

#define FDAPDE_DEFINE_MODEL_GETTER                                                                                     \
    inline const Model& model() const { return static_cast<const Model&>(*this); }                                     \
    inline Model& model() { return static_cast<Model&>(*this); }

// macro for runtime sanity checks on data, should be the first instruction in a solve() implementation
#define BLOCK_FRAME_SANITY_CHECKS                                                                                      \
    /* stop if incoming data has no observations */                                                                    \
    if (!data().has_block(OBSERVATIONS_BLK))                                                                           \
        throw std::logic_error("bad BlockFrame, model without observations is ill-formed");

// standardized definitions for model BlockFrame's blocks
#define OBSERVATIONS_BLK  "OBSERVATIONS"    // observations
#define INDEXES_BLK       "INDEXES"         // observation indices
#define DESIGN_MATRIX_BLK "DESIGN_MATRIX"   // covariates
#define WEIGHTS_BLK       "WEIGHTS"         // weights for heteroskedastic observations

#endif   // __MODEL_MACROS_H__
