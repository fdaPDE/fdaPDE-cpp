#ifndef __MODEL_MACROS_H__
#define __MODEL_MACROS_H__

  // macros for the import of common symbols to avoid long annoying lists of using declarations in model implemetations

  // this macro is intended to import all **common** symbols a model type can expect from its parent classes
#define IMPORT_MODEL_SYMBOLS                                                             \
  using Base::y;       /* vector of observations y = [y_1 ... y_n] */                    \
  using Base::n_obs;   /* number of observations n */                                    \
  using Base::n_basis; /* number of basis function for discretization in space N */      \
  using Base::n_locs;  /* number of locations p_1 ... p_n where data are observed */     \
  using Base::Psi;     /* n x N matrix of spatial basis evaluations at p_1 ... p_n */    \
  using Base::PsiTD;   /* block P^T*D, being D the matrix of subdomains' measure */      \
                       /* returns P^T if sampling is not areal */                        \
  using Base::R1;      /* discretization of differential operator L (tensorized for */   \
                       /* space-time problems) */                                        \
  using Base::R0;      /* mass matrix in space (tensorized for space-time problems) */   \
  using Base::u;       /* discretization of forcing term */                              \
  using Base::pde;     /* differential operator L (regularizing term) */                 \
  using Base::data;    /* BlockFrame object containing data */                           \
  
  // this macro is intended to import all **common** symbols a model can expect from a Regression base
  // symbols specific for the regularization type used need to be imported via dedicated using declaration
#define IMPORT_REGRESSION_SYMBOLS                                                        \
  IMPORT_MODEL_SYMBOLS;			              			                 \
  /* data access */                                                                      \
  using Base::W;             /* matrix of observation weights W_ = diag[W_1 ... W_n] */  \
  using Base::q;             /* number of covariates */                                  \
  using Base::X;             /* n x q design matrix X = [X_1 ... X_q] */                 \
  using Base::XtWX;          /* q x q dense matrix X^T*W*X */                            \
  using Base::invXtWX;       /* partialPivLU factorization of X^T*W*X */                 \
  /* utilities */                                                                        \
  using Base::hasCovariates; /* true if the model is semi-parametric */                  \
  using Base::hasWeights;    /* true if heteroscedastic observations are assumed */      \
  using Base::lmbQ;          /* efficient left multiplication by Q */                    \
  /* room for problem solution */                                                        \
  using Base::f_;            /* estimate of the nonparametric part of the model */       \
  using Base::g_;            /* PDE misfit */                                            \
  using Base::beta_;         /* estimate of coefficient vector for parametric part */    \
  using Base::U_;            /* woodbury matrix [\Psi^T*D*W*X, 0] */                     \
  using Base::V_;            /* woodbury matrix [X^T*W*\Psi,   0] */                     \


  // macro for the import of some common CRTP functionalities. Requires a Model
  // type to be in the scope of this macro
#define DEFINE_CRTP_MODEL_UTILS                                                  \
  inline const Model& model() const { return static_cast<const Model&>(*this); } \
  inline Model& model() { return static_cast<Model&>(*this); }                   \
  
  // standardized definitions for stat model BlockFrame. layers below will make heavy assumptions on
  // the layout of the BlockFrame, use these instead of manually typing the block name when accessing df_
#define OBSERVATIONS_BLK    "OBSERVATIONS"     // matrix of observations
#define INDEXES_BLK         "INDEXES"          // vector of observation indices
#define SPACE_LOCATIONS_BLK "SPACE_LOCATIONS"  // matrix of space-location
#define TIME_LOCATIONS_BLK  "TIME_LOCATIONS"   // vector of time-locations
#define SPACE_AREAL_BLK     "INCIDENCE_MATRIX" // incidence matrix of areal observations
#define DESIGN_MATRIX_BLK   "DESIGN_MATRIX"    // in regression is the design matrix
#define WEIGHTS_BLK         "WEIGHTS"          // in regression are the weights for heteroscedastic observations

#endif // __MODEL_MACROS_H__
