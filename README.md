<div align="center"> <h1> fdaPDE </h1>

<h5> Physics-Informed Spatial and Functional Data Analysis </h5> </div>

![test-linux-gcc](https://img.shields.io/github/actions/workflow/status/fdaPDE/fdaPDE-cpp/test-linux-gcc.yml?branch=stable&label=test-linux-gcc)
![test-linux-clang](https://img.shields.io/github/actions/workflow/status/fdaPDE/fdaPDE-cpp/test-linux-clang.yml?branch=stable&label=test-linux-clang)
![test-macos-clang](https://img.shields.io/github/actions/workflow/status/fdaPDE/fdaPDE-cpp/test-macos-clang.yml?branch=stable&label=test-macos-clang)

fdaPDE is a C++ library for the analysis of spatial and functional data observed over complex multidimensional domains, featuring a Partial Differential Equation regularization. 

It is built on top of the [fdaPDE Core Library](https://github.com/fdaPDE/fdaPDE-core).

## Documentation
Documentation can be found on our [documentation site](https://fdapde.github.io/)

## Functional regression models

Include `<fdaPDE/models.h>` to use `fPLS` or `fPCR`. Both take centered functional
predictors through a single GeoFrame layer and centered responses with subjects in
rows; scalar and multivariate responses retain every response column.

`fPLS` saves direction and loading GCV curves, their actual EDF estimates, candidate
tuples and zero-based selected candidate indices during calibration. The public
`direction_gcv_edfs()`, `loading_gcv_edfs()`, `direction_lambda_grid()`,
`loading_lambda_grid()`, `direction_selected_indices()` and
`loading_selected_indices()` accessors perform no solves or EDF estimates. Fixed
fits clear the calibration snapshots and record the fixed penalties; each repeated
fit replaces the previous snapshot.

`fPCR("X", Y, data, penalty)` reuses native power-iteration fPCA. Its `fit` takes a
component count, fixed penalty tuple or GCV candidate tuples, calibration/SVD flags,
iteration limit, tolerance, EDF probe count and seed. `X_latent_scores()` projects
final loadings sequentially; `score_coefficients(h)` refits OLS for every prefix.
`fitted(h)`, `predict(Xnew, h)`, `reconstructed(h)` and `Beta(h)` expose prefix
results. Beta maps the original centered predictors through the functional basis,
including all preceding deflations: `X * Psi * Beta(h) == fitted(h)`. Prefix zero
selects all fitted components. `gcv_values()`, `gcv_edf()`, `lambda_grid()` and
`selected_indices()` expose the evaluations used by its fPCA calibration.

Power-solver GCV now uses floating-point residual degrees of freedom (`m - edf`)
rather than truncating them to an integer. Explicit seeds and probe counts are
honored on repeated elliptic fits. These corrections can change GCV selections
relative to callers that previously used default random probes. Explicit EDF
controls are supported by the fPCA power policy; other policies reject nondefault
controls. This branch's pinned core still provides the legacy one-argument assertion
macro, so new permanent public validation uses exceptions pending the core assertion
API migration.
