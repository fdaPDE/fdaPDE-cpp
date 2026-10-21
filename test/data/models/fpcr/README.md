# fPCR regression references

These outputs freeze the fPCR implementation at library commit
`9864fa6eb94072ad2b79de989db8a46e4ae81b1e` with core commit
`472f9bd374a63fbf851106e5c45f8870e31e8189`.

The cases in `test/src/fpcr.cpp` read the existing `X.csv` and `Y.csv` from
`models/fpls/2D_test1` and `models/fpls/2D_test2`, respectively. They use the
same `unit_square_60` mesh, P1 Laplacian penalty, smooth predictor mean and
ordinary response means as the corresponding fPLS tests.

Both fits retain three components, use exact SVD initialization, at most
20 power updates and objective tolerance `1e-2`:

- `2D_test1`: fixed smoothing parameter `10` for centering and all components
- `2D_test2`: GCV grid `{1e-4, 1e-3, 1e-2, 1e-1, 1}`, 1000 EDF probes and seed
  `476813` for centering and component calibration; selected component penalties
  are `{1e-3, 1e-3, 1e-2}`

Each case stores the uncentered `Y_hat` (60 × 2), uncentered `X_hat`
(60 × 3600) and centered regression coefficients `B_hat` (3600 × 2).
Like the fPLS references, these `.csv` files contain MatrixMarket coordinate
data, written by `Eigen::saveMarket` with 17 significant digits. The tests
only read them and apply the existing `1e-7` comparison tolerance.

Run from `test/build` with `./fdapde_test --gtest_filter='fpcr.*'`.
