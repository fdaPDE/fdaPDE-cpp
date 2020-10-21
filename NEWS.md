# fdaPDE 1.1-1

## New features

1) Optimization methods (Newton's methods) to find best smoothing parameter through GCV minimization
2) smooth regression for non-gaussian data (GLM model)

# fdaPDE 1.1-0

## New features

1) smooth regression for space-time data
2) density estimation
3) faster search algorithm 

# fdaPDE 1.0-8

## Bug fixes
Compilation errors with clang fixed.

# fdaPDE 1.0-7

## Bug fixes
Compilation error in macOS fixed.

# fdaPDE 1.0-6

## Bug fixes
Compilation error in solaris fixed.


# fdaPDE 1.0-5

## Bug fixes
gcc-ASAN problem with RTriangle fixed.


# fdaPDE 1.0-1

## Bug fixes
Compilation errors with clang fixed.


# fdaPDE 1.0

## New features

1) smooth regression for manifold domains 
2) smooth regression with areal data 
3) functional smooth PCA (SF-PCA algorithm) : function FEM.FPCA
4) Stochastic GCV computation has been added: parameter 'GCVmethod' can be used both regression and FPCA, can be either 'Stochastic' or 'Exact'.
5) Kfold cross-validation is available in FPCA algorithm.

## Deprecated functions and name changes

1) smooth.FEM.basis, smooth.FEM.PDE.basis and smooth.FEM.PDE.sv.basis are deprecated, use smooth.FEM instead.
2) the R-only functions (whose names began with R_) have been deprecated, and will be removed soon.
2) The usage of 'MESH' in classes and function names has been deprecated, now use 'mesh'.
3) Old datasets have been removed. New datasets are provided.
4) The parameter CPP_CODE has been removed. 

## Bug fixes
A bug in R/C++ indexes conversion in the space-varying regression has been fixed.
