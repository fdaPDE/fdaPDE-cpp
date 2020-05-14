# fdaPDE_dev

This repository contains the development version of fdaPDE package (future 1.0 version).

New features wrt CRAN: smooth regression for manifold and volumetric domains, also with areal data. Smooth fPCA over 2D, 2.5D and 3D domains, also with areal data.

smooth.FEM.basis, smooth.PDE.FEM.basis, smooth.FEM.PDE.sv.basis are deprecated, smooth.FEM has to be used in all cases.

Image.FEM has been restored. Bugs in fPCA, boundary conditions and space-varying regression have been fixed. Issues of point location in 2.5D have been fixed. Areal data are still undergoing tests.

Compiled in both Win RStudio and Ubuntu 18.04 using g++ compiler. If using a Linux machine, it is advisable to install rgl, plot3D and plot3Drgl before fdaPDE.

## Subfolder structure:
/src contains all C++ code and a special file named Makevars necessary to build and install the R package,
/R contains the R functions that wrap the C++ calls,
/tests contains a script to test the package,
/data contains the data to run the tests in /tests.

## Remarks:

1) the shift of indexes from R to C++ is done within the R functions smooth.FEM and FPCA.FEM Do not use C++ scripts directly on the R mesh objects, unless you take care of shifing indexes by yourself.
