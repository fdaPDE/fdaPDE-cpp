<div align="center"> <h1> fdaPDE </h1>

<h5> Bring the expressiveness power of PDE into statistics </h5> </div>

<div align="center">

![GitHub Workflow Status](https://img.shields.io/github/workflow/status/AlePalu/fdaPDE/build-workflow?label=build&style=for-the-badge) ![GitHub Workflow Status](https://img.shields.io/github/workflow/status/AlePalu/fdaPDE/test-workflow?label=test-status&style=for-the-badge)

</div>


This is the rewriting repo of the fdaPDE project. Code here will be merged into mainstream when it is ready.

## Documentation
Documentation for both end-users and developers can be found on the [documentation site](https://alepalu.github.io/fdaPDE/)

## Dependencies
In order to be able to compile the C++ core of the library the following dependencies have to be installed on your system:
* a C++17 compliant compiler. (`gcc` versions higher than 7 should be fine)
* make
* CMake
* Eigen linear algebra library (version higher than 3.3)

## Installation (R development)
fdaPDE makes use of `Rcpp` and `RcppEigen` to interface the C++ core library. To install the package for development purposes execute the following from an R console 

```
# install dependencies (execute this only once, if you not have Rcpp installed yet)
install.packages("Rcpp")
install.packages("RcppEigen")

# load Rcpp library
library(Rcpp)

setwd("/fdaPDE_root_folder/wrappers/R/")
# update RCppExports.cpp
compileAttributes(".")
# install fdaPDE
install.packages(".", type="source", repos=NULL)

# you might be required to exit from the current R session to load the new installed version of the library in case you have executed library(fdaPDE) before
```

Compilation under `RcppEigen` produces a lot of annoying warnings. Those are kept hidden by `-Wno-ignored-attributes` flag in `wrappers/R/Makevars` file. Recall that CRAN policies doesn't allow to hide such kind of warnings.
