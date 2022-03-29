# RfdaPDE

This class is the main entry point to the core C++ library from R. 

!!! note
	fdaPDE uses Rcpp and RcppEigen to interface the core library to R. For more informations see

	* [**Rcpp: Seamless R and C++ Integration**](https://github.com/RcppCore/Rcpp)
	* [**RcppEigen: R and Eigen via Rcpp**](https://github.com/RcppCore/RcppEigen)

The fdaPDE interface can be loaded from R using the following instructions

```
library("Rcpp")
Rcpp::sourceCpp("RfdaPDE.cpp")

# create an fdaPDE interface accessible from R
fdaPDE_interface = new(RfdaPDE)

# you can now access methods using the $ notation
fdaPDE_interface$qualcosa ...

```

## Interface extension

To provide new functionalities to the outside of the C++ core library is required to add a method to the RfdaPDE class working as entry point.
