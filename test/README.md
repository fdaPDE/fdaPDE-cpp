<div align="center"> <h1> fdaPDE </h1>

<h5> Testing framework </h5> </div>

</div>

This folder containes the testing framework of the fdaPDE project.

## Dependencies

To run the tests the following dependencies have to be installed on your sistem:

* a C++17 compliant compiler. (`gcc` versions higher than 7 should be fine)
* make
* CMake
* Eigen linear algebra library (version higher than 3.3)
* Google testing framework ([Google Test](http://google.github.io/googletest/))

If you want to run the entire test suite locally execute from a shell:

```Shell
cmake CMakeLists.txt
make
./fdaPDE_test
```

---

Unit tests are written following the [Google Test](http://google.github.io/googletest/) API. Please refer to its documentation to understand more about the specific API used.

## Statistical models tests

### Regression

The following is a list of all the tests performed on regression models. The R equivalent version of each test is reported for reproducibility:

| model    | description                                                                                        | link                                         |
|:---------|:---------------------------------------------------------------------------------------------------|:---------------------------------------------|
| SRPDE 2D | Non parametric model with laplacian regularization. <br />Data sampled at mesh nodes               | [Data](data/models/SRPDE/2D_test1/README.md) |
| SRPDE 2D | Semi-parametric model with laplacian regularization. <br />Data sampled at general locations       | [Data](data/models/SRPDE/2D_test2/README.md) |
| SRPDE 2D | Non parametric model with costant coefficient PDE regularization. <br />Data sampled at mesh nodes | [Data](data/models/SRPDE/2D_test3/README.md) |
| SRPDE 2D | Non parametric model with space-varying PDE regularization. <br />Areal sampling                   | [Data](data/models/SRPDE/2D_test4/README.md) |
