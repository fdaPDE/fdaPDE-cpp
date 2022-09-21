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
ctest --output-on-failure
```

---

Unit tests are written following the [Google Test](http://google.github.io/googletest/) API. Please refer to its documentation to understand more about the specific API used.
