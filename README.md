<div align="center"> <h1> fdaPDE </h1>

<h5> Physics-Informed Spatial and Functional Data Analysis </h5> </div>

![test-linux-gcc](https://img.shields.io/github/actions/workflow/status/fdaPDE/fdaPDE-cpp/test-linux-gcc.yml?branch=stable&label=test-linux-gcc)
![test-linux-clang](https://img.shields.io/github/actions/workflow/status/fdaPDE/fdaPDE-cpp/test-linux-clang.yml?branch=stable&label=test-linux-clang)
![test-macos-clang](https://img.shields.io/github/actions/workflow/status/fdaPDE/fdaPDE-cpp/test-macos-clang.yml?branch=stable&label=test-macos-clang)

fdaPDE is a C++ library for the analysis of spatial and functional data observed over complex multidimensional domains, featuring a Partial Differential Equation regularization. 

It is built on top of the [fdaPDE Core Library](https://github.com/fdaPDE/fdaPDE-core).

## Documentation
Documentation can be found on our [documentation site](https://fdapde.github.io/)

## Installation
The source code of this project can be found at [https://github.com/RiccardoSena/fdaPDE-cpp.git](https://github.com/RiccardoSena/fdaPDE-cpp.git), which is a fork of the fdaPDE repository. The prerequisites to install and run test cases are:

- A C++17 compliant compiler
- Make
- CMake
- The **Eigen** library (at least version 3.3)

### Instructions to install the library:

1. Clone the repository:
   ```bash
   git clone https://github.com/RiccardoSena/fdaPDE-cpp.git
2. Then, navigate into the develop branch and update the core submodule:
   ```bash
    cd fdaPDE-cpp
    git submodule init
    git submodule update
3. To run the test cases, use the following commands inside the fdaPDE-cpp directory:
   ```bash
    cd test
    mkdir build
    cd build
    cmake ..
    make
    cd ..
    ./run_tests.sh
It should be noted that in the test/main.cpp file, one can choose to run all tests about different models and methods. The ones that are compliant to our development are test/src/inference_test.cpp and test/src/inferencetime_test.cpp

