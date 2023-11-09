#!/bin/sh

# cd into build directory
BUILD_DIR=build/
if [ -d "$BUILD_DIR" ];
then
    cd build/
else
    mkdir build/
    cd build/
fi

cmake -Wno-dev -H../ -B.
make
./fdapde_test
