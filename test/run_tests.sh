#!/bin/sh

# set defaults
SCRIPT_NAME=$(basename "$0")
BUILD_DIR=build/
MEMCHECK=false
COMPILER="gcc"

help()
{
    echo "Usage: .run_tests.sh [options]

       -m --memcheck         use valgrind to checks for memory errors
       -c --compiler         sets compiler (gcc/clang), default gcc
       -h --help             shows this message"
    exit 2
}


clean_build_dir()
{
    if [ -d "$BUILD_DIR" ];
    then
	rm -r build/CMakeCache.txt build/CMakeFiles/
    fi
}

## parse command line inputs
SHORT=m,c:,h
LONG=memcheck,compiler:,help
OPTS=$(getopt -a --n "$SCRIPT_NAME" --options $SHORT --longoptions $LONG -- "$@") 

eval set -- "$OPTS"

while :; do
    case "$1" in
	-m | --memcheck )
	    MEMCHECK=true
	    shift 1
	    ;;
	-c | --compiler )
	    COMPILER="$2"
	    shift 2
	    ;;
	-h | --help )
	    help
	    ;;
	--)
	    shift;
	    break
	    ;; 
	*)
	    echo "Unexpected option: $1"
	    help
	    ;;
    esac
done

## set CMake compiler
if [ "$COMPILER" = "gcc" ]; then
    # find GCC in the sytem
    export CC=$(which gcc)
    export CXX=$(which g++)
elif [ "$COMPILER" = "clang" ]; then
    # find Clang in the system
    export CC=$(which clang)
    export CXX=$(which clang++)
fi


if [ -d "$BUILD_DIR" ]; then
    # If the build directory exists, check if cmake has already been run
    if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
        echo "CMake not executed. Running configuration..."
        clean_build_dir
        cmake -Wno-dev ../CMakeLists.txt
    fi
    cd build/
else
    mkdir build/
    cd build/
    cmake -Wno-dev ../CMakeLists.txt
fi


# Check if the executable already exists, if so, do not run make
if [ ! -f "./fdapde_test" ]; then
    echo "Compilation needed. Running make..."
    make
fi


if [ "$MEMCHECK" = true ]; then
    valgrind --leak-check=full --track-origins=yes ./fdapde_test
else
    ./fdapde_test
    TEST_OUTPUT=$? ## get exit code
fi

rm fdapde_test
exit $TEST_OUTPUT
