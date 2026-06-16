#!/bin/sh

set -e

# set defaults
SCRIPT_NAME=$(basename "$0")
SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
BUILD_DIR="$SCRIPT_DIR/build"
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
	rm -rf "$BUILD_DIR/CMakeCache.txt" "$BUILD_DIR/CMakeFiles"
    fi
}

## parse command line inputs
while [ "$#" -gt 0 ]; do
    case "$1" in
	-m | --memcheck )
	    MEMCHECK=true
	    shift 1
	    ;;
	-c | --compiler )
	    if [ "$#" -lt 2 ]; then
		echo "Missing compiler after $1"
		help
	    fi
	    COMPILER="$2"
	    shift 2
	    ;;
	-h | --help )
	    help
	    ;;
	*)
	    echo "Unexpected option: $1"
	    help
	    ;;
    esac
done

## set CMake compiler
if [ "$COMPILER" = "gcc" ]; then
    if command -v gcc-15 >/dev/null 2>&1 && command -v g++-15 >/dev/null 2>&1; then
	export CC=$(command -v gcc-15)
	export CXX=$(command -v g++-15)
    elif [ -x /opt/homebrew/bin/gcc-15 ] && [ -x /opt/homebrew/bin/g++-15 ]; then
	export CC=/opt/homebrew/bin/gcc-15
	export CXX=/opt/homebrew/bin/g++-15
    elif command -v gcc >/dev/null 2>&1 && command -v g++ >/dev/null 2>&1; then
	export CC=gcc
	export CXX=g++
    else
	echo "Could not find gcc/g++"
	exit 1
    fi
elif [ "$COMPILER" = "clang" ]; then
    export CC=/usr/bin/clang
    export CXX=/usr/bin/clang++
else
    echo "Unknown compiler: $COMPILER"
    help
fi

# cd into build directory
if [ -d "$BUILD_DIR" ];
then
    clean_build_dir
else
    mkdir -p "$BUILD_DIR"
fi
cd "$BUILD_DIR"

cmake -Wno-dev -S "$SCRIPT_DIR" -B .
cmake --build .

if [ "$MEMCHECK" = true ]; then
    valgrind --leak-check=full --track-origins=yes ./fdapde_test
else
    ./fdapde_test
fi

# rm fdapde_test
