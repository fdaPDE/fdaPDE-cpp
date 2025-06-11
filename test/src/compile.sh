#!/bin/bash

# Check number of arguments
if [ "$#" -ne 2 ]; then
  echo "Usage: ./compile <executable_name> <source_file>"
  exit 1
fi

# Assign arguments
EXEC_NAME=$1
SOURCE_FILE=$2

# Choose compiler
COMPILER="/opt/homebrew/bin/g++-14"

rm $EXEC_NAME

# Run g++ with the specified options
$COMPILER -o "$EXEC_NAME" "$SOURCE_FILE" \
  -I../.. \
  -I../../fdaPDE/core \
  -I/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 \
  -std=c++20 -march=native -O2