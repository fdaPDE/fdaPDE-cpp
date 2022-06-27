#ifndef __BLOCK_MATRIX_H__
#define __BLOCK_MATRIX_H__

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <iostream>
#include <tuple>
#include <variant>

// this class is heavily based on C++17 std::variant<> capabilities.

// a class allowing to build NxM block matrices with both dense and sparse
// blocks. This allows some notational convenience to translate formal writing into software
// Observe that N and M are respectively the number of blocks by rows and columns.
template <typename... E> struct BlockMatrix {
  // data container, we store references to avoid expensive copies
  std::tuple<E&...> data_;
  
  BlockMatrix() = default;                         // initialize empty block matrix
  explicit BlockMatrix(std::tuple<E&...> data) :   // build block matrix from a tuple
    data_(data) {};

  // access block at position I. I must be known at compile time
  template <unsigned int I>
  decltype(std::get<I>(data_)) get() { return std::get<I>(data_); };
};

// comma initialization
template <typename B1, typename B2> // base case
BlockMatrix<B1, B2> operator,(B1 &op1, B2 &op2) {
  return BlockMatrix(std::tuple_cat(std::tie(op1), std::tie(op2)));
}

template<typename... B1, typename B2> // recursive step
BlockMatrix<B1..., B2> operator,(BlockMatrix<B1...>&& op1, B2& op2){
  return BlockMatrix(std::tuple_cat(op1.data_, std::tie(op2)));
}

#endif // __BLOCK_MATRIX_H__
