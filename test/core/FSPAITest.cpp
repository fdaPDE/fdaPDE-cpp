#include <gtest/gtest.h> // testing framework
#include <unsupported/Eigen/SparseExtra>

#include "../../fdaPDE/core/NLA/FSPAI.h"
using fdaPDE::core::NLA::FSPAI;    

// test FSPAI implementation using different matrices taken from matrix market with respect to output coming from
// fspai-1.1, Matous Sedlacek https://www5.in.tum.de/wiki/index.php/FSPAI.
//      output is obtained executing ./fspai-1.1 <matrix-file> -diag -ep 0.005 -mn 3 -ns 150

// the idea of the test is that even if matrices are different (we do not expect to get **exactly** the same values between the two
// implementations) they should in any way produce results approximatively of the same quality.
class FSPAITest : public ::testing::Test {
protected:
  Eigen::SparseMatrix<double> input_{}; // input matrix
  Eigen::SparseMatrix<double> reference_{}; // reference approximate inverse of the cholesky factor of input_

  Eigen::SparseMatrix<double> L_{}; // cholesky factor of input_

  // if L is the cholesky factor of the target matrix, a measure of the goodness for an approximate inverse to be considered as an inverse of L
  // is the norm of the quantity L*inv - I. \norm{L*inv - I} \approx 0 \iff L*inv \approx I \iff inv \approx L^{-1}. To measure the goodness of
  // our implementation we consider the ratio |\norm{L*inv - I}/\norm{L*reference_ - I}| and ask that this ratio is equal to 1 up to tolerance_
  double tolerance_ = std::pow(0.1, 4);
  
  FSPAITest() = default;
  // load matrices from files
  void loadData(const std::string& input, const std::string& reference){
    Eigen::loadMarket(input_, input);
    // matrix market format only stores the lower triangular part. Symmetrize the matrix here
    Eigen::SparseMatrix<double> tmp = input_.transpose().triangularView<Eigen::StrictlyUpper>();
    input_ = input_ + tmp;

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > choleskySolver;
    choleskySolver.analyzePattern(this->input_);
    choleskySolver.factorize(this->input_);
    // recover cholesky factor of input matrix
    L_ = choleskySolver.matrixL();
        
    Eigen::loadMarket(reference_, reference);
    return;
  }
};

TEST_F(FSPAITest, ex33) {
  // load matrix from files
  this->loadData("data/ex33.mtx", "data/ex33-fspai.mtx");

  // compute sparse approximate inverse
  FSPAI fspai(this->input_);
  fspai.compute(150, 3, 0.005);
  // recover inverse of cholesky factor of L  
  Eigen::SparseMatrix<double> invL = fspai.getL();

  // build identity matrix
  Eigen::SparseMatrix<double> I;
  I.resize(this->input_.rows(), this->input_.cols());
  I.setIdentity();

  EXPECT_NEAR(std::abs((this->L_*invL - I).norm()/(this->L_*this->reference_ - I).norm()), 1, this->tolerance_);
}

TEST_F(FSPAITest, msc01440) {
  // load matrix from files
  this->loadData("data/msc01440.mtx", "data/msc01440-fspai.mtx");
  // compute sparse approximate inverse
  FSPAI fspai(this->input_);
  fspai.compute(150, 3, 0.005);
  // recover inverse of cholesky factor of L  
  Eigen::SparseMatrix<double> invL = fspai.getL();

  // build identity matrix
  Eigen::SparseMatrix<double> I;
  I.resize(this->input_.rows(), this->input_.cols());
  I.setIdentity();

  EXPECT_NEAR(std::abs((this->L_*invL - I).norm()/(this->L_*this->reference_ - I).norm()), 1, this->tolerance_);
}
