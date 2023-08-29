#ifndef TauAnalysis_Entanglement_comp_concurrence_h
#define TauAnalysis_Entanglement_comp_concurrence_h

#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h" // math::Matrix3x3

#include <Eigen/Dense> // Eigen::

/**
 * @brief Compute Kronecker product (https://en.wikipedia.org/wiki/Kronecker_product)
 * @param A Matrix of dimensions (n, m) on the LHS of the product operator
 * @param B Matrix of dimensions (p, q) on the RHS of the product operator
 * @return Matrix of dimensions (np, mq)
 */
template <typename T, typename U>
Eigen::Matrix<typename T::Scalar,
              T::RowsAtCompileTime * U::RowsAtCompileTime,
              T::ColsAtCompileTime * U::ColsAtCompileTime>
kroneckerProduct(const Eigen::MatrixBase<T>& A,
                 const Eigen::MatrixBase<U>& B)
{
  Eigen::Matrix<typename T::Scalar,
                T::RowsAtCompileTime * U::RowsAtCompileTime,
                T::ColsAtCompileTime * U::ColsAtCompileTime> K;

  for (int i = 0; i < A.rows(); ++i)
  {
    for (int j = 0; j < A.cols(); ++j)
    {
      K.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i, j) * B;
    }
  }

  return K;
}

/**
 * @brief Computes concurrence according to Eq. 5 of the paper arXiv:2211.10513
 * @param Bplus Polarization vector of tau+ (B)
 * @param Bminus Polarization vector of tau- (Bbar)
 * @param C Spin correlation matrix
 * @param threshold Upper limit for the imaginary components of the eigenvalues
 * @param verbosity Verbosity level
 * @return Concurrence (or, in case there are errors, -1, which is an unphysical value)
 */
double
comp_concurrence(const math::Vector3& Bplus,
                 const math::Vector3& Bminus,
                 const math::Matrix3x3& C,
                 double threshold = 1e-6,
                 int verbosity = -1);

/**
 * @brief Computes concurrence according to Eq. 35 of the paper arXiv:2211.10513
 * @param C Spin correlation matrix
 * @return Concurrence, assuming only nonzero CP phase in Yukawa interactions
 */
double
comp_concurrence_simplified(const math::Matrix3x3& C);

#endif // TauAnalysis_Entanglement_comp_concurrence_h
