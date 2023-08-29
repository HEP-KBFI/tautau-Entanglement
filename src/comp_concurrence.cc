#include "TauAnalysis/Entanglement/interface/comp_concurrence.h"

#include "TauAnalysis/Entanglement/interface/square.h" // square()

#include <Eigen/Eigenvalues> // const Eigen::ComplexEigenSolver

namespace {
  /**
   * @brief Function for retrieving Pauli matrices
   * @return Vector of Pauli matrices
   */
  const std::vector<Eigen::Matrix2cd>
  getPauliMatrices()
  {
    static const std::vector<Eigen::Matrix2cd> sigma = []{
      std::vector<Eigen::Matrix2cd> matrices(3);
      matrices[0] << 0, 1,
                     1, 0;
      matrices[1] << 0,                          std::complex<double>(0, -1),
                     std::complex<double>(0, 1), 0;
      matrices[2] << 1,  0,
                     0, -1;
      return matrices;
    }();

    return sigma;
  }

  /**
   * @brief Density operator for biparticle system, computed according to Eq. 1 of the paper arXiv:2211.10513
   * @param B Polarization vector of one particle
   * @param Bbar Polarization vector of another particle
   * @param C Spin correlation matrix
   * @return Density matrix
   */
  Eigen::Matrix4cd
  density_op(const math::Vector3& B,
             const math::Vector3& Bbar,
             const math::Matrix3x3& C)
  {
    const auto& sigma = getPauliMatrices();
    Eigen::Matrix4cd result = Eigen::Matrix4cd::Identity();

    for(int i = 0; i < 3; ++i)
    {
      result += kroneckerProduct(sigma[i], Eigen::Matrix2cd::Identity()) * B(i);
      result += kroneckerProduct(Eigen::Matrix2cd::Identity(), sigma[i]) * Bbar(i);
      for(int j = 0; j < 3; ++j)
      {
        result += kroneckerProduct(sigma[i], sigma[j]) * C(i, j);
      }
    }

    return result / 4.0;
  }

  /**
   * @brief Complementary density operator (denoted as rho-tilde in the paper arXiv:2211.10513)
   * @param rho Initial desnity matrix
   * @return Complementary density operator
   */
  Eigen::Matrix4cd
  density_op_compl(const Eigen::Matrix4cd& rho)
  {
    const auto& sigma = getPauliMatrices();
    return kroneckerProduct(sigma[1], sigma[1]) * rho.conjugate() * kroneckerProduct(sigma[1], sigma[1]);
  }

  /**
   * @brief Finds eignevalues of a complex 4x4 matrix
   * @param R Input matrix
   * @param errorBit For detecting errors
   * @param threshold Upper limit for the imaginary components of the eigenvalues
   * @param verbosity Verbosity level
   * @return Square root of the eigenvalues of the input matrix, sorted in descending order
   */
  Eigen::Vector4d
  find_sqrtEigenValues(const Eigen::Matrix4cd& R,
                       int & errorBit,
                       double threshold = 1e-6,
                       int verbosity = -1)
  {

    const Eigen::ComplexEigenSolver<Eigen::Matrix4cd> es(R);
    Eigen::Vector4cd evals = es.eigenvalues();

    for(int i = 0; i < 4; ++i)
    {
      if(std::fabs(evals(i).imag()) < threshold)
      {
        evals(i) = evals(i).real();
      }

      if(evals(i).imag() != 0)
      {
        if(verbosity >= 0)
        {
          std::cerr << "Not all eigenvalues are real: " << evals.transpose() << '\n';
        }
        errorBit = 1;
        return {};
      }

      if(evals(i).real() < 0)
      {
        if(verbosity >= 0)
        {
          std::cerr << "Eigenvalues must be non-negative: " << evals.transpose() << '\n';
        }
        errorBit = 2;
        return {};
      }
    }

    Eigen::Vector4d realEvals = evals.real();
    std::sort(realEvals.data(), realEvals.data() + realEvals.size(), std::greater<double>());
    return realEvals.cwiseSqrt();
  }
}

double
comp_concurrence(const math::Vector3& Bplus,
                 const math::Vector3& Bminus,
                 const math::Matrix3x3& C,
                 double threshold,
                 int verbosity)
{
  const Eigen::Matrix4cd rho = ::density_op(Bplus, Bminus, C);
  const Eigen::Matrix4cd rho_tilde = ::density_op_compl(rho);
  const Eigen::Matrix4cd R2 = rho * rho_tilde;
  // Note that according to arXiv:quant-ph/9709029 (just after Eq. 10),
  // instead of finding the eigenvalues of matrix R, one can achieve the same by
  // finding the eigenvalues of matrix rho * rho-tilde and then take the square
  // root of them. That's what we're also doing here.

  int errorBit = 0;
  const Eigen::Vector4d eta = ::find_sqrtEigenValues(R2, errorBit, threshold, verbosity);
  return ! errorBit ? std::max(0.0, eta(0) - eta(1) - eta(2) - eta(3)) : -1;
}

double
comp_concurrence_simplified(const math::Matrix3x3& C)
{
  // Elements ordered as: n, r, k
  const double C_rn = C(1, 0);
  const double C_nr = C(0, 1);
  const double C_rr = C(1, 1);
  const double C_nn = C(0, 0);
  const double C_kk = C(2, 2);
  const double Dplus  = std::sqrt(square(C_rn + C_nr) + square(C_rr - C_nn));
  const double Dminus = std::sqrt(square(C_rn - C_nr) + square(C_rr + C_nn));
  return std::max(0., std::max(Dplus + C_kk, Dminus - C_kk) - 1.) / 2.;
}
