#define _USE_MATH_DEFINES // M_PI

#include "TauAnalysis/Entanglement/interface/comp_Rchsh.h"        // comp_Rchsh()
#include "TauAnalysis/Entanglement/interface/comp_Ek.h"           // comp_Ek()
#include "TauAnalysis/Entanglement/interface/comp_steerability.h" // comp_steerability()
#include "TauAnalysis/Entanglement/interface/comp_concurrence.h"  // comp_concurrence()

int main(int argc, char* argv[])
{
  std::cout <<
    "The following attempts to reproduce the results presented in the paper arXiv:2211.10513:\n"
    "1) CHSH inequality variable (Eq. 10)\n"
    "2) entanglement signature (Eq. 4)\n"
    "3) steerability (Eq. 21)\n"
    "4) concurrence (Eq. 5)\n\n"
    "Of the above variables it's only the concurrence that requires polarization vectors B & Bbar "
    "as input;\nthe rest hinges on the spin correlation matrix C. The variables are calculated for "
    "the following three scenarios:\n"
    "a) CP phase in H->tautau interactions, i.e., C is set to Eq. 19\n"
    "b) results extracted for the ILC case (Table III, 2nd column)\n"
    "c) results extracted for the FCC case (Table III, 3rd column)\n\n"
    "In all cases, the polarization vectors of tau leptons (Bplus, Bminus) are set to zero.\n\n"
  ;

  // Define the inputs:
  const math::Vector3 Bplus  { 0., 0., 0. }; // polarization vector of tau+
  const math::Vector3 Bminus { 0., 0., 0. }; // polarization vector of tau-

  const double delta_CP_degrees = 7; // set to zero for the SM
  const double delta_CP_radians = M_PI * delta_CP_degrees / 180.0;
  const std::vector<double> C_CP_elements {
     std::cos(2 * delta_CP_radians), std::sin(2 * delta_CP_radians),  0.0,
    -std::sin(2 * delta_CP_radians), std::cos(2 * delta_CP_radians),  0.0,
    0.0,                             0.0,                            -1.0,
  };
  const math::Matrix3x3 C_CP(C_CP_elements.cbegin(), C_CP_elements.cend());
  const std::vector<double> C_ILC_elements {
    +0.830, +0.020, -0.019,
    -0.034, +0.981, -0.029,
    -0.001, -0.021, -0.729,
  };
  const math::Matrix3x3 C_ILC(C_ILC_elements.cbegin(), C_ILC_elements.cend());
  const std::vector<double> C_FCC_elements {
    +0.925, -0.011, +0.038,
    -0.009, +0.929, -0.001,
    -0.026, -0.019, -0.879,
  };
  const math::Matrix3x3 C_FCC(C_FCC_elements.cbegin(), C_FCC_elements.cend());

  // Print 4 decimal places
  std::cout << std::fixed << std::setprecision(4);

  // Compute 1)
  const double Rchsh_CP  = comp_Rchsh(C_CP);
  const double Rchsh_ILC = comp_Rchsh(C_ILC);
  const double Rchsh_FCC = comp_Rchsh(C_FCC);

  std::cout << "1) RCHSH:\n"
            << "  CP phase: " << Rchsh_CP << " (expected: " << std::sqrt(2.) << ", see Eq. (22))\n"
            << "  ILC:      " << Rchsh_ILC << " (expected: 1.103+/-0.163, see Table III)\n"
            << "  FCC:      " << Rchsh_FCC << " (expected: 1.276+/-0.094, see Table III)\n\n";

  // Compute 2)
  const double Ek_CP_theory = 2 * std::fabs(std::cos(2 * delta_CP_radians)) + 1.; // Eq. 20
  const double Ek_CP  = comp_Ek(C_CP);
  const double Ek_ILC = comp_Ek(C_ILC);
  const double Ek_FCC = comp_Ek(C_FCC);

  std::cout << "2) Ek:\n"
            << "  CP phase: " << Ek_CP  << " (expected: " << Ek_CP_theory << ", see Eq. (20))\n"
            << "  ILC:      " << Ek_ILC << " (expected: 2.567+/-0.279, see Table III)\n"
            << "  FCC:      " << Ek_FCC << " (expected: 2.696+/-0.215, see Table III)\n\n";

  // Compute 3)
  const double S_CP  = comp_steerability(C_CP);
  const double S_ILC = comp_steerability(C_ILC);
  const double S_FCC = comp_steerability(C_FCC);

  std::cout << "3) S:\n"
            << "  CP phase: " << S_CP  << " (expected: 2)\n"
            << "  ILC:      " << S_ILC << " (expected: 1.760+/-0.161, see Table III)\n"
            << "  FCC:      " << S_FCC << " (expected: 1.851+/-0.111, see Table III)\n\n";

  // Compute 4)
  const double c_CP  = comp_concurrence(Bplus, Bminus, C_CP);
  const double c_ILC = comp_concurrence(Bplus, Bminus, C_ILC);
  const double c_FCC = comp_concurrence(Bplus, Bminus, C_FCC);

  const double cs_CP  = comp_concurrence_simplified(C_CP);
  const double cs_ILC = comp_concurrence_simplified(C_ILC);
  const double cs_FCC = comp_concurrence_simplified(C_FCC);

  std::cout << "4) C:\n"
            << "  CP phase: " << c_CP  << ", simplified " << cs_CP  << " (expected: 1)\n"
            << "  ILC:      " << c_ILC << ", simplified " << cs_ILC << " (expected: 0.778+/-0.126, see Table III)\n"
            << "  FCC:      " << c_FCC << ", simplified " << cs_FCC << " (expected: 0.871+/-0.084, see Table III)\n\n";
}
