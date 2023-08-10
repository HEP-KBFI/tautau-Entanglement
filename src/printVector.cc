#include "TauAnalysis/Entanglement/interface/printVector.h"

#include <cmath> // std::sqrt()

void
printVector(const std::string& label,
            const reco::Candidate::Vector& p3,
            bool cartesian)
{
  std::cout << label << ":";
  if ( cartesian )
  {
    std::cout << " Px = " << p3.x() << ", Py = " << p3.y() << ", Pz = " << p3.z() << "\n";
  }
  else
  {
    std::cout << " p = " << p3.r() << ", pT = " << p3.rho() << ", theta = " << p3.theta() << ", phi = " << p3.phi() << "\n";
  }
}
