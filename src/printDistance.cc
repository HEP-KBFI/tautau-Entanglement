#include "TauAnalysis/Entanglement/interface/printDistance.h"

#include <cmath> // std::sqrt()

void
printDistance(const std::string& label,
              const reco::Candidate::Vector& p3,
              bool cartesian)
{
  std::cout << label << ":";
  if ( cartesian )
  {
    std::cout << " dx = " << p3.x() << ", dy = " << p3.y() << ", dz = " << p3.z() << "\n";
  }
  else
  {
    std::cout << " d = " << std::sqrt(p3.mag2()) << ", theta = " << p3.theta() << ", phi = " << p3.phi() << "\n";
  }
}
