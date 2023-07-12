#include "TauAnalysis/Entanglement/interface/printVector.h"

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
    std::cout << " pT = " << p3.r() << ", eta = " << p3.eta() << ", phi = " << p3.phi() << "\n";
  }
}
