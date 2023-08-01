#include "TauAnalysis/Entanglement/interface/printLorentzVector.h"

void
printLorentzVector(const std::string& label,
                   const reco::Candidate::LorentzVector& p4,
                   bool cartesian)
{
  std::cout << label << ":";
  if ( cartesian )
  {
    std::cout << " E = " << p4.energy() << ", Px = " << p4.px() << ", Py = " << p4.py() << ", Pz = " << p4.pz() << "\n";
  }
  else
  {
    std::cout << " pT = " << p4.pt() << ", theta = " << p4.theta() << ", phi = " << p4.phi() << ", mass = " << p4.mass() << "\n";
  }
}
