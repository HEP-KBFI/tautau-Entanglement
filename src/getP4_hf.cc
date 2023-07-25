#include "TauAnalysis/Entanglement/interface/getP4_hf.h"

#include "TauAnalysis/Entanglement/interface/printVector.h" // ONLY FOR TESTING !!

reco::Candidate::LorentzVector
getP4_hf(const reco::Candidate::LorentzVector& p4,
         const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k)
{
  // CV: rotate given four-vector to helicity frame
  reco::Candidate::Vector p3 = p4.Vect();
printVector("p3", p3); 
  double Pr = p3.Dot(r);
  double Pn = p3.Dot(n);
  double Pk = p3.Dot(k);
std::cout << "Pr = " << Pr << ", Pn = " << Pn << ", Pk = " << Pk << "\n";
  reco::Candidate::LorentzVector p4_hf(Pr, Pn, Pk, p4.energy());
  return p4_hf;
}
