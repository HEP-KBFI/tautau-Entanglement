#include "TauAnalysis/Entanglement/interface/getP4_rf.h"

reco::Candidate::LorentzVector
getP4_rf(const reco::Candidate::LorentzVector& p4,
         const ROOT::Math::Boost& boost)
{
  // CV: boost given four-vector to restframe
  reco::Candidate::LorentzVector p4_rf = boost(p4);
  return p4_rf;
}
