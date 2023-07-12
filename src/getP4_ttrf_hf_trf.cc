#include "TauAnalysis/Entanglement/interface/getP4_ttrf_hf_trf.h"

#include "TauAnalysis/Entanglement/interface/getP4_rf.h" // getP4_rf()
#include "TauAnalysis/Entanglement/interface/getP4_hf.h" // getP4_hf()

reco::Candidate::LorentzVector
getP4_ttrf_hf_trf(const reco::Candidate::LorentzVector& p4,
                  const ROOT::Math::Boost& boost_ttrf,
                  const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k,
                  const ROOT::Math::Boost& boost_trf)
{
  // CV: boost given four-vector to restframe of tau pair,
  //     rotate to helicity frame,
  //     and finally boost to tau restframe
  reco::Candidate::LorentzVector p4_ttrf = getP4_rf(p4, boost_ttrf);
  reco::Candidate::LorentzVector p4_hf = getP4_hf(p4_ttrf, r, n, k);
  reco::Candidate::LorentzVector p4_trf = getP4_rf(p4_hf, boost_trf);
  return p4_trf;
}
