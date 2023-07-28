#include "TauAnalysis/Entanglement/interface/fixMass.h"

#include "TauAnalysis/Entanglement/interface/constants.h" // mHiggs, mTau

reco::Candidate::LorentzVector
fixMass(const reco::Candidate::LorentzVector& p4, double mass)
{
  double px     = p4.px();
  double py     = p4.py();
  double pz     = p4.pz();
  double energy = std::sqrt(px*px + py*py + pz*pz + mass*mass);
  reco::Candidate::LorentzVector p4_fixed(px, py, pz, energy);
  return p4_fixed;
}

reco::Candidate::LorentzVector
fixHiggsMass(const reco::Candidate::LorentzVector& higgsP4)
{
  return fixMass(higgsP4, mHiggs);
}

reco::Candidate::LorentzVector
fixTauMass(const reco::Candidate::LorentzVector& tauP4)
{
  return fixMass(tauP4, mTau);
}

reco::Candidate::LorentzVector
fixNuMass(const reco::Candidate::LorentzVector& nuP4)
{
  return fixMass(nuP4, 0.);
}
