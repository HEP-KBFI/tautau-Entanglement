#ifndef TauAnalysis_Entanglement_fixMass_h
#define TauAnalysis_Entanglement_fixMass_h

#include "DataFormats/Candidate/interface/Candidate.h"             // reco::Candidate::LorentzVector

reco::Candidate::LorentzVector
fixHiggsMass(const reco::Candidate::LorentzVector& higgsP4);

reco::Candidate::LorentzVector
fixTauMass(const reco::Candidate::LorentzVector& tauP4);

reco::Candidate::LorentzVector
fixNuMass(const reco::Candidate::LorentzVector& nuP4);

#endif // TauAnalysis_Entanglement_fixMass_h
