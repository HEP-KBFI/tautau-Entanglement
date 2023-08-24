#ifndef TauAnalysis_Entanglement_comp_cosThetaGJ_h
#define TauAnalysis_Entanglement_comp_cosThetaGJ_h

#include "DataFormats/Candidate/interface/Candidate.h" // reco::Candidate::LorentzVector

double
comp_cosThetaGJ(const reco::Candidate::LorentzVector& tauP4, const reco::Candidate::LorentzVector& visP4, bool* errorFlag = 0, int verbosity = -1);

double
comp_cosThetaGJ_solution(const reco::Candidate::LorentzVector& tauP4, const reco::Candidate::LorentzVector& visP4, bool* errorFlag = 0, int verbosity = -1);

#endif // TauAnalysis_Entanglement_comp_cosThetaGJ_h
