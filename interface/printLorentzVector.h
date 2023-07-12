#ifndef TauAnalysis_Entanglement_printLorentzVector_h
#define TauAnalysis_Entanglement_printLorentzVector_h

#include "DataFormats/Candidate/interface/Candidate.h" // reco::Candidate::LorentzVector

#include <string>                                      // std::string

void
printLorentzVector(const std::string& label,
                   const reco::Candidate::LorentzVector& p4,
                   bool cartesian = true);

#endif // TauAnalysis_Entanglement_printLorentzVector_h
