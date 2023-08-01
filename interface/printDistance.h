#ifndef TauAnalysis_Entanglement_printDistance_h
#define TauAnalysis_Entanglement_printDistance_h

#include "DataFormats/Candidate/interface/Candidate.h" // reco::Candidate::Vector

#include <string>                                      // std::string

void
printDistance(const std::string& label,
              const reco::Candidate::Vector& p3,
              bool cartesian = true);

#endif // TauAnalysis_Entanglement_printDistance_h
