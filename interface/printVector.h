#ifndef TauAnalysis_Entanglement_printVector_h
#define TauAnalysis_Entanglement_printVector_h

#include "DataFormats/Candidate/interface/Candidate.h" // reco::Candidate::Vector

#include <string>                                      // std::string

void
printVector(const std::string& label,
            const reco::Candidate::Vector& p3,
            bool cartesian = true);

#endif // TauAnalysis_Entanglement_printVector_h
