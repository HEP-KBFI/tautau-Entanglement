#ifndef TauAnalysis_Entanglement_printPoint_h
#define TauAnalysis_Entanglement_printPoint_h

#include "DataFormats/Candidate/interface/Candidate.h" // reco::Candidate::Point

#include <string>                                      // std::string

void
printPoint(const std::string& label,
           const reco::Candidate::Point& p);

#endif // TauAnalysis_Entanglement_printPoint_h
