#ifndef TauAnalysis_Entanglement_bookHistogram1d_h
#define TauAnalysis_Entanglement_bookHistogram1d_h

#include "PhysicsTools/FWLite/interface/TFileService.h" // fwlite::TFileService

#include <TH1.h>                                        // TH1

#include <string>                                       // std::string

TH1*
bookHistogram1d(fwlite::TFileService& fs, const std::string& name, int numBinsX, double xMin, double xMax);

#endif // TauAnalysis_Entanglement_bookHistogram1d_h
