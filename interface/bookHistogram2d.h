#ifndef TauAnalysis_Entanglement_bookHistogram2d_h
#define TauAnalysis_Entanglement_bookHistogram2d_h

#include "PhysicsTools/FWLite/interface/TFileService.h" // fwlite::TFileService

#include <TH2.h>                                        // TH2

#include <string>                                       // std::string

TH2*
bookHistogram2d(fwlite::TFileService& fs, const std::string& name, int numBinsX, double xMin, double xMax, int numBinsY, double yMin, double yMax);

#endif // TauAnalysis_Entanglement_bookHistogram2d_h
