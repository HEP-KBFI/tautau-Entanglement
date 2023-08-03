#ifndef TauAnalysis_Entanglement_showHistogram1d_h
#define TauAnalysis_Entanglement_showHistogram1d_h

#include <TH2.h>  // TH2

#include <string> // std::string

void showHistogram1d(TH1* histogram, 
                     const std::string& xAxisTitle, double xAxisOffset, 
                     bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
                     double avEvtWeight,
                     bool showStatsBox,
                     const std::string& outputFileName);

#endif // TauAnalysis_Entanglement_showHistogram1d_h
