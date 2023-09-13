#ifndef TauAnalysis_Entanglement_showHistogram1d_h
#define TauAnalysis_Entanglement_showHistogram1d_h

#include <TH2.h>  // TH2

#include <string> // std::string

void showHistogram1d(int canvasSizeX, int canvasSizeY,
                     TH1* histogram, 
                     const std::string& xAxisTitle, double xAxisOffset, 
                     bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
                     bool showStatsBox, const std::string& drawOption,
                     const std::string& outputFileName, bool addHistogramName = true,
                     int verbosity = -1);

#endif // TauAnalysis_Entanglement_showHistogram1d_h
