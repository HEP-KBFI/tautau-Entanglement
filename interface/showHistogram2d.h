#ifndef TauAnalysis_Entanglement_showHistogram2d_h
#define TauAnalysis_Entanglement_showHistogram2d_h

#include <TH2.h>  // TH2

#include <string> // std::string

void showHistogram2d(int canvasSizeX, int canvasSizeY,
                     TH2* histogram,
                     const std::string& xAxisTitle, double xAxisOffset, 
                     const std::string& yAxisTitle, double yAxisOffset,
                     bool showDiagonal, bool showStatsBox, const std::string& drawOption,
                     const std::string& outputFileName);

#endif // TauAnalysis_Entanglement_showHistogram2d_h
