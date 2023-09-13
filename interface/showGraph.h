#ifndef TauAnalysis_Entanglement_showGraph_h
#define TauAnalysis_Entanglement_showGraph_h

#include <TGraph.h> // TGraph

#include <string>   // std::string

void
showGraph(double canvasSizeX, double canvasSizeY,
          TGraph* graph,
          const std::string& xAxisTitle, double xMin, double xMax,
          const std::string& yAxisTitle, double yMin, double yMax,
          const std::string& drawingOption,
          const std::string& outputFileName,
          int verbosity = -1);

#endif // TauAnalysis_Entanglement_showGraph_h
