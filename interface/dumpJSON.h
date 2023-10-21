#ifndef TauAnalysis_Entanglement_dumpJSON_h
#define TauAnalysis_Entanglement_dumpJSON_h

#include "TauAnalysis/Entanglement/interface/Measurement.h" // Measurement

std::string
dumpJSON(const spin::Measurement & measurement,
         double absCosTheta_cut);

#endif // TauAnalysis_Entanglement_dumpJSON_h
