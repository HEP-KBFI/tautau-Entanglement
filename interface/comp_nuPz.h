#ifndef TauAnalysis_Entanglement_comp_nuPz_h
#define TauAnalysis_Entanglement_comp_nuPz_h

#include "DataFormats/Candidate/interface/Candidate.h" // reco::Candidate::LorentzVector

double
comp_nuPz(const reco::Candidate::LorentzVector& visP4, double nuPx, double nuPy, double sign, 
          double& nu_dPzdPx, double& nu_dPzdPy,
          bool& errorFlag,
          int verbosity = -1);

reco::Candidate::LorentzVector
build_nuP4(double nuPx, double nuPy, double nuPz);

#endif // TauAnalysis_Entanglement_comp_nuPz_h
