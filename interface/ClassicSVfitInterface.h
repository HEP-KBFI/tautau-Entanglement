#ifndef TauAnalysis_Entanglement_ClassicSVfitInterface_h
#define TauAnalysis_Entanglement_ClassicSVfitInterface_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"            // edm::ParameterSet
#include "DataFormats/Candidate/interface/Candidate.h"             // reco::Candidate::LorentzVector

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h"     // KinematicEvent
#include "TauAnalysis/Entanglement/interface/PolarimetricVector.h" // PolarimetricVector

class ClassicSVfitInterface
{
 public:
  ClassicSVfitInterface(const edm::ParameterSet& cfg);
  ~ClassicSVfitInterface();

  KinematicEvent
  operator()(const KinematicEvent& evt);

 private:
  int collider_;

  bool skip_;

  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_ClassicSVfitInterface_h
