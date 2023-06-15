#ifndef TauAnalysis_Entanglement_KinematicFit_h
#define TauAnalysis_Entanglement_KinematicFit_h

#include "DataFormats/Candidate/interface/Candidate.h"         // reco::Candidate::LorentzVector

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h" // KinematicEvent
#include "TauAnalysis/Entanglement/interface/Resolutions.h"    // Resolutions

class KinematicFit
{
 public:
  KinematicFit(const edm::ParameterSet& cfg);
  ~KinematicFit();

  KinematicEvent
  operator()(const KinematicEvent& evt);

 private:
  void
  findStartPosition(const KinematicEvent& evt);

  reco::Candidate::LorentzVector tauPlusP4_;
  reco::Candidate::LorentzVector tauMinusP4_;

  Resolutions* resolutions_;
};

#endif // TauAnalysis_Entanglement_KinematicFit_h
