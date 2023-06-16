#ifndef TauAnalysis_Entanglement_Smearing_h
#define TauAnalysis_Entanglement_Smearing_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"        // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h" // KinematicEvent
#include "TauAnalysis/Entanglement/interface/Resolutions.h"    // Resolutions

#include <TRandom3.h>                                          // TRandom3

class Smearing
{
 public:
  Smearing(const edm::ParameterSet& cfg);
  ~Smearing();

  KinematicEvent
  operator()(const KinematicEvent& evt);

 private:
  KinematicParticle
  smear_daughter(const KinematicParticle& daughter);

  reco::Candidate::Point
  smear_tipPCA(const std::vector<KinematicParticle>& daughters, const reco::Candidate::Point& tipPCA);

  reco::Candidate::Point
  smear_sv(const reco::Candidate::LorentzVector& p4, const reco::Candidate::Point& sv);

  TRandom3 rnd_;

  bool applySmearing_recoilPx_;
  bool applySmearing_recoilPy_;
  bool applySmearing_recoilPz_;
  bool applySmearing_recoilE_;
  bool applySmearing_pvXY_;
  bool applySmearing_pvZ_;
  bool applySmearing_svPerp_;
  bool applySmearing_svParl_;
  bool applySmearing_tipPerp_;

  Resolutions* resolutions_;
};

#endif // TauAnalysis_Entanglement_Smearing_h
