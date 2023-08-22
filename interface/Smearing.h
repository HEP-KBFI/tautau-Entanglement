#ifndef TauAnalysis_Entanglement_Smearing_h
#define TauAnalysis_Entanglement_Smearing_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"        // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h" // KinematicEvent
#include "TauAnalysis/Entanglement/interface/Resolutions.h"    // Resolutions

#include <TRandom.h>                                           // TRandom
#include <TRandom3.h>                                          // TRandom3

class Smearing
{
 public:
  Smearing(const edm::ParameterSet& cfg);
  ~Smearing();

  KinematicEvent
  operator()(const KinematicEvent& evt);

 private:
  reco::Candidate::Point
  smear_pv(const reco::Candidate::Point& pv);

  reco::Candidate::LorentzVector
  smear_recoil_p4(const reco::Candidate::LorentzVector& recoilP4);

  reco::Candidate::Point
  smear_sv(const reco::Candidate::LorentzVector& p4, const reco::Candidate::Point& sv);

  reco::Candidate::LorentzVector
  smear_daughter_p4(const KinematicParticle& daughter);

  reco::Candidate::Point
  smear_tipPCA(const std::vector<KinematicParticle>& daughters, const reco::Candidate::Point& tipPCA);

  //TRandom3 rnd_;
  // CV: switch to TRandom, because TRandom3::GetSeed() does not change 
  //     after each call to TRandom3::Uniform() and TRandom3::Gaus(),
  //     hampering reproducibility/debugging
  //    (TRandom is not affected by this problem)
  TRandom rnd_;
  ULong_t rndSeed_;

  bool applySmearing_recoil_px_;
  bool applySmearing_recoil_py_;
  bool applySmearing_recoil_pz_;
  bool applySmearing_recoil_mass_;

  bool applySmearing_pv_xy_;
  bool applySmearing_pv_z_;

  bool applySmearing_track_pt_;
  bool applySmearing_track_theta_;
  bool applySmearing_track_phi_;

  bool applySmearing_ecal_energy_;
  bool applySmearing_ecal_theta_;
  bool applySmearing_ecal_phi_;

  bool applySmearing_sv_perp_;
  bool applySmearing_sv_parl_;

  bool applySmearing_tip_perp_;

  Resolutions* resolutions_;

  int collider_;

  int verbosity_;
  bool cartesian_;
};

#endif // TauAnalysis_Entanglement_Smearing_h
