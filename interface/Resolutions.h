#ifndef TauAnalysis_Entanglement_Resolutions_h
#define TauAnalysis_Entanglement_Resolutions_h

#include "FWCore/ParameterSet/interface/ParameterSet.h" // edm::ParameterSet

#include "DataFormats/Candidate/interface/Candidate.h"  // reco::Candidate::LorentzVector, reco::Candidate::Point

#include <TFormula.h>                                   // TFormula

class Resolutions
{
 public:
  Resolutions(const edm::ParameterSet& cfg);
  ~Resolutions();

  double
  recoilResolution_px() const;

  double
  recoilResolution_py() const;
  
  double
  recoilResolution_pz() const;

  double
  recoilResolution_mass() const;

  double
  pvResolution_xy() const;
  
  double
  pvResolution_z() const;

  const TFormula*
  trackResolution_pt() const;

  double
  trackResolution_theta() const;

  double
  trackResolution_phi() const;

  const TFormula* 
  ecalResolution_energy() const;

  double
  ecalResolution_theta() const;

  double
  ecalResolution_phi() const;
  
  double
  svResolution_parl() const;
  
  double
  svResolution_perp() const;
  
  double
  tipResolution_perp() const;

 private:
  double recoilResolution_px_;      // [GeV]
  double recoilResolution_py_;      // [GeV]
  double recoilResolution_pz_;      // [GeV]
  double recoilResolution_mass_;    // [GeV]

  double pvResolution_xy_;          // [cm]
  double pvResolution_z_;           // [cm]

  TFormula* trackResolution_pt_;    // resolution on pT in units of GeV
  double trackResolution_theta_;    // [rad]
  double trackResolution_phi_;      // [rad]

  TFormula* ecalResolution_energy_; // resolution on energy in units of GeV
  double ecalResolution_theta_;     // [rad]
  double ecalResolution_phi_;       // [rad]

  double svResolution_parl_;        // [cm]
  double svResolution_perp_;        // [cm]

  double tipResolution_perp_;       // [cm]
};

double
get_trackResolution_pt(const reco::Candidate::LorentzVector& p4, const Resolutions& resolutions);

double
get_ecalResolution_pt(const reco::Candidate::LorentzVector& p4, const Resolutions& resolutions);

#endif // TauAnalysis_Entanglement_Resolutions_h
