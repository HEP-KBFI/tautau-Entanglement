#include "TauAnalysis/Entanglement/interface/Resolutions.h"

#include "TauAnalysis/Entanglement/interface/square.h" // square()

#include <cmath>                                       // std::sqrt()

Resolutions::Resolutions(const edm::ParameterSet& cfg)
  : recoilResolution_px_(cfg.getParameter<double>("recoilResolution_px"))
  , recoilResolution_py_(cfg.getParameter<double>("recoilResolution_py"))
  , recoilResolution_pz_(cfg.getParameter<double>("recoilResolution_pz"))
  , recoilResolution_energy_(cfg.getParameter<double>("recoilResolution_energy"))
  , pvResolution_xy_(cfg.getParameter<double>("pvResolution_xy"))
  , pvResolution_z_(cfg.getParameter<double>("pvResolution_z"))
  , trackResolution_pt_(cfg.getParameter<double>("trackResolution_pt"))
  , trackResolution_theta_(cfg.getParameter<double>("trackResolution_theta"))
  , trackResolution_phi_(cfg.getParameter<double>("trackResolution_phi"))
  , ecalResolution_energy_a_(cfg.getParameter<double>("ecalResolution_energy_a"))
  , ecalResolution_energy_b_(cfg.getParameter<double>("ecalResolution_energy_b"))
  , ecalResolution_theta_(cfg.getParameter<double>("ecalResolution_theta"))
  , ecalResolution_phi_(cfg.getParameter<double>("ecalResolution_phi"))
  , svResolution_parl_(cfg.getParameter<double>("svResolution_parl"))
  , svResolution_perp_(cfg.getParameter<double>("svResolution_perp"))
  , tipResolution_perp_(cfg.getParameter<double>("tipResolution_perp"))
{}

Resolutions::~Resolutions()
{}

double
Resolutions::get_recoilResolution_px() const
{
  return recoilResolution_px_;
}

double
Resolutions::get_recoilResolution_py() const
{
  return recoilResolution_py_;
}
  
double
Resolutions::get_recoilResolution_pz() const
{
  return recoilResolution_pz_;
}

double
Resolutions::get_recoilResolution_energy() const
{
  return recoilResolution_energy_;
}

double
Resolutions::get_pvResolution_xy() const
{
  return pvResolution_xy_;
}
  
double
Resolutions::get_pvResolution_z() const
{
  return pvResolution_z_;
}
  
double
Resolutions::get_trackResolution_pt() const
{
  return trackResolution_pt_;
}

double
Resolutions::get_trackResolution_theta() const
{
  return trackResolution_theta_;
}

double
Resolutions::get_trackResolution_phi() const
{
  return trackResolution_phi_;
}

double
Resolutions::get_ecalResolution_energy_a() const
{
  return ecalResolution_energy_a_;
}

double
Resolutions::get_ecalResolution_energy_b() const
{
  return ecalResolution_energy_b_;
}

double
Resolutions::get_ecalResolution_theta() const
{
  return ecalResolution_theta_;
}

double
Resolutions::get_ecalResolution_phi() const
{
  return ecalResolution_phi_;
}

double
Resolutions::get_svResolution_parl() const
{
  return svResolution_parl_;
}
  
double
Resolutions::get_svResolution_perp() const
{
  return svResolution_perp_;
}
  
double
Resolutions::get_tipResolution_perp() const
{
  return tipResolution_perp_;
}  

double
get_trackResolution_pt(const reco::Candidate::LorentzVector& p4, const Resolutions& resolutions)
{
  return square(p4.pt())*resolutions.get_trackResolution_pt();
}

double
get_ecalResolution_pt(const reco::Candidate::LorentzVector& p4, const Resolutions& resolutions)
{
  double ecalResolution2_energy = square(std::sqrt(p4.energy())*resolutions.get_ecalResolution_energy_a())
                                 + square(p4.energy()*resolutions.get_ecalResolution_energy_b());
  return (p4.pt()/p4.energy())*std::sqrt(ecalResolution2_energy);
}
