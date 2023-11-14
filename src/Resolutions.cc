#include "TauAnalysis/Entanglement/interface/Resolutions.h"

#include "TauAnalysis/Entanglement/interface/square.h" // square()

#include <TString.h>                                   // TString

#include <cmath>                                       // std::sqrt()

Resolutions::Resolutions(const edm::ParameterSet& cfg)
  : recoilResolution_px_(cfg.getParameter<double>("recoilResolution_px"))
  , recoilResolution_py_(cfg.getParameter<double>("recoilResolution_py"))
  , recoilResolution_pz_(cfg.getParameter<double>("recoilResolution_pz"))
  , recoilResolution_mass_(cfg.getParameter<double>("recoilResolution_mass"))
  , pvResolution_xy_(cfg.getParameter<double>("pvResolution_xy"))
  , pvResolution_z_(cfg.getParameter<double>("pvResolution_z"))
  , trackResolution_pt_(nullptr)
  , trackResolution_theta_(cfg.getParameter<double>("trackResolution_theta"))
  , trackResolution_phi_(cfg.getParameter<double>("trackResolution_phi"))
  , ecalResolution_energy_(nullptr)
  , ecalResolution_theta_(cfg.getParameter<double>("ecalResolution_theta"))
  , ecalResolution_phi_(cfg.getParameter<double>("ecalResolution_phi"))
  , svResolution_parl_(cfg.getParameter<double>("svResolution_parl"))
  , svResolution_perp_(cfg.getParameter<double>("svResolution_perp"))
  , tipResolution_perp_(cfg.getParameter<double>("tipResolution_perp"))
{
  TString trackResolution_pt = cfg.getParameter<std::string>("trackResolution_pt").c_str();
  trackResolution_pt = trackResolution_pt.ReplaceAll("pT", "x");
  trackResolution_pt = trackResolution_pt.ReplaceAll("beta", "y");
  trackResolution_pt_ = new TFormula("trackResolution_pt", trackResolution_pt.Data());
  
  TString ecalResolution_energy = cfg.getParameter<std::string>("ecalResolution_energy").c_str();
  ecalResolution_energy = ecalResolution_energy.ReplaceAll("E", "x");
  ecalResolution_energy_ = new TFormula("ecalResolution_energy", ecalResolution_energy.Data());
}

Resolutions::~Resolutions()
{
  delete trackResolution_pt_;
  delete ecalResolution_energy_;
}

double
Resolutions::recoilResolution_px() const
{
  return recoilResolution_px_;
}

double
Resolutions::recoilResolution_py() const
{
  return recoilResolution_py_;
}
  
double
Resolutions::recoilResolution_pz() const
{
  return recoilResolution_pz_;
}

double
Resolutions::recoilResolution_mass() const
{
  return recoilResolution_mass_;
}

double
Resolutions::pvResolution_xy() const
{
  return pvResolution_xy_;
}
  
double
Resolutions::pvResolution_z() const
{
  return pvResolution_z_;
}
  
const TFormula*
Resolutions::trackResolution_pt() const
{
  return trackResolution_pt_;
}

double
Resolutions::trackResolution_theta() const
{
  return trackResolution_theta_;
}

double
Resolutions::trackResolution_phi() const
{
  return trackResolution_phi_;
}

const TFormula*
Resolutions::ecalResolution_energy() const
{
  return ecalResolution_energy_;
}

double
Resolutions::ecalResolution_theta() const
{
  return ecalResolution_theta_;
}

double
Resolutions::ecalResolution_phi() const
{
  return ecalResolution_phi_;
}

double
Resolutions::svResolution_parl() const
{
  return svResolution_parl_;
}
  
double
Resolutions::svResolution_perp() const
{
  return svResolution_perp_;
}
  
double
Resolutions::tipResolution_perp() const
{
  return tipResolution_perp_;
}  

double
get_trackResolution_pt(const reco::Candidate::LorentzVector& p4, const Resolutions& resolutions)
{
  const TFormula* formula = resolutions.trackResolution_pt();
  assert(formula);
  double pT = p4.pt();
  double gamma = p4.energy()/p4.mass();
  double beta = std::sqrt((square(gamma) - 1)/square(gamma));
  double resolution_pT = formula->Eval(pT, beta);
  return resolution_pT;
}

double
get_ecalResolution_pt(const reco::Candidate::LorentzVector& p4, const Resolutions& resolutions)
{
  const TFormula* formula = resolutions.ecalResolution_energy();
  assert(formula);
  double E = p4.energy();
  double resolution_E = formula->Eval(E);
  double resolution_pT = (p4.pt()/p4.energy())*resolution_E;
  return resolution_pT;
}
