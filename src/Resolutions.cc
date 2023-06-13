#include "TauAnalysis/Entanglement/interface/Resolutions.h"

Resolutions::Resolutions(const edm::ParameterSet& cfg)
  : recoilResolutionPx_(cfg.getParameter<double>("recoilResolutionPx"))
  , recoilResolutionPy_(cfg.getParameter<double>("recoilResolutionPy"))
  , recoilResolutionPz_(cfg.getParameter<double>("recoilResolutionPz"))
  , recoilResolutionE_(cfg.getParameter<double>("recoilResolutionE"))
  , pvResolutionXY_(cfg.getParameter<double>("pvResolutionXY"))
  , pvResolutionZ_(cfg.getParameter<double>("pvResolutionZ"))
  , svResolutionParl_(cfg.getParameter<double>("svResolutionParl"))
  , svResolutionPerp_(cfg.getParameter<double>("svResolutionPerp"))
  , tipResolutionPerp_(cfg.getParameter<double>("tipResolutionPerp"))
{}

Resolutions::~Resolutions()
{}

double
Resolutions::get_recoilResolutionPx() const
{
  return recoilResolutionPx_;
}

double
Resolutions::get_recoilResolutionPy() const
{
  return recoilResolutionPy_;
}
  
double
Resolutions::get_recoilResolutionPz() const
{
  return recoilResolutionPz_;
}

double
Resolutions::get_recoilResolutionE() const
{
  return recoilResolutionE_;
}

double
Resolutions::get_pvResolutionXY() const
{
  return pvResolutionXY_;
}
  
double
Resolutions::get_pvResolutionZ() const
{
  return pvResolutionZ_;
}
  
double
Resolutions::get_svResolutionParl() const
{
  return svResolutionParl_;
}
  
double
Resolutions::get_svResolutionPerp() const
{
  return svResolutionPerp_;
}
  
double
Resolutions::get_tipResolutionPerp() const
{
  return tipResolutionPerp_;
}  
