#include "TauAnalysis/Entanglement/interface/StartPosAlgoBase.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h" // cmsException
#include "TauAnalysis/Entanglement/interface/constants.h"    // kLHC, kSuperKEKB

StartPosAlgoBase::StartPosAlgoBase(const edm::ParameterSet& cfg)
  : algo_(cfg.getParameter<int>("algo")) 
  , resolutions_(nullptr)
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{
  edm::ParameterSet cfg_resolutions = cfg.getParameterSet("resolutions");
  resolutions_ = new Resolutions(cfg_resolutions);

  std::string collider = cfg.getParameter<std::string>("collider");
  if      ( collider == "LHC"       ) collider_ = kLHC;
  else if ( collider == "SuperKEKB" ) collider_ = kSuperKEKB;
  else throw cmsException("StartPosAlgo2", __LINE__)
    << "Invalid Configuration parameter 'collider' = " << collider << " !!\n";
}

StartPosAlgoBase::~StartPosAlgoBase()
{
  delete resolutions_;
}

int
StartPosAlgoBase::get_algo() const
{
  return algo_;
}
