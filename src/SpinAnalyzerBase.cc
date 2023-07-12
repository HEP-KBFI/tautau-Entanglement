#include "TauAnalysis/Entanglement/interface/SpinAnalyzerBase.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h" // kBeam, kHiggs

SpinAnalyzerBase::SpinAnalyzerBase(const edm::ParameterSet& cfg)
  : verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{
  std::string hAxis = cfg.getParameter<std::string>("hAxis");
  if      ( hAxis == "beam"  ) hAxis_ = kBeam;
  else if ( hAxis == "higgs" ) hAxis_ = kHiggs;
  else throw cmsException("SpinAnalyzerBase", __LINE__)
    << "Invalid Configuration parameter 'hAxis' = " << hAxis << " !!\n";
}

SpinAnalyzerBase::~SpinAnalyzerBase()
{}
