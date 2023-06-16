#include "TauAnalysis/Entanglement/interface/SpinAnalyzer.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h" // kBeam, kHiggs

SpinAnalyzer::SpinAnalyzer(const edm::ParameterSet& cfg)
  : verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{
  std::string hAxis = cfg.getParameter<std::string>("hAxis");
  if      ( hAxis == "beam"  ) hAxis_ = kBeam;
  else if ( hAxis == "higgs" ) hAxis_ = kHiggs;
  else throw cmsException("SpinAnalyzer", __LINE__)
    << "Invalid Configuration parameter 'hAxis' = " << hAxis << " !!\n";
}

SpinAnalyzer::~SpinAnalyzer()
{}

