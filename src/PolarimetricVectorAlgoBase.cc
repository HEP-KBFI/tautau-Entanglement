#include "TauAnalysis/Entanglement/interface/PolarimetricVectorAlgoBase.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/constants.h"                 // kLHC, kSuperKEKB
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h" // kBeam, kHiggs

PolarimetricVectorAlgoBase::PolarimetricVectorAlgoBase(const edm::ParameterSet& cfg)
  : verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{
  std::string hAxis = cfg.getParameter<std::string>("hAxis");
  if      ( hAxis == "beam"  ) hAxis_ = kBeam;
  else if ( hAxis == "higgs" ) hAxis_ = kHiggs;
  else throw cmsException("PolarimetricVectorAlgoBase", __LINE__)
    << "Invalid Configuration parameter 'hAxis' = " << hAxis << " !!\n";

  std::string collider = cfg.getParameter<std::string>("collider");
  if      ( collider == "LHC"       ) collider_ = kLHC;
  else if ( collider == "SuperKEKB" ) collider_ = kSuperKEKB;
  else throw cmsException("PolarimetricVectorAlgoBase", __LINE__)
    << "Invalid Configuration parameter 'collider' = " << collider << " !!\n";
}

PolarimetricVectorAlgoBase::~PolarimetricVectorAlgoBase()
{}
