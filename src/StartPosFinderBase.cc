#include "TauAnalysis/Entanglement/interface/StartPosFinderBase.h"

StartPosFinderBase::StartPosFinderBase(const edm::ParameterSet& cfg)
  : verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{}

StartPosFinderBase::~StartPosFinderBase()
{}
