#include "TauAnalysis/Entanglement/interface/StartPosAlgoBase.h"

StartPosAlgoBase::StartPosAlgoBase(const edm::ParameterSet& cfg)
  : verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{}

StartPosAlgoBase::~StartPosAlgoBase()
{}
