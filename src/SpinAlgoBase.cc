#include "TauAnalysis/Entanglement/interface/SpinAlgoBase.h"

SpinAlgoBase::SpinAlgoBase(const edm::ParameterSet& cfg)
  : verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
{}

SpinAlgoBase::~SpinAlgoBase()
{}

void
SpinAlgoBase::set_verbosity(int verbosity)
{
  verbosity_ = verbosity;
}
