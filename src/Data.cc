#include "TauAnalysis/Entanglement/interface/Data.h"

using namespace spin;

Data::Data(float hPlus_r, float hPlus_n, float hPlus_k, 
           float hMinus_r, float hMinus_n, float hMinus_k,
           float evtWeight)
  : hPlus_r_(hPlus_r)
  , hPlus_n_(hPlus_n)
  , hPlus_k_(hPlus_k)
  , hMinus_r_(hMinus_r)
  , hMinus_n_(hMinus_n)
  , hMinus_k_(hMinus_k)
  , evtWeight_(evtWeight)
{}

Data::~Data()
{}
