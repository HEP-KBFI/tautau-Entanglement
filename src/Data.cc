#include "TauAnalysis/Entanglement/interface/Data.h"

using namespace spin;

Data::Data(float hPlus_n, float hPlus_r, float hPlus_k, 
           float hMinus_n, float hMinus_r, float hMinus_k,
           float evtWeight)
  : hPlus_n_(hPlus_n)
  , hPlus_r_(hPlus_r)
  , hPlus_k_(hPlus_k)
  , hMinus_n_(hMinus_n)
  , hMinus_r_(hMinus_r)
  , hMinus_k_(hMinus_k)
  , evtWeight_(evtWeight)
{}

Data::Data(const Data& data,
           float evtWeight)
  : hPlus_n_(data.hPlus_n_)
  , hPlus_r_(data.hPlus_r_)
  , hPlus_k_(data.hPlus_k_)
  , hMinus_n_(data.hMinus_n_)
  , hMinus_r_(data.hMinus_r_)
  , hMinus_k_(data.hMinus_k_)
  , evtWeight_(evtWeight)  
{}

Data::~Data()
{}
