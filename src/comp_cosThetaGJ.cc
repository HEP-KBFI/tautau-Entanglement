#include "TauAnalysis/Entanglement/interface/comp_cosThetaGJ.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h" // cmsException
#include "TauAnalysis/Entanglement/interface/square.h"       // square()

double
comp_cosThetaGJ(const reco::Candidate::LorentzVector& tauP4, const reco::Candidate::LorentzVector& visP4)
{
  double cosThetaGJ = tauP4.Vect().unit().Dot(visP4.Vect().unit());
  assert(cosThetaGJ >= -1. && cosThetaGJ <= +1.);
  return cosThetaGJ;
}

double
comp_cosThetaGJ_solution(const reco::Candidate::LorentzVector& tauP4, const reco::Candidate::LorentzVector& visP4, int sign)
{
  double signFactor = 0.;
  if      ( sign == kPlusSign  ) signFactor = +1.;
  else if ( sign == kMinusSign ) signFactor = -1.;
  else throw cmsException("comp_cosThetaGJ_solution", __LINE__) 
         << "Invalid parameter 'sign' !!";
  double cosThetaGJ_solution = (square(tauP4.P()) + square(visP4.P()) - square(tauP4.energy() + signFactor*visP4.energy()))/(2.*tauP4.P()*visP4.P());
  assert(cosThetaGJ_solution >= -1. && cosThetaGJ_solution <= +1.);
  return cosThetaGJ_solution;
}

