#include "TauAnalysis/Entanglement/interface/comp_cosThetaGJ.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h" // cmsException
#include "TauAnalysis/Entanglement/interface/square.h"       // square()

double
comp_cosThetaGJ(const reco::Candidate::LorentzVector& tauP4, const reco::Candidate::LorentzVector& visP4)
{
  //std::cout << "<comp_cosThetaGJ>:\n";
  double cosThetaGJ = tauP4.Vect().unit().Dot(visP4.Vect().unit());
  //std::cout << "cosThetaGJ = " << cosThetaGJ << "\n";
  assert(cosThetaGJ >= -1. && cosThetaGJ <= +1.);
  return cosThetaGJ;
}

double
comp_cosThetaGJ_solution(const reco::Candidate::LorentzVector& tauP4, const reco::Candidate::LorentzVector& visP4)
{
  //std::cout << "<comp_cosThetaGJ_solution>:\n";
  double cosThetaGJ_solution = (square(tauP4.P()) + square(visP4.P()) - square(tauP4.energy() - visP4.energy()))/(2.*tauP4.P()*visP4.P());
  //std::cout << "cosThetaGJ_solution = " << cosThetaGJ_solution << "\n";
  assert(cosThetaGJ_solution >= -1. && cosThetaGJ_solution <= +1.);
  return cosThetaGJ_solution;
}

