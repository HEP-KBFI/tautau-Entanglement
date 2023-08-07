#include "TauAnalysis/Entanglement/interface/comp_cosThetaGJ.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h" // cmsException
#include "TauAnalysis/Entanglement/interface/square.h"       // square()

double
comp_cosThetaGJ(const reco::Candidate::LorentzVector& tauP4, const reco::Candidate::LorentzVector& visP4)
{
  //std::cout << "<comp_cosThetaGJ>:\n";
  //std::cout << " tauP4: P = " << tauP4.P() << "\n";
  //std::cout << " visP4: P = " << visP4.P() << "\n";
  double cosThetaGJ = 0.;
  if ( tauP4.P() > 0. && visP4.P() > 0. )
  {
    cosThetaGJ = tauP4.Vect().unit().Dot(visP4.Vect().unit());
    // CV: allow for small rounding errors
    if ( cosThetaGJ < -1.01 || cosThetaGJ > +1.01 )
    {
      std::cerr << "Error in <comp_cosThetaGJ>: cosTheta = " << cosThetaGJ << " outside physical range !!\n";
      assert(0);
    }
    if ( cosThetaGJ < -1. ) cosThetaGJ = -1.;
    if ( cosThetaGJ > +1. ) cosThetaGJ = +1.;
  }
  //std::cout << "cosThetaGJ = " << cosThetaGJ << "\n";
  return cosThetaGJ;
}

double
comp_cosThetaGJ_solution(const reco::Candidate::LorentzVector& tauP4, const reco::Candidate::LorentzVector& visP4)
{
  //std::cout << "<comp_cosThetaGJ_solution>:\n";
  //std::cout << " tauP4: energy = " << tauP4.energy() << ", P = " << tauP4.P() << " (mass = " << tauP4.mass() << ")\n";
  //std::cout << " visP4: energy = " << visP4.energy() << ", P = " << visP4.P() << " (mass = " << visP4.mass() << ")\n";
  double cosThetaGJ_solution = 0.;
  if ( tauP4.P() > 0. && visP4.P() > 0. )
  {
    cosThetaGJ_solution = (square(tauP4.P()) + square(visP4.P()) - square(tauP4.energy() - visP4.energy()))/(2.*tauP4.P()*visP4.P());
    // CV: allow for small rounding errors
    if ( cosThetaGJ_solution < -1.01 || cosThetaGJ_solution > +1.01 )
    {
      std::cerr << "Error in <comp_cosThetaGJ_solution>: cosTheta = " << cosThetaGJ_solution << " outside physical range !!\n";
      assert(0);
    }
    if ( cosThetaGJ_solution < -1. ) cosThetaGJ_solution = -1.;
    if ( cosThetaGJ_solution > +1. ) cosThetaGJ_solution = +1.;
  }
  //std::cout << "cosThetaGJ_solution = " << cosThetaGJ_solution << "\n";
  return cosThetaGJ_solution;
}

