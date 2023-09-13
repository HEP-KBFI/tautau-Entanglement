#include "TauAnalysis/Entanglement/interface/comp_mT.h"

#include "TauAnalysis/Entanglement/interface/square.h" // square()

#include <cmath>                                       // std::sqrt()

double
comp_mT(double particle1_mass, double particle1_px, double particle1_py,
        double particle2_mass, double particle2_px, double particle2_py)
{
  // CV: formulas for transverse energy and transverse mass
  //     taken from https://en.wikipedia.org/wiki/Transverse_mass
  //    (using pT^2 = pX^2 + pY^2 and pTvec1*pTvec2 = pX1*pX2 + pY1*pY2)
  double particle1_mass2 = square(particle1_mass);
  double particle1_Et = std::sqrt(particle1_mass2 + square(particle1_px) + square(particle1_py));
  double particle2_mass2 = square(particle2_mass);
  double particle2_Et = std::sqrt(particle2_mass2 + square(particle2_px) + square(particle2_py));
  double mT2 = particle1_mass2 + particle2_mass2
              + 2.*((particle1_Et*particle2_Et) - (particle1_px*particle2_px + particle1_py*particle2_py));
  return std::sqrt(mT2);
}
