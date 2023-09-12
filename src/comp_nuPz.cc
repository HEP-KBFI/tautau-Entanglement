#include "TauAnalysis/Entanglement/interface/comp_nuPz.h"

#include "TauAnalysis/Entanglement/interface/constants.h" // mTau
#include "TauAnalysis/Entanglement/interface/square.h"    // square()

double
comp_nuPz(const reco::Candidate::LorentzVector& visP4, double nuPx, double nuPy, double sign, 
          double& nu_dPzdPx, double& nu_dPzdPy,
          int verbosity)
{
  if ( verbosity >= 3 )
  {
    std::cout << "<comp_nuPz>:\n";
  }

  double visPx    = visP4.px();
  double visPy    = visP4.py();
  double visPz    = visP4.pz();
  double visE     = visP4.energy();
  double visMass  = ( visMass < mTau ) ? visP4.mass() : mTau;
  double visMass2 = square(visMass);
  double mTau2    = square(mTau); 
  double term1    = visMass2 + square(visPx) + square(visPy);
  double term2    = (mTau2 - visMass2)*visPz + 2.*(nuPx*visPx + nuPy*visPy)*visPz;
  double term3    = square(mTau2) + square(visMass2) - 4.*square(visPx*nuPy - nuPx*visPy) 
                   + mTau2*(-2.*visMass2 + 4.*nuPx*visPx + 4.*nuPy*visPy) 
                   - 4.*visMass2*(nuPx*(nuPx + visPx) + nuPy*(nuPy + visPy));
  if ( verbosity >= 3 )
  {
    std::cout << "term1 = " << term1 << "\n";
    std::cout << "term2 = " << term2 << "\n";
    std::cout << "term3 = " << term3 << "\n";
  }

  double nuPz = (1./(2.*term1))*(term2 + sign*visE*sqrt(std::max(0., term3)));
  if ( verbosity >= 3 )
  {
    std::cout << "nuPz = " << nuPz << "\n";
  }

  nu_dPzdPx = (1./term1)*(visPx*visPz
             + sign*(visE/std::max(1.e-2, std::sqrt(term3)))*(mTau2*visPx - visMass2*(2.*nuPx + visPx) + 2.*visPx*nuPy*visPy - 2.*nuPx*square(visPy)));
  nu_dPzdPy = (1./term1)*(visPy*visPz
             + sign*(visE/std::max(1.e-2, std::sqrt(term3)))*(mTau2*visPy - visMass2*(2.*nuPy + visPy) + 2.*visPy*nuPx*visPx - 2.*nuPy*square(visPx)));
  if ( verbosity >= 3 )
  {
    std::cout << "nu_dPzdPx = " << nu_dPzdPx << "\n";
    std::cout << "nu_dPzdPy = " << nu_dPzdPy << "\n";
  }

  return nuPz;
}

reco::Candidate::LorentzVector
build_nuP4(double nuPx, double nuPy, double nuPz)
{
  double nuE = std::sqrt(square(nuPx) + square(nuPy) + square(nuPz));
  reco::Candidate::LorentzVector nuP4(nuPx, nuPy, nuPz, nuE);
  return nuP4;
}
