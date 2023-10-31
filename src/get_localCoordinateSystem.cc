#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h" // cmsException
#include "TauAnalysis/Entanglement/interface/constants.h"    // kLHC, kSuperKEKB, mElectron, mProton
#include "TauAnalysis/Entanglement/interface/getP4_rf.h"     // getP4_rf()
#include "TauAnalysis/Entanglement/interface/printVector.h"  // printVector()
#include "TauAnalysis/Entanglement/interface/square.h"       // square()

namespace
{
  reco::Candidate::Vector
  get_k(const reco::Candidate::LorentzVector& p4, const ROOT::Math::Boost* boost_ttrf, int verbosity, bool cartesian)
  {
    reco::Candidate::Vector k;
    if ( boost_ttrf )
    {
      reco::Candidate::LorentzVector p4_ttrf = getP4_rf(p4, *boost_ttrf);
      k = p4_ttrf.Vect().unit();
    }
    else
    {
      k = p4.Vect().unit();
    }
    if ( verbosity >= 3 )
    {
      printVector("k", k, cartesian);
    }
    return k;
  }

  reco::Candidate::Vector
  get_h_beamAxis(const ROOT::Math::Boost* boost_ttrf, int collider, int verbosity, bool cartesian)
  {
    double beamE, mBeamParticle;
    if ( collider == kLHC )
    {
      beamE = beamEnergy_LHC;
      mBeamParticle = mProton;
    }
    else if ( collider == kSuperKEKB )
    {
      // CV: electron beam defines +z direction
      beamE = beamEnergy_SuperKEKB_eMinus;
      mBeamParticle = mElectron;
    } else throw cmsException("get_h_beamAxis", __LINE__)
        << "Invalid function argument 'collider' = " << collider << " !!\n";
    const double beamPx = 0.;
    const double beamPy = 0.;
    const double beamPz = std::sqrt(square(beamE) - square(mBeamParticle));
    reco::Candidate::LorentzVector beamP4(beamPx, beamPy, beamPz, beamE);
    reco::Candidate::Vector h;
    if ( boost_ttrf )
    {
      reco::Candidate::LorentzVector beamP4_ttrf = getP4_rf(beamP4, *boost_ttrf);
      h = beamP4_ttrf.Vect().unit();
    }
    else
    {
      h = beamP4.Vect().unit();
    }
    if ( verbosity >= 3 )
    {
      printVector("h", h, cartesian);
    }
    return h;
  }

  reco::Candidate::Vector
  get_h_higgsAxis(const reco::Candidate::LorentzVector& recoilP4, const ROOT::Math::Boost& boost_ttrf, int verbosity, bool cartesian)
  {
    // CV: this code does not depend on the assumption that the tau pair originates from a Higgs boson decay;
    //     it also works for tau pairs originating from Z/gamma* -> tau+ tau- decays
    const double sf = 1.01;
    double higgsPx = sf*recoilP4.px();
    double higgsPy = sf*recoilP4.py();
    double higgsPz = sf*recoilP4.pz();
    double higgsE = std::sqrt(square(higgsPx) + square(higgsPy) + square(higgsPz) + square(recoilP4.mass()));
    reco::Candidate::LorentzVector higgsP4(higgsPx, higgsPy, higgsPz, higgsE);
    reco::Candidate::LorentzVector higgsP4_ttrf = getP4_rf(higgsP4, boost_ttrf);
    reco::Candidate::Vector h = higgsP4_ttrf.Vect().unit();
    if ( verbosity >= 3 )
    {
      printVector("h", h, cartesian);
    }
    return h;
  }
}

reco::Candidate::Vector
get_r(const reco::Candidate::Vector& k, const reco::Candidate::Vector& h, int verbosity, bool cartesian)
{
  double cosTheta = k.Dot(h);
  // CV: allow for small rounding errors
  if ( cosTheta < -1.01 || cosTheta > +1.01 )
  {
    std::cerr << "Error in <get_r>: cosTheta = " << cosTheta << " outside physical range !!\n";
    assert(0);
  }
  if ( cosTheta < -1. ) cosTheta = -1.;
  if ( cosTheta > +1. ) cosTheta = +1.;
  double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
  reco::Candidate::Vector r = (h - k*cosTheta)*(1./sinTheta);
  if ( verbosity >= 3 )
  {
    printVector("r", r, cartesian);
  }
  return r;
}

reco::Candidate::Vector
get_n(const reco::Candidate::Vector& k, const reco::Candidate::Vector& r, int verbosity, bool cartesian)
{
  // CV: The ordering of r and k in the cross product has been agreed with Luca on 06/09/2023.
  //     The definition n = r x k has been chosen for consistency with Eq. (2.5) in the paper arXiv:1508.05271,
  //     which Luca and Marco have used in their previous papers on Entanglement.
  //    (Whether one computes the vector n using n = r x k or using n = p x k makes no difference:
  //     in both cases, the vector n refers to the direction perpendicular to the scattering plane
  //     and the vectors { n, r, k } define a right-handed coordinate system)
  reco::Candidate::Vector n = r.Cross(k);
  if ( verbosity >= 3 )
  {
    printVector("n", n, cartesian);
  }
  return n;
}

void
get_localCoordinateSystem(const reco::Candidate::LorentzVector& p4,
                          const reco::Candidate::LorentzVector* recoilP4, const ROOT::Math::Boost* boost_ttrf,
                          int hAxis, int collider,
                          reco::Candidate::Vector& r, reco::Candidate::Vector& n, reco::Candidate::Vector& k,
                          int verbosity, bool cartesian)
{
  if ( verbosity >= 3 )
  {
    std::cout << "<get_localCoordinateSystem>:\n";
  }

  k = get_k(p4, boost_ttrf, verbosity, cartesian);

  reco::Candidate::Vector h;
  if ( hAxis == kBeam )
  {
    h = get_h_beamAxis(boost_ttrf, collider, verbosity, cartesian);
  }
  else if ( hAxis == kHiggs )
  {
    if ( !recoilP4 )
      throw cmsException("get_localCoordinateSystem", __LINE__)
        << "Configuration parameter 'hAxis' = " << hAxis << ", but no recoilP4 given !!\n";
    if ( !boost_ttrf )
      throw cmsException("get_localCoordinateSystem", __LINE__)
        << "Configuration parameter 'hAxis' = " << hAxis << ", but no boost_ttrf given !!\n";
    h = get_h_higgsAxis(*recoilP4, *boost_ttrf, verbosity, cartesian);
  }
  else assert(0);
  r = get_r(k, h, verbosity, cartesian);
  n = get_n(k, r, verbosity, cartesian);

  if ( verbosity >= 4 )
  {
    std::cout << "n*r = " << n.Dot(r) << "\n";
    std::cout << "n*k = " << n.Dot(k) << "\n";
    std::cout << "r*k = " << r.Dot(k) << "\n";
  }
}
