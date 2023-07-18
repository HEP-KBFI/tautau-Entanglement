#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h" // cmsException
#include "TauAnalysis/Entanglement/interface/constants.h"    // mProton
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
  get_h_beamAxis(const ROOT::Math::Boost* boost_ttrf, int verbosity, bool cartesian)
  {
    const double beamE  = 7.e+3;
    const double beamPx = 0.;
    const double beamPy = 0.;
    const double beamPz = std::sqrt(square(beamE) - square(mProton));
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

  reco::Candidate::Vector
  get_r(const reco::Candidate::Vector& k, const reco::Candidate::Vector& h, int verbosity, bool cartesian)
  {
    double cosTheta = k.Dot(h);
    assert(cosTheta >= -1. && cosTheta <= +1.);
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
    reco::Candidate::Vector n = k.Cross(r);
    if ( verbosity >= 3 )
    {
      printVector("n", n, cartesian);
    }
    return n;
  }
}

void
get_localCoordinateSystem(const reco::Candidate::LorentzVector& p4,
                          const reco::Candidate::LorentzVector* recoilP4, const ROOT::Math::Boost* boost_ttrf,
                          int hAxis,
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
    h = get_h_beamAxis(boost_ttrf, verbosity, cartesian);
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
    std::cout << "r*n = " << r.Dot(n) << "\n";
    std::cout << "r*k = " << r.Dot(k) << "\n";
    std::cout << "n*k = " << n.Dot(k) << "\n";
  }
}
