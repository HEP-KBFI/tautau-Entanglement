#include "TauAnalysis/Entanglement/interface/StartPosFinder1.h"

#include "DataFormats/Candidate/interface/Candidate.h"             // reco::Candidate::LorentzVector
#include "DataFormats/Math/interface/deltaR.h"                     // deltaR2()
#include "DataFormats/Math/interface/Matrix.h"                     // math::Matrix
#include "DataFormats/Math/interface/Vector.h"                     // math::Vector
#include "DataFormats/TauReco/interface/PFTau.h"                   // reco::PFTau::hadronicDecayMode

#include "TauAnalysis/Entanglement/interface/comp_PCA_line2line.h" // comp_PCA_line2line()
#include "TauAnalysis/Entanglement/interface/constants.h"          // mChargedPion, mHiggs, mTau
#include "TauAnalysis/Entanglement/interface/cmsException.h"       // cmsException
#include "TauAnalysis/Entanglement/interface/cube.h"               // cube()
#include "TauAnalysis/Entanglement/interface/fixMass.h"            // fixNuMass(), fixTauMass()
#include "TauAnalysis/Entanglement/interface/get_decayMode.h"      // is1Prong()
#include "TauAnalysis/Entanglement/interface/get_leadTrack.h"      // get_leadTrack()
#include "TauAnalysis/Entanglement/interface/printLorentzVector.h" // printLorentzVector()
#include "TauAnalysis/Entanglement/interface/square.h"             // square()

#include <cmath>                                                   // std::sqrt()
#include <utility>                                                 // std::make_pair()
#include <iostream>                                                // std::cout

namespace math
{
  typedef Matrix<3,3>::type Matrix3x3;
  typedef Vector<3>::type   Vector3;
}

StartPosFinder1::StartPosFinder1(const edm::ParameterSet& cfg)
  : StartPosFinderBase(cfg)
  , resolutions_(nullptr)
  , applyHiggsMassConstraint_(cfg.getParameter<bool>("applyHiggsMassConstraint"))
{
  edm::ParameterSet cfg_resolutions = cfg.getParameterSet("resolutions");
  resolutions_ = new Resolutions(cfg_resolutions);
}

StartPosFinder1::~StartPosFinder1()
{
  delete resolutions_;
}

namespace
{
  reco::Candidate::LorentzVector
  buildHiggsP4_corr(const reco::Candidate::LorentzVector& higgsP4, double sf)
  {
    double higgsPx_corr = sf*higgsP4.px();
    double higgsPy_corr = sf*higgsP4.py();
    double higgsPz_corr = sf*higgsP4.pz();
    double higgsE_corr  = std::sqrt(square(higgsPx_corr) + square(higgsPy_corr) + square(higgsPz_corr) + square(mHiggs));
    reco::Candidate::LorentzVector higgsP4_corr(higgsPx_corr, higgsPy_corr, higgsPz_corr, higgsE_corr);
    return higgsP4_corr;
  }

  double
  get_component(const reco::Candidate::LorentzVector& p4, int idx)
  {
    if      ( idx == 0 ) return p4.energy();
    else if ( idx == 1 ) return p4.px();
    else if ( idx == 2 ) return p4.py();
    else if ( idx == 3 ) return p4.pz();
    else assert(0);  
  }

  int
  sgn(int x)
  {
    if ( x > 0 ) return +1;
    if ( x < 0 ) return -1;
    return 0;
  }

  double
  get_epsilon(int mu, int nu, int rho, int sigma)
  {
    //std::cout << "<get_epsilon>:\n";
    // CV: formula for computation of four-dimensional Levi-Civita symbol taken from
    //       https://en.wikipedia.org/wiki/Levi-Civita_symbol
    //    (Section "Definition" -> "Generalization to n dimensions")
    int a[4];
    a[0] = mu;
    a[1] = nu;
    a[2] = rho;
    a[3] = sigma;
    double epsilon = 1.;
    for ( int i = 0; i < 4; ++i )
    {
      for ( int j = i + 1; j < 4; ++j )
      {
        epsilon *= sgn(a[j] - a[i]); 
      }
    }
    //std::cout << " mu = " << mu << ", nu = " << nu << ", rho = " << rho << ", sigma = " << sigma << ":" 
    //          << " epsilon = " << epsilon << "\n";
    return epsilon;
  }

  std::vector<double>
  comp_d(double a, double b, double c, double x, double y, double z,
         const reco::Candidate::LorentzVector& higgsP4,
         const reco::Candidate::LorentzVector& visTauPlusP4, double mVisTauPlus2,
         const reco::Candidate::LorentzVector& visTauMinusP4, double mVisTauMinus2,
         reco::Candidate::LorentzVector& qP4,
         int verbosity, bool cartesian = true)
  {
    double a2 = square(a);
    double b2 = square(b);
    double c2 = square(c);

    double mH2 = square(higgsP4.mass());

    double qPx = 0.;
    double qPy = 0.;
    double qPz = 0.;
    double qE  = 0.;
    for ( int mu = 0; mu < 4; ++mu )
    {
      for ( int nu = 0; nu < 4; ++nu )
      {
        double higgsP4_nu = get_component(higgsP4, nu);
        for ( int rho = 0; rho < 4; ++rho )
        {
          double visTauPlusP4_rho = get_component(visTauPlusP4, rho);
          for ( int sigma = 0; sigma < 4; ++sigma )
          {
            double visTauMinusP4_sigma = get_component(visTauMinusP4, sigma);
            // CV: minus sign on four-dimensional Levi-Civita symbol 
            //     due to raising indices from covariant from contravariant from
            //     by multiplication with metric tensor (+---),
            //    cf. https://en.wikipedia.org/wiki/Levi-Civita_symbol
            //    (Section "Levi-Civita tensors")
            double epsilon = -get_epsilon(mu, nu, rho, sigma);
            double product = (1./mH2)*epsilon*higgsP4_nu*visTauPlusP4_rho*visTauMinusP4_sigma;
            if      ( mu == 0 ) qE  += product;
            else if ( mu == 1 ) qPx += product;
            else if ( mu == 2 ) qPy += product;
            else if ( mu == 3 ) qPz += product;
            else assert(0);
          }
        }
      }
    }
    qP4 = reco::Candidate::LorentzVector(qPx, qPy, qPz, qE);
    // CV: compute q2 "by hand" instead of calling qP4.mass() function,
    //     as the four-vector qP4 may be space-like, i.e. have negative "mass"
    double q2 = qE*qE - (qPx*qPx + qPy*qPy + qPz*qPz);
    if ( verbosity >= 1 )
    {
      printLorentzVector("qP4", qP4, cartesian);
      std::cout << "q2 = " << q2 << "\n";
    }

    const double mTau2 = square(mTau);

    double d2 = (-0.25/q2)*((1. + a2)*mH2 + 0.5*(b2 + c2)*(mVisTauPlus2 + mVisTauMinus2) - 4.*mTau2 + 2.*(a*c*y - a*b*x - b*c*z));
    if ( verbosity >= 1 )
    {
      std::cout << "d2 = " << d2 << "\n";
    }

    std::vector<double> d;
    if ( d2 > 0. )
    {
      d.push_back(+std::sqrt(d2));
      d.push_back(-std::sqrt(d2));
    }
    else
    {
      d.push_back(0.);
    }
    return d;
  }
  
  double
  get_tipPerp(const reco::Candidate::LorentzVector& tauP4, 
              const reco::Candidate::LorentzVector& visTauP4, 
              const reco::Candidate::Point& pv, const reco::Candidate::Point& tipPCA)
  {
    reco::Candidate::Vector e_tau  = tauP4.Vect().unit();
    reco::Candidate::Vector e_vis  = visTauP4.Vect().unit();
    reco::Candidate::Vector e_perp = e_tau.Cross(e_vis).unit();
    reco::Candidate::Vector flightlength = tipPCA - pv;
    double tipPerp = std::fabs(e_perp.Dot(flightlength));
    return tipPerp;
  }

  std::pair<reco::Candidate::LorentzVector, reco::Candidate::LorentzVector>
  comp_tauP4(double a, double b, double c, const std::vector<double>& d,
             const reco::Candidate::LorentzVector& higgsP4,
             const reco::Candidate::LorentzVector& visTauPlusP4,
             const reco::Candidate::LorentzVector& visTauMinusP4,
             const reco::Candidate::LorentzVector& qP4,
             const KinematicEvent& evt,
             double& tauMassFix,
             int verbosity, bool cartesian = true)
  {
    reco::Candidate::LorentzVector tauPlusP4_min_tipPerp;
    reco::Candidate::LorentzVector tauMinusP4_min_tipPerp;
    double min_tipPerp = 1.e+3;
    for ( size_t idx = 0; idx < d.size(); ++idx )
    {
      reco::Candidate::LorentzVector tauPlusP4  = 0.5*(1. - a)*higgsP4 + 0.5*b*visTauPlusP4 - 0.5*c*visTauMinusP4 + d[idx]*qP4;
      reco::Candidate::LorentzVector tauMinusP4 = 0.5*(1. + a)*higgsP4 - 0.5*b*visTauPlusP4 + 0.5*c*visTauMinusP4 - d[idx]*qP4;
      double tipPerp = get_tipPerp(tauPlusP4, visTauPlusP4, evt.pv(), evt.tipPCATauPlus())
                      + get_tipPerp(tauMinusP4, visTauMinusP4, evt.pv(), evt.tipPCATauMinus());
      if ( verbosity >= 1 )
      {
        std::cout << "solution #" << idx << ": tipPerp = " << tipPerp << "\n";
        printLorentzVector("tauPlusP4", tauPlusP4, cartesian);
        if ( cartesian )
        {
          std::cout << " mass = " << tauPlusP4.mass() << "\n";
        }
        printLorentzVector("tauMinusP4", tauMinusP4, cartesian);
        {
          std::cout << " mass = " << tauMinusP4.mass() << "\n";
        }
      }
      if ( idx == 0 || tipPerp < min_tipPerp )
      {        
        // CV: the reconstruction of the tau+ and tau- four-vectors using Eq. (C3)
        //     of the paper arXiv:2211.10513 yields mass values that a few GeV off,
        //     presumably due to rounding errors.
        //     Adjust the energy component of the tau+ and tau- four-vectors 
        //     such that their mass values match the correct (PDG) value,
        //     while keeping the Px, Py, Pz momentum components fixed
        tauPlusP4_min_tipPerp = fixTauMass(tauPlusP4);
        tauMinusP4_min_tipPerp = fixTauMass(tauMinusP4);
        min_tipPerp = tipPerp;
        tauMassFix = square(tauPlusP4.mass() - mTau) + square(tauMinusP4.mass() - mTau);
      }
    }
    if ( verbosity >= 1 )
    {
      std::cout << "best solution:\n";
      printLorentzVector("tauPlusP4", tauPlusP4_min_tipPerp, cartesian);
      if ( cartesian )
      {
        std::cout << " mass = " << tauPlusP4_min_tipPerp.mass() << "\n";
      }
      printLorentzVector("tauMinusP4", tauMinusP4_min_tipPerp, cartesian);
      if ( cartesian )
      {
        std::cout << " mass = " << tauMinusP4_min_tipPerp.mass() << "\n";
      }
    }

    return std::make_pair(tauPlusP4_min_tipPerp, tauMinusP4_min_tipPerp);
  }
}

KinematicEvent
StartPosFinder1::operator()(const KinematicEvent& kineEvt)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<StartPosFinder1::operator()>:\n";
  }

  reco::Candidate::LorentzVector higgsP4 = kineEvt.recoilP4();
  if ( applyHiggsMassConstraint_ )
  {
    //printLorentzVector("higgsP4", higgsP4);
    double higgsP     = higgsP4.P();
    double higgsP_2   = square(higgsP);
    double higgsP_3   = higgsP_2*higgsP;
    double higgsP_4   = higgsP_3*higgsP;
    double sigmaP_2   = square(higgsP4.px()/higgsP)*square(resolutions_->recoilResolution_px())
                       + square(higgsP4.py()/higgsP)*square(resolutions_->recoilResolution_py())
                       + square(higgsP4.pz()/higgsP)*square(resolutions_->recoilResolution_pz());
    double sigmaP_4   = square(sigmaP_2);
    double higgsE     = higgsP4.energy();
    double sigmaE_2   = square(higgsP4.px()/higgsE)*square(resolutions_->recoilResolution_px())
                       + square(higgsP4.py()/higgsE)*square(resolutions_->recoilResolution_py())
                       + square(higgsP4.pz()/higgsE)*square(resolutions_->recoilResolution_pz())
                       + square(higgsP4.mass()/higgsE)*square(resolutions_->recoilResolution_mass());
    double sigmaE_4   = square(sigmaE_2);
    double mHiggs_2   = square(mHiggs);
    double mHiggs_4   = square(mHiggs_2);
    double higgsEc    = std::sqrt(higgsP_2 + mHiggs_2);
    double higgsEc_3  = cube(higgsEc);
    double higgsEc_5  = higgsEc_3*square(higgsEc);
    double higgsEc_10 = square(higgsEc_5);
    double term1 = higgsEc_5*sigmaE_2/(3.*higgsE*mHiggs_2*higgsP_2);
    double term2 = -1./sigmaE_2 + 3.*higgsE*mHiggs_2*higgsP_2/(higgsEc_5*sigmaE_2) - higgsE*higgsP_2/(higgsEc_3*sigmaE_2) + higgsE/(higgsEc*sigmaE_2) - 1./sigmaP_2;
    double term3 = 3.*higgsE*mHiggs_2*higgsP_2*sigmaP_2*(2.*higgsEc_5*sigmaE_2 - higgsE*mHiggs_2*higgsP_2*sigmaP_2 + 2.*higgsE*higgsP_4*sigmaP_2)
                  + square(higgsEc_5*sigmaE_2 - higgsE*mHiggs_4*sigmaP_2 - 4.*higgsE*mHiggs_2*higgsP_2*sigmaP_2 + higgsEc_5*sigmaP_2);
    double term4 = higgsEc_10*sigmaE_4*sigmaP_4;
    double sf1 = term1*(term2 + std::sqrt(term3/term4));
    double sf2 = term1*(term2 - std::sqrt(term3/term4));
    double sf = -1.;
    if ( sf1 > 0. && sf2 > 0. )
    {
      if ( std::fabs(sf1 - 1.) < std::fabs(sf2 - 1.) ) sf = sf1;
      else                                             sf = sf2;
    }
    else if ( sf1 > 0. )
    {
      sf = sf1;
    }
    else if ( sf2 > 0. )
    {
      sf = sf2;
    }
    else 
    {
      std::cerr << "WARNING: Failed to correct four-vector of recoil system !!" << std::endl;
      sf = 1.;
    }
    reco::Candidate::LorentzVector higgsP4_corr = buildHiggsP4_corr(higgsP4, sf);
    if ( verbosity_ >= 1 )
    {
      std::cout << "correcting four-vector of recoil system:\n";
      std::cout << "sf = " << sf << "\n";
      printLorentzVector("higgsP4_corr", higgsP4_corr);
      std::cout << " mass = " << higgsP4_corr.mass() << "\n";
    }
    higgsP4 = higgsP4_corr;
  }

  const reco::Candidate::LorentzVector& visTauPlusP4  = kineEvt.visTauPlusP4();
  const reco::Candidate::LorentzVector& visTauMinusP4 = kineEvt.visTauMinusP4();

  // reconstruct the four-vectors of tau+ and tau- in the laboratory (detector) system.
  // The formulas are taken from Appendix C of the paper arXiv:2211.10513
  double mVisTauPlus  = kineEvt.tauPlus_decayMode()  == reco::PFTau::kOneProng0PiZero ?
    mChargedPion : visTauPlusP4.mass();
  double mVisTauMinus = kineEvt.tauMinus_decayMode() == reco::PFTau::kOneProng0PiZero ?
    mChargedPion : visTauMinusP4.mass();

  double x = higgsP4.Dot(visTauPlusP4);
  double y = higgsP4.Dot(visTauMinusP4);
  double z = visTauPlusP4.Dot(visTauMinusP4);

  double mH2 = square(higgsP4.mass());
  double mVisTauPlus2 = square(mVisTauPlus);
  double mVisTauMinus2 = square(mVisTauMinus);

  math::Matrix3x3 M;
  M(0,0) = -x;
  M(0,1) = mVisTauPlus2;
  M(0,2) = -z;
  M(1,0) =  y;
  M(1,1) = -z;
  M(1,2) = mVisTauMinus2;
  M(2,0) = mH2;
  M(2,1) = -x;
  M(2,2) =  y;
  if ( verbosity_ >= 1 )
  {
    std::cout << "M0:\n";
    std::cout << M << "\n";
  }

  const double mTau2 = square(mTau);

  math::Vector3 lambda;
  lambda(0) = mTau2 + mVisTauPlus2  - x;
  lambda(1) = mTau2 + mVisTauMinus2 - y;
  lambda(2) = 0.;
  if ( verbosity_ >= 1 )
  {
    std::cout << "lambda0:\n";
    std::cout << lambda << "\n";
  }

  // CV: invert matrix M;
  //     see Section "Linear algebra functions" of the ROOT documentation https://root.cern.ch/doc/v608/SMatrixDoc.html for the syntax
  int errorFlag = 0;
  math::Matrix3x3 Minv = M.Inverse(errorFlag);
  if ( errorFlag != 0 )
    throw cmsException("StartPosFinder1::operator()", __LINE__)
      << "Failed to invert matrix M !!\n";

  math::Vector3 v = Minv*lambda;

  double a0 = v(0);
  double b0 = v(1);
  double c0 = v(2);
  if ( verbosity_ >= 1 )
  {
    std::cout << "a0 = " << a0 << "\n";
    std::cout << "b0 = " << b0 << "\n";
    std::cout << "c0 = " << c0 << "\n";
  }

  reco::Candidate::LorentzVector q0P4;
  std::vector<double> d0 = comp_d(a0, b0, c0, x, y, z, higgsP4, visTauPlusP4, mVisTauPlus2, visTauMinusP4, mVisTauMinus2, q0P4, verbosity_, cartesian_);
  double tauMassFix0 = 0.;
  auto tau0P4 = comp_tauP4(a0, b0, c0, d0, higgsP4, visTauPlusP4, visTauMinusP4, q0P4, kineEvt, tauMassFix0, verbosity_, cartesian_);
  reco::Candidate::LorentzVector tauPlusP4 = tau0P4.first;
  reco::Candidate::LorentzVector tauMinusP4 = tau0P4.second;

  int iteration = 0;
  const int max_iterations = 10;
  bool hasConverged = false;
  while ( !hasConverged && iteration < max_iterations )
  {
    if ( verbosity_ >= 1 )
    {
      std::cout << "iteration #" << iteration << ":\n";
    }

    double dm2 = mVisTauMinus2 - mVisTauPlus2;

    M(0,0) = -x;
    M(0,1) = mVisTauPlus2;
    M(0,2) = -z;
    M(1,0) =  y;
    M(1,1) = -z;
    M(1,2) = mVisTauMinus2;
    M(2,0) = mH2;
    M(2,1) = -x + 0.5*b0*dm2;
    M(2,2) =  y + 0.5*c0*dm2;
    if ( verbosity_ >= 1 )
    {
      std::cout << "M:\n";
      std::cout << M << "\n";
    }
  
    lambda(0) = mTau2 + mVisTauPlus2  - x;
    lambda(1) = mTau2 + mVisTauMinus2 - y;
    lambda(2) = 0.25*(square(b0) + square(c0))*dm2;
    if ( verbosity_ >= 1 )
    {
      std::cout << "lambda:\n";
      std::cout << lambda << "\n";
    }
 
    v = Minv*lambda;

    double a1 = v(0);
    double b1 = v(1);
    double c1 = v(2);
    if ( verbosity_ >= 1 )
    {
      std::cout << "a = " << a1 << "\n";
      std::cout << "b = " << b1 << "\n";
      std::cout << "c = " << c1 << "\n";
    }

    reco::Candidate::LorentzVector q1P4;
    std::vector<double> d1 = comp_d(a1, b1, c1, x, y, z, higgsP4, visTauPlusP4, mVisTauPlus2, visTauMinusP4, mVisTauMinus2, q1P4, verbosity_, cartesian_);
    double tauMassFix1 = 0.;
    auto tau1P4 = comp_tauP4(a1, b1, c1, d1, higgsP4, visTauPlusP4, visTauMinusP4, q1P4, kineEvt, tauMassFix1, verbosity_, cartesian_);
    if ( verbosity_ >= 1 )
    {
      std::cout << "tauMassFix1 = " << tauMassFix1 << " (tauMassFix0 = " << tauMassFix0 << ")\n";
    }

    if ( tauMassFix1 < tauMassFix0 )
    {
      a0 = a1;
      b0 = b1;
      c0 = c1;
      tau0P4 = tau1P4;
      tauPlusP4 = tau1P4.first;
      tauMinusP4 = tau1P4.second;
      tauMassFix0 = tauMassFix1;
    }
    else
    {
      hasConverged = true;
    }
    ++iteration;
  }
  if ( !hasConverged )
  {
    std::cerr << "WARNING: StartPosFinder1 failed to converge !!" << std::endl;
  }
  if ( verbosity_ >= 1 )
  {
    std::cout << "#iterations = " << iteration << " (max_iterations = " << max_iterations << ")\n";
  }

  if ( verbosity_ >= 1 )
  {
    std::cout << "final start position for kinematic fit:\n";
    printLorentzVector("tauPlusP4", tauPlusP4, cartesian_);
    if ( cartesian_ )
    {
      std::cout << " mass = " << tauPlusP4.mass() << "\n";
    }
    printLorentzVector("tauMinusP4", tauMinusP4, cartesian_);
    {
      std::cout << " mass = " << tauMinusP4.mass() << "\n";
    }
  }

  KinematicEvent kineEvt_startpos = kineEvt;

  kineEvt_startpos.tauPlusP4_ = tauPlusP4;
  kineEvt_startpos.tauPlusP4_isValid_ = true;
  //kineEvt_startpos.nuTauPlusP4_ = fixNuMass(tauPlusP4 - visTauPlusP4);
  kineEvt_startpos.nuTauPlusP4_ = tauPlusP4 - visTauPlusP4;
  kineEvt_startpos.nuTauPlusP4_isValid_ = true;
  if ( is1Prong(kineEvt_startpos.tauPlus_decayMode()) )
  {
    const KinematicParticle* tauPlus_leadTrack = get_leadTrack(kineEvt_startpos.daughtersTauPlus());
    assert(tauPlus_leadTrack);
    const reco::Candidate::Point& tauPlus_tipPCA = kineEvt_startpos.tipPCATauPlus();
    kineEvt_startpos.svTauPlus_ = comp_PCA_line2line(kineEvt_startpos.pv(), tauPlusP4, tauPlus_tipPCA, tauPlus_leadTrack->p4(), verbosity_);
    kineEvt_startpos.svTauPlus_isValid_ = true;
  }

  kineEvt_startpos.tauMinusP4_ = tauMinusP4;
  kineEvt_startpos.tauMinusP4_isValid_ = true;
  //kineEvt_startpos.nuTauMinusP4_ = fixNuMass(tauMinusP4 - visTauMinusP4);
  kineEvt_startpos.nuTauMinusP4_ = tauMinusP4 - visTauMinusP4;
  kineEvt_startpos.nuTauMinusP4_isValid_ = true;
  if ( is1Prong(kineEvt_startpos.tauMinus_decayMode()) )
  {
    const KinematicParticle* tauMinus_leadTrack = get_leadTrack(kineEvt_startpos.daughtersTauMinus());
    assert(tauMinus_leadTrack);
    const reco::Candidate::Point& tauMinus_tipPCA = kineEvt_startpos.tipPCATauMinus();
    kineEvt_startpos.svTauMinus_ = comp_PCA_line2line(kineEvt_startpos.pv(), tauMinusP4, tauMinus_tipPCA, tauMinus_leadTrack->p4(), verbosity_);
    kineEvt_startpos.svTauMinus_isValid_ = true;
  }

  return kineEvt_startpos;
}
