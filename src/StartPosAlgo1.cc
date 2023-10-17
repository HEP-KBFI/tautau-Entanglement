#include "TauAnalysis/Entanglement/interface/StartPosAlgo1.h"

#include "DataFormats/Candidate/interface/Candidate.h"                    // reco::Candidate::LorentzVector, reco::Candidate::Point
#include "DataFormats/Math/interface/deltaR.h"                            // deltaR2()
#include "DataFormats/Math/interface/Matrix.h"                            // math::Matrix
#include "DataFormats/Math/interface/Vector.h"                            // math::Vector
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wswitch-enum"
#include "DataFormats/TauReco/interface/PFTau.h"                          // reco::PFTau::hadronicDecayMode
#pragma GCC diagnostic pop

#include "TauAnalysis/Entanglement/interface/comp_PCA_line2line.h"        // comp_PCA_line2line()
#include "TauAnalysis/Entanglement/interface/constants.h"                 // ct, mChargedPion, mHiggs, mTau
#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/cube.h"                      // cube()
#include "TauAnalysis/Entanglement/interface/fixMass.h"                   // fixNuMass(), fixTauMass()
#include "TauAnalysis/Entanglement/interface/get_decayMode.h"             // is1Prong()
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h" // get_localCoordinateSystem()
#include "TauAnalysis/Entanglement/interface/get_leadTrack.h"             // get_leadTrack()
#include "TauAnalysis/Entanglement/interface/get_tipPerp.h"               // get_tipPerp()
#include "TauAnalysis/Entanglement/interface/printLorentzVector.h"        // printLorentzVector()
#include "TauAnalysis/Entanglement/interface/square.h"                    // square()

#include <algorithm>                                                      // std::sort()
#include <cmath>                                                          // std::sqrt()
#include <utility>                                                        // std::make_pair()
#include <iostream>                                                       // std::cout

namespace math
{
  typedef Matrix<3,3>::type Matrix3x3;
  typedef Vector<3>::type   Vector3;
}

StartPosAlgo1::StartPosAlgo1(const edm::ParameterSet& cfg)
  : StartPosAlgoBase(cfg)
  , applyHiggsMassConstraint_(cfg.getParameter<bool>("applyHiggsMassConstraint"))
{}

StartPosAlgo1::~StartPosAlgo1()
{}

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
  get_epsilon(int mu, int nu, int rho, int sigma, int verbosity)
  {
    //if ( verbosity >= 3 )
    //{
    //  std::cout << "<get_epsilon>:\n";
    //}
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
    //if ( verbosity >= 3 )
    //{
    //  std::cout << " mu = " << mu << ", nu = " << nu << ", rho = " << rho << ", sigma = " << sigma << ":" 
    //            << " epsilon = " << epsilon << "\n";
    //}
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
    if ( verbosity >= 3 )
    {
      std::cout << "<comp_d>:\n";
      printLorentzVector("p_H", higgsP4);
      printLorentzVector("p_vis+", visTauPlusP4);
      printLorentzVector("p_vis-", visTauMinusP4);
    }
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
            // CV: compute Levi-Civita symbol in Eq. (C3) of the paper arXiv:2211.10513,
            //    cf. https://en.wikipedia.org/wiki/Levi-Civita_symbol
            //    (Section "Levi-Civita tensors")
            double epsilon = get_epsilon(mu, nu, rho, sigma, verbosity);
            double product = (1./mH2)*epsilon*higgsP4_nu*visTauPlusP4_rho*visTauMinusP4_sigma;
            // CV: minus sign for spatial dimensions 
            //     due to raising indices from covariant from contravariant from
            //     by multiplication with metric tensor (+---)
            double sign = ( mu == 0 ) ? +1. : -1.;            
            if      ( mu == 0 ) qE  += sign*product;
            else if ( mu == 1 ) qPx += sign*product;
            else if ( mu == 2 ) qPy += sign*product;
            else if ( mu == 3 ) qPz += sign*product;
            else assert(0);
          }
        }
      }
    }
    qP4 = reco::Candidate::LorentzVector(qPx, qPy, qPz, qE);
    if ( verbosity >= 3 )
    {
      printLorentzVector("q", qP4);
      std::cout << "q*p_H = " << qP4.Dot(higgsP4) << "\n";
      std::cout << "q*p_vis+ = " << qP4.Dot(visTauPlusP4) << "\n";
      std::cout << "q*p_vis- = " << qP4.Dot(visTauMinusP4) << "\n";
    }
    // CV: compute q2 "by hand" instead of calling qP4.mass() function,
    //     as the four-vector qP4 may be space-like, i.e. have negative "mass"
    double q2 = qE*qE - (qPx*qPx + qPy*qPy + qPz*qPz);
    if ( verbosity >= 1 )
    {
      printLorentzVector("qP4", qP4, cartesian);
      std::cout << "q2 = " << q2 << "\n";
    }

    const double mTau2 = square(mTau);

    double d2 = (-0.25/q2)*((1. + a2)*mH2 + b2*mVisTauPlus2 + c2*mVisTauMinus2 - 4.*mTau2 + 2.*(a*c*y - a*b*x - b*c*z));
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
}

std::vector<KinematicEvent>
StartPosAlgo1::operator()(const KinematicEvent& kineEvt)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<StartPosAlgo1::operator()>:\n";
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
    std::cout << "M:\n";
    std::cout << M << "\n";
  }

  const double mTau2 = square(mTau);

  math::Vector3 lambda;
  lambda(0) = mTau2 + mVisTauPlus2  - x;
  lambda(1) = mTau2 + mVisTauMinus2 - y;
  lambda(2) = 0.;
  if ( verbosity_ >= 1 )
  {
    std::cout << "lambda:\n";
    std::cout << lambda << "\n";
  }

  // CV: invert matrix M;
  //     see Section "Linear algebra functions" of the ROOT documentation https://root.cern.ch/doc/v608/SMatrixDoc.html for the syntax
  int errorFlag = 0;
  math::Matrix3x3 Minv = M.Inverse(errorFlag);
  if ( errorFlag != 0 )
    throw cmsException("StartPosAlgo1::operator()", __LINE__)
      << "Failed to invert matrix M !!\n";

  math::Vector3 v = Minv*lambda;

  double a = v(0);
  double b = v(1);
  double c = v(2);
  if ( verbosity_ >= 1 )
  {
    std::cout << "a = " << a << "\n";
    std::cout << "b = " << b << "\n";
    std::cout << "c = " << c << "\n";
  }

  std::vector<KinematicEvent> kineEvts_startpos;

  reco::Candidate::LorentzVector qP4;
  std::vector<double> d = comp_d(a, b, c, x, y, z, higgsP4, visTauPlusP4, mVisTauPlus2, visTauMinusP4, mVisTauMinus2, qP4, verbosity_, cartesian_);
  for ( size_t idx = 0; idx < d.size(); ++idx )
  {
    reco::Candidate::LorentzVector tauPlusP4  = 0.5*(1. - a)*higgsP4 + 0.5*b*visTauPlusP4 - 0.5*c*visTauMinusP4 + d[idx]*qP4;
    reco::Candidate::LorentzVector tauMinusP4 = 0.5*(1. + a)*higgsP4 - 0.5*b*visTauPlusP4 + 0.5*c*visTauMinusP4 - d[idx]*qP4;
    if ( verbosity_ >= 1 )
    {
      std::cout << "solution #" << idx << ":\n";
      printLorentzVector("tauPlusP4", tauPlusP4, cartesian_);
      if ( cartesian_ )
      {
        std::cout << " mass = " << tauPlusP4.mass() << "\n";
      }
      printLorentzVector("tauMinusP4", tauMinusP4, cartesian_);
      {
        std::cout << " mass = " << tauMinusP4.mass() << "\n";
      }
      double tipPerp = get_tipPerp(tauPlusP4, visTauPlusP4, kineEvt.pv(), kineEvt.tipPCATauPlus())
                      + get_tipPerp(tauMinusP4, visTauMinusP4, kineEvt.pv(), kineEvt.tipPCATauMinus());
      std::cout << "tipPerp = " << tipPerp << "\n";
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

    kineEvts_startpos.push_back(kineEvt_startpos);
  }

  return kineEvts_startpos;
}
