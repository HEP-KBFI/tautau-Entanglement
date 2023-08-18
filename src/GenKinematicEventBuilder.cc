#include "TauAnalysis/Entanglement/interface/GenKinematicEventBuilder.h"

#include "DataFormats/Candidate/interface/Candidate.h"                            // reco::Candidate::LorentzVector, reco::Candidate::Vector
#include "DataFormats/TauReco/interface/PFTau.h"                                  // reco::PFTau::hadronicDecayMode

#include "TauAnalysis/Entanglement/interface/comp_nuP4.h"                         // comp_nuP4()
#include "TauAnalysis/Entanglement/interface/comp_visP4.h"                        // comp_visP4()
#include "TauAnalysis/Entanglement/interface/findDecayProducts.h"                 // findDecayProducts()
#include "TauAnalysis/Entanglement/interface/findLastTau.h"                       // findLastTau()
#include "TauAnalysis/Entanglement/interface/getCov_hf.h"                         // getCov_hf()
#include "TauAnalysis/Entanglement/interface/get_decayMode.h"                     // get_decayMode()
#include "TauAnalysis/Entanglement/interface/get_leadTrack.h"                     // get_leadTrack()
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h"         // get_localCoordinateSystem()
#include "TauAnalysis/Entanglement/interface/get_particles_of_type.h"             // get_chargedHadrons(), get_neutralPions(), get_neutrinos()
#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h"                 // math::Matrix*, math::Vector*
#include "TauAnalysis/Entanglement/interface/KinematicParticle.h"                 // KinematicParticle
#include "TauAnalysis/Entanglement/interface/printCovMatrix.h"                    // printCovMatrix()
#include "TauAnalysis/Entanglement/interface/printDistance.h"                     // printDistance()
#include "TauAnalysis/Entanglement/interface/printEigenVectors_and_EigenValues.h" // printEigenVectors_and_EigenValues()
#include "TauAnalysis/Entanglement/interface/printGenParticles.h"                 // printGenParticles()
#include "TauAnalysis/Entanglement/interface/printInverseCovMatrix.h"             // printInverseCovMatrix()
#include "TauAnalysis/Entanglement/interface/printLorentzVector.h"                // printLorentzVector()
#include "TauAnalysis/Entanglement/interface/printPoint.h"                        // printPoint()
#include "TauAnalysis/Entanglement/interface/printVector.h"                       // printVector()
#include "TauAnalysis/Entanglement/interface/PolarimetricVectorAlgoBase.h"        // ePolarimetricVector::kTauPlus, PolarimetricVector::kTauMinus
#include "TauAnalysis/Entanglement/interface/square.h"                            // square()

#include <TMath.h>                                                                // TMath::Pi()

#include <iostream>                                                               // std::cout
#include <cmath>                                                                  // std::fabs(), std::sqrt()

GenKinematicEventBuilder::GenKinematicEventBuilder(const edm::ParameterSet& cfg)
  : resolutions_(nullptr)
  , smearing_(cfg)
  , applySmearing_(cfg.getParameter<bool>("applySmearing"))
  , polarimetricVector_(cfg)
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
  , cartesian_(cfg.getUntrackedParameter<bool>("cartesian"))
{
  edm::ParameterSet cfg_resolutions = cfg.getParameterSet("resolutions");
  resolutions_ = new Resolutions(cfg_resolutions);
}

GenKinematicEventBuilder::~GenKinematicEventBuilder()
{
  delete resolutions_;
}

namespace
{
  math::Matrix3x3
  build_pvCov(const Resolutions& resolutions)
  {
    math::Matrix3x3 pvCov;
    pvCov(0,0) = square(resolutions.pvResolution_xy());
    pvCov(1,1) = square(resolutions.pvResolution_xy());
    pvCov(2,2) = square(resolutions.pvResolution_z());
    return pvCov;
  }

  math::Matrix4x4
  build_recoilCov(const reco::Candidate::LorentzVector& recoilP4, const Resolutions& resolutions)
  {
    double sigma2_px   = square(resolutions.recoilResolution_px());
    double sigma2_py   = square(resolutions.recoilResolution_py());
    double sigma2_pz   = square(resolutions.recoilResolution_pz());
    double sigma2_mass = square(resolutions.recoilResolution_mass());
    math::Matrix4x4 recoilCov;
    recoilCov(0,0) = sigma2_px;
    recoilCov(0,3) = (recoilP4.px()/recoilP4.energy())*sigma2_px;
    recoilCov(1,1) = sigma2_py;
    recoilCov(1,3) = (recoilP4.py()/recoilP4.energy())*sigma2_py;
    recoilCov(2,2) = sigma2_pz;
    recoilCov(2,3) = (recoilP4.pz()/recoilP4.energy())*sigma2_pz;
    recoilCov(3,0) = recoilCov(0,3);
    recoilCov(3,1) = recoilCov(1,3);
    recoilCov(3,2) = recoilCov(2,3);
    recoilCov(3,3) = square(recoilP4.px()/recoilP4.energy())*sigma2_px
                    + square(recoilP4.py()/recoilP4.energy())*sigma2_py
                    + square(recoilP4.pz()/recoilP4.energy())*sigma2_pz
                    + square(recoilP4.mass()/recoilP4.energy())*sigma2_mass;
    return recoilCov;
  }

  math::Matrix4x4
  comp_visTauCov(const reco::Candidate::LorentzVector& visTauP4, const std::vector<KinematicParticle>& daughters)
  {
    math::Matrix4x4 visTauCov;
    for ( const KinematicParticle& daughter : daughters )
    {
      int daughter_absPdgId = std::abs(daughter.pdgId());
      if ( daughter_absPdgId == 12 || daughter_absPdgId == 14 || daughter_absPdgId == 16 ) continue;
      const math::Matrix7x7& cov7x7 = daughter.cov7x7();
      // CV: compute 4x4 covariance matrix of the tauh four-vector by summing the covariance matrices of all visible tau decay products;
      //     for the syntax of retrieving a small covariance matrix from a larger one, 
      //     see Section "Accessing and setting methods" of the ROOT documentation https://root.cern.ch/doc/v608/SMatrixDoc.html
      visTauCov += cov7x7.Sub<math::Matrix4x4>(0,0);
    }
    // CV: add uncertainty on mass of visible tau decay products.
    //     This uncertainty needs to be added even if tau is reconstructed in OneProng0Pi0 decay mode,
    //     to avoid that covariance matrix is rank-deficient and cannot be inverted
    const double visTau_sigma2_mass = square(1.e-1);
    visTauCov(3,3) += square(visTauP4.mass()/visTauP4.energy())*visTau_sigma2_mass;
    return visTauCov;
  }

  math::Matrix3x3
  comp_svTauCov(const std::vector<KinematicParticle>& daughters)
  {
    math::Matrix3x3 svTauCov;
    bool isFirst = true;
    for ( const KinematicParticle& daughter : daughters )
    {
      int daughter_absPdgId = std::abs(daughter.pdgId());
      if ( daughter_absPdgId == 12 || daughter_absPdgId == 14 || daughter_absPdgId == 16 ) continue;
      const math::Matrix7x7& cov7x7 = daughter.cov7x7();
      // CV: compute 3x3 covariance matrix of the tauh vertex position, simply take covariance matrix of first tau decay product
      //    (assuming the vertex positions of all tau decay products to be equal)
      //     For the syntax of retrieving a small covariance matrix from a larger one, 
      //     see Section "Accessing and setting methods" of the ROOT documentation https://root.cern.ch/doc/v608/SMatrixDoc.html
      if ( isFirst )
      {
        svTauCov = cov7x7.Sub<math::Matrix3x3>(4,4);
        isFirst = false;
      }
    }
    return svTauCov;
  }

  reco::Candidate::Point
  comp_PCA_line2point(const reco::Candidate::Point& pv, const reco::GenParticle* leadTrack, int verbosity)
  {
    // CV: compute point of closest approach (PCA) between primary event vertex and track
    if ( verbosity >= 2 )
    {
      std::cout << "<comp_PCA_line2point>:" << std::endl;
    }
    reco::Candidate::Point sv = leadTrack->vertex();
    auto flightlength = sv - pv;
    if ( verbosity >= 2 )
    {
      printDistance("flightlength", flightlength, true);
      printDistance("flightlength", flightlength, false);
    }
    auto e_trk = leadTrack->momentum().unit();
    reco::Candidate::Point pca = pv + flightlength - flightlength.Dot(e_trk)*e_trk;
    if ( verbosity >= 2 )
    {
      printPoint("pca", pca);
    }
    return pca;
  }

  std::vector<KinematicParticle>
  build_kineDaughters(const reco::Candidate::LorentzVector& visTauP4, int tau_decayMode, const std::vector<const reco::GenParticle*>& tau_decayProducts,
                      const Resolutions& resolutions,
                      int verbosity)
  {
    if ( verbosity >= 3 )
    {
      std::cout << "<build_kineDaughters>:" << std::endl;
    }

    reco::Candidate::Vector r, n, k; 
    get_localCoordinateSystem(visTauP4, nullptr, nullptr, kBeam, r, n, k);
    if ( verbosity >= 3 )
    {
      printVector("r", r);
      printVector("n", n);
      printVector("k", k);
    }

    std::vector<KinematicParticle> kineDaughters;
    for ( const reco::GenParticle* decayProduct : tau_decayProducts )
    {
      KinematicParticle kineDaughter(decayProduct->pdgId());

      math::Vector7 params7;
      params7(0) = decayProduct->px();
      params7(1) = decayProduct->py();
      params7(2) = decayProduct->pz();
      params7(3) = decayProduct->energy();
      params7(4) = decayProduct->vertex().x();
      params7(5) = decayProduct->vertex().y();
      params7(6) = decayProduct->vertex().z();

      math::Matrix7x7 cov7x7;
      double sigma_pt     = 0.;
      double sigma_theta  = 0.;
      double sigma_phi    = 0.;
      if ( std::fabs(decayProduct->charge()) > 0.5 )
      {
        sigma_pt          = get_trackResolution_pt(decayProduct->p4(), resolutions);
        sigma_theta       = TMath::Pi()*resolutions.trackResolution_theta();
        sigma_phi         = TMath::Pi()*resolutions.trackResolution_phi();
      }
      else if ( decayProduct->pdgId() == 111 )
      {
        sigma_pt          = get_ecalResolution_pt(decayProduct->p4(), resolutions);
        sigma_theta       = TMath::Pi()*resolutions.ecalResolution_theta();
        sigma_phi         = TMath::Pi()*resolutions.ecalResolution_phi();
      }
      double sigma2_pt    = square(sigma_pt);
      double sigma2_theta = square(sigma_theta);
      double sigma2_phi   = square(sigma_phi);
      double cot_theta    = cos(decayProduct->theta())/sin(decayProduct->theta());
      double cot2_theta   = square(cot_theta);
      double csc_theta    = 1./sin(decayProduct->theta());
      double csc2_theta   = square(csc_theta);
      double csc4_theta   = square(csc2_theta);
      double cos_phi      = cos(decayProduct->phi());
      double cos2_phi     = square(cos_phi);
      double sin_phi      = sin(decayProduct->phi());
      double sin2_phi     = square(sin_phi);
      double energy       = decayProduct->energy();
      double energy2      = square(energy);
      double pt           = decayProduct->pt();
      double pt2          = square(pt);
      double pt3          = pt*pt2;
      double pt4          = pt*pt3;
      cov7x7(0,0) = cos2_phi*sigma2_pt + pt2*sin2_phi*sigma2_phi;
      cov7x7(0,1) = cos_phi*sin_phi*sigma2_pt - pt2*cos_phi*sin_phi*sigma2_phi;
      cov7x7(0,2) = cos_phi*cot_theta*sigma2_pt;
      cov7x7(0,3) = (pt*cos_phi*csc2_theta/energy)*sigma2_pt;
      cov7x7(1,0) = cov7x7(0,1);
      cov7x7(1,1) = sin2_phi*sigma2_pt + pt2*cos2_phi*sigma2_phi;
      cov7x7(1,2) = sin_phi*cot_theta*sigma2_pt;
      cov7x7(1,3) = (pt*sin_phi*csc2_theta/energy)*sigma2_pt;
      cov7x7(2,0) = cov7x7(0,2);
      cov7x7(2,1) = cov7x7(1,2);
      cov7x7(2,2) = cot2_theta*sigma2_pt + pt2*csc4_theta*sigma2_theta;
      cov7x7(2,3) = (pt*cot_theta*csc2_theta/energy)*sigma2_pt + (pt3*cot_theta*csc4_theta/energy)*sigma2_theta;
      cov7x7(3,0) = cov7x7(0,3);
      cov7x7(3,1) = cov7x7(1,3);
      cov7x7(3,2) = cov7x7(2,3);
      cov7x7(3,3) = (pt2*csc4_theta/energy2)*sigma2_pt + (pt4*cot2_theta*csc4_theta/energy2)*sigma2_theta ;
      double dk = 0.;
      double dr = 0.;
      double dn = 0.;
      if ( tau_decayMode == reco::PFTau::kThreeProng0PiZero )
      {
        // CV: set uncertainty on position of tau decay vertex 
        //     in direction parallel and perpendicular to the momentum vector of the visible decay products
        //     for three-prongs
        dk = resolutions.svResolution_parl();
        dr = resolutions.svResolution_perp();
        dn = resolutions.svResolution_perp();
      }
      else
      {
        // CV: for one-prongs, the position of the tau decay vertex in direction perpendicular to the momentum vector of the visible decay products
        //     is constrained by the transverse impact parameter of the leading track,
        //     while the position in parallel direction is not known.
        //     Following the "huge error method" described in Section 6 of https://www.phys.ufl.edu/~avery/fitting/fitting1.pdf
        //     we set the uncertainty in parallel direction to a large (but not very large, to avoid numerical instabilities) value,
        //     so that the position in parallel direction is practically unconstrained.
        dk = 1.e+1;
        dr = resolutions.tipResolution_perp();
        dn = resolutions.tipResolution_perp();
      }
      math::Matrix3x3 cov3x3 = getCov_hf(dk, dr, dn, r, n, k);
      cov7x7.Place_at(cov3x3, 4, 4);
      kineDaughter.set_params7(params7, cov7x7);

      kineDaughters.push_back(kineDaughter);
    }
    return kineDaughters;
  }
 
  math::Matrix3x3
  comp_nuCov(const reco::Candidate::LorentzVector& visTauP4, int verbosity)
  {
    // CV: the four-vectors of tau+ and tau- are not really "measured";
    //     we acccount for this by setting their covariance matrix to an ellipsoid.
    //     The size of the ellipsoid in direction of the visible tau decay products (k) is set to a large value,
    //     following the "huge error method" described in Section 6 of https://www.phys.ufl.edu/~avery/fitting/fitting1.pdf
    //     In the direction transverse to the visible tau decay products (r, n),
    //     the size of the ellipsoid is set to be 1 GeV each.
    //     The motivation of using a smaller uncertainty in transverse direction is to guide the kinematic fit 
    //     that the tau lepton momentum in this direction is limited by the tau lepton mass.
    //     Note that the neutrino momentum in transverse direction is limited to half the tau lepton mass in the tau restframe
    //     and is not affected by the Lorentz boost in tau direction.
    //     We enlarge the uncertainty somewhat to allow for more flexibility in the fit.
    reco::Candidate::Vector r, n, k; 
    get_localCoordinateSystem(visTauP4, nullptr, nullptr, kBeam, r, n, k);
    double dk = 2.5e+2;
    double dr = 2.e0;
    double dn = dr;
    math::Matrix3x3 nuCov = getCov_hf(dk, dr, dn, r, n, k);
    return nuCov;
  }
}

KinematicEvent
GenKinematicEventBuilder::operator()(const reco::GenParticleCollection& genParticles)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<GenKinematicEventBuilder::operator()>:\n";
    std::cout << " applySmearing = " << applySmearing_ << "\n";
  }

  const reco::GenParticle* tauPlus  = nullptr;
  const reco::GenParticle* tauMinus = nullptr;
  for ( const reco::GenParticle& genParticle : genParticles )
  {
    if ( genParticle.pdgId() == -15 && !tauPlus  ) tauPlus  = findLastTau(&genParticle);
    if ( genParticle.pdgId() == +15 && !tauMinus ) tauMinus = findLastTau(&genParticle);
  }
  if ( !(tauPlus && tauMinus) ) 
  {
    std::cerr << "WARNING: Failed to find tau+ tau- pair --> skipping the event !!" << std::endl;
    return KinematicEvent();
  }
  reco::Candidate::LorentzVector tauPlusP4 = tauPlus->p4();
  reco::Candidate::LorentzVector tauMinusP4 = tauMinus->p4();

  std::vector<const reco::GenParticle*> tauPlus_daughters;
  findDecayProducts(tauPlus, tauPlus_daughters);
  reco::Candidate::LorentzVector visTauPlusP4 = comp_visP4(tauPlus_daughters);
  reco::Candidate::LorentzVector nuTauPlusP4 = comp_nuP4(tauPlus_daughters);  
  math::Matrix3x3 nuTauPlusCov = comp_nuCov(visTauPlusP4, verbosity_);
  std::vector<const reco::GenParticle*> tauPlus_ch = get_chargedHadrons(tauPlus_daughters);
  std::vector<const reco::GenParticle*> tauPlus_pi0 = get_neutralPions(tauPlus_daughters);
  std::vector<const reco::GenParticle*> tauPlus_nu = get_neutrinos(tauPlus_daughters);
  int tauPlus_decayMode = get_decayMode(tauPlus_ch, tauPlus_pi0, tauPlus_nu);
  if ( verbosity_ >= 1 )
  {
    printLorentzVector("tauPlusP4", tauPlusP4, cartesian_);
    if ( cartesian_ )
    {
      std::cout << " mass = " << tauPlusP4.mass() << "\n";
    }
    std::cout << "tau+ decay products:" << "\n";
    printGenParticles(tauPlus_daughters, cartesian_);
    printLorentzVector("visTauPlusP4", visTauPlusP4, cartesian_);
    if ( cartesian_ )
    {
      std::cout << " mass = " << visTauPlusP4.mass() << "\n";
    }
    std::cout << "tauPlus_decayMode = " << tauPlus_decayMode << "\n";
    printCovMatrix("nuTauPlusCov", nuTauPlusCov);
    if ( verbosity_ >= 3 )
    {
      printEigenVectors_and_EigenValues(nuTauPlusCov);
    }
    printInverseCovMatrix("nuTauPlusCov", nuTauPlusCov);
  }
  std::vector<KinematicParticle> daughtersTauPlus = build_kineDaughters(visTauPlusP4, tauPlus_decayMode, tauPlus_daughters, *resolutions_, verbosity_);

  std::vector<const reco::GenParticle*> tauMinus_daughters;
  findDecayProducts(tauMinus, tauMinus_daughters);
  reco::Candidate::LorentzVector visTauMinusP4 = comp_visP4(tauMinus_daughters);
  reco::Candidate::LorentzVector nuTauMinusP4 = comp_nuP4(tauMinus_daughters);
  math::Matrix3x3 nuTauMinusCov = comp_nuCov(visTauMinusP4, verbosity_);
  std::vector<const reco::GenParticle*> tauMinus_ch = get_chargedHadrons(tauMinus_daughters);
  std::vector<const reco::GenParticle*> tauMinus_pi0 = get_neutralPions(tauMinus_daughters);
  std::vector<const reco::GenParticle*> tauMinus_nu = get_neutrinos(tauMinus_daughters);
  int tauMinus_decayMode = get_decayMode(tauMinus_ch, tauMinus_pi0, tauMinus_nu);
  if ( verbosity_ >= 1 )
  {
    printLorentzVector("tauMinusP4", tauMinusP4, cartesian_);
    if ( cartesian_ )
    {
      std::cout << " mass = " << tauPlusP4.mass() << "\n";
    }
    std::cout << "tau- decay products:" << "\n";
    printGenParticles(tauMinus_daughters, cartesian_);
    printLorentzVector("visTauMinusP4", visTauMinusP4, cartesian_);
    if ( cartesian_ )
    {
      std::cout << " mass = " << visTauMinusP4.mass() << "\n";
    }
    std::cout << "tauMinus_decayMode = " << tauMinus_decayMode << "\n";
    printCovMatrix("nuTauMinusCov", nuTauMinusCov);
    if ( verbosity_ >= 3 )
    {
      printEigenVectors_and_EigenValues(nuTauMinusCov);
    }
    printInverseCovMatrix("nuTauMinusCov", nuTauMinusCov);
  }
  std::vector<KinematicParticle> daughtersTauMinus = build_kineDaughters(visTauMinusP4, tauMinus_decayMode, tauMinus_daughters, *resolutions_, verbosity_);
 
  if ( !(std::fabs(tauPlus->vertex().x() - tauMinus->vertex().x()) < 1.e-3 &&
         std::fabs(tauPlus->vertex().y() - tauMinus->vertex().y()) < 1.e-3 &&
         std::fabs(tauPlus->vertex().z() - tauMinus->vertex().z()) < 1.e-3) )
  {
    std::cerr << "WARNING: Failed to find vertex of tau+ tau- pair --> skipping the event !!" << std::endl;
    return KinematicEvent();
  }
  // CV: set primary event vertex (PV),
  //     using generator-level information
  reco::Candidate::Point pv = tauMinus->vertex();
  math::Matrix3x3 pvCov = build_pvCov(*resolutions_);
  if ( verbosity_ >= 1 )
  {
    printPoint("pv", pv);
    printCovMatrix("pvCov", pvCov);
    printInverseCovMatrix("pvCov", pvCov);
  }

  reco::Candidate::LorentzVector recoilP4 = tauPlusP4 + tauMinusP4;
  math::Matrix4x4 recoilCov = build_recoilCov(recoilP4, *resolutions_);
  if ( verbosity_ >= 1 )
  {
    printLorentzVector("recoilP4", recoilP4, cartesian_);
    printCovMatrix("recoilCov", recoilCov);
    printInverseCovMatrix("recoilCov", recoilCov);
  }

  const reco::GenParticle* tauPlus_leadTrack = get_leadTrack(tauPlus_daughters);
  const reco::GenParticle* tauMinus_leadTrack = get_leadTrack(tauMinus_daughters);
  if ( !(tauPlus_leadTrack && tauMinus_leadTrack) )
  {
    std::cerr << "WARNING: Failed to find leading tracks of tau+ tau- pair --> skipping the event !!" << std::endl;
    return KinematicEvent();
  }
  reco::Candidate::Point svTauPlus = tauPlus_leadTrack->vertex();
  math::Matrix3x3 svTauPlusCov = comp_svTauCov(daughtersTauPlus);
  reco::Candidate::Point svTauMinus = tauMinus_leadTrack->vertex();
  math::Matrix3x3 svTauMinusCov = comp_svTauCov(daughtersTauMinus);
  if ( verbosity_ >= 1 )
  {
    printPoint("svTauPlus", svTauPlus);
    printDistance("svTauPlus - pv", svTauPlus - pv, cartesian_);
    printDistance("svTauPlus - pv", svTauPlus - pv, false);
    printCovMatrix("svTauPlusCov", svTauPlusCov);
    if ( verbosity_ >= 3 )
    {
      printEigenVectors_and_EigenValues(svTauPlusCov);
    }
    printInverseCovMatrix("svTauPlusCov", svTauPlusCov);
    printPoint("svTauMinus", svTauMinus);
    printDistance("svTauMinus - pv", svTauMinus - pv, cartesian_);
    printDistance("svTauMinus - pv", svTauMinus - pv, false);
    printCovMatrix("svTauMinusCov", svTauMinusCov);
    if ( verbosity_ >= 3 )
    {
      printEigenVectors_and_EigenValues(svTauMinusCov);
    }
    printInverseCovMatrix("svTauMinusCov", svTauMinusCov);
  }

  KinematicEvent kineEvt;

  kineEvt.pv_ = pv;
  kineEvt.pvCov_ = pvCov;

  kineEvt.recoilP4_ = recoilP4;
  kineEvt.recoilCov_ = recoilCov;

  if ( !applySmearing_ )
  {
    kineEvt.tauPlusP4_ = tauPlusP4;
    kineEvt.tauPlusP4_isValid_ = true;
    kineEvt.nuTauPlusP4_ = nuTauPlusP4;
    kineEvt.nuTauPlusP4_isValid_ = true;
  }
  kineEvt.nuTauPlusCov_ = nuTauPlusCov;
  kineEvt.nuTauPlusCov_isValid_ = true;
  kineEvt.visTauPlusP4_ = visTauPlusP4;
  kineEvt.visTauPlusCov_ = comp_visTauCov(visTauPlusP4, daughtersTauPlus);
  kineEvt.tauPlus_decayMode_ = tauPlus_decayMode;
  kineEvt.daughtersTauPlus_ = daughtersTauPlus;
  kineEvt.tipPCATauPlus_ = comp_PCA_line2point(pv, tauPlus_leadTrack, verbosity_);
  kineEvt.svTauPlus_ = svTauPlus;
  kineEvt.svTauPlus_isValid_ = true;
  kineEvt.svTauPlusCov_ = svTauPlusCov;

  if ( !applySmearing_ )
  {
    kineEvt.tauMinusP4_ = tauMinusP4;
    kineEvt.tauMinusP4_isValid_ = true;
    kineEvt.nuTauMinusP4_ = nuTauMinusP4;
    kineEvt.nuTauMinusP4_isValid_ = true;
  }
  kineEvt.nuTauMinusCov_ = nuTauMinusCov;
  kineEvt.nuTauMinusCov_isValid_ = true;
  kineEvt.visTauMinusP4_ = visTauMinusP4;
  kineEvt.visTauMinusCov_ = comp_visTauCov(visTauMinusP4, daughtersTauMinus);
  kineEvt.tauMinus_decayMode_ = tauMinus_decayMode;
  kineEvt.daughtersTauMinus_ = daughtersTauMinus;
  kineEvt.tipPCATauMinus_ = comp_PCA_line2point(pv, tauMinus_leadTrack, verbosity_);
  kineEvt.svTauMinus_ = tauMinus_leadTrack->vertex();
  kineEvt.svTauMinus_isValid_ = true;
  kineEvt.svTauMinusCov_ = svTauMinusCov;

  if ( applySmearing_ )
  {
    kineEvt = smearing_(kineEvt);
  }

  reco::Candidate::Vector hPlus = polarimetricVector_(kineEvt, pol::kTauPlus);
  kineEvt.hPlus_ = hPlus;
  kineEvt.hPlus_isValid_ = true;
  reco::Candidate::Vector hMinus = polarimetricVector_(kineEvt, pol::kTauMinus);
  kineEvt.hMinus_ = hMinus;
  kineEvt.hMinus_isValid_ = true;

  return kineEvt;
}
