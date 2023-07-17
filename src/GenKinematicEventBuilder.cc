#include "TauAnalysis/Entanglement/interface/GenKinematicEventBuilder.h"

#include "DataFormats/Candidate/interface/Candidate.h"                    // reco::Candidate::LorentzVector, reco::Candidate::Vector
#include "DataFormats/TauReco/interface/PFTau.h"                          // reco::PFTau::hadronicDecayMode

#include "TauAnalysis/Entanglement/interface/compVisP4.h"                 // compVisP4()
#include "TauAnalysis/Entanglement/interface/findDecayProducts.h"         // findDecayProducts()
#include "TauAnalysis/Entanglement/interface/findLastTau.h"               // findLastTau()
#include "TauAnalysis/Entanglement/interface/get_decayMode.h"             // get_decayMode()
#include "TauAnalysis/Entanglement/interface/get_leadTrack.h"             // get_leadTrack()
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h" // get_localCoordinateSystem()
#include "TauAnalysis/Entanglement/interface/get_particles_of_type.h"     // get_chargedHadrons(), get_neutralPions(), get_neutrinos()
#include "TauAnalysis/Entanglement/interface/KinematicParticle.h"         // math::Matrix7x7, math::Vector7
#include "TauAnalysis/Entanglement/interface/printGenParticles.h"         // printGenParticles()
#include "TauAnalysis/Entanglement/interface/printLorentzVector.h"        // printLorentzVector()
#include "TauAnalysis/Entanglement/interface/printPoint.h"                // printPoint()
#include "TauAnalysis/Entanglement/interface/printVector.h"               // printVector()
#include "TauAnalysis/Entanglement/interface/SpinAnalyzerBase.h"          // SpinAnalyzerBase::kTauPlus, SpinAnalyzerBase::kTauMinus
#include "TauAnalysis/Entanglement/interface/square.h"                    // square()

#include <Math/Boost.h>                                                   // ROOT::Math::Boost
#include <TMath.h>                                                        // TMath::Pi()
#include <TMatrixD.h>                                                     // TMatrixD
#include <TVectorD.h>                                                     // TVectorD

#include <iostream>                                                       // std::cout
#include <cmath>                                                          // std::fabs(), std::sqrt()

GenKinematicEventBuilder::GenKinematicEventBuilder(const edm::ParameterSet& cfg)
  : resolutions_(nullptr)
  , smearing_(cfg)
  , applySmearing_(cfg.getParameter<bool>("applySmearing"))
  , spinAnalyzer_(cfg)
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
  reco::Candidate::Point
  get_tipPCA(const reco::Candidate::Point& pv, const reco::GenParticle* leadTrack)
  {
    // CV: compute point of closest approach (PCA) between primary event vertex and track
    reco::Candidate::Point sv = leadTrack->vertex();
    auto flightlength = sv - pv;
    auto e_trk = leadTrack->momentum().unit();
    reco::Candidate::Point pca = pv + flightlength - flightlength.Dot(e_trk)*e_trk;
    return pca;
  }

  std::vector<KinematicParticle>
  get_kineDaughters(const reco::Candidate::LorentzVector& visTauP4, int tau_decayMode, const std::vector<const reco::GenParticle*>& tau_decayProducts,
                    const Resolutions& resolutions,
                    int verbosity)
  {
    if ( verbosity >= 3 )
    {
      std::cout << "<get_kineDaughters>:" << std::endl;
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
        sigma_theta       = TMath::Pi()*resolutions.get_trackResolution_theta();
        sigma_phi         = TMath::Pi()*resolutions.get_trackResolution_phi();
      }
      else if ( decayProduct->pdgId() == 111 )
      {
        sigma_pt          = get_ecalResolution_pt(decayProduct->p4(), resolutions);
        sigma_theta       = TMath::Pi()*resolutions.get_ecalResolution_theta();
        sigma_phi         = TMath::Pi()*resolutions.get_ecalResolution_phi();
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
      cov7x7(2,2) = cot2_theta*sigma2_pt + pt2*csc4_theta*sigma2_phi;
      cov7x7(2,3) = (pt*cot_theta*csc2_theta/energy)*sigma2_pt - (pt3*cot_theta*csc4_theta/energy)*sigma2_theta;
      cov7x7(3,0) = cov7x7(0,3);
      cov7x7(3,1) = cov7x7(1,3);
      cov7x7(3,2) = cov7x7(2,3);
      cov7x7(3,3) = (pt2*csc4_theta/energy2)*sigma2_pt + (pt4*cot2_theta*csc4_theta/energy2)*sigma2_theta;
      double dk = 0.;
      double dr = 0.;
      double dn = 0.;
      if ( tau_decayMode == reco::PFTau::kThreeProng0PiZero )
      {
        // CV: set uncertainty on position of tau decay vertex 
        //     in direction parallel and perpendicular to the momentum vector of the visible decay products
        //     for three-prongs
        dk = resolutions.get_svResolution_parl();
        dr = resolutions.get_svResolution_perp();
        dn = resolutions.get_svResolution_perp();
      }
      else
      {
        // CV: for one-prongs, the position of the tau decay vertex in direction perpendicular to the momentum vector of the visible decay products
        //     is constrained by the transverse impact parameter of the leading track,
        //     while the position in parallel direction is not known.
        //     Following the "huge error method" described in Section 6 of https://www.phys.ufl.edu/~avery/fitting/fitting1.pdf
        //     we set the uncertainty in parallel direction to a large (but not very large, to avoid numerical instabilities) value,
        //     so that the position in parallel direction is practically unconstrained.
        dk = 1.e+2;
        dr = resolutions.get_tipResolution_perp();
        dn = resolutions.get_tipResolution_perp();
      }
      double dk2 = square(dk);
      double dr2 = square(dr);
      double dn2 = square(dn);
      reco::Candidate::Vector x(1., 0., 0.);
      reco::Candidate::Vector y(0., 1., 0.);
      reco::Candidate::Vector z(0., 0., 1.);
      double k_x = k.Dot(x);
      double k_y = k.Dot(y);
      double k_z = k.Dot(z);
      double r_x = r.Dot(x);
      double r_y = r.Dot(y);
      double r_z = r.Dot(z);
      double n_x = n.Dot(x);
      double n_y = n.Dot(y);
      double n_z = n.Dot(z);
      // CV: computation of covariance matrix in rotated coordinates taken from the paper
      //       "On transformation of covariance matrices between local Cartesian cooridinate systems and commutative diagrams",
      //       T. Soler and M. Chin; 
      //     also see Appendix I of https://www.phys.ufl.edu/~avery/fitting/kinematic.pdf
      cov7x7(4,4) = dk2*k_x*k_x + dr2*r_x*r_x + dn2*n_x*n_x;
      cov7x7(4,5) = dk2*k_x*k_y + dr2*r_x*r_y + dn2*n_x*n_y;
      cov7x7(4,6) = dk2*k_x*k_z + dr2*r_x*r_z + dn2*n_x*n_z;
      cov7x7(5,4) = dk2*k_y*k_x + dr2*r_y*r_x + dn2*n_y*n_x;
      cov7x7(5,5) = dk2*k_y*k_y + dr2*r_y*r_y + dn2*n_y*n_y;
      cov7x7(5,6) = dk2*k_y*k_z + dr2*r_y*r_z + dn2*n_y*n_z;
      cov7x7(6,4) = dk2*k_z*k_x + dr2*r_z*r_x + dn2*n_z*n_x;
      cov7x7(6,5) = dk2*k_z*k_y + dr2*r_z*r_y + dn2*n_z*n_y;
      cov7x7(6,6) = dk2*k_z*k_z + dr2*r_z*r_z + dn2*n_z*n_z;
      if ( verbosity >= 3 )
      {
        TMatrixD cov3x3;
        for ( int idxRow = 0; idxRow < 3; ++idxRow )
        {
          for ( int idxColumn = 0; idxColumn < 3; ++idxColumn )
          {
            cov3x3(idxRow,idxColumn) = cov7x7(idxRow + 4,idxColumn + 4);
          }
        }
        std::cout << "cov3x3:\n";
        cov3x3.Print();

        TVectorD eigenValues(3);
        TMatrixD eigenVectors = cov3x3.EigenVectors(eigenValues);
        for ( int idx = 0; idx < 3; ++idx )
        {
          TVectorD eigenVector(3);
          eigenVector(0) = eigenVectors(0,idx);
          eigenVector(1) = eigenVectors(1,idx);
          eigenVector(2) = eigenVectors(2,idx);
          std::cout << "EigenVector #" << idx << " (EigenValue = " << eigenValues(idx) << "):\n";
          eigenVector.Print();
          reco::Candidate::Vector e(eigenVector(0), eigenVector(1), eigenVector(2));
          std::cout << "e*k = " << e.Dot(k) << "\n";
          std::cout << "e*r = " << e.Dot(r) << "\n";
          std::cout << "e*n = " << e.Dot(n) << "\n";
        }
      }
      kineDaughter.set_params7(params7, cov7x7);

      kineDaughters.push_back(kineDaughter);
    }
    return kineDaughters;
  }
}

KinematicEvent
GenKinematicEventBuilder::operator()(const reco::GenParticleCollection& genParticles)
{
  KinematicEvent kineEvt;

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
    return kineEvt;
  }
  reco::Candidate::LorentzVector tauPlusP4 = tauPlus->p4();
  reco::Candidate::LorentzVector tauMinusP4 = tauMinus->p4();

  std::vector<const reco::GenParticle*> tauPlus_daughters;
  findDecayProducts(tauPlus, tauPlus_daughters);
  reco::Candidate::LorentzVector visTauPlusP4 = compVisP4(tauPlus_daughters);
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
  }

  std::vector<const reco::GenParticle*> tauMinus_daughters;
  findDecayProducts(tauMinus, tauMinus_daughters);
  reco::Candidate::LorentzVector visTauMinusP4 = compVisP4(tauMinus_daughters);
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
  }
 
  if ( !(std::fabs(tauPlus->vertex().x() - tauMinus->vertex().x()) < 1.e-3 &&
         std::fabs(tauPlus->vertex().y() - tauMinus->vertex().y()) < 1.e-3 &&
         std::fabs(tauPlus->vertex().z() - tauMinus->vertex().z()) < 1.e-3) )
  {
    std::cerr << "WARNING: Failed to find vertex of tau+ tau- pair --> skipping the event !!" << std::endl;
    return kineEvt;
  }
  // CV: set primary event vertex (PV),
  //     using generator-level information
  reco::Candidate::Point pv = tauMinus->vertex();
  if ( verbosity_ >= 1 )
  {
    printPoint("GenPV", pv);
  }

  const reco::GenParticle* tauPlus_leadTrack = get_leadTrack(tauPlus_daughters);
  const reco::GenParticle* tauMinus_leadTrack = get_leadTrack(tauMinus_daughters);
  if ( !(tauPlus_leadTrack && tauMinus_leadTrack) )
  {
    std::cerr << "WARNING: Failed to find leading tracks of tau+ tau- pair --> skipping the event !!" << std::endl;
    return kineEvt;
  }
  if ( verbosity_ >= 1 )
  {
    printPoint("GenSV(tau+)", tauPlus_leadTrack->vertex());
    printPoint("GenSV(tau-)", tauMinus_leadTrack->vertex());
  }

  kineEvt.pv_ = pv;
  reco::Candidate::LorentzVector recoilP4 = tauPlusP4 + tauMinusP4;
  kineEvt.recoilP4_ = recoilP4;

  if ( !applySmearing_ )
  {
    kineEvt.tauPlusP4_ = tauPlusP4;
    kineEvt.tauPlusP4_isValid_ = true;
  }
  kineEvt.visTauPlusP4_ = visTauPlusP4;
  kineEvt.tauPlus_decayMode_ = tauPlus_decayMode;
  kineEvt.daughtersTauPlus_ = get_kineDaughters(visTauPlusP4, tauPlus_decayMode, tauPlus_daughters, *resolutions_, verbosity_);
  kineEvt.tipPCATauPlus_ = get_tipPCA(pv, tauPlus_leadTrack);
  // CV: set tau decay vertex (SV) for three-prongs
  if ( tauPlus_ch.size() >= 3 )
  {
    kineEvt.svTauPlus_ = tauPlus_leadTrack->vertex();
    kineEvt.svTauPlus_isValid_ = true;
  }

  if ( !applySmearing_ )
  {
    kineEvt.tauMinusP4_ = tauMinusP4;
    kineEvt.tauMinusP4_isValid_ = true;
  }
  kineEvt.visTauMinusP4_ = visTauMinusP4;
  kineEvt.tauMinus_decayMode_ = tauMinus_decayMode;
  kineEvt.daughtersTauMinus_ = get_kineDaughters(visTauMinusP4, tauMinus_decayMode, tauMinus_daughters, *resolutions_, verbosity_);
  kineEvt.tipPCATauMinus_ = get_tipPCA(pv, tauMinus_leadTrack);
  // CV: set tau decay vertex (SV) for three-prongs
  if ( tauMinus_ch.size() >= 3 )
  {
    kineEvt.svTauMinus_ = tauMinus_leadTrack->vertex();
    kineEvt.svTauMinus_isValid_ = true;
  }

  if ( applySmearing_ )
  {
    kineEvt = smearing_(kineEvt);
  }

  reco::Candidate::Vector hPlus = spinAnalyzer_(kineEvt, SpinAnalyzerBase::kTauPlus);
  kineEvt.hPlus_ = hPlus;
  reco::Candidate::Vector hMinus = spinAnalyzer_(kineEvt, SpinAnalyzerBase::kTauMinus);
  kineEvt.hMinus_ = hMinus;

  return kineEvt;
}
