#include "TauAnalysis/Entanglement/interface/KinematicParticle.h"

#include "TauAnalysis/Entanglement/interface/auxFunctions.h" // square()
#include "TauAnalysis/Entanglement/interface/cmsException.h" // cmsException
#include "TauAnalysis/Entanglement/interface/constants.h"    // Bfield; xr, yr, zr

#include <cstdlib>                                           // std::getenv()

TDatabasePDG* KinematicParticle::pdg_ = nullptr;

KinematicParticle::KinematicParticle(int pdgId)
  : pdgId_(pdgId)
  , params5_isValid_(false)
  , params7_isValid_(false)
{
  if ( !pdg_ )
  {
    pdg_ = TDatabasePDG::Instance();
    std::string rootsys = std::getenv("ROOTSYS");
    if ( rootsys == "" )
      throw cmsException("KinematicParticle", __LINE__)
        << "Environment variable '$ROOTSYS' not initialized !!\n";
    std::string filename = rootsys.append("/etc/pdg_table.txt");
    pdg_->ReadPDGTable(filename.c_str());
  }

  TParticlePDG* pdgEntry =  pdg_->GetParticle(pdgId);
  if ( !pdgEntry )
    throw cmsException("KinematicParticle", __LINE__)
      << "Failed to find entry in PDG table for particle with pdgId = " << pdgId << " !!\n";
  mass_ = pdgEntry->Mass();
  charge_ = round(pdgEntry->Charge());
}

KinematicParticle::~KinematicParticle()
{}

void 
KinematicParticle::set_params5(const TVectorD& params5, const TMatrixD& cov5x5)
{
  params5_ = params5;
  cov5x5_ = cov5x5;
  params5_isValid_ = true;

  // unpack "C" parameters
  double c      = params5(0);
  double phi0   = params5(1);
  double d0     = params5(2);
  double lambda = params5(3);
  double z0     = params5(4);

  // convert from "C" to "W" parameters, 
  // using Eq. (43) of http://www.phys.ufl.edu/~avery/fitting/kinematic.pdf
  double a = -0.2998*Bfield*charge_;
  double a_over_2c = a/(2.*c);

  double px = a_over_2c*cos(phi0);
  double py = a_over_2c*sin(phi0);
  double pz = a_over_2c*lambda;
  double E  = sqrt(square(a_over_2c)*(1. + square(lambda)) + square(mass_));
  p4_ = reco::Candidate::LorentzVector(px, py, pz, E);

  double x = xr - d0*sin(phi0);
  double y = yr + d0*cos(phi0);
  double z = xr + z0;
  vertex_ = reco::Candidate::Point(x, y, z);

  // pack "W" paramaters
  params7_.ResizeTo(7);
  params7_(0) = px;
  params7_(1) = py;
  params7_(2) = pz;
  params7_(3) = E;
  params7_(4) = x;
  params7_(5) = y;
  params7_(6) = z;

  // compute covariance matrix of "W" paramaters using error propagation:
  //   cov(params7) = B*cov(params5)*B^T,
  // where B refers to the partial derivatives of the "new coordinates" (7 params)
  // with respect to the "old coordinates" (5 params), cf. Eq. (4.50) in 
  //   V. Blobel and E. Lohrmann "Statistische und numerische Methoden der Datenanalyse".
  // The partial derivatives are taken from appendix I of
  //   http://www.phys.ufl.edu/~avery/fitting/kinematic.pdf

  double pt = sqrt(square(px) + square(py));
  double p = sqrt(square(px) + square(py) + square(pz));

  TMatrixD B(7,7);
  B(0,0) = -px/c;
  B(0,1) = -py;
  B(1,0) = -py/c;
  B(1,1) =  px;
  B(2,0) = -pz/c;
  B(2,3) =  pt;
  B(3,0) = -square(p)/(E*c);
  B(3,3) =  lambda*square(pt)/E;
  B(4,1) = -y;
  B(4,2) = -py/pt;
  B(5,1) =  x;
  B(5,2) =  py/pt;
  B(6,4) =  1.;

  TMatrixD BT = TMatrixD(B);
  BT.Transpose(B);

  cov7x7_.ResizeTo(7,7);
  cov7x7_ = B*cov5x5_*BT;

  params7_isValid_ = true;
}

void 
KinematicParticle::set_params7(const TVectorD& params7, const TMatrixD& cov7x7)
{
  params7_ = params7;
  cov7x7_ = cov7x7;
  params7_isValid_ = true;
}

const reco::Candidate::LorentzVector&
KinematicParticle::get_p4() const
{
  return p4_;
}

const reco::Candidate::Point&
KinematicParticle::get_vertex() const
{
  return vertex_;
}

int
KinematicParticle::get_pdgId() const
{
  return pdgId_;
}

double
KinematicParticle::get_mass() const
{
  return mass_;
}

int
KinematicParticle::get_charge() const
{
  return charge_;
}

const TVectorD&
KinematicParticle::get_params5() const
{
  throw cmsException("KinematicParticle", __LINE__)
    << "Parametrization 'C' not initialized !!\n";
  return params5_;
}

const TMatrixD&
KinematicParticle::get_cov5x5() const
{
  throw cmsException("KinematicParticle", __LINE__)
    << "Parametrization 'C' not initialized !!\n";
  return cov5x5_;
}

const TVectorD&
KinematicParticle::get_params7() const
{
  throw cmsException("KinematicParticle", __LINE__)
    << "Parametrization 'W' not initialized !!\n";
  return params7_;
}

const TMatrixD&
KinematicParticle::get_cov7x7() const
{
  throw cmsException("KinematicParticle", __LINE__)
    << "Parametrization 'W' not initialized !!\n";
  return cov7x7_;
}
