#include "TauAnalysis/Entanglement/interface/auxFunctions.h"

double
square(double x)
{
  return x*x;
}
  
double
cube(double x)
{
  return x*x*x;
}

reco::Candidate::LorentzVector
getP4_rf(const reco::Candidate::LorentzVector& p4,
         const ROOT::Math::Boost& boost)
{
  // CV: boost given four-vector to restframe
  reco::Candidate::LorentzVector p4_rf = boost(p4);
  return p4_rf;
}

reco::Candidate::LorentzVector
getP4_hf(const reco::Candidate::LorentzVector& p4,
         const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k)
{
  // CV: rotate given four-vector to helicity frame
  reco::Candidate::Vector p3 = p4.Vect();
  double Pr = p3.Dot(r);
  double Pn = p3.Dot(n);
  double Pk = p3.Dot(k);
  reco::Candidate::LorentzVector p4_hf(Pr, Pn, Pk, p4.energy());
  return p4_hf;
}

reco::Candidate::LorentzVector
getP4_ttrf_hf_trf(const reco::Candidate::LorentzVector& p4,
                  const ROOT::Math::Boost& boost_ttrf,
                  const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k,
                  const ROOT::Math::Boost& boost_trf)
{
  // CV: boost given four-vector to restframe of tau pair,
  //     rotate to helicity frame,
  //     and finally boost to tau restframe
  reco::Candidate::LorentzVector p4_ttrf = getP4_rf(p4, boost_ttrf);
  reco::Candidate::LorentzVector p4_hf = getP4_hf(p4_ttrf, r, n, k);
  reco::Candidate::LorentzVector p4_trf = getP4_rf(p4_hf, boost_trf);
  return p4_trf;
}

void
printLorentzVector(const std::string& label,
                   const reco::Candidate::LorentzVector& p4,
                   bool cartesian)
{
  std::cout << label << ":";
  if ( cartesian )
  {
    std::cout << " E = " << p4.energy() << ", Px = " << p4.px() << ", Py = " << p4.py() << ", Pz = " << p4.pz() << "\n";
  }
  else
  {
    std::cout << " pT = " << p4.pt() << ", eta = " << p4.eta() << ", phi = " << p4.phi() << ", mass = " << p4.mass() << "\n";
  }
}

void
printVector(const std::string& label,
            const reco::Candidate::Vector& p3,
            bool cartesian)
{
  std::cout << label << ":";
  if ( cartesian )
  {
    std::cout << " Px = " << p3.x() << ", Py = " << p3.y() << ", Pz = " << p3.z() << "\n";
  }
  else
  {
    std::cout << " pT = " << p3.r() << ", eta = " << p3.eta() << ", phi = " << p3.phi() << "\n";
  }
}
