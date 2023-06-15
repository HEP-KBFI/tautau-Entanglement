#ifndef TauAnalysis_Entanglement_auxFunctions_h
#define TauAnalysis_Entanglement_auxFunctions_h

#include "DataFormats/Candidate/interface/Candidate.h" // Candidate::LorentzVector, Candidate::Vector

#include <Math/Boost.h>                                // Boost

#include <string>                                      // std::string

double
square(double x);
  
double
cube(double x);

reco::Candidate::LorentzVector
getP4_rf(const reco::Candidate::LorentzVector& p4,
         const ROOT::Math::Boost& boost);

reco::Candidate::LorentzVector
getP4_hf(const reco::Candidate::LorentzVector& p4,
         const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k);

reco::Candidate::LorentzVector
getP4_ttrf_hf_trf(const reco::Candidate::LorentzVector& p4,
                  const ROOT::Math::Boost& boost_ttrf,
                  const reco::Candidate::Vector& r, const reco::Candidate::Vector& n, const reco::Candidate::Vector& k,
                  const ROOT::Math::Boost& boost_trf);

void
printLorentzVector(const std::string& label,
                   const reco::Candidate::LorentzVector& p4,
                   bool cartesian = true);

void
printVector(const std::string& label,
            const reco::Candidate::Vector& p3,
            bool cartesian = true);

#endif // TauAnalysis_Entanglement_auxFunctions_h
