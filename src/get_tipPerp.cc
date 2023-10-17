#include "TauAnalysis/Entanglement/interface/get_tipPerp.h"

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
