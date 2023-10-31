
#include "FWCore/ParameterSet/interface/ParameterSet.h"             // edm::ParameterSet

#include "DataFormats/Candidate/interface/Candidate.h"              // reco::Candidate::LorentzVector, reco::Candidate::Point

#include "TauAnalysis/Entanglement/interface/cmsException.h"        // cmsException
#include "TauAnalysis/Entanglement/interface/comp_PCA_line2line.h"  // comp_PCA_line2line()
#include "TauAnalysis/Entanglement/interface/printDistance.h"       // printDistance()
#include "TauAnalysis/Entanglement/interface/printLorentzVector.h"  // printLorentzVector()
#include "TauAnalysis/Entanglement/interface/printPoint.h"          // printPoint()
#include "TauAnalysis/Entanglement/interface/printVector.h"         // printVector()

#include <TBenchmark.h>                                                    // TBenchmark

#include <cmath>                                                    // std::sqrt()
#include <string>                                                   // std::string
#include <utility>                                                  // std::pair
#include <vector>                                                   // std::vector

int main(int argc, char* argv[])
{
//--- throw an exception in case ROOT encounters an error
  gErrorAbortLevel = kError;

  std::cout << "<testTIPCompatibility>:\n";

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("testTIPCompatibility");  

  // CV: event 1:1:37 taken from file /local/karl/belle_eeToTauTau/aod/unwgt_pythia_extended/aodsim_1.root
  reco::Candidate::Point pv(0.105206, 0.168923, -0.0264207);
  printPoint("pv", pv);

  reco::Candidate::LorentzVector tauPlusP4(-1.53789, 1.59619, -5.83331, 6.48831);
  reco::Candidate::LorentzVector visTauPlusP4(-0.271275, 0.349832, -2.85522, 3.02035);
  reco::Candidate::LorentzVector tauPlus_leadTrackP4(-0.145479, 0.620737, -1.75896, 1.87614);
  reco::Candidate::Point svTauPlus(0.0898842, 0.184825, -0.0845357);
  reco::Candidate::Point tipPCATauPlus(0.0946355, 0.164552, -0.0270889);
  printLorentzVector("tauPlusP4", tauPlusP4, true);
  printLorentzVector("visTauPlusP4", visTauPlusP4, true);
  printLorentzVector("tauPlus_leadTrackP4", tauPlus_leadTrackP4, true);
  printPoint("svTauPlus", svTauPlus);

  printDistance("svTauPlus - pv", svTauPlus - pv, true);
  printDistance("svTauPlus - pv", svTauPlus - pv, false);
  printLorentzVector("tauPlusP4", tauPlusP4, false);

  printDistance("svTauPlus - tipPCATauPlus", svTauPlus - tipPCATauPlus, false);
  printLorentzVector("tauPlus_leadTrackP4", tauPlus_leadTrackP4, false);

  reco::Candidate::LorentzVector tauMinusP4(1.53789, -1.59619, 2.94331, 4.09069);
  reco::Candidate::LorentzVector visTauMinusP4(1.60931, -1.00572, 2.73146, 3.45932);
  reco::Candidate::LorentzVector tauMinus_leadTrackP4(-0.0772096, -0.0106584, 0.200282, 0.256257);
  reco::Candidate::Point svTauMinus(0.106321, 0.167765, -0.024285);
  reco::Candidate::Point tipPCATauMinus(0.106913, 0.167847, -0.0258197);
  printLorentzVector("tauMinusP4", tauMinusP4, true);
  printLorentzVector("visTauMinusP4", visTauMinusP4, true);
  printLorentzVector("tauMinus_leadTrackP4", tauMinus_leadTrackP4, true);
  printPoint("svTauMinus", svTauMinus);

  printDistance("svTauMinus - pv", svTauMinus - pv, true);
  printDistance("svTauMinus - pv", svTauMinus - pv, false);
  printLorentzVector("tauMinusP4", tauMinusP4, false);

  printDistance("svTauMinus - tipPCATauMinus", svTauMinus - tipPCATauMinus, false);
  printLorentzVector("tauMinus_leadTrackP4", tauMinus_leadTrackP4, false);

  std::pair<reco::Candidate::Point, reco::Candidate::Point> Qs = comp_PCA_line2line(svTauPlus, tauPlus_leadTrackP4.Vect(), svTauMinus, tauMinus_leadTrackP4.Vect());
  const reco::Candidate::Point& Q1_1 = Qs.first;
  const reco::Candidate::Point& Q2_1 = Qs.second;
  printPoint("Q1@1", Q1_1);
  printPoint("Q2@1", Q2_1);
  std::cout << "(dmin = " << std::sqrt((Q2_1 - Q1_1).mag2()) << ")\n";

  reco::Candidate::Vector e_tauPlus = tauPlus_leadTrackP4.Vect().unit();
  reco::Candidate::Vector e_tauMinus = tauMinus_leadTrackP4.Vect().unit();

  reco::Candidate::Point Q1_2;
  reco::Candidate::Point Q2_2;
  double dmin = 1.e+6;
  const int numSteps = 30000;
  for ( int idxStep1 = 0; idxStep1 < numSteps; ++idxStep1 )
  {
    auto p1 = tipPCATauPlus + 1.e-1*((idxStep1 - 0.5*numSteps)/numSteps)*e_tauPlus;
    for ( int idxStep2 = 0; idxStep2 < numSteps; ++idxStep2 )
    {
      auto p2 = tipPCATauMinus + 1.e-1*((idxStep2 - 0.5*numSteps)/numSteps)*e_tauMinus;
      double d = std::sqrt((p1 - p2).mag2());
      if ( d < dmin )
      {
        Q1_2 = p1;
        Q2_2 = p2;
        dmin = d; 
      }
    }
  }
  printPoint("Q1@2", Q1_2);
  printPoint("Q2@2", Q2_2);
  std::cout << "(dmin = " << std::sqrt((Q2_2 - Q1_2).mag2()) << ")\n";

  clock.Show("testTIPCompatibility");

  return EXIT_SUCCESS;
}
