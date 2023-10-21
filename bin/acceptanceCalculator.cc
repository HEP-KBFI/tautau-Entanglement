#include "TauAnalysis/Entanglement/interface/constants.h" // beamEnergy_SuperKEKB_ePlus, beamEnergy_SuperKEKB_eMinus
#include "DataFormats/Candidate/interface/Candidate.h"    // reco::Candidate::LorentzVector
#include <Math/Boost.h>                                   // ROOT::Math::Boost()
#include <TMath.h>                                        // TMath::DegToRad()

reco::Candidate::LorentzVector
get_vector(double angle)
{
  reco::Candidate::LorentzVector v;
  const double radians = TMath::DegToRad() * angle;
  v.SetXYZT(0., std::sin(radians), std::cos(radians), 1.);
  return v;
}

reco::Candidate::LorentzVector
get_beam(double pz)
{
  reco::Candidate::LorentzVector v;
  v.SetPxPyPzE(0., 0, pz, std::fabs(pz)); // ignore electron mass
  return v;
}

int main(int argc, char* argv[])
{
  // +z = the direction of e-, which has higher energy (HER).
  // Polar angles are defined wrt +z direction. If the HER beam goes in the opposite -z
  // direction (like happened in our simulation), then the direction of beams should be
  // changed as well as acceptance angles should be switched around from (17, 150) to
  // (30, 163) so that the smaller opening angle is associated with the HER beam because
  // otherwise we won't undo the boost when going from lab frame to c.o.m. frame but
  // instead amplify the boost. The acceptance region should be roughly symmetric in
  // c.o.m. frame, which is why they look asymmetric in the lab frame.

  if(argc < 2)
  {
    std::cerr << "Please provide a boolean argument (e.g., true or false)\n";
    return EXIT_FAILURE;
  }
  std::string arg = argv[1];
  std::transform(arg.begin(), arg.end(), arg.begin(), [](unsigned char c) { return std::tolower(c); });

  bool switch_directions = false; // set to true to match the simulation
  if(arg == "true" || arg == "1")
  {
    switch_directions = true;
  }
  else if(arg == "false" || arg == "0")
  {
    switch_directions = false;
  }
  else
  {
    std::cerr << "Invalid boolean argument provided. Use true/false or 1/0\n";
    return EXIT_FAILURE;
  }

  const double angle1 = 17;
  const double angle2 = 30;
  const double angle1_before = switch_directions ? angle2 : angle1;
  const double angle2_before = 180 - (switch_directions ? angle1 : angle2);

  const int sign = switch_directions ? -1 : +1;
  const auto beam1 = get_beam(sign * beamEnergy_SuperKEKB_eMinus); // e- (HER) in +z direction
  const auto beam2 = get_beam(-sign * beamEnergy_SuperKEKB_ePlus); // e+ (LER) in -z direction
  const auto boost = ROOT::Math::Boost((beam1 + beam2).BoostToCM());
  const auto vec1_post = boost(get_vector(angle1_before));
  const auto vec2_post = boost(get_vector(angle2_before));

  const double angle1_post = vec1_post.Theta() * TMath::RadToDeg();
  const double angle2_post = vec2_post.Theta() * TMath::RadToDeg();

  std::cout
    << std::fixed << std::setprecision(1)
    << "Acceptance in lab frame changes from [" << angle1_before << ", " << angle2_before << "] "
       "to [" << angle1_post << ", " << angle2_post << "] when boosted into c.o.m. frame\n"
  ;
  const double cosThetaCut = std::min(
    std::fabs(std::cos(vec1_post.Theta())),
    std::fabs(std::cos(vec2_post.Theta()))
  );
  std::cout << "Use |cos(theta)| < " << std::setprecision(3) << cosThetaCut << " as the acceptance cut\n";

  return EXIT_SUCCESS;
}
