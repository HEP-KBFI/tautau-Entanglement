#include "TauAnalysis/Entanglement/interface/printPoint.h"

void
printPoint(const std::string& label,
           const reco::Candidate::Point& p)
{
  std::cout << label << ":";
  std::cout << " x = " << p.x() << ", y = " << p.y() << ", z = " << p.z() << "\n";
}
