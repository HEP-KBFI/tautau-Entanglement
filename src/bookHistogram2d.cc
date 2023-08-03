#include "TauAnalysis/Entanglement/interface/bookHistogram2d.h"

TH2*
bookHistogram2d(fwlite::TFileService& fs, const std::string& name, int numBinsX, double xMin, double xMax, int numBinsY, double yMin, double yMax)
{
  TH2* histogram = fs.make<TH2D>(name.c_str(), name.c_str(), numBinsX, xMin, xMax, numBinsY, yMin, yMax);
  histogram->Sumw2();
  return histogram;
}
