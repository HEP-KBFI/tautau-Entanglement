#include "TauAnalysis/Entanglement/interface/bookHistogram1d.h"

TH1*
bookHistogram1d(fwlite::TFileService& fs, const std::string& name, int numBinsX, double xMin, double xMax)
{
  TH1* histogram = fs.make<TH1D>(name.c_str(), name.c_str(), numBinsX, xMin, xMax);
  histogram->Sumw2();
  return histogram;
}
