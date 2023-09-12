#include "TauAnalysis/Entanglement/interface/scaleHistogram.h"

void
scaleHistogram(TH1* histogram, double avEvtWeight)
{
  if ( avEvtWeight > 0. )
  {
    histogram->Scale(1./avEvtWeight);
  }
}
