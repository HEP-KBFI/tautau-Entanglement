#include "TauAnalysis/Entanglement/interface/fillWithOverFlow.h"

#include <TH1.h> // TH1
#include <TH2.h> // TH2

const double epsilon = 1.e-6;

void
fill1D(TH1 * histogram,
       double x,
       double evtWeight)
{
  if ( !histogram ) return;
  const TAxis * const xAxis = histogram->GetXaxis();
  if ( x < (xAxis->GetXmin() + epsilon) || x > (xAxis->GetXmax() - epsilon) ) return;
  histogram->Fill(x, evtWeight);
}

void
fillWithOverFlow1D(TH1 * histogram,
                   double x,
                   double evtWeight)
{
  if ( !histogram ) return;
  const TAxis * const xAxis = histogram->GetXaxis();
  double x_bounded = x;
  if ( x_bounded < (xAxis->GetXmin() + epsilon) ) x_bounded = xAxis->GetXmin() + epsilon;
  if ( x_bounded > (xAxis->GetXmax() - epsilon) ) x_bounded = xAxis->GetXmax() - epsilon;
  histogram->Fill(x_bounded, evtWeight);
}

void
fill2D(TH2 * histogram,
       double x,
       double y,
       double evtWeight)
{
  if( !histogram ) return;
  const TAxis * const xAxis = histogram->GetXaxis();
  if ( x < (xAxis->GetXmin() + epsilon) || x > (xAxis->GetXmax() - epsilon) ) return;
  const TAxis * const yAxis = histogram->GetYaxis();
  if ( y < (yAxis->GetXmin() + epsilon) || y > (yAxis->GetXmax() - epsilon) ) return;
  histogram->Fill(x, y, evtWeight);
}

void
fillWithOverFlow2D(TH2 * histogram,
                   double x,
                   double y,
                   double evtWeight)
{
  if( !histogram ) return;
  const TAxis * const xAxis = histogram->GetXaxis();
  double x_bounded = x;
  if ( x_bounded < (xAxis->GetXmin() + epsilon) ) x_bounded = xAxis->GetXmin() + epsilon;
  if ( x_bounded > (xAxis->GetXmax() - epsilon) ) x_bounded = xAxis->GetXmax() - epsilon;
  const TAxis * const yAxis = histogram->GetYaxis();
  double y_bounded = y;
  if ( y_bounded < (yAxis->GetXmin() + epsilon) ) y_bounded = yAxis->GetXmin() + epsilon;
  if ( y_bounded > (yAxis->GetXmax() - epsilon) ) y_bounded = yAxis->GetXmax() - epsilon;
  histogram->Fill(x_bounded, y_bounded, evtWeight);
}
