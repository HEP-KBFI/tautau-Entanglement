#include "TauAnalysis/Entanglement/interface/showHistogram2d.h"

#include <TAxis.h>   // TAxis
#include <TCanvas.h> // TCanvas

#include <stdlib.h>  // sleep()

void showHistogram2d(TH2* histogram, 
                     const std::string& xAxisTitle, double xAxisOffset, 
                     const std::string& yAxisTitle, double yAxisOffset,
                     double avEvtWeight,
                     bool showStatsBox,
                     const std::string& outputFileName)
{
  //double integral = histogram->Integral();
  //if ( integral > 0. ) histogram->Scale(1./integral);
  if ( avEvtWeight > 0. )
  {
    histogram->Scale(1./avEvtWeight);
  }

  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2); 
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.12);
  
  histogram->SetTitle("");
  histogram->SetStats(showStatsBox);

  TAxis* xAxis = histogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.c_str());
  xAxis->SetTitleSize(0.045);
  xAxis->SetTitleOffset(xAxisOffset);  

  TAxis* yAxis = histogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.c_str());
  yAxis->SetTitleSize(0.045);
  yAxis->SetTitleOffset(yAxisOffset);

  histogram->Draw("BOX");

  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  outputFileName_plot.append("_");
  outputFileName_plot.append(histogram->GetName());
  canvas->Print(std::string(outputFileName_plot).append(".png").c_str());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").c_str());
  
  // CV: pause for 1 second to reduce load on file system
  sleep(1);

  delete canvas;  
}

