#include "TauAnalysis/Entanglement/interface/showHistogram1d.h"

#include <TAxis.h>   // TAxis
#include <TCanvas.h> // TCanvas

#include <stdlib.h>  // sleep()

void showHistogram1d(int canvasSizeX, int canvasSizeY,
                     TH1* histogram, 
                     const std::string& xAxisTitle, double xAxisOffset, 
                     bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
                     bool showStatsBox, const std::string& drawOption,
                     const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2); 
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.12);
  canvas->SetLogy(useLogScale);
  
  histogram->SetTitle("");
  histogram->SetStats(showStatsBox);
  if ( yMax > yMin )
  {
    histogram->SetMinimum(yMin*histogram->Integral());
    histogram->SetMaximum(yMax*histogram->Integral());
  }

  TAxis* xAxis = histogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleSize(0.045);
  xAxis->SetTitleOffset(xAxisOffset);  

  TAxis* yAxis = histogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleSize(0.045);
  yAxis->SetTitleOffset(yAxisOffset);

  histogram->SetLineColor(1);
  histogram->SetLineWidth(1);
  histogram->SetLineStyle(1);
  histogram->SetMarkerColor(1);
  histogram->SetMarkerSize(1);
  histogram->SetMarkerStyle(8);
  histogram->Draw(drawOption.c_str());

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
