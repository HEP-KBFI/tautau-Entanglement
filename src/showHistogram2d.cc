#include "TauAnalysis/Entanglement/interface/showHistogram2d.h"

#include <TAxis.h>   // TAxis
#include <TCanvas.h> // TCanvas
#include <TGraph.h>  // TGraph

#include <iostream>  // std::cout
#include <stdlib.h>  // sleep()

void showHistogram2d(int canvasSizeX, int canvasSizeY,
                     TH2* histogram,
                     const std::string& xAxisTitle, double xAxisOffset, 
                     const std::string& yAxisTitle, double yAxisOffset,
                     bool showDiagonal, bool showStatsBox, const std::string& drawOption,
                     const std::string& outputFileName, bool addHistogramName,
                     int verbosity)
{
  if ( verbosity >= 1 )
  {
    std::cout << "<showHistogram2d>:\n";
    std::cout << " outputFileName = '" << outputFileName << "'\n";
  }

  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetTopMargin(0.06);
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.14);
  canvas->SetRightMargin(0.06);
  
  TAxis* xAxis = histogram->GetXaxis();
  double xMin = xAxis->GetXmin();
  double xMax = xAxis->GetXmax();

  TAxis* yAxis = histogram->GetYaxis();
  double yMin = yAxis->GetXmin();
  double yMax = yAxis->GetXmax();

  histogram->SetTitle("");
  histogram->SetStats(showStatsBox);
  xAxis->SetTitle(xAxisTitle.c_str());
  xAxis->SetTitleSize(0.045);
  xAxis->SetTitleOffset(xAxisOffset);  
  yAxis->SetTitle(yAxisTitle.c_str());
  yAxis->SetTitleSize(0.045);
  yAxis->SetTitleOffset(yAxisOffset);

  histogram->Draw(drawOption.c_str());

  TGraph* diagonal = nullptr;
  if ( showDiagonal )
  {
    diagonal = new TGraph(2);
    diagonal->SetPoint(0, xMin, yMin);
    diagonal->SetPoint(1, xMax, yMax);
    diagonal->SetLineColor(2);
    diagonal->SetLineWidth(1);
    diagonal->Draw("L");
  }

  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  if ( addHistogramName )
  {
    outputFileName_plot.append("_");
    outputFileName_plot.append(histogram->GetName());
  }
  canvas->Print(std::string(outputFileName_plot).append(".png").c_str());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").c_str());
  
  // CV: pause for 1 second to reduce load on file system
  sleep(1);

  delete diagonal;
  delete canvas;  
}

