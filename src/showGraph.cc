#include "TauAnalysis/Entanglement/interface/showGraph.h"

#include <TCanvas.h> // TCanvas
#include <TH1.h>     // TH1D

#include <iostream>  // std::cout
#include <stdlib.h>  // sleep()

void
showGraph(double canvasSizeX, double canvasSizeY,
          TGraph* graph,
          const std::string& xAxisTitle, double xMin, double xMax,
          const std::string& yAxisTitle, double yMin, double yMax,
          const std::string& drawingOption,
          const std::string& outputFileName,
          int verbosity)
{
  if ( verbosity >= 1 )
  {
    std::cout << "<showGraph>:\n";
    std::cout << " outputFileName = '" << outputFileName << "'\n";
  }

  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.12);

  graph->SetMarkerStyle(8);
  graph->SetMarkerSize(1);
  graph->SetMarkerColor(8);
  graph->SetLineStyle(8);
  graph->SetLineWidth(1);
  graph->SetLineColor(8);

  TH1* dummyHistogram = new TH1D("dummyHistogram", "dummyHistogram", graph->GetN(), xMin, xMax);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  dummyHistogram->SetMinimum(yMin);
  dummyHistogram->SetMaximum(yMax);
  dummyHistogram->GetXaxis()->SetTitle("Iteration");
  dummyHistogram->GetYaxis()->SetTitle("#chi^{2}");
  dummyHistogram->Draw("axis");

  graph->Draw(drawingOption.c_str());

  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  canvas->Print(std::string(outputFileName_plot).append(".png").c_str());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").c_str());

  // CV: pause for 1 second to reduce load on file system
  sleep(1);

  delete dummyHistogram;
  delete canvas;
}
