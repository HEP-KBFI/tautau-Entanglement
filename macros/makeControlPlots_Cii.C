
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <math.h>
#include <limits>

TFile*
openFile(const std::string& inputFilePath, const std::string& inputFileName)
{
  std::string inputFileName_full = inputFilePath;
  if ( inputFileName_full.find_last_of("/") != (inputFileName_full.size() - 1) ) inputFileName_full.append("/");
  inputFileName_full.append(inputFileName);
  TFile* inputFile = new TFile(inputFileName_full.data());
  if ( !inputFile ) {
    std::cerr << "Failed to open input file = '" << inputFileName_full << "' !!" << std::endl;
    assert(0);
  }
  return inputFile;
}

TH1*
loadHistogram(TFile* inputFile, const std::string& histogramName)
{
  TH1* histogram = dynamic_cast<TH1*>(inputFile->Get(histogramName.data()));
  if ( !histogram ) {
    std::cerr << "Failed to load histogram = " << histogramName << " from file = " << inputFile->GetName() << " !!" << std::endl;
    assert(0);
  }
  if ( !histogram->GetSumw2N() ) histogram->Sumw2();
  return histogram;
}

void
showHistograms(double canvasSizeX, double canvasSizeY,
               TH1* histogram1, const std::string& legendEntry1,
               TH1* histogram2, const std::string& legendEntry2,
               TH1* histogram3, const std::string& legendEntry3,
               TH1* histogram4, const std::string& legendEntry4,
               TH1* histogram5, const std::string& legendEntry5,
               TH1* histogram6, const std::string& legendEntry6,
               int colors[], int lineStyles[], 
               double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, 
               double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
               bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
               const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetTopMargin(0.05);
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.14);
  canvas->SetRightMargin(0.05);
  canvas->SetLogy(useLogScale);
  
  canvas->SetGridx(1);
  canvas->SetGridy(1);

  if ( !histogram1 ) {
    std::cerr << "<showHistograms>: histogram1 = NULL --> skipping !!" << std::endl;
    return;
  }

  histogram1->SetTitle("");
  histogram1->SetStats(false);
  histogram1->SetMinimum(yMin);
  histogram1->SetMaximum(yMax);

  TAxis* xAxis = histogram1->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleSize(0.045);
  xAxis->SetTitleOffset(xAxisOffset);  
  if ( xMax > xMin ) {
    std::cout << "limiting x-axis range to " << xMin << ".." << xMax << std::endl;
    xAxis->SetRangeUser(xMin, xMax);
  }

  TAxis* yAxis = histogram1->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleSize(0.045);
  yAxis->SetTitleOffset(yAxisOffset);

  histogram1->SetLineColor(colors[0]);
  histogram1->SetLineWidth(2);
  histogram1->SetLineStyle(lineStyles[0]);
  histogram1->Draw("hist");

  if ( histogram2 ) {
    histogram2->SetLineColor(colors[1]);
    histogram2->SetLineWidth(2);
    histogram2->SetLineStyle(lineStyles[1]);
    histogram2->Draw("histsame");
  }

  if ( histogram3 ) {
    histogram3->SetLineColor(colors[2]);
    histogram3->SetLineWidth(2);
    histogram3->SetLineStyle(lineStyles[2]);
    histogram3->Draw("histsame");
  }

  if ( histogram4 ) {
    histogram4->SetLineColor(colors[3]);
    histogram4->SetLineWidth(2);
    histogram4->SetLineStyle(lineStyles[3]);
    histogram4->Draw("histsame");
  }

  if ( histogram5 ) {
    histogram5->SetLineColor(colors[4]);
    histogram5->SetLineWidth(2);
    histogram5->SetLineStyle(lineStyles[4]);
    histogram5->Draw("histsame");
  }

  if ( histogram6 ) {
    histogram6->SetLineColor(colors[5]);
    histogram6->SetLineWidth(2);
    histogram6->SetLineStyle(lineStyles[5]);
    histogram6->Draw("histsame");
  }

  TLegend* legend = 0;
  if ( legendEntry1 != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetTextSize(legendTextSize);
    legend->AddEntry(histogram1, legendEntry1.data(), "l");
    if ( histogram2 ) legend->AddEntry(histogram2, legendEntry2.data(), "l");
    if ( histogram3 ) legend->AddEntry(histogram3, legendEntry3.data(), "l");
    if ( histogram4 ) legend->AddEntry(histogram4, legendEntry4.data(), "l");
    if ( histogram5 ) legend->AddEntry(histogram5, legendEntry5.data(), "l");
    if ( histogram6 ) legend->AddEntry(histogram6, legendEntry6.data(), "l");
    legend->Draw();
  }

  histogram1->Draw("axissame");

  canvas->Update();
  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  //if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  
  delete legend;
  delete canvas;  
}

void
makeControlPlots_Cii()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = "/home/veelken/Entanglement/analysis/SuperKEKB/2023Oct15_wSmearing/plots/";

  std::vector<std::string> modes;
  modes.push_back("gen");
  modes.push_back("startPos");
  modes.push_back("kinFit");

  std::map<std::string, TFile*> inputFiles;
  for ( std::vector<std::string>::const_iterator mode = modes.begin(); mode != modes.end(); ++mode )
  {
    std::string inputFileName = Form("makeControlPlots_dy_lo_pythia8_ext_%sMode_beamAxis_pi_piDecayMode.root", mode->c_str());
    TFile* inputFile = openFile(inputFilePath, inputFileName);
    inputFiles[*mode] = inputFile;
  }

  std::vector<std::string> observables;
  observables.push_back("C_rr");
  observables.push_back("C_nn");
  observables.push_back("C_kk");

  std::map<std::string, int> rebin; // key = observable
  rebin["C_rr"]                =   1;
  rebin["C_nn"]                =   1;
  rebin["C_kk"]                =   1;

  std::map<std::string, double> xMin; // key = observable
  xMin["C_rr"]                 =  -9.;
  xMin["C_nn"]                 =  -9.;
  xMin["C_kk"]                 =  -9.;

  std::map<std::string, double> xMax; // key = observable
  xMax["C_rr"]                 =  +9.;
  xMax["C_nn"]                 =  +9.;
  xMax["C_kk"]                 =  +9.;
  
  std::map<std::string, std::string> xAxisTitles; // key = observable
  xAxisTitles["C_rr"]          = "C_{rr}";
  xAxisTitles["C_nn"]          = "C_{nn}";
  xAxisTitles["C_kk"]          = "C_{kk}";
  
  int colors[6]       = {  1,  2,  8,  4,  6,  7 };
  int lineStyles[6]   = {  1,  1,  1,  1,  1,  1 };
  int markerStyles[6] = { 22, 32, 20, 24, 21, 25 };

  std::vector<std::string> histograms;
  histograms.push_back("Crr");
  histograms.push_back("Cnn");
  histograms.push_back("Ckk");

  for ( std::vector<std::string>::const_iterator observable = observables.begin(); observable != observables.end(); ++observable )
  {
    std::map<std::string, TH1*> histograms;
    for ( std::vector<std::string>::const_iterator mode = modes.begin(); mode != modes.end(); ++mode )
    {
      TFile* inputFile = inputFiles[*mode];
      assert(inputFile);

      TH1* histogram = loadHistogram(inputFile, *observable);
      if ( rebin[*observable] > 1 )
      {
        histogram->Rebin(rebin[*observable]);
      }
      histogram->Scale(1./histogram->Integral());

      histograms[*mode] = histogram;
    }

    std::string outputFileName = Form("makeControlPlots_%s.png", observable->c_str());
    showHistograms(1150, 950,
                   histograms["gen"],      "gen",
                   histograms["startPos"], "startPos",
                   histograms["kinFit"],   "kinFit",
		   0, "",
		   0, "",
		   0, "",
		   colors, lineStyles, 
		   0.040, 0.68, 0.64, 0.22, 0.30, 
		   xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.2, 
		   true, 1.e-4, 1.99, "Events", 1.4, 
		   outputFileName);
  }

  for ( std::vector<std::string>::const_iterator mode = modes.begin(); mode != modes.end(); ++mode )
  {
    delete inputFiles[*mode];
  }
}

