
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

double
compIntegral(const TH1* histogram, bool includeUnderflowBin = false, bool includeOverflowBin = false)
{
  int firstBin = 1;
  if ( includeUnderflowBin ) firstBin -= 1;
  int lastBin = histogram->GetNbinsX();
  if ( includeOverflowBin ) lastBin += 1;
  double integral = 0.;
  for ( int idxBin = firstBin; idxBin <= lastBin; ++idxBin ) {
    integral += histogram->GetBinContent(idxBin);
  }
  return integral;
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
makeRecoEffectPlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = "/scratch/persistent/veelken/Entanglement/studyRecoEffects/2023Jun02/";
  std::string inputFileName_baseline = "makeControlPlots_ggH_htt_pythia8_beamAxis_maxSumPhotonEnLt5.root";
  TFile* inputFile_baseline = openFile(inputFilePath, inputFileName_baseline);
  
  std::string inputFileName_minVisTauPtGt10 = "makeControlPlots_ggH_htt_pythia8_beamAxis_minVisTauPtGt10.root";
  TFile* inputFile_minVisTauPtGt10 = openFile(inputFilePath, inputFileName_minVisTauPtGt10);
  std::string inputFileName_minVisTauPtGt20 = "makeControlPlots_ggH_htt_pythia8_beamAxis_minVisTauPtGt20.root";
  TFile* inputFile_minVisTauPtGt20 = openFile(inputFilePath, inputFileName_minVisTauPtGt20);
  std::string inputFileName_minVisTauPtGt30 = "makeControlPlots_ggH_htt_pythia8_beamAxis_minVisTauPtGt30.root";
  TFile* inputFile_minVisTauPtGt30 = openFile(inputFilePath, inputFileName_minVisTauPtGt30);
  std::string inputFileName_minVisTauPtGt40 = "makeControlPlots_ggH_htt_pythia8_beamAxis_minVisTauPtGt40.root";
  TFile* inputFile_minVisTauPtGt40 = openFile(inputFilePath, inputFileName_minVisTauPtGt40);

  std::string inputFileName_noNeutralKaonCut = "makeControlPlots_ggH_htt_pythia8_beamAxis_noNeutralKaonCut.root";
  TFile* inputFile_noNeutralKaonCut = openFile(inputFilePath, inputFileName_noNeutralKaonCut);

  std::vector<std::string> observables;
  observables.push_back("Bp_n");
  observables.push_back("Bp_r");
  observables.push_back("Bp_k");
  observables.push_back("Bm_n");
  observables.push_back("Bm_r");
  observables.push_back("Bm_k");
  observables.push_back("C_rr");
  observables.push_back("C_nn");
  observables.push_back("C_kk");

  std::map<std::string, int> rebin; // key = observable
  rebin["Bp_n"]                =   1;
  rebin["Bp_r"]                =   1;
  rebin["Bp_k"]                =   1;
  rebin["Bm_n"]                =   1;
  rebin["Bm_r"]                =   1;
  rebin["Bm_k"]                =   1;
  rebin["C_rr"]                =   1;
  rebin["C_nn"]                =   1;
  rebin["C_kk"]                =   1;

  std::map<std::string, double> xMin; // key = observable
  xMin["Bp_n"]                 =  -3.;
  xMin["Bp_r"]                 =  -3.;
  xMin["Bp_k"]                 =  -3.;
  xMin["Bm_n"]                 =  -3.;
  xMin["Bm_r"]                 =  -3.;
  xMin["Bm_k"]                 =  -3.;
  xMin["C_rr"]                 =  -9.;
  xMin["C_nn"]                 =  -9.;
  xMin["C_kk"]                 =  -9.;

  std::map<std::string, double> xMax; // key = observable
  xMax["Bp_n"]                 =  +3.;
  xMax["Bp_r"]                 =  +3.;
  xMax["Bp_k"]                 =  +3.;
  xMax["Bm_n"]                 =  +3.;
  xMax["Bm_r"]                 =  +3.;
  xMax["Bm_k"]                 =  +3.;
  xMax["C_rr"]                 =  +9.;
  xMax["C_nn"]                 =  +9.;
  xMax["C_kk"]                 =  +9.;
  
  std::map<std::string, std::string> xAxisTitles; // key = observable
  xAxisTitles["Bp_n"]          = "B^{+}_{n}";
  xAxisTitles["Bp_r"]          = "B^{+}_{r}";
  xAxisTitles["Bp_k"]          = "B^{+}_{k}";
  xAxisTitles["Bm_n"]          = "B^{-}_{n}";
  xAxisTitles["Bm_r"]          = "B^{-}_{r}";
  xAxisTitles["Bm_k"]          = "B^{-}_{k}";
  xAxisTitles["C_rr"]          = "C_{rr}";
  xAxisTitles["C_nn"]          = "C_{nn}";
  xAxisTitles["C_kk"]          = "C_{kk}";
  
  int colors[6]       = {  1,  2,  8,  4,  6,  7 };
  int lineStyles[6]   = {  1,  1,  1,  1,  1,  1 };
  int markerStyles[6] = { 22, 32, 20, 24, 21, 25 };

  for ( std::vector<std::string>::const_iterator observable = observables.begin();
	observable != observables.end(); ++observable ) {
    std::string histogramName = observable->c_str();

    TH1* histogram_baseline = loadHistogram(inputFile_baseline, histogramName);

    TH1* histogram_minVisTauPtGt10 = loadHistogram(inputFile_minVisTauPtGt10, histogramName);
    TH1* histogram_minVisTauPtGt20 = loadHistogram(inputFile_minVisTauPtGt20, histogramName);
    TH1* histogram_minVisTauPtGt30 = loadHistogram(inputFile_minVisTauPtGt30, histogramName);
    TH1* histogram_minVisTauPtGt40 = loadHistogram(inputFile_minVisTauPtGt40, histogramName);

    TH1* histogram_noNeutralKaonCut = loadHistogram(inputFile_noNeutralKaonCut, histogramName);

    double normFactor = 1./histogram_baseline->Integral();
    histogram_baseline->Scale(normFactor);

    histogram_minVisTauPtGt10->Scale(normFactor);
    histogram_minVisTauPtGt20->Scale(normFactor);
    histogram_minVisTauPtGt30->Scale(normFactor);
    histogram_minVisTauPtGt40->Scale(normFactor);

    histogram_noNeutralKaonCut->Scale(normFactor);

    if ( rebin[*observable] > 1 )
    {
      histogram_baseline->Rebin();

      histogram_minVisTauPtGt10->Rebin();
      histogram_minVisTauPtGt20->Rebin();
      histogram_minVisTauPtGt30->Rebin();
      histogram_minVisTauPtGt40->Rebin();

      histogram_noNeutralKaonCut->Rebin();
    }

    std::string outputFileName_minVisTauPt = Form("makeRecoEffectPlots_%s_minVisTauPt.png", observable->c_str());
    showHistograms(1150, 950,
                   histogram_baseline,        "No cut",
                   histogram_minVisTauPtGt10, "p_{T}^{vis} > 10 GeV",
                   histogram_minVisTauPtGt20, "p_{T}^{vis} > 20 GeV",
                   histogram_minVisTauPtGt30, "p_{T}^{vis} > 30 GeV",
                   histogram_minVisTauPtGt40, "p_{T}^{vis} > 40 GeV",
		   0, "",
		   colors, lineStyles, 
		   0.040, 0.68, 0.64, 0.22, 0.30, 
		   xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.2, 
		   true, 1.e-4, 1.99, "Events", 1.4, 
		   outputFileName_minVisTauPt);

    std::string outputFileName_noNeutralKaonCut = Form("makeRecoEffectPlots_%s_noNeutralKaonCut.png", observable->c_str());
    showHistograms(1150, 950,
                   histogram_baseline,         "Num. K^{0} = 0",
                   histogram_noNeutralKaonCut, "No cut",
		   0, "",
		   0, "",
		   0, "",
		   0, "",
		   colors, lineStyles, 
		   0.040, 0.68, 0.78, 0.22, 0.14, 
		   xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.2, 
		   true, 1.e-4, 1.99, "Events", 1.4, 
		   outputFileName_noNeutralKaonCut);
  }

  delete inputFile_baseline;

  delete inputFile_minVisTauPtGt10;
  delete inputFile_minVisTauPtGt20;
  delete inputFile_minVisTauPtGt30;
  delete inputFile_minVisTauPtGt40;

  delete inputFile_noNeutralKaonCut;
}

