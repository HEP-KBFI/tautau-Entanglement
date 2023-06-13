
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TVectorD.h>
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

TTree*
loadTree(TFile* inputFile, const std::string& treeName)
{
  TTree* tree = dynamic_cast<TTree*>(inputFile->Get(treeName.data()));
  if ( !tree ) {
    std::cerr << "Failed to load tree = " << treeName << " from file = " << inputFile->GetName() << " !!" << std::endl;
    assert(0);
  }
  return tree;
}

double
getVectorElement(TTree* tree, const std::string& branchName, size_t idx)
{
  TVectorD* B = new TVectorD(3);
  tree->SetBranchAddress(branchName.c_str(), &B);
  tree->GetEntry(0);
  double retVal = (*B)[idx];
  delete B;
  return retVal;
}

double
getMatrixElement(TTree* tree, const std::string& branchName, size_t idxRow, size_t idxColumn)
{
  TMatrixD* C = new TMatrixD(3,3);
  tree->SetBranchAddress(branchName.c_str(), &C);
  tree->GetEntry(0);
  double retVal = (*C)[idxRow][idxColumn];
  delete C;
  return retVal;
}

void 
showGraphs(double canvasSizeX, double canvasSizeY,
           TGraph* graph1, const std::string& legendEntry1,
           TGraph* graph2, const std::string& legendEntry2,
           TGraph* graph3, const std::string& legendEntry3,
           TGraph* graph4, const std::string& legendEntry4,
           TGraph* graph5, const std::string& legendEntry5,
           TGraph* graph6, const std::string& legendEntry6,
           int colors[], int markerStyles[], int lineStyles[], 
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
  canvas->SetBottomMargin(0.12);
  canvas->SetRightMargin(0.05);
  canvas->SetLogy(useLogScale);
  
  canvas->SetGridx(1);
  canvas->SetGridy(1);

  if ( !graph1 ) {
    std::cerr << "<showGraphs>: graph1 = NULL --> skipping !!" << std::endl;
    return;
  }

  TH1* dummyHistogram = new TH1F("dummyHistogram", "dummyHistogram", 10, xMin, xMax);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  dummyHistogram->SetMinimum(yMin);
  dummyHistogram->SetMaximum(yMax);

  TAxis* xAxis = dummyHistogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleSize(0.045);
  xAxis->SetTitleOffset(xAxisOffset);  
  if ( xMax > xMin ) {
    std::cout << "limiting x-axis range to " << xMin << ".." << xMax << std::endl;
    xAxis->SetRangeUser(xMin, xMax);
  }

  TAxis* yAxis = dummyHistogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleSize(0.045);
  yAxis->SetTitleOffset(yAxisOffset);

  dummyHistogram->Draw();
  //dummyHistogram->Draw("axis");

  graph1->SetMarkerColor(colors[0]);
  graph1->SetMarkerSize(2);
  graph1->SetMarkerStyle(markerStyles[0]);
  graph1->SetLineColor(colors[0]);
  graph1->SetLineWidth(2);
  graph1->SetLineStyle(lineStyles[0]);
  graph1->Draw("PL");

  if ( graph2 ) {
    graph2->SetMarkerColor(colors[1]);
    graph2->SetMarkerSize(2);
    graph2->SetMarkerStyle(markerStyles[1]);
    graph2->SetLineColor(colors[1]);
    graph2->SetLineWidth(2);
    graph2->SetLineStyle(lineStyles[1]);
    graph2->Draw("PL");
  }

  if ( graph3 ) {
    graph3->SetMarkerColor(colors[2]);
    graph3->SetMarkerSize(2);
    graph3->SetMarkerStyle(markerStyles[2]);
    graph3->SetLineColor(colors[2]);
    graph3->SetLineWidth(2);
    graph3->SetLineStyle(lineStyles[2]);
    graph3->Draw("PL");
  }

  if ( graph4 ) {
    graph4->SetMarkerColor(colors[3]);
    graph4->SetMarkerSize(2);
    graph4->SetMarkerStyle(markerStyles[3]);
    graph4->SetLineColor(colors[3]);
    graph4->SetLineWidth(2);
    graph4->SetLineStyle(lineStyles[3]);
    graph4->Draw("PL");
  }

  if ( graph5 ) {
    graph5->SetMarkerColor(colors[4]);
    graph5->SetMarkerSize(2);
    graph5->SetMarkerStyle(markerStyles[4]);
    graph5->SetLineColor(colors[4]);
    graph5->SetLineWidth(2);
    graph5->SetLineStyle(lineStyles[4]);
    graph5->Draw("PL");
  }

  if ( graph6 ) {
    graph6->SetMarkerColor(colors[5]);
    graph6->SetMarkerSize(2);
    graph6->SetMarkerStyle(markerStyles[5]);
    graph6->SetLineColor(colors[5]);
    graph6->SetLineWidth(2);
    graph6->SetLineStyle(lineStyles[5]);
    graph6->Draw("PL");
  }

  TLegend* legend = 0;
  if ( legendEntry1 != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetTextSize(legendTextSize);
    legend->AddEntry(graph1, legendEntry1.data(), "p");
    if ( graph2 ) legend->AddEntry(graph2, legendEntry2.data(), "L");
    if ( graph3 ) legend->AddEntry(graph3, legendEntry3.data(), "L");
    if ( graph4 ) legend->AddEntry(graph4, legendEntry4.data(), "L");
    if ( graph5 ) legend->AddEntry(graph5, legendEntry5.data(), "L");
    if ( graph6 ) legend->AddEntry(graph6, legendEntry6.data(), "L");
    legend->Draw();
  }

  dummyHistogram->Draw("axissame");

  canvas->RedrawAxis();

  canvas->Update();
  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  //if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());

  delete dummyHistogram;
  delete legend;
  delete canvas;  
}

void
makeUncertaintyPlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = "/scratch/persistent/veelken/Entanglement/studyRecoEffects/2023Jun02/";
  std::string inputFileName_baseline = "analyzeEntanglementNtuple2_ggH_htt_pythia8_beamAxis_maxSumPhotonEnLt5.root";
  TFile* inputFile_baseline = openFile(inputFilePath, inputFileName_baseline);
  
  std::string inputFileName_minVisTauPtGt10 = "analyzeEntanglementNtuple2_ggH_htt_pythia8_beamAxis_minVisTauPtGt10.root";
  TFile* inputFile_minVisTauPtGt10 = openFile(inputFilePath, inputFileName_minVisTauPtGt10);
  std::string inputFileName_minVisTauPtGt20 = "analyzeEntanglementNtuple2_ggH_htt_pythia8_beamAxis_minVisTauPtGt20.root";
  TFile* inputFile_minVisTauPtGt20 = openFile(inputFilePath, inputFileName_minVisTauPtGt20);
  std::string inputFileName_minVisTauPtGt30 = "analyzeEntanglementNtuple2_ggH_htt_pythia8_beamAxis_minVisTauPtGt30.root";
  TFile* inputFile_minVisTauPtGt30 = openFile(inputFilePath, inputFileName_minVisTauPtGt30);
  std::string inputFileName_minVisTauPtGt40 = "analyzeEntanglementNtuple2_ggH_htt_pythia8_beamAxis_minVisTauPtGt40.root";
  TFile* inputFile_minVisTauPtGt40 = openFile(inputFilePath, inputFileName_minVisTauPtGt40);

  const std::string treeName = "fitResult";
  
  int colors[6]       = {  1,  2,  8,  4,  6,  7 };
  int lineStyles[6]   = {  1,  1,  1,  1,  1,  1 };
  int markerStyles[6] = { 22, 32, 20, 24, 21, 25 };

  TTree* tree_baseline = loadTree(inputFile_baseline, treeName);

  TTree* tree_minVisTauPtGt10 = loadTree(inputFile_minVisTauPtGt10, treeName);
  TTree* tree_minVisTauPtGt20 = loadTree(inputFile_minVisTauPtGt20, treeName);
  TTree* tree_minVisTauPtGt30 = loadTree(inputFile_minVisTauPtGt30, treeName);
  TTree* tree_minVisTauPtGt40 = loadTree(inputFile_minVisTauPtGt40, treeName);

  const size_t idxC_rr = 0;
  const size_t idxC_nn = 1;
  const size_t idxC_kk = 2;
    
  TGraph* graph_C_rr_fitValue = new TGraph(5);
  graph_C_rr_fitValue->SetPoint(0,  0., getMatrixElement(tree_baseline,        "C", idxC_rr, idxC_rr));
  graph_C_rr_fitValue->SetPoint(1, 10., getMatrixElement(tree_minVisTauPtGt10, "C", idxC_rr, idxC_rr));
  graph_C_rr_fitValue->SetPoint(2, 20., getMatrixElement(tree_minVisTauPtGt20, "C", idxC_rr, idxC_rr));
  graph_C_rr_fitValue->SetPoint(3, 30., getMatrixElement(tree_minVisTauPtGt30, "C", idxC_rr, idxC_rr));
  graph_C_rr_fitValue->SetPoint(4, 40., getMatrixElement(tree_minVisTauPtGt40, "C", idxC_rr, idxC_rr));

  TGraph* graph_C_nn_fitValue = new TGraph(5);
  graph_C_nn_fitValue->SetPoint(0,  0., getMatrixElement(tree_baseline,        "C", idxC_nn, idxC_nn));
  graph_C_nn_fitValue->SetPoint(1, 10., getMatrixElement(tree_minVisTauPtGt10, "C", idxC_nn, idxC_nn));
  graph_C_nn_fitValue->SetPoint(2, 20., getMatrixElement(tree_minVisTauPtGt20, "C", idxC_nn, idxC_nn));
  graph_C_nn_fitValue->SetPoint(3, 30., getMatrixElement(tree_minVisTauPtGt30, "C", idxC_nn, idxC_nn));
  graph_C_nn_fitValue->SetPoint(4, 40., getMatrixElement(tree_minVisTauPtGt40, "C", idxC_nn, idxC_nn));

  TGraph* graph_C_kk_fitValue = new TGraph(5);
  graph_C_kk_fitValue->SetPoint(0,  0., getMatrixElement(tree_baseline,        "C", idxC_kk, idxC_kk));
  graph_C_kk_fitValue->SetPoint(1, 10., getMatrixElement(tree_minVisTauPtGt10, "C", idxC_kk, idxC_kk));
  graph_C_kk_fitValue->SetPoint(2, 20., getMatrixElement(tree_minVisTauPtGt20, "C", idxC_kk, idxC_kk));
  graph_C_kk_fitValue->SetPoint(3, 30., getMatrixElement(tree_minVisTauPtGt30, "C", idxC_kk, idxC_kk));
  graph_C_kk_fitValue->SetPoint(4, 40., getMatrixElement(tree_minVisTauPtGt40, "C", idxC_kk, idxC_kk));

  std::string outputFileName_C = Form("makeUncertaintyPlots_C.png");
  showGraphs(1150, 850,
             graph_C_rr_fitValue, "C_{rr}",
             graph_C_nn_fitValue, "C_{nn}",
             graph_C_kk_fitValue, "C_{kk}",
             nullptr, "",
             nullptr, "",
             nullptr, "",
             colors, markerStyles, lineStyles, 
             0.040, 0.76, 0.72, 0.16, 0.21, 
             -5., 45., "p_{T} Threshold [GeV]", 1.2, 
             false, -1.2, +2.8, "C", 1.4, 
             outputFileName_C);

  TGraph* graph_C_rr_fitError = new TGraph(5);
  graph_C_rr_fitError->SetPoint(0,  0., getMatrixElement(tree_baseline,        "CErr", idxC_rr, idxC_rr));
  graph_C_rr_fitError->SetPoint(1, 10., getMatrixElement(tree_minVisTauPtGt10, "CErr", idxC_rr, idxC_rr));
  graph_C_rr_fitError->SetPoint(2, 20., getMatrixElement(tree_minVisTauPtGt20, "CErr", idxC_rr, idxC_rr));
  graph_C_rr_fitError->SetPoint(3, 30., getMatrixElement(tree_minVisTauPtGt30, "CErr", idxC_rr, idxC_rr));
  graph_C_rr_fitError->SetPoint(4, 40., getMatrixElement(tree_minVisTauPtGt40, "CErr", idxC_rr, idxC_rr));

  TGraph* graph_C_nn_fitError = new TGraph(5);
  graph_C_nn_fitError->SetPoint(0,  0., getMatrixElement(tree_baseline,        "CErr", idxC_nn, idxC_nn));
  graph_C_nn_fitError->SetPoint(1, 10., getMatrixElement(tree_minVisTauPtGt10, "CErr", idxC_nn, idxC_nn));
  graph_C_nn_fitError->SetPoint(2, 20., getMatrixElement(tree_minVisTauPtGt20, "CErr", idxC_nn, idxC_nn));
  graph_C_nn_fitError->SetPoint(3, 30., getMatrixElement(tree_minVisTauPtGt30, "CErr", idxC_nn, idxC_nn));
  graph_C_nn_fitError->SetPoint(4, 40., getMatrixElement(tree_minVisTauPtGt40, "CErr", idxC_nn, idxC_nn));

  TGraph* graph_C_kk_fitError = new TGraph(5);
  graph_C_kk_fitError->SetPoint(0,  0., getMatrixElement(tree_baseline,        "CErr", idxC_kk, idxC_kk));
  graph_C_kk_fitError->SetPoint(1, 10., getMatrixElement(tree_minVisTauPtGt10, "CErr", idxC_kk, idxC_kk));
  graph_C_kk_fitError->SetPoint(2, 20., getMatrixElement(tree_minVisTauPtGt20, "CErr", idxC_kk, idxC_kk));
  graph_C_kk_fitError->SetPoint(3, 30., getMatrixElement(tree_minVisTauPtGt30, "CErr", idxC_kk, idxC_kk));
  graph_C_kk_fitError->SetPoint(4, 40., getMatrixElement(tree_minVisTauPtGt40, "CErr", idxC_kk, idxC_kk));

  std::string outputFileName_CErr = Form("makeUncertaintyPlots_CErr.png");
  showGraphs(1150, 850,
             graph_C_rr_fitError, "C_{rr}",
             graph_C_nn_fitError, "C_{nn}",
             graph_C_kk_fitError, "C_{kk}",
             nullptr, "",
             nullptr, "",
             nullptr, "",
             colors, markerStyles, lineStyles, 
             0.040, 0.22, 0.72, 0.16, 0.21, 
             -5., 45., "p_{T} Threshold [GeV]", 1.2, 
             false, 0., 0.15, "#sigma_{C}", 1.4, 
             outputFileName_CErr);

  delete inputFile_baseline;

  delete inputFile_minVisTauPtGt10;
  delete inputFile_minVisTauPtGt20;
  delete inputFile_minVisTauPtGt30;
  delete inputFile_minVisTauPtGt40;
}

