
#include "FWCore/ParameterSet/interface/ParameterSet.h"                    // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/cmsException.h"               // cmsException
#include "TauAnalysis/Entanglement/interface/fillWithOverFlow.h"           // fillWithOverFlow2D()
#include "TauAnalysis/Entanglement/interface/KinFitAlgo.h"                 // KinFitAlgo, KinFitSummary
#include "TauAnalysis/Entanglement/interface/KinFitConstraintBase.h"       // KinFitConstraintBase
#include "TauAnalysis/Entanglement/interface/KinFitConstraintTesterBase.h" // KinFitConstraintTesterBase
#include "TauAnalysis/Entanglement/interface/printMatrix.h"                // printMatrix()
#include "TauAnalysis/Entanglement/interface/printVector.h"                // printVector()
#include "TauAnalysis/Entanglement/interface/showGraph.h"                  // showGraph()
#include "TauAnalysis/Entanglement/interface/square.h"                     // square()

#include <TAxis.h>                                                         // TAxis
#include <TBenchmark.h>                                                    // TBenchmark
#include <TCanvas.h>                                                       // TCanvas
#include <TF1.h>                                                           // TF1
#include <TF2.h>                                                           // TF2
#include <TGraph.h>                                                        // TGraph
#include <TH1.h>                                                           // TH1D
#include <TH2.h>                                                           // TH2D
#include <TMath.h>                                                         // TMath::Pi()
#include <TMatrixD.h>                                                      // TMatrixD
#include <TRandom3.h>                                                      // TRandom3
#include <TString.h>                                                       // TString
#include <TVectorD.h>                                                      // TVectorD

#include <cmath>                                                           // std::isnan(), std::sqrt()
#include <sstream>                                                         // std::ostringstream
#include <string>                                                          // std::string
#include <utility>                                                         // std::pair
#include <vector>                                                          // std::vector

typedef std::pair<double, double> point2d;

template <unsigned int rank>
class PolynomialConstraint : public KinFitConstraintBase
{
 public:
  PolynomialConstraint(const std::vector<point2d>& points, int verbosity = -1)
    : KinFitConstraintBase(2, 1, 0, verbosity)
    , points_(points)
    , coeff_(rank + 1)
  {
    if ( verbosity_ >= 1 )
    {
      std::cout << "<PolynomialConstraint::PolynomialConstraint>:\n";
      for ( size_t idxPoint = 0; idxPoint < points.size(); ++idxPoint )
      {
        const point2d& point = points.at(idxPoint);
        std::cout << "point #" << idxPoint << ": x = " << point.first << ", y = " << point.second << "\n";
      }
    }

    if ( points.size() != (rank + 1) )
      throw cmsException("PolynomialConstraint", __LINE__) 
        << "Given number of points = " << points.size() << " does not match rank + 1 = " << rank + 1 << " !!\n";
    
    TMatrixD M(rank + 1, rank + 1);
    TVectorD y(rank + 1);
    for ( unsigned int idxRow = 0; idxRow < (rank + 1); ++idxRow )
    {
      const point2d& point = points.at(idxRow);
      for ( unsigned int idxColumn = 0; idxColumn < (rank + 1); ++idxColumn )
      {
        M(idxRow,idxColumn) = pow(point.first, (rank + 1) - (idxColumn + 1));
      }
      y(idxRow) = point.second;
    }

    if ( verbosity_ >= 1 )
    {
      printMatrix("M", M);
    }

    if ( M.Determinant() == 0. )
    {
      throw cmsException("PolynomialConstraint", __LINE__) 
        << "Failed to invert matrix M !!\n";
    }
    TMatrixD Minv(TMatrixD::kInverted, M);
    if ( verbosity_ >= 1 )
    {
      printMatrix("Minv", Minv);
    }

    coeff_ = Minv*y;
    if ( verbosity_ >= 1 )
    {
      printVector("coeff", coeff_);
    }
  }
  ~PolynomialConstraint()
  {}

  operator std::string() const
  {
    std::ostringstream os;
    for ( unsigned int idx = 0; idx < (rank + 1); ++idx )
    {
      if ( idx > 0 ) os << " + ";
      os << coeff_(idx) << "*x^" << (rank - idx);
    }
    return os.str();
  }

  void
  set_alphaA(const TVectorD& alphaA)
  {
    assert(alphaA.GetNrows() == 2);
    double x = alphaA(0);
    double y = alphaA(1);
    
    D_eq_(0,0) = get_derrivative(x);
    D_eq_(0,1) = -1.;
    d_eq_(0) = get_polynomial(x) - y;
    d_eq_metric_(0,0) = 1.;
  }

 private:
  double 
  get_derrivative(double x) const
  {
    // CV: for a polynomial of rank 3,
    //     this function returns 3*coeff[0]*x^2 + 2*coeff[1]*x + coeff[2]
    double retVal = 0.;
    for ( unsigned int idx = 0; idx < rank; ++idx )
    {
      retVal += (rank - idx)*coeff_(idx)*pow(x, rank - (idx + 1));
    }
    return retVal;
  }

  double 
  get_polynomial(double x) const
  {
    // CV: for a polynomial of rank 3,
    //     this function returns coeff[0]*x^3 + coeff[1]*x^2 + coeff[2]*x + coeff[3]
    double retVal = 0.;
    for ( unsigned int idx = 0; idx < (rank + 1); ++idx )
    {
      retVal += coeff_(idx)*pow(x, rank - idx);
    }
    return retVal;
  }

  unsigned int rank_;

  std::vector<point2d> points_;

  TVectorD coeff_;
};

class PolynomialConstraintTester : public KinFitConstraintTesterBase
{
 public:
  PolynomialConstraintTester(const TVectorD& alpha0, int verbosity = -1)
    : KinFitConstraintTesterBase(alpha0, verbosity)
  {}

  void
  operator()(KinFitConstraintBase& constraint, const std::string& outputFileName)
  {
    const unsigned int Np = constraint.get_Np();
    const unsigned int Nc_eq = constraint.get_Nc_eq();

    rndMean_.ResizeTo(Np); 
    rndMean_ = alpha0_;
    
    rndWidth_.ResizeTo(Np);
    for ( unsigned int idxParameter = 0; idxParameter < Np; ++idxParameter )
    {
      rndWidth_(idxParameter) = 2.;
    }

    plotRange_.ResizeTo(Nc_eq,Np);
    for ( unsigned int idxConstraint = 0; idxConstraint < Nc_eq; ++idxConstraint )
    {
      for ( unsigned int idxParameter = 0; idxParameter < Np; ++idxParameter )
      {
        plotRange_(idxConstraint, idxParameter) = 5.;
      }
    }

    KinFitConstraintTesterBase::operator()(constraint, outputFileName);
  }
};

TVectorD
get_startPos(const point2d& point)
{
  TVectorD startPos(2);
  startPos(0) = point.first;
  startPos(1) = point.second;
  return startPos;
}

template <unsigned int rank>
void 
showFit(const TVectorD& alpha0, const TMatrixD& V_alpha0, const PolynomialConstraint<rank>& constraint, 
        const std::vector<KinFitSummary>& fitHistory,
        const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", 900, 800);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2); 
  canvas->SetTopMargin(0.05);
  canvas->SetLeftMargin(0.08);
  canvas->SetBottomMargin(0.07);
  canvas->SetRightMargin(0.14);

  if ( V_alpha0.Determinant() == 0. )
  {
    throw cmsException("showFit", __LINE__) 
      << "Failed to invert matrix V_alpha0 !!\n";
  }
  TMatrixD Vinv_alpha0(TMatrixD::kInverted, V_alpha0);
  printMatrix("Vinv_alpha0", Vinv_alpha0);

  double V_alpha0_det = V_alpha0.Determinant();
  std::cout << "det(V_alpha0) = " << V_alpha0_det << "\n";

  const double xMin = -10.;
  const double xMax = +10.;
  const double yMin = -10.;
  const double yMax = +10.;

  std::string chi2_formula = "(x - [1])*[3]*(x - [1]) + (x - [1])*[4]*(y - [2]) + (y - [2])*[5]*(x - [1]) + (y - [2])*[6]*(y - [2])";
  TF2* chi2_function = new TF2("chi2", chi2_formula.c_str(), xMin, xMax, yMin, yMax);
  chi2_function->SetParameter(0, 1./(2.*TMath::Pi()*V_alpha0_det));
  chi2_function->SetParameter(1, alpha0(0));
  chi2_function->SetParameter(2, alpha0(1));
  chi2_function->SetParameter(3, Vinv_alpha0(0,0));
  chi2_function->SetParameter(4, Vinv_alpha0(0,1));
  chi2_function->SetParameter(5, Vinv_alpha0(1,0));
  chi2_function->SetParameter(6, Vinv_alpha0(1,1));

  chi2_function->SetTitle("");
  chi2_function->Draw("Colz");

  TF1* constraint_function = new TF1("constraint", ((std::string)constraint).c_str(), xMin, xMax);
  constraint_function->SetLineStyle(1);
  constraint_function->SetLineWidth(2);
  constraint_function->SetLineColor(1);

  constraint_function->Draw("same");

  int numPoints = fitHistory.size();
  TGraph* fit_graph = new TGraph(numPoints);
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint )
  {
    const TVectorD& alpha = fitHistory.at(idxPoint).get_alpha();
    fit_graph->SetPoint(idxPoint, alpha(0), alpha(1));
  }
  fit_graph->SetMarkerStyle(8);
  fit_graph->SetMarkerSize(1);
  fit_graph->SetMarkerColor(8);
  fit_graph->SetLineStyle(7);
  fit_graph->SetLineWidth(1);
  fit_graph->SetLineColor(8);

  fit_graph->Draw("LPsame");

  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  canvas->Print(std::string(outputFileName_plot).append(".png").c_str());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").c_str());

  delete fit_graph;
  delete chi2_function;
  delete constraint_function;
  delete canvas;
}

void
showChi2(const std::vector<KinFitSummary>& fitHistory,
         const std::string& outputFileName)
{
  // CV: delete startPos (iteration = -1), 
  //     as it does not have a valid chi^2
  std::vector<KinFitSummary> fitHistory_posIterations;
  for ( const KinFitSummary& it : fitHistory )
  {
    if ( it.get_iteration() >= 0 )
    {
      fitHistory_posIterations.push_back(it);
    }
  }

  int numPoints = fitHistory_posIterations.size();
  TGraph* graph = new TGraph(numPoints);
  double yMin = 0.;
  double yMax = 1.;
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint )
  {
    double chi2 = fitHistory_posIterations.at(idxPoint).get_chi2();
    graph->SetPoint(idxPoint, fitHistory_posIterations.at(idxPoint).get_iteration(), chi2);
    if ( chi2 > yMax )
    {
      yMax = chi2;
    }
  }
  
  showGraph(800, 600, graph, "Iteration", -0.5, numPoints - 0.5, "#chi^{2}", yMin, yMax, "LPsame", outputFileName);

  delete graph;
}

void
showDistance(const TVectorD& mean, 
             const std::vector<KinFitSummary>& fitHistory,
             const std::string& outputFileName)
{
  int numPoints = fitHistory.size();
  TGraph* graph = new TGraph(numPoints);
  double yMin = 0.;
  double yMax = 1.;
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint )
  {
    const TVectorD& alpha = fitHistory.at(idxPoint).get_alpha();
    double distance = std::sqrt(square(alpha(0) - mean(0)) + square(alpha(1) - mean(1)));
    graph->SetPoint(idxPoint, fitHistory.at(idxPoint).get_iteration(), distance);
    if ( distance > yMax )
    {
      yMax = distance;
    }
  }

  showGraph(800, 600, graph, "Iteration", -1.5, numPoints - 1.5, "Distance", yMin, yMax, "LPsame", outputFileName);

  delete graph;
}

template <unsigned int rank>
void 
showResults(const TVectorD& mean, const TMatrixD& cov, const PolynomialConstraint<rank>& constraint, 
            const std::vector<KinFitSummary>& fitHistory,
            const std::string& outputFileName)
{
  showFit<rank>(mean, cov, constraint, fitHistory, outputFileName);
  showChi2(fitHistory, TString(outputFileName.c_str()).ReplaceAll(".png", "_chi2.png").Data());
  showDistance(mean, fitHistory, TString(outputFileName.c_str()).ReplaceAll(".png", "_distance.png").Data());
}

int main(int argc, char* argv[])
{
//--- throw an exception in case ROOT encounters an error
  gErrorAbortLevel = kError;

  std::cout << "<testKinematicFit>:\n";

  int verbosity = 1;

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("testKinematicFit");  

  TVectorD alpha0(2);
  alpha0(0) = 0.;
  alpha0(1) = 0.;

  TMatrixD V_alpha0_uncorr(2,2);
  V_alpha0_uncorr(0,0) = 1.;
  V_alpha0_uncorr(0,1) = 0.;
  V_alpha0_uncorr(1,0) = V_alpha0_uncorr(0,1);
  V_alpha0_uncorr(1,1) = 4.;

  TMatrixD V_alpha0_corr_wide(2,2);
  V_alpha0_corr_wide(0,0) = 1.;
  V_alpha0_corr_wide(0,1) = 1.6;
  V_alpha0_corr_wide(1,0) = V_alpha0_corr_wide(0,1);
  V_alpha0_corr_wide(1,1) = 4.;

  TMatrixD V_alpha0_corr_narrow(2,2);
  V_alpha0_corr_narrow(0,0) = 1./16;
  V_alpha0_corr_narrow(0,1) = 0.4;
  V_alpha0_corr_narrow(1,0) = V_alpha0_corr_narrow(0,1);
  V_alpha0_corr_narrow(1,1) = 4.;
  
  std::vector<point2d> points_pol1;
  points_pol1.push_back(point2d(-3., -2.));
  points_pol1.push_back(point2d(+3., +4.));  
  PolynomialConstraint<1> constraint_pol1(points_pol1, verbosity);
  PolynomialConstraintTester constraintTester_pol1(alpha0, verbosity);
  constraintTester_pol1(constraint_pol1, "testConstraint_pol1.png");

  std::vector<point2d> points_pol3;
  points_pol3.push_back(point2d(-3., -2.));
  points_pol3.push_back(point2d(-1., +2.));
  points_pol3.push_back(point2d(+1.,  0.));
  points_pol3.push_back(point2d(+3., +4.));
  PolynomialConstraint<3> constraint_pol3(points_pol3, verbosity);
  PolynomialConstraintTester constraintTester_pol3(alpha0, verbosity);
  constraintTester_pol3(constraint_pol3, "testConstraint_pol3.png");
  
  edm::ParameterSet cfg_kinFit;
  cfg_kinFit.addParameter<int>("verbosity", -1);
  KinFitAlgo kinFit(cfg_kinFit);

  TVectorD startPos_pol1 = get_startPos(points_pol1.at(0));

  TVectorD startPos_pol3 = get_startPos(points_pol3.at(0));

  std::cout << "Fitting covariance matrix:\n";
  V_alpha0_uncorr.Print();
  std::cout << "with constraint = '" << (std::string)constraint_pol1 << "'...\n";  
  std::vector<KinFitSummary> fitHistory_uncorr_pol1;
  kinFit(alpha0, V_alpha0_uncorr, constraint_pol1, startPos_pol1, &fitHistory_uncorr_pol1);
  std::string outputFileName_uncorr_pol1 = "testKinematicFit_uncorr_pol1.png";
  showResults(alpha0, V_alpha0_uncorr, constraint_pol1, fitHistory_uncorr_pol1, outputFileName_uncorr_pol1);
  std::cout << " Done.\n";

  std::cout << "Fitting covariance matrix:\n";
  V_alpha0_corr_wide.Print();
  std::cout << "with constraint = '" << (std::string)constraint_pol1 << "'...\n";
  std::vector<KinFitSummary> fitHistory_corr_wide_pol1;
  kinFit(alpha0, V_alpha0_corr_wide, constraint_pol1, startPos_pol1, &fitHistory_corr_wide_pol1);
  std::string outputFileName_corr_wide_pol1 = "testKinematicFit_corr_wide_pol1.png";
  showResults(alpha0, V_alpha0_corr_wide, constraint_pol1, fitHistory_corr_wide_pol1, outputFileName_corr_wide_pol1);
  std::cout << " Done.\n";

  std::cout << "Fitting covariance matrix:\n";
  V_alpha0_corr_narrow.Print();
  std::cout << "with constraint = '" << (std::string)constraint_pol1 << "'...\n";
  std::vector<KinFitSummary> fitHistory_corr_narrow_pol1;
  kinFit(alpha0, V_alpha0_corr_narrow, constraint_pol1, startPos_pol1, &fitHistory_corr_narrow_pol1);
  std::string outputFileName_corr_narrow_pol1 = "testKinematicFit_corr_narrow_pol1.png";
  showResults(alpha0, V_alpha0_corr_narrow, constraint_pol1, fitHistory_corr_narrow_pol1, outputFileName_corr_narrow_pol1);

  std::cout << "Fitting covariance matrix:\n";
  V_alpha0_uncorr.Print();
  std::cout << "with constraint = '" << (std::string)constraint_pol3 << "'...\n";  
  std::vector<KinFitSummary> fitHistory_uncorr_pol3;
  kinFit(alpha0, V_alpha0_uncorr, constraint_pol3, startPos_pol3, &fitHistory_uncorr_pol3);
  std::string outputFileName_uncorr_pol3 = "testKinematicFit_uncorr_pol3.png";
  showResults(alpha0, V_alpha0_uncorr, constraint_pol3, fitHistory_uncorr_pol3, outputFileName_uncorr_pol3);
  std::cout << " Done.\n";

  std::cout << "Fitting covariance matrix:\n";
  V_alpha0_corr_wide.Print();
  std::cout << "with constraint = '" << (std::string)constraint_pol3<< "'...\n";
  std::vector<KinFitSummary> fitHistory_corr_wide_pol3;
  kinFit(alpha0, V_alpha0_corr_wide, constraint_pol3, startPos_pol3, &fitHistory_corr_wide_pol3);
  std::string outputFileName_corr_wide_pol3 = "testKinematicFit_corr_wide_pol3.png";
  showResults(alpha0, V_alpha0_corr_wide, constraint_pol3, fitHistory_corr_wide_pol3, outputFileName_corr_wide_pol3);
  std::cout << " Done.\n";

  std::cout << "Fitting covariance matrix:\n";
  V_alpha0_corr_narrow.Print();
  std::cout << "with constraint = '" << (std::string)constraint_pol3 << "'...\n";
  std::vector<KinFitSummary> fitHistory_corr_narrow_pol3;
  kinFit(alpha0, V_alpha0_corr_narrow, constraint_pol3, startPos_pol3, &fitHistory_corr_narrow_pol3);
  std::string outputFileName_corr_narrow_pol3 = "testKinematicFit_corr_narrow_pol3.png";
  showResults(alpha0, V_alpha0_corr_narrow, constraint_pol3, fitHistory_corr_narrow_pol3, outputFileName_corr_narrow_pol3);

  clock.Show("testKinematicFit");

  return EXIT_SUCCESS;
}
