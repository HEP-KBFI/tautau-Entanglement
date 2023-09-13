
#include "DataFormats/Math/interface/Matrix.h"                             // math::Matrix
#include "DataFormats/Math/interface/Vector.h"                             // math::Vector
#include "FWCore/ParameterSet/interface/ParameterSet.h"                    // edm::ParameterSet

#include "TauAnalysis/Entanglement/interface/cmsException.h"               // cmsException
#include "TauAnalysis/Entanglement/interface/fillWithOverFlow.h"           // fillWithOverFlow2D()
#include "TauAnalysis/Entanglement/interface/KinFitAlgo.h"                 // KinFitAlgo<>, KinFitSummary<>
#include "TauAnalysis/Entanglement/interface/KinFitConstraintBase.h"       // KinFitConstraintBase<>
#include "TauAnalysis/Entanglement/interface/KinFitConstraintTesterBase.h" // KinFitConstraintTesterBase<>
#include "TauAnalysis/Entanglement/interface/printCovMatrix.h"             // printCovMatrix()
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
#include <TRandom3.h>                                                      // TRandom3
#include <TString.h>                                                       // TString

#include <cmath>                                                           // std::isnan(), std::sqrt()
#include <sstream>                                                         // std::ostringstream
#include <string>                                                          // std::string
#include <utility>                                                         // std::pair
#include <vector>                                                          // std::vector

typedef std::pair<double, double> point2d;

typedef typename math::Matrix<1,1>::type Matrix1x1;
typedef typename math::Matrix<1,2>::type Matrix1x2;
typedef typename math::Matrix<2,1>::type Matrix2x1;
typedef typename math::Matrix<2,2>::type Matrix2x2;
typedef typename math::Vector<1>::type Vector1;
typedef typename math::Vector<2>::type Vector2;

template <unsigned int rank>
class PolynomialConstraint : public KinFitConstraintBase<2,1>
{
  typedef typename math::Matrix<rank + 1,rank + 1>::type MatrixRxR;
  typedef typename math::Vector<rank + 1>::type VectorR;

 public:
  PolynomialConstraint(const std::vector<point2d>& points, int verbosity = -1)
    : KinFitConstraintBase<2,1>(verbosity)
    , points_(points)
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
    
    MatrixRxR M;
    VectorR y;
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
      std::cout << "M:\n";
      std::cout << M << "\n";
    }

    int Minv_errorFlag = 0;
    MatrixRxR Minv = M.Inverse(Minv_errorFlag);
    if ( Minv_errorFlag != 0 )
      throw cmsException("PolynomialConstraint", __LINE__) 
        << "Failed to invert matrix M !!\n";
    if ( verbosity_ >= 1 )
    {
      std::cout << "Minv:\n";
      std::cout << Minv << "\n";
    }

    coeff_ = Minv*y;
    if ( verbosity_ >= 1 )
    {
      std::cout << "coeff:\n";
      std::cout << coeff_ << "\n";
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
  set_alphaA(const typename KinFitConstraintBase<2,1>::VectorP& alphaA)
  {
    double x = alphaA(0);
    double y = alphaA(1);
    
    D_(0,0) = get_derrivative(x);
    D_(0,1) = -1.;
    d_(0) = get_polynomial(x) - y;
    d_metric_(0,0) = 1.;
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

  using KinFitConstraintBase<2,1>::D_;
  using KinFitConstraintBase<2,1>::d_;
  using KinFitConstraintBase<2,1>::d_metric_;
  using KinFitConstraintBase<2,1>::verbosity_;

  unsigned int rank_;

  std::vector<point2d> points_;

  VectorR coeff_;
};

class PolynomialConstraintTester : public KinFitConstraintTesterBase<2,1>
{
 public:
  PolynomialConstraintTester(const Vector2& alpha0, int verbosity = -1)
    : KinFitConstraintTesterBase<2,1>(alpha0, verbosity)
  {
    rndMean_ = alpha0_;
    
    for ( unsigned int idxParameter = 0; idxParameter < 2; ++idxParameter )
    {
      rndWidth_(idxParameter) = 2.;
    }

    for ( unsigned int idxConstraint = 0; idxConstraint < 1; ++idxConstraint )
    {
      for ( unsigned int idxParameter = 0; idxParameter < 2; ++idxParameter )
      {
        plotRange_(idxConstraint, idxParameter) = 5.;
      }
    }
  }
 
 private:
  using KinFitConstraintTesterBase<2,1>::alpha0_;
  using KinFitConstraintTesterBase<2,1>::rndMean_;
  using KinFitConstraintTesterBase<2,1>::rndWidth_;
  using KinFitConstraintTesterBase<2,1>::plotRange_;
};

Vector2 
get_startPos(const point2d& point)
{
  Vector2 startPos;
  startPos(0) = point.first;
  startPos(1) = point.second;
  return startPos;
}

template <unsigned int rank>
void 
showFit(const Vector2& alpha0, const Matrix2x2& V_alpha0, const PolynomialConstraint<rank>& constraint, 
        const std::vector<KinFitSummary<2>>& fitHistory,
        const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", 900, 800);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2); 
  canvas->SetTopMargin(0.05);
  canvas->SetLeftMargin(0.08);
  canvas->SetBottomMargin(0.07);
  canvas->SetRightMargin(0.14);

  int Vinv_alpha0_errorFlag = 0;
  Matrix2x2 Vinv_alpha0 = V_alpha0.Inverse(Vinv_alpha0_errorFlag);
  if ( Vinv_alpha0_errorFlag != 0 )
    throw cmsException("showFit", __LINE__) 
      << "Failed to invert matrix V_alpha0 !!\n";
  std::cout << "Vinv_alpha0:\n";
  std::cout << Vinv_alpha0 << "\n";

  double V_alpha0_det = -1.;
  V_alpha0.Det2(V_alpha0_det);
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
    const KinFitSummary<2>::VectorP& alpha = fitHistory.at(idxPoint).get_alpha();
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
showChi2(const std::vector<KinFitSummary<2>>& fitHistory,
         const std::string& outputFileName)
{
  // CV: delete startPos (iteration = -1), 
  //     as it does not have a valid chi^2
  std::vector<KinFitSummary<2>> fitHistory_posIterations;
  for ( const KinFitSummary<2>& it : fitHistory )
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
showDistance(const Vector2& mean, 
             const std::vector<KinFitSummary<2>>& fitHistory,
             const std::string& outputFileName)
{
  int numPoints = fitHistory.size();
  TGraph* graph = new TGraph(numPoints);
  double yMin = 0.;
  double yMax = 1.;
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint )
  {
    const KinFitSummary<2>::VectorP& alpha = fitHistory.at(idxPoint).get_alpha();
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
showResults(const Vector2& mean, const Matrix2x2& cov, const PolynomialConstraint<rank>& constraint, 
            const std::vector<KinFitSummary<2>>& fitHistory,
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

  Vector2 alpha0;
  alpha0(0) = 0.;
  alpha0(1) = 0.;

  Matrix2x2 V_alpha0_uncorr;
  V_alpha0_uncorr(0,0) = 1.;
  V_alpha0_uncorr(0,1) = 0.;
  V_alpha0_uncorr(1,0) = V_alpha0_uncorr(0,1);
  V_alpha0_uncorr(1,1) = 4.;

  Matrix2x2 V_alpha0_corr_wide;
  V_alpha0_corr_wide(0,0) = 1.;
  V_alpha0_corr_wide(0,1) = 1.6;
  V_alpha0_corr_wide(1,0) = V_alpha0_corr_wide(0,1);
  V_alpha0_corr_wide(1,1) = 4.;

  Matrix2x2 V_alpha0_corr_narrow;
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
  KinFitAlgo<2,1> kinFit(cfg_kinFit);

  Vector2 startPos_pol1 = get_startPos(points_pol1.at(0));

  Vector2 startPos_pol3 = get_startPos(points_pol3.at(0));

  std::cout << "Fitting covariance matrix:\n";
  std::cout << V_alpha0_uncorr << "\n";
  std::cout << "with constraint = '" << (std::string)constraint_pol1 << "'...\n";  
  std::vector<KinFitSummary<2>> fitHistory_uncorr_pol1;
  kinFit(alpha0, V_alpha0_uncorr, constraint_pol1, startPos_pol1, &fitHistory_uncorr_pol1);
  std::string outputFileName_uncorr_pol1 = "testKinematicFit_uncorr_pol1.png";
  showResults(alpha0, V_alpha0_uncorr, constraint_pol1, fitHistory_uncorr_pol1, outputFileName_uncorr_pol1);
  std::cout << " Done.\n";

  std::cout << "Fitting covariance matrix:\n";
  std::cout << V_alpha0_corr_wide << "\n";
  std::cout << "with constraint = '" << (std::string)constraint_pol1 << "'...\n";
  std::vector<KinFitSummary<2>> fitHistory_corr_wide_pol1;
  kinFit(alpha0, V_alpha0_corr_wide, constraint_pol1, startPos_pol1, &fitHistory_corr_wide_pol1);
  std::string outputFileName_corr_wide_pol1 = "testKinematicFit_corr_wide_pol1.png";
  showResults(alpha0, V_alpha0_corr_wide, constraint_pol1, fitHistory_corr_wide_pol1, outputFileName_corr_wide_pol1);
  std::cout << " Done.\n";

  std::cout << "Fitting covariance matrix:\n";
  std::cout << V_alpha0_corr_narrow << "\n";
  std::cout << "with constraint = '" << (std::string)constraint_pol1 << "'...\n";
  std::vector<KinFitSummary<2>> fitHistory_corr_narrow_pol1;
  kinFit(alpha0, V_alpha0_corr_narrow, constraint_pol1, startPos_pol1, &fitHistory_corr_narrow_pol1);
  std::string outputFileName_corr_narrow_pol1 = "testKinematicFit_corr_narrow_pol1.png";
  showResults(alpha0, V_alpha0_corr_narrow, constraint_pol1, fitHistory_corr_narrow_pol1, outputFileName_corr_narrow_pol1);

  std::cout << "Fitting covariance matrix:\n";
  std::cout << V_alpha0_uncorr << "\n";
  std::cout << "with constraint = '" << (std::string)constraint_pol3 << "'...\n";  
  std::vector<KinFitSummary<2>> fitHistory_uncorr_pol3;
  kinFit(alpha0, V_alpha0_uncorr, constraint_pol3, startPos_pol3, &fitHistory_uncorr_pol3);
  std::string outputFileName_uncorr_pol3 = "testKinematicFit_uncorr_pol3.png";
  showResults(alpha0, V_alpha0_uncorr, constraint_pol3, fitHistory_uncorr_pol3, outputFileName_uncorr_pol3);
  std::cout << " Done.\n";

  std::cout << "Fitting covariance matrix:\n";
  std::cout << V_alpha0_corr_wide << "\n";
  std::cout << "with constraint = '" << (std::string)constraint_pol3<< "'...\n";
  std::vector<KinFitSummary<2>> fitHistory_corr_wide_pol3;
  kinFit(alpha0, V_alpha0_corr_wide, constraint_pol3, startPos_pol3, &fitHistory_corr_wide_pol3);
  std::string outputFileName_corr_wide_pol3 = "testKinematicFit_corr_wide_pol3.png";
  showResults(alpha0, V_alpha0_corr_wide, constraint_pol3, fitHistory_corr_wide_pol3, outputFileName_corr_wide_pol3);
  std::cout << " Done.\n";

  std::cout << "Fitting covariance matrix:\n";
  std::cout << V_alpha0_corr_narrow << "\n";
  std::cout << "with constraint = '" << (std::string)constraint_pol3 << "'...\n";
  std::vector<KinFitSummary<2>> fitHistory_corr_narrow_pol3;
  kinFit(alpha0, V_alpha0_corr_narrow, constraint_pol3, startPos_pol3, &fitHistory_corr_narrow_pol3);
  std::string outputFileName_corr_narrow_pol3 = "testKinematicFit_corr_narrow_pol3.png";
  showResults(alpha0, V_alpha0_corr_narrow, constraint_pol3, fitHistory_corr_narrow_pol3, outputFileName_corr_narrow_pol3);

  clock.Show("testKinematicFit");

  return EXIT_SUCCESS;
}
