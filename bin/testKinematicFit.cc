
#include "DataFormats/Math/interface/Matrix.h"                 // math::Matrix
#include "DataFormats/Math/interface/Vector.h"                 // math::Vector

#include "TauAnalysis/Entanglement/interface/cmsException.h"   // cmsException
#include "TauAnalysis/Entanglement/interface/printCovMatrix.h" // printCovMatrix()
#include "TauAnalysis/Entanglement/interface/square.h"         // square()

#include <TBenchmark.h>                                        // TBenchmark
#include <TCanvas.h>                                           // TCanvas
#include <TF1.h>                                               // TF1
#include <TF2.h>                                               // TF2
#include <TGraph.h>                                            // TGraph
#include <TH1.h>                                               // TH1D
#include <TMath.h>                                             // TMath::Pi()
#include <TString.h>                                           // TString

#include <cmath>                                               // std::sqrt()
#include <sstream>                                             // std::ostringstream
#include <string>                                              // std::string
#include <utility>                                             // std::pair
#include <vector>                                              // std::vector

typedef std::pair<double, double> point2d;

typedef math::Matrix<1,1>::type Matrix1x1;
typedef math::Matrix<1,2>::type Matrix1x2;
typedef math::Matrix<2,1>::type Matrix2x1;
typedef math::Matrix<2,2>::type Matrix2x2;
typedef math::Vector<1>::type   Vector1;
typedef math::Vector<2>::type   Vector2;

template <unsigned int rank>
class Polynomial
{
 public:
  typedef typename math::Matrix<rank + 1,rank + 1>::type MatrixPxP;
  typedef typename math::Vector<rank + 1>::type VectorP;

  Polynomial(const std::vector<point2d>& points, int verbosity = -1)
    : rank_(rank)
    , points_(points)
    , verbosity_(verbosity)
  {
    if ( verbosity_ >= 1 )
    {
      std::cout << "<Polynomial::Polynomial>:\n";
      for ( size_t idxPoint = 0; idxPoint < points.size(); ++idxPoint )
      {
        const point2d& point = points.at(idxPoint);
        std::cout << "point #" << idxPoint << ": x = " << point.first << ", y = " << point.second << "\n";
      }
    }

    if ( points.size() != (rank + 1) )
      throw cmsException("Polynomial", __LINE__) 
        << "Given number of points = " << points.size() << " does not match rank + 1 = " << (rank + 1) << " !!\n";
    
    const int P = rank + 1;
    MatrixPxP M;
    VectorP y;
    for ( int idxRow = 0; idxRow < P; ++idxRow )
    {
      const point2d& point = points.at(idxRow);
      for ( int idxColumn = 0; idxColumn < P; ++idxColumn )
      {
        M(idxRow,idxColumn) = pow(point.first, P - (idxColumn + 1));
      }
      y(idxRow) = point.second;
    }

    if ( verbosity_ >= 1 )
    {
      std::cout << "M:\n";
      std::cout << M << "\n";
    }

    int Minv_errorFlag = 0;
    MatrixPxP Minv = M.Inverse(Minv_errorFlag);
    if ( Minv_errorFlag != 0 )
      throw cmsException("Polynomial", __LINE__) 
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
  ~Polynomial()
  {}

  std::string
  get_formula() const
  {
    std::ostringstream os;
    for ( unsigned int idx = 0; idx <= rank_; ++idx )
    {
      if ( idx > 0 ) os << " + ";
      os << coeff_(idx) << "*x^" << (rank_ - idx);
    }
    return os.str();
  }

  double 
  get_derrivative(double x) const
  {
    // CV: for a polynomial of rank 3,
    //     this function returns 3*coeff[0]*x^2 + 2*coeff[1]*x + coeff[2]
    double retVal = 0.;
    for ( unsigned int idx = 0; idx < rank_; ++idx )
    {
      retVal += (rank - idx)*coeff_(idx)*pow(x, rank_ - (idx + 1));
    }
    return retVal;
  }

  Matrix1x2
  get_D(const point2d& point) const
  {
    double x = point.first;
    Matrix1x2 D;
    D(0,0) = get_derrivative(x);
    D(0,1) = -1.;
    return D;
  }

  double 
  get_polynomial(double x) const
  {
    // CV: for a polynomial of rank 3,
    //     this function returns coeff[0]*x^3 + coeff[1]*x^2 + coeff[2]*x + coeff[3]
    double retVal = 0.;
    for ( unsigned int idx = 0; idx <= rank_; ++idx )
    {
      retVal += coeff_(idx)*pow(x, rank_ - idx);
    }
    return retVal;
  }

  Vector1
  get_d(const point2d& point) const
  {
    double x = point.first;
    double y = point.second;
    Vector1 d;
    d(0) = get_polynomial(x) - y;
    return d;
  }

 private:
  unsigned int rank_;

  std::vector<point2d> points_;

  VectorP coeff_;

  int verbosity_;
};

struct FitResult
{
  FitResult(int iteration, const point2d& point, double chi2)
    : iteration_(iteration)
    , point_(point)
    , chi2_(chi2)
  {}
  ~FitResult()
  {}

  int iteration_;
  point2d point_;
  double chi2_;
};

template <unsigned int rank>
std::vector<FitResult>
fit(const Vector2& alpha0, const Matrix2x2& V_alpha0, const Polynomial<rank>& constraint, const point2d& startPos, int verbosity)
{
  int Vinv_alpha0_errorFlag = 0;
  Matrix2x2 Vinv_alpha0 = V_alpha0.Inverse(Vinv_alpha0_errorFlag);
  if ( Vinv_alpha0_errorFlag != 0 ) {
    printCovMatrix("V_alpha0", V_alpha0);
    throw cmsException("fit", __LINE__) 
      << "Failed to invert matrix V_alpha0 !!\n";
  }

  const int numParameters = 2;
  const int numConstraints = 1;

  std::vector<FitResult> fitResult;
  fitResult.push_back(FitResult(-1, startPos, -1.));

  Vector2 alpha_bestfit;
  double min_chi2 = -1.;
  int status = -1;
  int iteration_bestfit = -1;
  const int max_iterations = 25;
  bool hasConverged = false;
  int iteration = 0;
  Vector2 alphaA;
  alphaA(0) = startPos.first;
  alphaA(1) = startPos.second;
  while ( !hasConverged && iteration < max_iterations )
  {
    if ( verbosity >= 1 )
    {
      std::cout << "iteration #" << iteration << ":\n";
      std::cout << "alphaA:\n";
      std::cout << alphaA << "\n";
    }

    Matrix1x2 D = constraint.get_D(point2d(alphaA(0), alphaA(1)));
    Vector1   d = constraint.get_d(point2d(alphaA(0), alphaA(1)));
    if ( verbosity >= 1 )
    {
      std::cout << "D:\n";
      std::cout << D << "\n";
      std::cout << "d:\n";
      std::cout << d << "\n";
    }

    Matrix2x1 DT = ROOT::Math::Transpose(D);

    Matrix1x1 Vinv_D = D*V_alpha0*DT;
    int V_D_errorFlag = 0;
    Matrix1x1 V_D = Vinv_D.Inverse(V_D_errorFlag);
    if ( V_D_errorFlag != 0 )
    {
      printCovMatrix("Vinv_D", Vinv_D);
      throw cmsException("fit", __LINE__) 
        << "Failed to invert matrix Vinv_D !!\n";
    }

    Vector2 dalpha0 = alpha0 - alphaA;
    Vector1 lambda = V_D*(D*dalpha0 + d);

    //const double sfStepSize = +1.;
    double sfStepSize = std::min(1., 2.*iteration/max_iterations);
    Vector2 alpha = alpha0 - sfStepSize*V_alpha0*DT*lambda;

    Matrix2x2 V_alpha = V_alpha0 - V_alpha0*DT*V_D*D*V_alpha0;
    if ( verbosity >= 1 )
    {
      printCovMatrix("V_alpha", V_alpha);
    }

    Vector2 dalpha = alpha - alphaA;
    double dalpha_mag = std::sqrt(ROOT::Math::Dot(dalpha, dalpha));
    if ( verbosity >= 1 )
    {
      std::cout << "dalpha:\n";
      std::cout << dalpha << "\n";
      std::cout << "|alpha - alphaA| = " << dalpha_mag << "\n";
    }

    Vector2 alpha_minus_alpha0 = alpha - alpha0;
    double chi2 = ROOT::Math::Dot(alpha_minus_alpha0, Vinv_alpha0*alpha_minus_alpha0) + ROOT::Math::Dot(lambda, D*dalpha + d);
    chi2 /= numParameters;
    if ( verbosity >= 1 )
    {
      std::cout << "chi^2/DoF = " << chi2 << "\n";
    }

    Vector1 residuals = D*dalpha + d;
    double residuals_sum = 0.;
    for ( int idx = 0; idx < numConstraints; ++idx )
    {
      residuals_sum += std::fabs(residuals(idx));
    }
    if ( verbosity >= 1 )
    {
      std::cout << "residuals of constraint equations:\n";
      std::cout << residuals << "\n";
      std::cout << "sum_i |residual[i]| = " << residuals_sum << "\n";
    }

    if ( verbosity >= 1 )
    {
      Vector2 pulls;
      for ( int idx = 0; idx < numParameters; ++idx )
      {      
        pulls(idx) = alpha_minus_alpha0(idx)/std::sqrt(V_alpha0(idx,idx));
      }
      std::cout << "pulls:\n";
      std::cout << pulls << "\n";
    }

    double d_mag = std::sqrt(ROOT::Math::Dot(d, d));
    if ( d_mag < 1.e-3 )
    {
      if ( status == -1 || chi2 < min_chi2 )
      {
        alpha_bestfit = alphaA;
        min_chi2 = chi2;
        status = 0;
        iteration_bestfit = iteration;
      }
      if ( status == 0 && dalpha_mag < 1.e-3 )
      {
        status = 1;
        hasConverged = true;
      }
    }

    fitResult.push_back(FitResult(iteration, point2d(alpha(0), alpha(1)), chi2));

    alphaA = alpha;

    ++iteration;
  }

  if ( verbosity >= 1 )
  {
    std::cout << "best fit:\n";
    std::cout << "iteration = " << iteration_bestfit << " (max_iterations = " << max_iterations << ")\n";
    std::cout << "alpha:\n";
    std::cout << alpha_bestfit << "\n";
    std::cout << "status = " << status << "\n";
    std::cout << "min(chi^2/DoF) = " << min_chi2 << "\n";
  }

  return fitResult;
}

template <unsigned int rank>
void 
showFit(const Vector2& mean, const Matrix2x2& cov, const Polynomial<rank>& constraint, 
        const std::vector<FitResult>& fitResult,
        const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", 900, 800);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2); 
  canvas->SetTopMargin(0.05);
  canvas->SetLeftMargin(0.08);
  canvas->SetBottomMargin(0.07);
  canvas->SetRightMargin(0.14);

  int covInv_errorFlag = 0;
  Matrix2x2 covInv = cov.Inverse(covInv_errorFlag);
  if ( covInv_errorFlag != 0 )
    throw cmsException("showFit", __LINE__) 
      << "Failed to invert matrix cov !!\n";
  std::cout << "covInv:\n";
  std::cout << covInv << "\n";

  double det = -1.;
  cov.Det2(det);
  std::cout << "det = " << det << "\n";

  const double xMin = -10.;
  const double xMax = +10.;
  const double yMin = -10.;
  const double yMax = +10.;

  std::string chi2_formula = "(x - [1])*[3]*(x - [1]) + (x - [1])*[4]*(y - [2]) + (y - [2])*[5]*(x - [1]) + (y - [2])*[6]*(y - [2])";
  TF2* chi2_function = new TF2("chi2", chi2_formula.c_str(), xMin, xMax, yMin, yMax);
  chi2_function->SetParameter(0, 1./(2.*TMath::Pi()*det));
  chi2_function->SetParameter(1, mean(0));
  chi2_function->SetParameter(2, mean(1));
  chi2_function->SetParameter(3, covInv(0,0));
  chi2_function->SetParameter(4, covInv(0,1));
  chi2_function->SetParameter(5, covInv(1,0));
  chi2_function->SetParameter(6, covInv(1,1));

  chi2_function->SetTitle("");
  chi2_function->Draw("Colz");

  TF1* constraint_function = new TF1("constraint", constraint.get_formula().c_str(), xMin, xMax);
  constraint_function->SetLineStyle(1);
  constraint_function->SetLineWidth(2);
  constraint_function->SetLineColor(1);

  constraint_function->Draw("same");

  int numPoints = fitResult.size();
  TGraph* fit_graph = new TGraph(numPoints);
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint )
  {
    const point2d& point = fitResult.at(idxPoint).point_;
    fit_graph->SetPoint(idxPoint, point.first, point.second);
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
showChi2(const std::vector<FitResult>& fitResult, 
         const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.12);

  double yMin = 0.;
  double yMax = 1.;

  // CV: delete startPos (iteration = -1), 
  //     as it does not have a valid chi^2
  std::vector<FitResult> fitResult_posIterations;
  for ( const FitResult& it : fitResult )
  {
    if ( it.iteration_ >= 0 )
    {
      fitResult_posIterations.push_back(it);
    }
  }

  int numPoints = fitResult_posIterations.size();
  TGraph* chi2_graph = new TGraph(numPoints);
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint )
  {
    double chi2 = fitResult_posIterations.at(idxPoint).chi2_;
    chi2_graph->SetPoint(idxPoint, fitResult_posIterations.at(idxPoint).iteration_, chi2);
    if ( chi2 > yMax )
    {
      yMax = chi2;
    }
  }
  chi2_graph->SetMarkerStyle(8);
  chi2_graph->SetMarkerSize(1);
  chi2_graph->SetMarkerColor(8);
  chi2_graph->SetLineStyle(8);
  chi2_graph->SetLineWidth(1);
  chi2_graph->SetLineColor(8);

  TH1* dummyHistogram = new TH1D("dummyHistogram", "dummyHistogram", numPoints, -0.5, numPoints - 0.5);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  dummyHistogram->SetMinimum(yMin);
  dummyHistogram->SetMaximum(yMax);
  dummyHistogram->GetXaxis()->SetTitle("Iteration");
  dummyHistogram->GetYaxis()->SetTitle("#chi^{2}");
  dummyHistogram->Draw("axis");

  chi2_graph->Draw("LPsame");

  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  canvas->Print(std::string(outputFileName_plot).append(".png").c_str());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").c_str());

  delete chi2_graph;
  delete dummyHistogram;
  delete canvas;
}

void
showDistance(const Vector2& mean, 
             const std::vector<FitResult>& fitResult,
             const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.12);

  double yMin = 0.;
  double yMax = 1.;

  int numPoints = fitResult.size();
  TGraph* distance_graph = new TGraph(numPoints);
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint )
  {
    const point2d& point = fitResult.at(idxPoint).point_;
    double distance = std::sqrt(square(point.first - mean(0)) + square(point.second - mean(1)));
    distance_graph->SetPoint(idxPoint, fitResult.at(idxPoint).iteration_, distance);
    if ( distance > yMax )
    {
      yMax = distance;
    }
  }
  distance_graph->SetMarkerStyle(8);
  distance_graph->SetMarkerSize(1);
  distance_graph->SetMarkerColor(8);
  distance_graph->SetLineStyle(8);
  distance_graph->SetLineWidth(1);
  distance_graph->SetLineColor(8);

  TH1* dummyHistogram = new TH1D("dummyHistogram", "dummyHistogram", numPoints, -1.5, numPoints - 1.5);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  dummyHistogram->SetMinimum(yMin);
  dummyHistogram->SetMaximum(yMax);
  dummyHistogram->GetXaxis()->SetTitle("Iteration");
  dummyHistogram->GetYaxis()->SetTitle("Distance");
  dummyHistogram->Draw("axis");

  distance_graph->Draw("LPsame");

  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  canvas->Print(std::string(outputFileName_plot).append(".png").c_str());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").c_str());

  delete distance_graph;
  delete dummyHistogram;
  delete canvas;
}

template <unsigned int rank>
void 
showResults(const Vector2& mean, const Matrix2x2& cov, const Polynomial<rank>& constraint, 
            const std::vector<FitResult>& fitResult, 
            const std::string& outputFileName)
{
  showFit(mean, cov, constraint, fitResult, outputFileName);
  showChi2(fitResult, TString(outputFileName.c_str()).ReplaceAll(".png", "_chi2.png").Data());
  showDistance(mean, fitResult, TString(outputFileName.c_str()).ReplaceAll(".png", "_distance.png").Data());
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

  std::vector<point2d> points_pol1;
  points_pol1.push_back(point2d(-3., -2.));
  points_pol1.push_back(point2d(+3., +4.));
  Polynomial<1> constraint_pol1(points_pol1, verbosity);

  std::vector<point2d> points_pol3;
  points_pol3.push_back(point2d(-3., -2.));
  points_pol3.push_back(point2d(-1., +2.));
  points_pol3.push_back(point2d(+1.,  0.));
  points_pol3.push_back(point2d(+3., +4.));
  Polynomial<3> constraint_pol3(points_pol3, verbosity);

  Vector2 mean;
  mean(0) = 0.;
  mean(1) = 0.;

  Matrix2x2 cov_uncorr;
  cov_uncorr(0,0) = 1.;
  cov_uncorr(0,1) = 0.;
  cov_uncorr(1,0) = cov_uncorr(0,1);
  cov_uncorr(1,1) = 4.;

  Matrix2x2 cov_corr_wide;
  cov_corr_wide(0,0) = 1.;
  cov_corr_wide(0,1) = 1.6;
  cov_corr_wide(1,0) = cov_corr_wide(0,1);
  cov_corr_wide(1,1) = 4.;

  Matrix2x2 cov_corr_narrow;
  cov_corr_narrow(0,0) = 1./16;
  cov_corr_narrow(0,1) = 0.4;
  cov_corr_narrow(1,0) = cov_corr_narrow(0,1);
  cov_corr_narrow(1,1) = 4.;

  std::cout << "Fitting covariance matrix:\n";
  std::cout << cov_uncorr << "\n";
  std::cout << "with constraint = '" << constraint_pol1.get_formula() << "'...\n";
  std::vector<FitResult> fitResult_uncorr_pol1 = fit(mean, cov_uncorr, constraint_pol1, points_pol1.at(0), verbosity);
  std::string outputFileName_uncorr_pol1 = "testKinematicFit_uncorr_pol1.png";
  showResults(mean, cov_uncorr, constraint_pol1, fitResult_uncorr_pol1, outputFileName_uncorr_pol1);
  std::cout << " Done.\n";

  std::cout << "Fitting covariance matrix:\n";
  std::cout << cov_corr_wide << "\n";
  std::cout << "with constraint = '" << constraint_pol1.get_formula() << "'...\n";
  std::vector<FitResult> fitResult_corr_wide_pol1 = fit(mean, cov_corr_wide, constraint_pol1, points_pol1.at(0), verbosity);
  std::string outputFileName_corr_wide_pol1 = "testKinematicFit_corr_wide_pol1.png";
  showResults(mean, cov_corr_wide, constraint_pol1, fitResult_corr_wide_pol1, outputFileName_corr_wide_pol1);
  std::cout << " Done.\n";

  std::cout << "Fitting covariance matrix:\n";
  std::cout << cov_corr_narrow << "\n";
  std::cout << "with constraint = '" << constraint_pol1.get_formula() << "'...\n";
  std::vector<FitResult> fitResult_corr_narrow_pol1 = fit(mean, cov_corr_narrow, constraint_pol1, points_pol1.at(0), verbosity);
  std::string outputFileName_corr_narrow_pol1 = "testKinematicFit_corr_narrow_pol1.png";
  showResults(mean, cov_corr_narrow, constraint_pol1, fitResult_corr_narrow_pol1, outputFileName_corr_narrow_pol1);
  std::cout << " Done.\n";

  std::cout << "Fitting covariance matrix:\n";
  std::cout << cov_uncorr << "\n";
  std::cout << "with constraint = '" << constraint_pol3.get_formula() << "'...\n";
  std::vector<FitResult> fitResult_uncorr_pol3 = fit(mean, cov_uncorr, constraint_pol3, points_pol3.at(0), verbosity);
  std::string outputFileName_uncorr_pol3 = "testKinematicFit_uncorr_pol3.png";
  showResults(mean, cov_uncorr, constraint_pol3, fitResult_uncorr_pol3, outputFileName_uncorr_pol3);
  std::cout << " Done.\n";

  std::cout << "Fitting covariance matrix:\n";
  std::cout << cov_corr_wide << "\n";
  std::cout << "with constraint = '" << constraint_pol3.get_formula() << "'...\n";
  std::vector<FitResult> fitResult_corr_wide_pol3 = fit(mean, cov_corr_wide, constraint_pol3, points_pol3.at(0), verbosity);
  std::string outputFileName_corr_wide_pol3 = "testKinematicFit_corr_wide_pol3.png";
  showResults(mean, cov_corr_wide, constraint_pol3, fitResult_corr_wide_pol3, outputFileName_corr_wide_pol3);
  std::cout << " Done.\n";

  std::cout << "Fitting covariance matrix:\n";
  std::cout << cov_corr_narrow << "\n";
  std::cout << "with constraint = '" << constraint_pol3.get_formula() << "'...\n";
  std::vector<FitResult> fitResult_corr_narrow_pol3 = fit(mean, cov_corr_narrow, constraint_pol3, points_pol3.at(0), verbosity);
  std::string outputFileName_corr_narrow_pol3 = "testKinematicFit_corr_narrow_pol3.png";
  showResults(mean, cov_corr_narrow, constraint_pol3, fitResult_corr_narrow_pol3, outputFileName_corr_narrow_pol3);
  std::cout << " Done.\n";

  clock.Show("testKinematicFit");

  return EXIT_SUCCESS;
}
