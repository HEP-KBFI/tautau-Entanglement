
#include "DataFormats/Math/interface/Matrix.h"               // math::Matrix
#include "DataFormats/Math/interface/Vector.h"               // math::Vector

#include "TauAnalysis/Entanglement/interface/cmsException.h" // cmsException

#include <TCanvas.h>                                         // TCanvas
#include <TF2.h>                                             // TF2
#include <TGraph.h>                                          // TGraph
#include <TMath.h>                                           // TMath::Pi()

#include <cmath>                                             // std::sqrt()
#include <sstream>                                           // std::ostringstream
#include <string>                                            // std::string
#include <utility>                                           // std::pair
#include <vector>                                            // std::vector

typdef std::pair<double, double> point2d;

typedef Matrix<2,2>::type Matrix2x2;
typedef Vector<2>::type   Vector2;

template <unsigned int rank>
class Polynomial
{
 public:
  Polynomial(const std::vector<point2d>& points, int verbosity)
    : rank_(rank)
    , points_(points)
    , verbosity_(verbosity)
  {
    if ( verbosity_ )
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
    typedef Matrix<P,P>::type MatrixPxP;
    MatrixPxP M;
    typedef Vector<P>::type VectorP;
    VectorP y;
    for ( int idxRow = 0; idxRow < P; ++idxRow )
    {
      const point2d& point = points.at(idxRow);
      double x = point.first;
      double y = point.second;
      for ( int idxColumn = 0; idxColumn < P; ++idxColumn )
      {
        M(idxRow,idxColumn) = pow(x, P - (idxColumn + 1));
      }
      y(idxRow) = point.second;
    }

    if ( verbosity_ )
    {
      std::cout << "M:\n";
      std::cout << M << "\n";
    }

    int Minv_errorFlag = 0;
    MatrixPxP Minv = M.Inverse(Minv_errorFlag);
    if ( Minv_errorFlag != 0 )
      throw cmsException("Polynomial", __LINE__) 
        << "Failed to invert matrix M !!\n";
    if ( verbosity >= 1 )
    {
      std::cout << "Minv:\n";
      std::cout << Minv << "\n";
    }

    coeff_ = Minv*y;
    if ( coeff_ >= 1 )
    {
      std::cout << "coeff:\n";
      std::cout << coeff << "\n";
    }
  }
  ~Polynomial()
  {}

  std::string
  get_formula()
  {
    std::ostringstream os;
    for ( int idx = 0; idx <= rank_; ++idx )
    {
      if ( idx > 0 ) retVal << " + ";
      os << coeff_(idx) << "*x^" << (rank_ - idx);
    }
    return os.str();
  }

  double 
  get_derrivative(double x)
  {
    // CV: for a polynomial of rank 3,
    //     this function returns 3*coeff[0]*x^2 + 2*coeff[1]*x + coeff[2]
    double retVal = 0.;
    for ( int idx = 0; idx < rank_; ++idx )
    {
      retVal += (rank - idx)*coeff_(idx)*pow(x, rank_ - (idx + 1));
    }
    return retVal;
  }

  typedef Matrix<1,2>::type Matrix1x2;
  Matrix1x2
  get_D(const point2d& point)
  {
    double x = point.first;
    Matrix1x2 D;
    D(0,0) = get_derrivative(x);
    D(0,1) = -1.;
    return D;
  }

  double 
  get_polynomial(double x)
  {
    // CV: for a polynomial of rank 3,
    //     this function returns coeff[0]*x^3 + coeff[1]*x^2 + coeff[2]*x + coeff[3]
    double retVal = 0.;
    for ( int idx = 0; idx <= rank_; ++idx )
    {
      retVal += coeff_(idx)*pow(x, rank_ - idx);
    }
    return retVal;
  }

  typedef Vector<1>::type Vector1;
  Vector1
  get_d(const point2d& point)
  {
    double x = point.first;
    double y = point.second;
    Vector1 d;
    d(0) = y - get_polynomial(x); 
  }

 private:
  unsigned int rank_;

  std::vector<point2d> points_;

  VectorP coeff_;

  int verbosity_;
};

template <unsigned int rank>
std::vector<point2d>&
fit(const Vector2& mean, const Matrix2x2& cov, const Polynomial<rank>& constraint, const point2d& startPos, int verbosity)
{
  int covInv_errorFlag = 0;
  MatrixPxP covInv = cov.Inverse(covInv_errorFlag);
  if ( covInv_errorFlag != 0 )
    throw cmsException("fit", __LINE__) 
      << "Failed to invert matrix cov !!\n";

  std::vector<point2d> fit_points;

  double min_chi2 = -1.;
  int status = -1;
  int iteration_bestfit = -1;
  const int max_iterations = 10;
  bool hasConverged = false;
  int iteration = 0;
  Vector2 alpha0 = mean;
  Matrix2x2 V_alpha0 = cov;
  Vector2 alphaA;
  while ( !hasConverged && iteration < max_iterations )
  {
    if ( verbosity >= 1 )
    {
      std::cout << "iteration #" << iteration << ":\n";
    }

    if ( iteration == 0 )
    {
      alphaA(0) = startPos.first;
      alphaA(1) = startPos.second;
    }

    Vector2 dy = alphaA - mean;

    dy(0) = alphaA(0) - mean(0)constraint.get_polynomial(double x)

    TMatrix2x2 Vinv_A = AT*V_alpha0*A;
  }

  return fit_points; 
}

void 
showFit(const Vector2& mean, const Matrix2x2& cov, const Polynomial<rank>& constraint, const std::vector<point2d>& fit_points, const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2); 
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.12);

  int covInv_errorFlag = 0;
  MatrixPxP covInv = cov.Inverse(covInv_errorFlag);
  if ( covInv_errorFlag != 0 )
    throw cmsException("showFit", __LINE__) 
      << "Failed to invert matrix cov !!\n";

  double det = -1.;
  cov.Det2(det);

  const double xMin = -3.;
  const double xMax = +3.;
  const double yMin = -3.;
  const double yMax = +3.;

  std::string gaussian_formula = "[0]*exp(-0.5*((x - [1])*[3]*(x - [1]) + (x - [1])*[4]*(y - [2]) + (y - [2])*[5]*(x - [1]) + (y - [2])*[6]*(y - [2])))";
  TF2* gaussian = new TF2("gaussian", gaussian_formula.c_str(), xMin, xMax, yMin, yMax);
  gaussian->SetParameter(0, 1./(2.*TMath::Pi()*det));
  gaussian->SetParameter(1, mean(0));
  gaussian->SetParameter(2, mean(1));
  gaussian->SetParameter(3, covInv(0,0));
  gaussian->SetParameter(4, covInv(0,1));
  gaussian->SetParameter(5, covInv(1,0));
  gaussian->SetParameter(6, covInv(1,1));

  gaussian->Draw("C");

  int numPoints = fit_points.size();
  TGraph* fit_graph = new TGraph(numPoints);
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint )
  {
    const point2d& point = fit_points.at(idxPoint);
    fit_graph->SetPoint(idxPoint, point.first, point.second);
  }
  fit_graph->SetMarkerStyle(8);
  fit_graph->SetMarkerSize(1);
  fit_graph->SetMarkerColor(1);
  fit_graph->SetLineStyle(8);
  fit_graph->SetLineWidth(1);
  fit_graph->SetLineColor(1);

  fit_graph->Draw("LPsame");

  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  outputFileName_plot.append("_");
  outputFileName_plot.append(histogram->GetName());
  canvas->Print(std::string(outputFileName_plot).append(".png").c_str());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").c_str());

  delete fit_graph;
  delete gaussian;
  delete canvas;
}

int main(int argc, char* argv[])
{
//--- throw an exception in case ROOT encounters an error
  gErrorAbortLevel = kError;

  std::cout << "<testKinematicFit>:\n";

  int verbosity = -1;

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("testKinematicFit");  

  std::vector<point2d> points_pol1;
  points_pol1.push_back(point2d(-3., -2.));
  points_pol1.push_back(point2d(+3., +4.));
  Polynomial<1> constraint_pol1(points_pol1);

  std::vector<point2d> points_pol3;
  points_pol3.push_back(point2d(-3., -2.));
  points_pol3.push_back(point2d(-1., +2.));
  points_pol3.push_back(point2d(+1.,  0.));
  points_pol3.push_back(point2d(+3., +4.));
  Polynomial<3> constraint_pol3(points_pol3);

  Vector2 mean;
  mean(0) = 0.;
  mean(1) = 0.;

  Matrix2x2 cov_uncorr;
  cov_uncorr(0,0) = 1.;
  cov_uncorr(0,1) = 0.;
  cov_uncorr(1,0) = 0.;
  cov_uncorr(1,1) = 4.;

  Matrix2x2 cov_corr_wide;
  cov_corr_wide(0,0) = 1.;
  cov_corr_wide(0,1) = 2.;
  cov_corr_wide(1,0) = 2.;
  cov_corr_wide(1,1) = 4.;

  Matrix2x2 cov_corr_narrow;
  cov_corr_narrow(0,0) = 1./16;
  cov_corr_narrow(0,1) = 1./std::sqrt(2.);
  cov_corr_narrow(1,0) = 1./std::sqrt(2.);
  cov_corr_narrow(1,1) = 4.;

  std::cout << "Fitting Gaussian with covariance matrix:\n";
  std::cout << cov_uncorr << "\n";
  std::cout << "with constraint = '" << constraint_pol1->get_formula() << "'...\n";
  std::vector<point2d> fit_points_uncorr_linear = fit(mean, cov_uncorr, constraint_pol1, points_pol1.at(0), verbosity);
  std::string outputFileName_uncorr_pol1 = "testKinematicFit_uncorr_pol1.png";
  showFit(mean, cov_uncorr, constraint_pol1, fit_points_uncorr_pol1, outputFileName_uncorr_pol1);
  std::cout << " Done.\n";

  std::cout << "Fitting Gaussian with covariance matrix:\n";
  std::cout << cov_corr_wide << "\n";
  std::cout << "with constraint = '" << constraint_pol1->get_formula() << "'...\n";
  std::vector<point2d> fit_points_corr_wide_linear = fit(mean, cov_corr_wide, constraint_pol1, points_pol1.at(0), verbosity);
  std::string outputFileName_corr_wide_pol1 = "testKinematicFit_corr_wide_pol1.png";
  showFit(mean, cov_corr_wide, constraint_pol1, fit_points_corr_wide_pol1, outputFileName_corr_wide_pol1);
  std::cout << " Done.\n";

  std::cout << "Fitting Gaussian with covariance matrix:\n";
  std::cout << cov_corr_narrow << "\n";
  std::cout << "with constraint = '" << constraint_pol1->get_formula() << "'...\n";
  std::vector<point2d> fit_points_corr_narrow_linear = fit(mean, cov_corr_narrow, constraint_pol1, points_pol1.at(0), verbosity);
  std::string outputFileName_corr_narrow_pol1 = "testKinematicFit_corr_narrow_pol1.png";
  showFit(mean, cov_corr_narrow, constraint_pol1, fit_points_corr_narrow_pol1, outputFileName_corr_narrow_pol1);
  std::cout << " Done.\n";

  std::cout << "Fitting Gaussian with covariance matrix:\n";
  std::cout << cov_uncorr << "\n";
  std::cout << "with constraint = '" << constraint_pol3->get_formula() << "'...\n";
  std::vector<point2d> fit_points_uncorr_linear = fit(mean, cov_uncorr, constraint_pol3, points_pol3.at(0), verbosity);
  std::string outputFileName_uncorr_pol3 = "testKinematicFit_uncorr_pol3.png";
  showFit(mean, cov_uncorr, constraint_pol3, fit_points_uncorr_pol3, outputFileName_uncorr_pol3);
  std::cout << " Done.\n";

  std::cout << "Fitting Gaussian with covariance matrix:\n";
  std::cout << cov_corr_wide << "\n";
  std::cout << "with constraint = '" << constraint_pol3->get_formula() << "'...\n";
  std::vector<point2d> fit_points_corr_wide_linear = fit(mean, cov_corr_wide, constraint_pol3, points_pol3.at(0), verbosity);
  std::string outputFileName_corr_wide_pol3 = "testKinematicFit_corr_wide_pol3.png";
  showFit(mean, cov_corr_wide, constraint_pol3, fit_points_corr_wide_pol3, outputFileName_corr_wide_pol3);
  std::cout << " Done.\n";

  std::cout << "Fitting Gaussian with covariance matrix:\n";
  std::cout << cov_corr_narrow << "\n";
  std::cout << "with constraint = '" << constraint_pol3->get_formula() << "'...\n";
  std::vector<point2d> fit_points_corr_narrow_linear = fit(mean, cov_corr_narrow, constraint_pol3, points_pol3.at(0), verbosity);
  std::string outputFileName_corr_narrow_pol3 = "testKinematicFit_corr_narrow_pol3.png";
  showFit(mean, cov_corr_narrow, constraint_pol3, fit_points_corr_narrow_pol3, outputFileName_corr_narrow_pol3);
  std::cout << " Done.\n";

  clock.Show("testKinematicFit");

  return EXIT_SUCCESS;
}
