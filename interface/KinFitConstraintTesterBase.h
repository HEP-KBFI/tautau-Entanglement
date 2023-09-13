#ifndef TauAnalysis_Entanglement_KinFitConstraintTesterBase_h
#define TauAnalysis_Entanglement_KinFitConstraintTesterBase_h

#include "DataFormats/Math/interface/Matrix.h"                   // math::Matrix
#include "DataFormats/Math/interface/Vector.h"                   // math::Vector

#include "TauAnalysis/Entanglement/interface/fillWithOverFlow.h" // fillWithOverFlow2D()
#include "TauAnalysis/Entanglement/interface/showHistogram2d.h"  // showHistogram2d()

#include <TH2.h>                                                 // TH2D
#include <TRandom3.h>                                            // TRandom3
#include <TString.h>                                             // Form()

template <unsigned int P, unsigned int C>
class KinFitConstraintTesterBase
{
 public:
  KinFitConstraintTesterBase(const typename KinFitConstraintBase<P,C>::VectorP& alpha0, int verbosity = -1)
    : alpha0_(alpha0)
    , verbosity_(verbosity)
  {}
  ~KinFitConstraintTesterBase()
  {}

  void
  operator()(KinFitConstraintBase<P,C>& constraint, const std::string& outputFileName)
  {
    if ( verbosity_ >= 1 )
    {
      std::cout << "<KinFitConstraintTesterBase::operator()>:\n";
    }

    std::map<unsigned, std::map<unsigned, TH2*>> histograms; // keys = idxConstraint, idxParameter
    for ( unsigned idxConstraint = 0; idxConstraint < C; ++idxConstraint )
    {
      for ( unsigned idxParameter = 0; idxParameter < P; ++idxParameter )
      {
        std::string histogramName = Form("histogram_constr%u_par%u", idxConstraint, idxParameter);
        double xMin = -plotRange_(idxConstraint,idxParameter);
        double xMax = +plotRange_(idxConstraint,idxParameter);
        double yMin = -plotRange_(idxConstraint,idxParameter);
        double yMax = +plotRange_(idxConstraint,idxParameter);
        TH2* histogram = new TH2D(histogramName.c_str(), histogramName.c_str(), 100, xMin, xMax, 100, yMin, yMax);
        histograms[idxConstraint][idxParameter] = histogram;
      }
    }

    const int numToys = 10000;
    int numToys_selected = 0;
    for ( int idxToy = 0; idxToy < numToys; ++idxToy )
    {
      if ( ((idxToy % 1000) == 0 && verbosity_ >= 1) || verbosity_ >= 3 )
      {
        std::cout << " processing toy " << idxToy << "/" << numToys << " (selected = " << numToys_selected << ")\n";
      }

      typename KinFitConstraintBase<P,C>::VectorP alphaA;
      for ( unsigned idxParameter = 0; idxParameter < P; ++idxParameter )
      {
        alphaA(idxParameter) = rnd_.Gaus(rndMean_(idxParameter), rndWidth_(idxParameter));
      }

      constraint.set_alphaA(alphaA);
      if ( constraint.get_errorFlag() )
      {
	if ( verbosity_ >= 1 )
	{
	  std::cout << "skipping toy, because constraint returned an error !!\n";
	}
        continue;
      }	  
      
      typename KinFitConstraintBase<P,C>::MatrixCxP D = constraint.get_D();

      const double epsilon = 1.e-2;
      for ( unsigned int idxParameter = 0; idxParameter < P; ++idxParameter )
      {
        typename KinFitConstraintBase<P,C>::VectorP alphaAp = alphaA;
        alphaAp(idxParameter) = alphaA(idxParameter) + 0.5*epsilon;
        constraint.set_alphaA(alphaAp);
        typename KinFitConstraintBase<P,C>::VectorC d_alphaAp = constraint.get_d();

        typename KinFitConstraintBase<P,C>::VectorP alphaAm = alphaA;
        alphaAm(idxParameter) = alphaA(idxParameter) - 0.5*epsilon;
        constraint.set_alphaA(alphaAm);
        typename KinFitConstraintBase<P,C>::VectorC d_alphaAm = constraint.get_d();
            
        typename KinFitConstraintBase<P,C>::VectorC dd_dalpha = (1./epsilon)*(d_alphaAp - d_alphaAm);
        for ( unsigned int idxConstraint = 0; idxConstraint < C; ++idxConstraint )
        {
          fillWithOverFlow2D(histograms[idxConstraint][idxParameter], dd_dalpha(idxConstraint), D(idxConstraint,idxParameter));
        }
      }

      ++numToys_selected;
    }

    if ( verbosity_ >= 1 )
    {
      std::cout << "#toys: analyzed = " << numToys << ", selected = " << numToys_selected << "\n";
    }

    for ( unsigned int idxConstraint = 0; idxConstraint < C; ++idxConstraint )
    {
      for ( unsigned int idxParameter = 0; idxParameter < P; ++idxParameter )
      {
        std::string xAxisTitle = Form("#frac{#partiald}{#partialpar%u}", idxParameter);
        std::string yAxisTitle = Form("D(%u,%u)", idxConstraint, idxParameter);
        size_t idx = outputFileName.find_last_of('.');
        std::string outputFileName_plot = std::string(outputFileName, 0, idx);
        outputFileName_plot.append("_");
        outputFileName_plot.append(Form("constr%u_par%u", idxConstraint, idxParameter));
        outputFileName_plot.append(std::string(outputFileName, idx));
        showHistogram2d(800, 800, histograms[idxConstraint][idxParameter], xAxisTitle, 1.3, yAxisTitle, 1.3, true, false, "BOX", outputFileName_plot, false, verbosity_);
        delete histograms[idxConstraint][idxParameter];
      }
    }
  }

 protected:
  typename KinFitConstraintBase<P,C>::VectorP alpha0_;
  
  TRandom3 rnd_;
  typename KinFitConstraintBase<P,C>::VectorP rndMean_;
  typename KinFitConstraintBase<P,C>::VectorP rndWidth_;

  typename KinFitConstraintBase<P,C>::MatrixCxP plotRange_;

  int verbosity_;
};

#endif // TauAnalysis_Entanglement_KinFitConstraintTesterBase_h
